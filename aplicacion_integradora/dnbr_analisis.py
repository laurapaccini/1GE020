"""
dnbr_analisis.py

Funciones para pasar del dNBR mensual de Sentinel-2 (salida de dnbr_mensual.py,
leída con dnbr_vis.leer_dnbr) a una fracción mensual de área quemada en la
grilla de 0.25 grados del modelo Sparky, y para leer y alinear el modelo.

Secuencia típica
----------------
anom   = anomalia_espacial(da)               # quita la señal estacional común
flag   = marcar_quema(anom, umbral=0.25)     # 1 quemado, 0 no quemado, NaN sin dato
nueva  = quema_nueva(flag, ventana=3)        # cuenta cada cicatriz una sola vez
obs    = a_grilla_modelo(nueva)              # fracción por celda de 0.25 grados
modelo = leer_modelo("datos_modelo", [2019, 2020, 2021], bbox)
obs, modelo = alinear(obs["BAF_obs"], modelo)

Dependencias: numpy, xarray y un lector de NetCDF (netCDF4 o h5netcdf).
"""

import glob
import os

import numpy as np
import xarray as xr

CELDA_MODELO = 0.25   # tamaño de celda del modelo, en grados


def anomalia_espacial(da):
    """Resta a cada mes la mediana espacial de la caja.

    El dNBR de todo el paisaje sube y baja con la estación seca y la lluviosa
    (la línea base de tres meses no es estacional). Esa señal común se estima
    con la mediana de los píxeles de cada mes y se resta, de modo que quedan
    resaltados los parches que se desvían de su entorno. Supone que menos de la
    mitad de la caja está quemada en un mes.

    Parámetros
    ----------
    da : xarray.DataArray
        dNBR con dimensiones (tiempo, lat, lon), salida de dnbr_vis.leer_dnbr.

    Devuelve
    --------
    xarray.DataArray
        Anomalía del dNBR, con la misma forma que `da`. Los píxeles sin dato
        siguen siendo NaN.
    """
    fondo = da.median(dim=("lat", "lon"), skipna=True)
    return (da - fondo).rename("anomalia_dNBR")


def marcar_quema(anomalia, umbral):
    """Clasifica cada píxel como quemado o no quemado con un umbral.

    Parámetros
    ----------
    anomalia : xarray.DataArray
        Anomalía del dNBR (salida de anomalia_espacial).
    umbral : float
        Valor desde el cual un píxel se considera quemado. Debe calibrarse
        mirando el histograma y los mapas; no hay un valor universal.

    Devuelve
    --------
    xarray.DataArray
        Arreglo de tipo flotante con 1.0 donde la anomalía es mayor o igual al
        umbral, 0.0 donde es menor, y NaN donde no hay dato.
    """
    bandera = (anomalia >= umbral).astype("float32")
    return bandera.where(anomalia.notnull()).rename("quemado")


def quema_nueva(bandera, ventana=3):
    """Conserva solo el primer mes en que un píxel aparece quemado.

    Una cicatriz sigue marcada en los meses siguientes al incendio. Para
    obtener área quemada por mes, un píxel se cuenta solo si no estuvo marcado
    en ninguno de los `ventana` meses anteriores. Los primeros meses de la serie
    no tienen historia previa, por lo que pueden incluir cicatrices anteriores.

    Parámetros
    ----------
    bandera : xarray.DataArray
        Salida de marcar_quema, con dimensiones (tiempo, lat, lon).
    ventana : int, opcional
        Número de meses previos que se revisan (3 por defecto).

    Devuelve
    --------
    xarray.DataArray
        Arreglo con 1.0 en el primer mes de cada evento, 0.0 en los demás meses
        y NaN donde no hay dato.
    """
    valores = bandera.values
    marcado = np.nan_to_num(valores, nan=0.0) > 0
    nuevo = np.zeros_like(marcado)
    for t in range(marcado.shape[0]):
        previo = marcado[max(0, t - ventana):t].any(axis=0) if t > 0 else False
        nuevo[t] = marcado[t] & ~previo
    salida = np.where(np.isnan(valores), np.nan, nuevo.astype("float32"))
    return bandera.copy(data=salida).rename("quemado_nuevo")


def a_grilla_modelo(campo, minimo_valido=0.5):
    """Promedia un campo de píxeles a la grilla de 0.25 grados del modelo.

    Cada celda recibe la fracción de píxeles marcados, calculada sobre los
    píxeles con dato. Si la fracción de píxeles con dato de una celda es menor
    que `minimo_valido`, la celda queda en NaN.

    Parámetros
    ----------
    campo : xarray.DataArray
        Campo de 0.0, 1.0 y NaN (por ejemplo, la salida de quema_nueva) con
        dimensiones (tiempo, lat, lon). Los bordes de la caja deben ser
        múltiplos de 0.25 grados.
    minimo_valido : float, opcional
        Fracción mínima de píxeles con dato para aceptar una celda (0.5 por
        defecto).

    Devuelve
    --------
    xarray.Dataset
        Variables "BAF_obs" (fracción quemada observada, entre 0 y 1) y
        "frac_valida" (fracción de píxeles con dato), ambas con dimensiones
        (tiempo, lat, lon) en la grilla de 0.25 grados.
    """
    res = abs(float(campo.lon[1] - campo.lon[0]))
    n = int(round(CELDA_MODELO / res))
    if abs(n * res - CELDA_MODELO) > 1e-9:
        raise ValueError("La resolución del campo no divide a 0.25 grados")
    bloques = campo.coarsen(lat=n, lon=n, boundary="exact")
    fraccion = bloques.mean()
    valida = campo.notnull().astype("float32").coarsen(lat=n, lon=n, boundary="exact").mean()
    fraccion = fraccion.where(valida >= minimo_valido)
    return xr.Dataset({"BAF_obs": fraccion, "frac_valida": valida})


def leer_modelo(carpeta, anios, bbox, variable="BAF_pred"):
    """Lee los archivos anuales del modelo y recorta la caja de interés.

    Se buscan archivos cuyo nombre empiece por BAF_pred_<año> dentro de la
    carpeta (por ejemplo, BAF_pred_2019.area-subset.0.-65.-40.-85.nc).

    Parámetros
    ----------
    carpeta : str
        Carpeta que contiene los archivos NetCDF del modelo.
    anios : list de int
        Años a leer.
    bbox : tuple de 4 float
        (oeste, sur, este, norte) en grados. Se conservan las celdas cuyo
        centro cae dentro de la caja.
    variable : str, opcional
        Nombre de la variable del modelo ("BAF_pred" por defecto).

    Devuelve
    --------
    xarray.DataArray
        Fracción quemada del modelo con dimensiones (tiempo, lat, lon), con la
        latitud en orden descendente y el atributo "product_type" con los tipos
        de producto de los archivos leídos. Lanza FileNotFoundError si falta
        algún año.
    """
    oeste, sur, este, norte = bbox
    piezas, tipos = [], set()
    for anio in anios:
        archivos = sorted(glob.glob(os.path.join(carpeta, f"BAF_pred_{anio}*.nc")))
        if not archivos:
            raise FileNotFoundError(f"No hay archivo del modelo para {anio} en {carpeta}")
        with xr.open_dataset(archivos[0]) as ds:
            tipos.add(str(ds.attrs.get("product_type", "desconocido")))
            sub = ds[variable].sel(lat=slice(norte, sur), lon=slice(oeste, este)).load()
        piezas.append(sub)
    da = xr.concat(piezas, dim="time").rename({"time": "tiempo"})
    da = da.assign_coords(tiempo=da.tiempo.astype("datetime64[ns]"))
    da.attrs["product_type"] = ", ".join(sorted(tipos))
    return da.rename("BAF_modelo")


def alinear(obs, modelo):
    """Lleva la observación y el modelo a las mismas celdas y los mismos meses.

    Parámetros
    ----------
    obs : xarray.DataArray
        Fracción quemada observada en la grilla de 0.25 grados (BAF_obs).
    modelo : xarray.DataArray
        Fracción quemada del modelo (salida de leer_modelo).

    Devuelve
    --------
    tuple (xarray.DataArray, xarray.DataArray)
        (obs, modelo) con coordenadas idénticas en tiempo, lat y lon. Se usa el
        centro de celda más cercano (tolerancia de 0.001 grados) y la
        intersección de los meses disponibles. Lanza ValueError si no hay
        meses o celdas en común.
    """
    obs = obs.assign_coords(tiempo=obs.tiempo.astype("datetime64[ns]"))
    modelo = modelo.sel(lat=obs.lat, lon=obs.lon, method="nearest", tolerance=1e-3)
    modelo = modelo.assign_coords(lat=obs.lat, lon=obs.lon)
    comunes = np.intersect1d(obs.tiempo.values, modelo.tiempo.values)
    if comunes.size == 0:
        raise ValueError("Observación y modelo no tienen meses en común")
    return obs.sel(tiempo=comunes), modelo.sel(tiempo=comunes)


def leer_mapbiomas(ruta, bbox, anio, res=CELDA_MODELO):
    """Calcula la fracción quemada mensual de MapBiomas Fuego en la grilla del modelo.

    Se espera el mapa anual de MapBiomas Fuego "monthly burned", de una sola
    banda, en el que cada píxel guarda el mes de la quema (1 a 12) y 0 indica
    que no se detectó quema. Para cada mes se promedia el indicador de quema
    de los píxeles de 30 m dentro de cada celda de `res` grados.

    Parámetros
    ----------
    ruta : str
        Ruta del GeoTIFF anual de MapBiomas Fuego (EPSG:4326). Debe cubrir la
        caja completa.
    bbox : tuple de 4 float
        (oeste, sur, este, norte) en grados, con bordes múltiplos de `res`.
    anio : int
        Año del mapa, usado para construir la coordenada de tiempo.
    res : float, opcional
        Tamaño de celda de salida en grados (0.25 por defecto, el del modelo).

    Devuelve
    --------
    xarray.DataArray
        Variable "BAF_mapbiomas" con dimensiones (tiempo, lat, lon): fracción
        de cada celda quemada en cada mes del año (entre 0 y 1). La grilla y la
        coordenada de tiempo coinciden con las de a_grilla_modelo y leer_modelo.
    """
    import rasterio
    from rasterio.transform import from_origin
    from rasterio.warp import Resampling, reproject
    from rasterio.windows import Window

    oeste, sur, este, norte = bbox
    with rasterio.open(ruta) as src:
        f0, c0 = src.index(oeste, norte)
        f1, c1 = src.index(este, sur)
        if f0 < 0 or c0 < 0 or f1 > src.height or c1 > src.width:
            raise ValueError("La caja no está contenida en el mapa de MapBiomas")
        ventana = Window(c0, f0, c1 - c0, f1 - f0)
        mes_quema = src.read(1, window=ventana)
        origen = src.window_transform(ventana)

    n_col = int(round((este - oeste) / res))
    n_fil = int(round((norte - sur) / res))
    destino = from_origin(oeste, norte, res, res)
    meses = []
    for mes in range(1, 13):
        quemado = (mes_quema == mes).astype("float32")
        fraccion = np.zeros((n_fil, n_col), dtype="float32")
        reproject(quemado, fraccion, src_transform=origen, src_crs="EPSG:4326",
                  dst_transform=destino, dst_crs="EPSG:4326",
                  resampling=Resampling.average)
        meses.append(fraccion)

    lat = norte - res * (np.arange(n_fil) + 0.5)
    lon = oeste + res * (np.arange(n_col) + 0.5)
    tiempo = np.array([np.datetime64(f"{anio}-{m:02d}-01", "ns") for m in range(1, 13)])
    return xr.DataArray(np.stack(meses), dims=("tiempo", "lat", "lon"),
                        coords={"tiempo": tiempo, "lat": lat, "lon": lon},
                        name="BAF_mapbiomas")
