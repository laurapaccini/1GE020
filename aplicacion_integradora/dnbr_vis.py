"""
dnbr_vis.py

Lectura y visualización de los GeoTIFF mensuales de dNBR generados por
dnbr_mensual.py: panel de mapas por mes y animación (GIF) de la evolución
temporal.

Dependencias: numpy, rasterio, xarray, matplotlib (incluye pillow para el GIF).

Uso típico en un notebook
-------------------------
import dnbr_vis as dv
da = dv.leer_dnbr("salidas_dnbr/dnbr_bolivia_beni_2019.tif")
fig = dv.panel_mensual(da)                    # 12 mapas con escala común
ruta = dv.animar(da, "bolivia_beni_2019.gif")  # GIF con un cuadro por mes

# Varios años en una sola serie temporal
rutas = [f"salidas_dnbr/dnbr_bolivia_beni_{y}.tif" for y in (2019, 2020, 2021)]
da = dv.leer_dnbr(rutas)

# Mostrar el GIF en el notebook
from IPython.display import Image
Image(filename=ruta)
"""

import os
import re

import numpy as np

NODATA = -32768     # mismo valor de "sin dato" que usa dnbr_mensual.py


def _tiempos(src, ruta):
    """Obtiene la fecha (mes) de cada banda de un GeoTIFF abierto con rasterio.

    Parámetros
    ----------
    src : rasterio.io.DatasetReader
        Archivo abierto. Si las bandas tienen descripción dNBR_AAAAMM, se usa
        esa información.
    ruta : str
        Ruta del archivo. Si las bandas no tienen descripción, se busca el año
        al final del nombre (por ejemplo, dnbr_ucayali_2020.tif) y se asume que
        las bandas son enero a diciembre.

    Devuelve
    --------
    numpy.ndarray de datetime64[ns]
        Una fecha (primer día del mes) por banda. Lanza ValueError si no se puede determinar.
    """
    fechas = []
    for desc in src.descriptions:
        m = re.fullmatch(r"dNBR_(\d{4})(\d{2})", desc or "")
        if not m:
            fechas = []
            break
        fechas.append(np.datetime64(f"{m.group(1)}-{m.group(2)}", "M"))
    if len(fechas) == src.count:
        return np.array(fechas).astype("datetime64[ns]")
    m = re.search(r"_(\d{4})\.tif$", os.fspath(ruta))
    if m and src.count == 12:
        return np.array([np.datetime64(f"{m.group(1)}-{i:02d}", "M")
                         for i in range(1, 13)]).astype("datetime64[ns]")
    raise ValueError(f"No se pudo determinar el mes de las bandas de {ruta}")


def leer_dnbr(rutas, escala=1000.0):
    """Lee uno o varios GeoTIFF mensuales de dNBR como un DataArray de xarray.

    Parámetros
    ----------
    rutas : str o lista de str
        Ruta de un GeoTIFF, o lista de rutas (por ejemplo, un archivo por año
        de la misma caja). Todos deben tener la misma grilla.
    escala : float, opcional
        Factor por el que se dividen los valores enteros almacenados. 1000 por
        defecto, que devuelve el dNBR en sus unidades originales.

    Devuelve
    --------
    xarray.DataArray
        Variable "dNBR" con dimensiones (tiempo, lat, lon), tipo flotante.
        Los píxeles sin dato son NaN. La coordenada "tiempo" es mensual y está
        ordenada. Lanza ValueError si las grillas de los archivos no coinciden
        o si no se puede determinar el mes de las bandas.
    """
    import rasterio
    import xarray as xr

    if isinstance(rutas, (str, os.PathLike)):
        rutas = [rutas]
    piezas = []
    for ruta in rutas:
        with rasterio.open(ruta) as src:
            datos = src.read().astype("float32")
            nodata = src.nodata if src.nodata is not None else NODATA
            datos[datos == nodata] = np.nan
            datos /= escala
            t = src.transform
            lon = t.c + t.a * (np.arange(src.width) + 0.5)
            lat = t.f + t.e * (np.arange(src.height) + 0.5)
            tiempo = _tiempos(src, ruta)
        piezas.append(xr.DataArray(
            datos, dims=("tiempo", "lat", "lon"), name="dNBR",
            coords={"tiempo": tiempo, "lat": lat, "lon": lon}))
    if len(piezas) == 1:
        da = piezas[0]
    else:
        try:
            da = xr.concat(piezas, dim="tiempo", join="exact")
        except ValueError as error:
            raise ValueError("Las grillas de los archivos no coinciden") from error
    return da.sortby("tiempo")


def _extension(da):
    """Calcula la extensión (en grados) de un DataArray para usarla con imshow.

    Parámetros
    ----------
    da : xarray.DataArray
        Arreglo con coordenadas lat (descendente) y lon (ascendente), tal como
        lo devuelve leer_dnbr.

    Devuelve
    --------
    tuple de 4 float
        (oeste, este, sur, norte), contando el borde de los píxeles extremos.
    """
    dx = abs(float(da.lon[1] - da.lon[0]))
    dy = abs(float(da.lat[1] - da.lat[0]))
    return (float(da.lon[0]) - dx / 2, float(da.lon[-1]) + dx / 2,
            float(da.lat[-1]) - dy / 2, float(da.lat[0]) + dy / 2)


def _limite_color(da, vlim):
    """Define el valor máximo (en valor absoluto) de la escala de color.

    Parámetros
    ----------
    da : xarray.DataArray
        Datos de dNBR.
    vlim : float o None
        Límite de la escala. Si es None, se usa el percentil 99 del valor
        absoluto de los datos válidos.

    Devuelve
    --------
    float
        Límite positivo; la escala va de -límite a +límite.
    """
    if vlim is not None:
        return float(vlim)
    return float(np.nanpercentile(np.abs(da.values), 99))


def _mapa_color(cmap):
    """Crea un mapa de color cuyo valor sin dato (NaN) se dibuja en gris claro.

    Parámetros
    ----------
    cmap : str
        Nombre de un mapa de color de matplotlib.

    Devuelve
    --------
    matplotlib.colors.Colormap
        Copia del mapa de color con el color "bad" en gris claro.
    """
    import matplotlib

    mapa = matplotlib.colormaps[cmap].copy()
    mapa.set_bad("0.85")
    return mapa


def panel_mensual(da, vlim=None, cmap="RdBu_r", ncols=4, titulo=None, vmin=None):
    """Dibuja un panel con un mapa por mes y una escala de color común.

    Con el mapa de color por defecto, los valores positivos (caída del NBR
    respecto a la línea base, posible quema) salen en rojo y los negativos en
    azul. Los píxeles sin dato salen en gris.

    Parámetros
    ----------
    da : xarray.DataArray
        Datos devueltos por leer_dnbr, con dimensiones (tiempo, lat, lon).
    vlim : float, opcional
        Límite de la escala de color: va de -vlim a +vlim. Si no se indica, se
        calcula con el percentil 99 de los valores absolutos.
    cmap : str, opcional
        Nombre del mapa de color de matplotlib ("RdBu_r" por defecto).
    ncols : int, opcional
        Número de columnas del panel (4 por defecto).
    titulo : str, opcional
        Título general de la figura.
    vmin : float, opcional
        Límite inferior de la escala. Si no se indica, la escala es simétrica
        (-vlim a +vlim). Para campos que no toman valores negativos, como una
        fracción quemada, conviene vmin=0.

    Devuelve
    --------
    matplotlib.figure.Figure
        Figura con un mapa por cada elemento de la dimensión tiempo.
    """
    import matplotlib.pyplot as plt

    n = da.sizes["tiempo"]
    nrows = int(np.ceil(n / ncols))
    lim = _limite_color(da, vlim)
    inf = -lim if vmin is None else vmin
    mapa = _mapa_color(cmap)
    fig, axes = plt.subplots(nrows, ncols, sharex=True, sharey=True,
                             figsize=(3.0 * ncols, 2.8 * nrows),
                             constrained_layout=True, squeeze=False)
    ext = _extension(da)
    im = None
    for k, ax in enumerate(axes.ravel()):
        if k >= n:
            ax.set_visible(False)
            continue
        im = ax.imshow(da.isel(tiempo=k).values, extent=ext, origin="upper",
                       cmap=mapa, vmin=inf, vmax=lim, interpolation="nearest")
        ax.set_title(str(da.tiempo.values[k])[:7], fontsize=9)
        ax.tick_params(labelsize=7)
    fig.colorbar(im, ax=axes, shrink=0.8, label=da.name or "valor")
    fig.supxlabel("Longitud (°)", fontsize=9)
    fig.supylabel("Latitud (°)", fontsize=9)
    if titulo:
        fig.suptitle(titulo)
    return fig


def animar(da, ruta_gif, vlim=None, cmap="RdBu_r", fps=2, titulo=None, vmin=None):
    """Genera una animación (GIF) con un cuadro por cada mes.

    Parámetros
    ----------
    da : xarray.DataArray
        Datos devueltos por leer_dnbr, con dimensiones (tiempo, lat, lon).
        Para varios años, leer los archivos juntos con leer_dnbr(lista_de_rutas).
    ruta_gif : str
        Ruta del archivo .gif de salida.
    vlim : float, opcional
        Límite de la escala de color: va de -vlim a +vlim. Si no se indica, se
        calcula con el percentil 99 de los valores absolutos.
    cmap : str, opcional
        Nombre del mapa de color de matplotlib ("RdBu_r" por defecto).
    fps : float, opcional
        Cuadros por segundo (2 por defecto, adecuado para series mensuales).
    titulo : str, opcional
        Texto que antecede a la fecha en el título de cada cuadro.
    vmin : float, opcional
        Límite inferior de la escala (simétrica por defecto); ver panel_mensual.

    Devuelve
    --------
    str
        Ruta del GIF generado (la misma de `ruta_gif`).
    """
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    lim = _limite_color(da, vlim)
    inf = -lim if vmin is None else vmin
    mapa = _mapa_color(cmap)
    datos = da.values
    fechas = [str(t)[:7] for t in da.tiempo.values]
    prefijo = f"{titulo}  " if titulo else ""

    fig, ax = plt.subplots(figsize=(5.5, 5))
    im = ax.imshow(datos[0], extent=_extension(da), origin="upper", cmap=mapa,
                   vmin=inf, vmax=lim, interpolation="nearest")
    fig.colorbar(im, ax=ax, shrink=0.8, label=da.name or "valor")
    ax.set_xlabel("Longitud (°)")
    ax.set_ylabel("Latitud (°)")
    ax.set_title(prefijo + fechas[0])

    def actualizar(i):
        im.set_data(datos[i])
        ax.set_title(prefijo + fechas[i])
        return (im,)

    animacion = FuncAnimation(fig, actualizar, frames=len(fechas), blit=False)
    animacion.save(ruta_gif, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return ruta_gif
