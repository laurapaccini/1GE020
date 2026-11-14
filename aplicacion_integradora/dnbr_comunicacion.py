"""
dnbr_comunicacion.py

Productos para comunicar los resultados de la aplicación integradora:
mapas interactivos con folium (ubicación, dNBR sobre imagen satelital y
fracción quemada por celda de 0.25 grados) y gráficos interactivos con plotly
(series regionales y mapa animado).

Entradas
--------
- dNBR o anomalía de dNBR a escala de píxel: xarray.DataArray con dimensiones
  (tiempo, lat, lon), como lo devuelven dnbr_vis.leer_dnbr y
  dnbr_analisis.anomalia_espacial.
- Fracción quemada en la grilla de 0.25 grados: xarray.Dataset con variables
  como BAF_obs, BAF_modelo y BAF_mapbiomas, por ejemplo el archivo
  datos_areas_quemadas_pc6.nc que guarda integrador_dnbr_modelo.ipynb.

Uso típico en un notebook
-------------------------
import xarray as xr
import dnbr_vis as dv
import dnbr_comunicacion as dc

ds = xr.open_dataset("datos_areas_quemadas_pc6.nc")
m = dc.mapa_ubicacion((-67.0, -14.0, -66.0, -13.0), nombre="Beni")
m = dc.mapa_celdas(ds, "BAF_obs", "2019-09")
fig = dc.grafico_series(dc.serie_regional(ds))
fig = dc.mapa_animado(ds["BAF_obs"], vmin=0, cmap="YlOrRd")

# Guardar para compartir o presentar
m.save("mapa_celdas.html")
fig.write_html("series.html")

Dependencias: numpy, pandas, xarray, matplotlib, folium y plotly.
"""

import numpy as np

# Imagen satelital de Esri, usada como mapa base opcional en folium
TESELAS_SATELITE = ("https://server.arcgisonline.com/ArcGIS/rest/services/"
                    "World_Imagery/MapServer/tile/{z}/{y}/{x}")
ATRIBUCION_SATELITE = "Esri, Maxar, Earthstar Geographics"

# Nombres legibles de las variables, para leyendas y etiquetas
ETIQUETAS = {
    "BAF_obs": "Observación (dNBR)",
    "BAF_modelo": "Modelo Sparky",
    "BAF_mapbiomas": "MapBiomas Fuego",
}


def _mapa_base(centro, zoom, satelite=True):
    """Crea un mapa de folium con OpenStreetMap y, opcionalmente, imagen satelital.

    Parámetros
    ----------
    centro : tuple de 2 float
        (lat, lon) del centro del mapa.
    zoom : int
        Nivel de acercamiento inicial.
    satelite : bool, opcional
        Si es True, agrega la imagen satelital como capa base seleccionable.

    Devuelve
    --------
    folium.Map
        Mapa vacío con las capas base.
    """
    import folium

    m = folium.Map(location=centro, zoom_start=zoom, tiles="OpenStreetMap")
    if satelite:
        folium.TileLayer(tiles=TESELAS_SATELITE, attr=ATRIBUCION_SATELITE,
                         name="Imagen satelital").add_to(m)
    return m


def _centro(bbox):
    """Devuelve el centro (lat, lon) de una caja (oeste, sur, este, norte)."""
    oeste, sur, este, norte = bbox
    return ((sur + norte) / 2, (oeste + este) / 2)


def _caja_de(da):
    """Devuelve la caja (oeste, sur, este, norte) que cubren los píxeles de `da`.

    Usa los centros de píxel de las coordenadas lat y lon y les suma medio
    píxel a cada lado.
    """
    dlon = abs(float(da.lon[1] - da.lon[0]))
    dlat = abs(float(da.lat[1] - da.lat[0]))
    return (float(da.lon.min()) - dlon / 2, float(da.lat.min()) - dlat / 2,
            float(da.lon.max()) + dlon / 2, float(da.lat.max()) + dlat / 2)


def mapa_ubicacion(bbox, nombre=None, zoom=6, satelite=True):
    """Mapa interactivo que ubica la caja de estudio.

    Sirve para dar contexto: dónde está la región antes de mostrar resultados.

    Parámetros
    ----------
    bbox : tuple de 4 float
        (oeste, sur, este, norte) en grados.
    nombre : str, opcional
        Texto que aparece al pasar el cursor sobre la caja.
    zoom : int, opcional
        Nivel de acercamiento inicial (6 por defecto).
    satelite : bool, opcional
        Si es True (por defecto), agrega la imagen satelital como capa base.

    Devuelve
    --------
    folium.Map
        Mapa con la caja dibujada y un selector de capas. Se muestra en el
        notebook al dejarlo como última línea de la celda, y se guarda con
        m.save("archivo.html").
    """
    import folium

    oeste, sur, este, norte = bbox
    m = _mapa_base(_centro(bbox), zoom, satelite)
    folium.Rectangle(bounds=[[sur, oeste], [norte, este]], color="#1F2D4A",
                     weight=2, fill=False, tooltip=nombre).add_to(m)
    folium.LayerControl().add_to(m)
    return m


def mapa_raster(da, mes, vlim=None, cmap="RdBu_r", vmin=None, opacidad=0.7,
                zoom=9, titulo=None):
    """Superpone un mes de un campo de píxeles (dNBR o anomalía) sobre un mapa.

    Permite comparar las zonas de dNBR alto con la imagen satelital del lugar.
    Los píxeles sin dato quedan transparentes.

    Parámetros
    ----------
    da : xarray.DataArray
        Campo con dimensiones (tiempo, lat, lon), como la salida de
        dnbr_vis.leer_dnbr o dnbr_analisis.anomalia_espacial.
    mes : str
        Mes a mostrar, en formato "AAAA-MM".
    vlim : float, opcional
        Límite superior de la escala de color. Si no se indica, se usa el
        percentil 99 de los valores absolutos del mes.
    cmap : str, opcional
        Mapa de color de matplotlib ("RdBu_r" por defecto).
    vmin : float, opcional
        Límite inferior de la escala. Si no se indica, la escala es simétrica
        (-vlim a vlim).
    opacidad : float, opcional
        Opacidad de la capa, entre 0 y 1 (0.7 por defecto).
    zoom : int, opcional
        Nivel de acercamiento inicial (9 por defecto).
    titulo : str, opcional
        Texto de la leyenda. Si no se indica, se usa el nombre de `da` y el mes.

    Devuelve
    --------
    folium.Map
        Mapa con la capa del mes, su leyenda y un selector de capas.
    """
    import branca.colormap as bcm
    import folium
    import matplotlib
    from matplotlib import colors

    campo = da.sel(tiempo=mes).squeeze()
    if "tiempo" in campo.dims:
        raise ValueError(f"El mes {mes} corresponde a más de un paso de tiempo")
    valores = campo.values
    if float(campo.lat[0]) < float(campo.lat[-1]):      # norte arriba
        valores = valores[::-1]

    sup = float(np.nanpercentile(np.abs(valores), 99)) if vlim is None else vlim
    inf = -sup if vmin is None else vmin
    mapa_color = matplotlib.colormaps[cmap]
    norma = colors.Normalize(vmin=inf, vmax=sup)
    rgba = mapa_color(norma(valores))
    rgba[np.isnan(valores), 3] = 0.0                      # sin dato, transparente

    oeste, sur, este, norte = _caja_de(campo)
    m = _mapa_base(_centro((oeste, sur, este, norte)), zoom, satelite=True)
    nombre_capa = f"{da.name or 'valor'} {mes}"
    folium.raster_layers.ImageOverlay(
        image=rgba, bounds=[[sur, oeste], [norte, este]], opacity=opacidad,
        mercator_project=True, name=nombre_capa).add_to(m)

    muestras = [colors.to_hex(mapa_color(x)) for x in np.linspace(0, 1, 9)]
    leyenda = bcm.LinearColormap(muestras, vmin=inf, vmax=sup)
    leyenda.caption = titulo or nombre_capa
    leyenda.add_to(m)
    folium.LayerControl().add_to(m)
    return m


def celdas_geojson(ds, mes, variables=None, decimales=3):
    """Convierte las celdas de 0.25 grados de un mes en una colección GeoJSON.

    Cada celda es un polígono cuadrado con los valores de las variables como
    propiedades. Las celdas sin dato guardan None.

    Parámetros
    ----------
    ds : xarray.Dataset
        Datos en la grilla de 0.25 grados, con dimensiones (tiempo, lat, lon).
    mes : str
        Mes a convertir, en formato "AAAA-MM".
    variables : list de str, opcional
        Variables a incluir. Por defecto, todas las de ETIQUETAS presentes en `ds`.
    decimales : int, opcional
        Decimales con que se redondean los valores (3 por defecto).

    Devuelve
    --------
    dict
        Colección GeoJSON (FeatureCollection) lista para folium.GeoJson, con las
        propiedades "lat", "lon" y una por variable.
    """
    if variables is None:
        variables = [v for v in ETIQUETAS if v in ds]
    sub = ds[variables].sel(tiempo=mes).squeeze()
    res_lat = abs(float(sub.lat[1] - sub.lat[0]))
    res_lon = abs(float(sub.lon[1] - sub.lon[0]))

    elementos = []
    for lat in sub.lat.values:
        for lon in sub.lon.values:
            o, e = lon - res_lon / 2, lon + res_lon / 2
            s, n = lat - res_lat / 2, lat + res_lat / 2
            props = {"lat": round(float(lat), 3), "lon": round(float(lon), 3)}
            for v in variables:
                x = float(sub[v].sel(lat=lat, lon=lon))
                props[v] = None if np.isnan(x) else round(x, decimales)
            elementos.append({
                "type": "Feature",
                "geometry": {"type": "Polygon",
                             "coordinates": [[[o, s], [e, s], [e, n], [o, n], [o, s]]]},
                "properties": props,
            })
    return {"type": "FeatureCollection", "features": elementos}


def mapa_celdas(ds, variable, mes, vmax=None, cmap="YlOrRd", zoom=8):
    """Mapa interactivo de la fracción quemada por celda de 0.25 grados.

    Las celdas se colorean según `variable`. Al pasar el cursor se muestran los
    valores de todas las fuentes de esa celda (observación, modelo, MapBiomas),
    lo que permite comparar celda a celda.

    Parámetros
    ----------
    ds : xarray.Dataset
        Datos en la grilla de 0.25 grados, con dimensiones (tiempo, lat, lon).
    variable : str
        Variable que define el color (por ejemplo "BAF_obs").
    mes : str
        Mes a mostrar, en formato "AAAA-MM".
    vmax : float, opcional
        Límite superior de la escala de color. Si no se indica, se usa el
        máximo del mes entre todas las variables de ETIQUETAS presentes, para
        que los mapas de distintas fuentes sean comparables.
    cmap : str, opcional
        Mapa de color de matplotlib ("YlOrRd" por defecto).
    zoom : int, opcional
        Nivel de acercamiento inicial (8 por defecto).

    Devuelve
    --------
    folium.Map
        Mapa con las celdas, su leyenda, la información al pasar el cursor y un
        selector de capas.
    """
    import branca.colormap as bcm
    import folium
    import matplotlib
    from matplotlib import colors

    variables = [v for v in ETIQUETAS if v in ds]
    if variable not in variables:
        variables.append(variable)
    geo = celdas_geojson(ds, mes, variables)

    if vmax is None:
        vmax = float(ds[variables].sel(tiempo=mes).to_array().max())
        vmax = vmax if vmax > 0 else 1.0
    mapa_color = matplotlib.colormaps[cmap]
    muestras = [colors.to_hex(mapa_color(x)) for x in np.linspace(0, 1, 9)]
    escala = bcm.LinearColormap(muestras, vmin=0, vmax=vmax)
    escala.caption = f"{ETIQUETAS.get(variable, variable)}, {mes} (fracción de la celda)"

    def estilo(elemento):
        valor = elemento["properties"][variable]
        if valor is None:
            return {"fillOpacity": 0.0, "color": "#888888", "weight": 0.5}
        return {"fillColor": escala(min(valor, vmax)), "fillOpacity": 0.75,
                "color": "#555555", "weight": 0.5}

    oeste = min(c[0] for f in geo["features"] for c in f["geometry"]["coordinates"][0])
    este = max(c[0] for f in geo["features"] for c in f["geometry"]["coordinates"][0])
    sur = min(c[1] for f in geo["features"] for c in f["geometry"]["coordinates"][0])
    norte = max(c[1] for f in geo["features"] for c in f["geometry"]["coordinates"][0])

    m = _mapa_base(_centro((oeste, sur, este, norte)), zoom, satelite=True)
    campos = ["lat", "lon"] + variables
    alias = ["Lat", "Lon"] + [ETIQUETAS.get(v, v) for v in variables]
    folium.GeoJson(geo, name=ETIQUETAS.get(variable, variable), style_function=estilo,
                   tooltip=folium.GeoJsonTooltip(fields=campos, aliases=alias)).add_to(m)
    escala.add_to(m)
    folium.LayerControl().add_to(m)
    return m


def serie_regional(ds, variables=None):
    """Promedia cada variable sobre las celdas de la caja, mes a mes.

    Es un promedio simple de las celdas con dato. En una caja de pocos grados
    cerca del ecuador las celdas de 0.25 grados tienen áreas parecidas; en
    cajas grandes o de latitudes altas convendría ponderar por el área.

    Parámetros
    ----------
    ds : xarray.Dataset
        Datos en la grilla de 0.25 grados, con dimensiones (tiempo, lat, lon).
    variables : list de str, opcional
        Variables a promediar. Por defecto, todas las de ETIQUETAS presentes en `ds`.

    Devuelve
    --------
    pandas.DataFrame
        Una fila por mes (índice "tiempo") y una columna por variable, con la
        fracción quemada media de la caja.
    """
    if variables is None:
        variables = [v for v in ETIQUETAS if v in ds]
    return ds[variables].mean(dim=("lat", "lon"), skipna=True).to_dataframe()[variables]


def grafico_series(df, titulo=None, etiqueta_y="Fracción quemada media de la caja"):
    """Gráfico interactivo de líneas, una por columna de `df`.

    Al pasar el cursor se ven los valores de todas las series en ese mes. Las
    series se pueden ocultar o mostrar haciendo clic en la leyenda.

    Parámetros
    ----------
    df : pandas.DataFrame
        Series con el tiempo como índice (por ejemplo, la salida de
        serie_regional).
    titulo : str, opcional
        Título del gráfico.
    etiqueta_y : str, opcional
        Texto del eje vertical.

    Devuelve
    --------
    plotly.graph_objects.Figure
        Figura que se muestra en el notebook y se guarda con
        fig.write_html("archivo.html").
    """
    import plotly.graph_objects as go

    fig = go.Figure()
    for col in df.columns:
        fig.add_trace(go.Scatter(x=df.index, y=df[col], mode="lines+markers",
                                 name=ETIQUETAS.get(col, col)))
    fig.update_layout(title=titulo, xaxis_title="Mes", yaxis_title=etiqueta_y,
                      hovermode="x unified", template="simple_white",
                      legend=dict(orientation="h", y=-0.2))
    return fig


def mapa_animado(da, vlim=None, vmin=None, cmap="RdBu_r", paso=1, titulo=None,
                 ancho=700, alto=650):
    """Mapa interactivo con un deslizador de meses y botón de reproducción.

    Es la versión interactiva de dnbr_vis.animar: además de reproducir la
    secuencia, permite detenerse en un mes y leer el valor de cada píxel.

    Parámetros
    ----------
    da : xarray.DataArray
        Campo con dimensiones (tiempo, lat, lon).
    vlim : float, opcional
        Límite superior de la escala de color. Si no se indica, se usa el
        percentil 99 de los valores absolutos.
    vmin : float, opcional
        Límite inferior de la escala. Si no se indica, la escala es simétrica
        (-vlim a vlim); para una fracción quemada conviene vmin=0.
    cmap : str, opcional
        Nombre del mapa de color ("RdBu_r" por defecto). Se aceptan los nombres
        de plotly, que incluyen los habituales de matplotlib.
    paso : int, opcional
        Factor de agregación espacial (1 por defecto, sin agregar). Para el
        dNBR a escala de píxel conviene paso=10 o más: cada píxel de salida es
        el promedio de paso x paso píxeles, y el archivo resultante pesa mucho
        menos.
    titulo : str, opcional
        Título del gráfico.
    ancho, alto : int, opcional
        Tamaño de la figura en píxeles (700 x 650 por defecto). Los grados de
        latitud y longitud se dibujan con la misma escala, por lo que el mapa
        se ajusta dentro de ese tamaño sin deformarse.

    Devuelve
    --------
    plotly.graph_objects.Figure
        Figura con un cuadro por mes, con el norte arriba.
    """
    import plotly.express as px

    if paso > 1:
        da = da.coarsen(lat=paso, lon=paso, boundary="trim").mean()
    sup = float(np.nanpercentile(np.abs(da.values), 99)) if vlim is None else vlim
    inf = -sup if vmin is None else vmin
    da = da.assign_coords(tiempo=[str(t)[:7] for t in da.tiempo.values])

    # Con coordenadas numéricas, plotly ubica cada fila según su latitud;
    # origin="lower" deja el eje creciendo hacia arriba (norte arriba).
    fig = px.imshow(da, animation_frame="tiempo", color_continuous_scale=cmap,
                    zmin=inf, zmax=sup, origin="lower", aspect="equal",
                    labels={"color": da.name or "valor", "lat": "Latitud (°)",
                            "lon": "Longitud (°)", "animation_frame": "Mes"})
    fig.update_layout(title=titulo, template="simple_white", width=ancho, height=alto)
    return fig
