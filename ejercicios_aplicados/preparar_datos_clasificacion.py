"""
Prepara los datos para la clasificación de cobertura de suelo (semana 11).

Para cada región genera, dentro de datos_clasificacion/:
  - s2_bandas_<año>_<region>.tif          6 bandas de Sentinel-2 L2A (B02, B03, B04, B08, B11, B12) en
                                          reflectancia, con nubes y sombras enmascaradas (NaN) a partir de la
                                          banda SCL, en la misma grilla que MapBiomas (EPSG:4326, ~30 m)
  - mapbiomas_col4_<año>_<region>.tif     Recorte del mapa de cobertura de MapBiomas Perú, Colección 4
y registra en metadata_escenas_clasificacion.csv qué escena se usó en cada región.

Elección de la escena: entre las candidatas de la época seca que cubren toda la caja, se toman las N_EVALUAR
de menor nubosidad de escena y se queda la que tiene menor porcentaje de píxeles enmascarados dentro de
la caja (según SCL). A igualdad, la de menor nubosidad de escena.

Las bandas de 10 m (B02, B03, B04, B08) y de 20 m (B11, B12) se llevan a la grilla de MapBiomas con el
promedio de los píxeles que caen en cada celda de 30 m.

Uso, desde la raíz del repo y con el environment 1geo20:
    python scripts/preparar_datos_clasificacion.py
"""
import csv
import os
from pathlib import Path

import numpy as np
import planetary_computer
import pystac_client
import rasterio
import rioxarray
import xarray as xr
from rasterio.enums import Resampling

# Configuración de GDAL para lectura remota de COG
os.environ["GDAL_DISABLE_READDIR_ON_OPEN"] = "EMPTY_DIR"
os.environ["GDAL_HTTP_MERGE_CONSECUTIVE_RANGES"] = "YES"
os.environ["VSI_CACHE"] = "TRUE"

# ----- Parámetros -----
# Cada región es [lon_min, lat_min, lon_max, lat_max]
AOIS = {
    # "puerto_maldonado": [-69.30, -12.70, -69.10, -12.50],   # demo del lunes
    "la_pampa": [-70.50, -13.05, -70.30, -12.85],           # práctica del martes
    "puerto_inca": [-74.98, -9.48, -74.78, -9.28]
}
ANIO = 2025
EPOCA = ("07-01", "09-30")          # época seca (mes-día de inicio y fin)
NUBES_MAX = 20                      # nubosidad máxima de la escena (%)
N_EVALUAR = 8                       # cuántas candidatas se evalúan con SCL dentro de la caja
BANDAS = ["B02", "B03", "B04", "B08", "B11", "B12"]

MAPBIOMAS = "/Users/lpaccini/Documents/PUCP_clases/geocomputacion/peru_coverage-col4_{anio}2.tif"
MAPBIOMAS = (
    "https://storage.googleapis.com/mapbiomas-public/initiatives/peru/"
    "collection4/lulc/coverage/peru_coverage/peru_coverage-col4_{anio}.tif"
)                                   # también puede ser una ruta local con {anio}

STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
SALIDA = Path("datos_clasificacion")
# Clases de la banda SCL que se enmascaran: 0 sin dato, 1 saturado o defectuoso,
# 3 sombra de nube, 8 y 9 nube (probabilidad media y alta), 10 cirro delgado
CLASES_ENMASCARADAS = [0, 1, 3, 8, 9, 10]


# ----- Selección de escena -----
def cubre(item, bbox):
    """True si la huella de la escena cubre toda la caja."""
    return (item.bbox[0] <= bbox[0] and item.bbox[1] <= bbox[1]
            and item.bbox[2] >= bbox[2] and item.bbox[3] >= bbox[3])


def buscar_candidatas(catalogo, bbox, anio):
    """Escenas de la época seca del año, con poca nubosidad, que cubren toda la caja."""
    items = catalogo.search(
        collections=["sentinel-2-l2a"],
        bbox=bbox,
        datetime=f"{anio}-{EPOCA[0]}/{anio}-{EPOCA[1]}",
        query={"eo:cloud_cover": {"lt": NUBES_MAX}},
    ).item_collection()
    return [i for i in items if cubre(i, bbox)]


def porcentaje_enmascarado(scl):
    """Porcentaje de píxeles de SCL que se enmascaran (nubes, sombras, sin dato)."""
    return 100 * float(scl.isin(CLASES_ENMASCARADAS).mean())


def evaluar_candidatas(candidatas, bbox):
    """
    Evalúa las N_EVALUAR candidatas de menor nubosidad de escena y devuelve una lista de
    (porcentaje enmascarado en la caja, nubosidad de escena, escena firmada), de mejor a peor.
    """
    if not candidatas:
        raise RuntimeError("Sin escenas candidatas; amplía EPOCA o sube NUBES_MAX")
    mejores = sorted(candidatas, key=lambda i: i.properties["eo:cloud_cover"])[:N_EVALUAR]
    evaluadas = []
    for item in mejores:
        firmada = planetary_computer.sign(item)
        evaluadas.append((porcentaje_enmascarado(leer_scl(firmada, bbox)),
                          item.properties["eo:cloud_cover"], firmada))
    return sorted(evaluadas, key=lambda t: (t[0], t[1]))


# ----- Lectura -----
def leer_banda(escena, clave, bbox):
    """Lee un asset recortado a la caja y lo convierte a reflectancia (escala 0.0001 y offset según la línea base)."""
    da = (
        rioxarray.open_rasterio(escena.assets[clave].href, masked=True)
        .squeeze("band", drop=True)
        .rio.clip_box(*bbox, crs="EPSG:4326")
        .load()
    )
    linea_base = escena.properties.get("s2:processing_baseline")
    if linea_base is None:
        raise RuntimeError(f"Sin línea base en las propiedades de {escena.id}")
    offset = -0.1 if float(linea_base) >= 4.0 else 0.0
    return (da * 0.0001 + offset).clip(min=0)


def leer_scl(escena, bbox):
    """Lee la banda de clasificación de escena (SCL) recortada a la caja."""
    return (
        rioxarray.open_rasterio(escena.assets["SCL"].href)
        .squeeze("band", drop=True)
        .rio.clip_box(*bbox, crs="EPSG:4326")
        .load()
    )


def recortar_mapbiomas(fuente, bbox):
    """Recorte de MapBiomas a la caja. Se conserva el tipo de dato y el nodata del original."""
    return (
        rioxarray.open_rasterio(fuente)
        .squeeze("band", drop=True)
        .rio.clip_box(*bbox, crs="EPSG:4326")
        .load()
    )


# ----- Alineación -----
def alinear_y_apilar(bandas, scl, referencia):
    """
    Lleva cada banda y la máscara de SCL a la grilla de referencia (MapBiomas) y apila las bandas
    en el orden de BANDAS. Los píxeles de nubes, sombras o fuera de la escena quedan en NaN.
    """
    mascara = scl.rio.reproject_match(referencia, resampling=Resampling.nearest).isin(CLASES_ENMASCARADAS)
    alineadas = []
    for nombre in BANDAS:
        banda = bandas[nombre].rio.write_nodata(np.nan)
        banda = banda.rio.reproject_match(referencia, resampling=Resampling.average)
        alineadas.append(banda.where(~mascara).astype("float32"))
    pila = xr.concat(alineadas, dim="band").assign_coords(band=list(BANDAS))
    pila = pila.rio.write_crs(referencia.rio.crs)
    pila.rio.write_nodata(np.nan, inplace=True)
    return pila, 100 * float(mascara.mean())


# ----- Comprobaciones y escritura -----
def misma_grilla(a, b):
    return (a.shape[-2:] == b.shape[-2:]
            and np.allclose(a.x.values, b.x.values)
            and np.allclose(a.y.values, b.y.values))


def resumen_cobertura(mb, anio):
    """Imprime las clases presentes y su número de píxeles (para completar la agrupación del notebook)."""
    valores, cuentas = np.unique(mb.values, return_counts=True)
    total = cuentas.sum()
    print(f"  clases de MapBiomas presentes en {anio} (código: píxeles, %):")
    for v, c in zip(valores, cuentas):
        print(f"    {int(v):>3}: {int(c):>8}  {100 * c / total:5.1f} %")


def guardar(da, ruta, descripciones=None):
    da.rio.to_raster(ruta, compress="LZW", tiled=True)
    if descripciones:
        with rasterio.open(ruta, "r+") as dst:
            dst.descriptions = tuple(descripciones)
    print(f"  guardado: {ruta}  forma {da.shape}  CRS {da.rio.crs}  tipo {da.dtype}")


def main():
    SALIDA.mkdir(parents=True, exist_ok=True)
    catalogo = pystac_client.Client.open(STAC_URL)

    filas = []
    for nombre, bbox in AOIS.items():
        print(f"\n{nombre}")
        candidatas = buscar_candidatas(catalogo, bbox, ANIO)
        print(f"  escenas candidatas en {ANIO}: {len(candidatas)}")
        evaluadas = evaluar_candidatas(candidatas, bbox)

        print("  mejores escenas (% enmascarado en la caja, nubosidad de escena):")
        for pct_caja, nubes, escena in evaluadas[:3]:
            print(f"    {escena.datetime.date()}  {pct_caja:.1f} %  {nubes:.1f} %  {escena.id}")
        pct_caja, nubes, escena = evaluadas[0]

        mb = recortar_mapbiomas(MAPBIOMAS.format(anio=ANIO), bbox)
        bandas = {b: leer_banda(escena, b, bbox) for b in BANDAS}
        scl = leer_scl(escena, bbox)
        pila, pct_mascara = alinear_y_apilar(bandas, scl, mb)

        if not misma_grilla(pila, mb):
            raise RuntimeError(f"Las bandas y MapBiomas no comparten grilla en {nombre}")
        pct_nan = 100 * float(np.isnan(pila.values).any(axis=0).mean())
        print(f"  escena elegida: {escena.id}  {escena.datetime.date()}")
        print(f"    enmascarado en la caja (SCL a la grilla de MapBiomas): {pct_mascara:.1f} %   "
              f"píxeles con NaN en alguna banda: {pct_nan:.1f} %")

        # guardar(pila, SALIDA / f"s2_bandas_{ANIO}_{nombre}.tif", descripciones=BANDAS)
        # guardar(mb, SALIDA / f"mapbiomas_col4_{ANIO}_{nombre}.tif")
        resumen_cobertura(mb, ANIO)

        filas.append({
            "aoi": nombre,
            "anio": ANIO,
            "bbox": bbox,
            "escena": escena.id,
            "fecha": str(escena.datetime.date()),
            "nubosidad_escena_pct": escena.properties["eo:cloud_cover"],
            "pixeles_enmascarados_caja_pct": round(pct_caja, 1),
            "linea_base": escena.properties.get("s2:processing_baseline"),
            "bandas": " ".join(BANDAS),
        })

    with open(SALIDA / "metadata_escenas_clasificacion.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(filas[0]))
        writer.writeheader()
        writer.writerows(filas)
    print(f"\nListo. Archivos en {SALIDA}")


if __name__ == "__main__":
    main()
