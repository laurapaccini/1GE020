#!/usr/bin/env python
"""
dnbr_mensual.py

Calcula, con Google Earth Engine, el dNBR mensual de un año sobre una caja
(bbox) y descarga el resultado como un GeoTIFF multibanda (una banda por mes).

Criterio (heurístico, no corresponde a un método publicado)
-----------------------------------------------------------
1. Sentinel-2 L2A armonizado (COPERNICUS/S2_SR_HARMONIZED) con máscara de
   nubes, sombras, cirros, nieve y agua a partir de la banda SCL.
2. NBR por imagen = (B8A - B12) / (B8A + B12) y mediana mensual.
3. dNBR del mes m = mediana de los NBR mensuales de los `--meses-previos`
   meses anteriores menos el NBR del mes m. Un valor positivo indica caída
   del NBR (posible quema nueva).
4. El resultado se agrega por promedio a una grilla regular en grados, cuyos
   bordes coinciden con la grilla de 0.25 grados del modelo. No se aplica
   ningún umbral aquí; el umbral se calibra después, en Python.

Salida
------
<salida>/dnbr_<nombre>_<año>.tif
  - 12 bandas (enero a diciembre), enteros de 16 bits.
  - Valor almacenado = dNBR x 1000.
  - Sin dato = -32768 (mes sin observaciones válidas, o sin línea base),
    declarado como nodata en el archivo.
  - CRS EPSG:4326.
  - Cada banda lleva como descripción su nombre (dNBR_AAAAMM).

Descarga por teselas
--------------------
Earth Engine limita la memoria de las descargas directas, de modo que una caja
grande en una sola petición puede fallar con "User memory limit exceeded". Por
eso el script divide la caja en teselas de `--tesela` grados (0.25 por
defecto), descarga cada una por separado y las ensambla en el GeoTIFF final.
Las teselas ya descargadas se conservan si el proceso se interrumpe y se
reutilizan al repetir la llamada. Para ensamblar se requiere rasterio.

Uso (ejemplo, una caja de 1 x 1 grado alineada a la grilla de 0.25 grados)
--------------------------------------------------------------------------
python dnbr_mensual.py --year 2020 --nombre ucayali \
    --bbox -74.5 -9.5 -73.5 -8.5 --res 0.001 --project MI-PROYECTO-EE

Bucle por años (Colab o Git Bash)
---------------------------------
for y in 2019 2020 2021; do
    python dnbr_mensual.py --year $y --nombre ucayali \
        --bbox -74.5 -9.5 -73.5 -8.5 --project MI-PROYECTO-EE
done

Uso en Colab (o en cualquier notebook)
--------------------------------------
1. Autenticar en una celda del notebook:
       import ee
       ee.Authenticate(auth_mode='notebook')
2. Disponer del script en la sesión: celda con %%writefile dnbr_mensual.py,
   clonar el repositorio del curso o copiarlo desde Drive.
3. Ejecutar desde el notebook, en una sola línea:
       !python dnbr_mensual.py --year 2020 --nombre ucayali --bbox -74.5 -9.5 -73.5 -8.5 --project MI-PROYECTO-EE
   o con el bucle por años en una celda %%bash (ver arriba).
4. Alternativa sin línea de comandos, importando el módulo:
       import dnbr_mensual as dn
       dn.procesar(2020, 'ucayali', (-74.5, -9.5, -73.5, -8.5), project='MI-PROYECTO-EE')
5. Los archivos quedan en el disco temporal de Colab; copiarlos a Drive o
   descargarlos antes de cerrar la sesión.

Notas
-----
- Requiere haber autenticado Earth Engine antes y un proyecto de Google Cloud
  habilitado para Earth Engine.
- Script no probado contra el servicio real. Primera comprobación: que el
  aviso final confirme que la extensión del tif coincide con la caja. Si no
  coincide, revisar el orden de crs_transform.
- Si getDownloadURL rechaza el tamaño de una tesela, aumentar --res.
- Si una tesela falla por memoria, reducir --tesela (debe ser múltiplo de --res).
"""

import argparse
import math
import os
import shutil
import sys

import ee
import requests

CELDA_MODELO = 0.25          # tamaño de celda del modelo, en grados
# Valor de "sin dato" del GeoTIFF. Es el mínimo de un entero de 16 bits con signo
# (numpy.iinfo(numpy.int16).min). Es una convención propia de este script, no un
# estándar de Earth Engine. Queda fuera del rango posible del dato: el NBR va de
# -1 a 1, el dNBR de -2 a 2 y el valor almacenado (x1000) de -2000 a 2000.
# El GeoTIFF final lo declara como nodata; las teselas intermedias no.
NODATA = -32768
LIMITE_BYTES = 45e6          # margen bajo el límite de getDownloadURL (~48 MB)
M_POR_GRADO = 111_320.0      # aproximación en el ecuador
RES_NATIVA_M = 20.0          # B8A y B12 tienen 20 m


def parsear_argumentos():
    """Lee los argumentos de la línea de comandos.

    Parámetros
    ----------
    Ninguno. Los valores se toman de sys.argv.

    Devuelve
    --------
    argparse.Namespace
        Objeto con los atributos year, nombre, bbox, res, meses_previos,
        nubes_max, tesela, salida y project, equivalentes a los parámetros de
        procesar().
    """
    p = argparse.ArgumentParser(description="dNBR mensual con Earth Engine")
    p.add_argument("--year", type=int, required=True, help="Año a procesar")
    p.add_argument("--nombre", required=True, help="Nombre corto de la región")
    p.add_argument("--bbox", type=float, nargs=4, required=True,
                   metavar=("OESTE", "SUR", "ESTE", "NORTE"),
                   help="Caja en grados; bordes en múltiplos de 0.25")
    p.add_argument("--res", type=float, default=0.001,
                   help="Resolución de salida en grados (0.25 debe ser múltiplo)")
    p.add_argument("--meses-previos", type=int, default=3,
                   help="Meses usados para la línea base")
    p.add_argument("--nubes-max", type=float, default=70.0,
                   help="Porcentaje máximo de nubes por imagen (prefiltro)")
    p.add_argument("--tesela", type=float, default=0.25,
                   help="Lado de las teselas de descarga en grados (múltiplo de --res)")
    p.add_argument("--salida", default="salidas_dnbr", help="Carpeta de salida")
    p.add_argument("--project", default=os.environ.get("EE_PROJECT"),
                   help="Proyecto de Earth Engine (o variable EE_PROJECT)")
    return p.parse_args()


def validar_caja(bbox, res):
    """Comprueba que la caja y la resolución encajen con la grilla de 0.25 grados.

    Parámetros
    ----------
    bbox : tuple de 4 float
        (oeste, sur, este, norte) en grados. Cada borde debe ser múltiplo de 0.25.
    res : float
        Tamaño de píxel de salida en grados. 0.25 debe ser múltiplo de res.

    Devuelve
    --------
    None
        No devuelve nada. Lanza ValueError si algún borde no es múltiplo de 0.25,
        si la caja está invertida o si res no divide exactamente a 0.25.
    """
    oeste, sur, este, norte = bbox
    for v in bbox:
        if abs(v / CELDA_MODELO - round(v / CELDA_MODELO)) > 1e-9:
            raise ValueError(f"El borde {v} no es múltiplo de {CELDA_MODELO}")
    if oeste >= este or sur >= norte:
        raise ValueError("bbox debe cumplir OESTE < ESTE y SUR < NORTE")
    pix_por_celda = CELDA_MODELO / res
    if abs(pix_por_celda - round(pix_por_celda)) > 1e-9:
        raise ValueError(f"{CELDA_MODELO} no es múltiplo de --res={res}")


def estimar_bytes(bbox, res, n_bandas=12):
    """Estima el tamaño de la descarga y las dimensiones de la imagen.

    Parámetros
    ----------
    bbox : tuple de 4 float
        (oeste, sur, este, norte) en grados.
    res : float
        Tamaño de píxel de salida en grados.
    n_bandas : int, opcional
        Número de bandas de la imagen (12 meses por defecto).

    Devuelve
    --------
    tuple (int, int, int)
        (bytes estimados, ancho en píxeles, alto en píxeles). Se asumen 2 bytes
        por valor (entero de 16 bits) y no se considera la compresión.
    """
    oeste, sur, este, norte = bbox
    ancho = round((este - oeste) / res)
    alto = round((norte - sur) / res)
    return ancho * alto * n_bandas * 2, ancho, alto


def mascara_scl(img):
    """Enmascara con la banda SCL los píxeles no utilizables de una imagen Sentinel-2.

    Se enmascaran las clases 0 (sin dato), 1 (saturado), 3 (sombra de nube),
    6 (agua), 8 y 9 (nube de probabilidad media y alta), 10 (cirro) y 11 (nieve).

    Parámetros
    ----------
    img : ee.Image
        Imagen de COPERNICUS/S2_SR_HARMONIZED, que incluye la banda SCL.

    Devuelve
    --------
    ee.Image
        La misma imagen, con máscara en los píxeles descartados.
    """
    scl = img.select("SCL")
    malos = (scl.eq(0).Or(scl.eq(1)).Or(scl.eq(3)).Or(scl.eq(6))
             .Or(scl.eq(8)).Or(scl.eq(9)).Or(scl.eq(10)).Or(scl.eq(11)))
    return img.updateMask(malos.Not())


def nbr(img):
    """Calcula el NBR = (B8A - B12) / (B8A + B12).

    La colección armonizada ya corrige el desfase de los datos posteriores a
    enero de 2022, por lo que no se aplica ningún offset adicional.

    Parámetros
    ----------
    img : ee.Image
        Imagen Sentinel-2 con las bandas B8A (infrarrojo cercano) y B12 (SWIR).

    Devuelve
    --------
    ee.Image
        Imagen de una banda llamada "NBR", con valores entre -1 y 1.
    """
    return img.normalizedDifference(["B8A", "B12"]).rename("NBR")


def compuesto_mensual(col, anio, mes):
    """Calcula la mediana del NBR de todas las imágenes de un mes.

    Parámetros
    ----------
    col : ee.ImageCollection
        Colección Sentinel-2 ya filtrada y con máscara de nubes.
    anio : int
        Año del mes a compuestar.
    mes : int
        Mes a compuestar, de 1 a 12.

    Devuelve
    --------
    ee.Image
        Imagen de una banda "NBR" de tipo flotante. Si el mes no tiene imágenes,
        devuelve una imagen totalmente enmascarada del mismo tipo, para que los
        pasos siguientes no fallen.
    """
    ini = ee.Date.fromYMD(anio, mes, 1)
    sub = col.filterDate(ini, ini.advance(1, "month")).map(nbr)
    vacio = (ee.Image.constant(0).toFloat().rename("NBR")
             .updateMask(ee.Image.constant(0)))
    # El tipo de banda debe ser idéntico en todos los meses (flotante); sin la
    # conversión explícita, un mes sin imágenes produce un tipo distinto y
    # Earth Engine rechaza la colección de compuestos.
    return ee.Image(ee.Algorithms.If(sub.size().gt(0), sub.median(), vacio)).toFloat()


def mes_relativo(anio, i):
    """Convierte un índice de mes relativo a un par (año, mes).

    El índice 0 es enero de `anio`; el 11 es diciembre; los negativos son los
    meses previos (-1 es diciembre del año anterior).

    Parámetros
    ----------
    anio : int
        Año de referencia.
    i : int
        Índice de mes relativo a enero de `anio`.

    Devuelve
    --------
    tuple (int, int)
        (año, mes) con el mes de 1 a 12.
    """
    return anio + (i // 12), (i % 12) + 1


def construir_dnbr(anio, region, meses_previos, nubes_max):
    """Construye la imagen de dNBR mensual de un año.

    El dNBR del mes m es la mediana de los NBR mensuales de los `meses_previos`
    meses anteriores menos el NBR del mes m. Los valores positivos indican caída
    del NBR respecto a la línea base.

    Parámetros
    ----------
    anio : int
        Año a procesar.
    region : ee.Geometry
        Área de interés, usada para filtrar las imágenes.
    meses_previos : int
        Número de meses anteriores que forman la línea base.
    nubes_max : float
        Porcentaje máximo de nubes por imagen (CLOUDY_PIXEL_PERCENTAGE).

    Devuelve
    --------
    ee.Image
        Imagen de 12 bandas flotantes llamadas dNBR_AAAAMM (enero a diciembre),
        con máscara donde el mes o la línea base no tienen observaciones válidas.
    """
    a_ini, m_ini = mes_relativo(anio, -meses_previos)
    inicio = ee.Date.fromYMD(a_ini, m_ini, 1)
    fin = ee.Date.fromYMD(anio, 12, 31).advance(1, "day")

    col = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
           .filterBounds(region)
           .filterDate(inicio, fin)
           .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", nubes_max))
           .map(mascara_scl))

    proj = ee.Image(col.first()).select("B8A").projection()

    comp = {}
    for i in range(-meses_previos, 12):
        a, m = mes_relativo(anio, i)
        comp[i] = compuesto_mensual(col, a, m)

    bandas = []
    for i in range(0, 12):
        base = ee.ImageCollection.fromImages(
            [comp[j] for j in range(i - meses_previos, i)]).median()
        dnbr = base.subtract(comp[i]).rename(f"dNBR_{anio}{i + 1:02d}")
        bandas.append(dnbr)

    img = ee.Image.cat(bandas).setDefaultProjection(proj)
    return img


def agregar_a_grilla(img, bbox, res):
    """Promedia la imagen a una grilla regular en grados y la convierte a entero.

    Parámetros
    ----------
    img : ee.Image
        Imagen de dNBR mensual a su resolución nativa (salida de construir_dnbr).
    bbox : tuple de 4 float
        (oeste, sur, este, norte) en grados; define la esquina superior izquierda
        de la grilla.
    res : float
        Tamaño de píxel de salida en grados. Debe dividir exactamente a 0.25.

    Devuelve
    --------
    tuple (ee.Image, list)
        La imagen agregada, como entero de 16 bits con el dNBR multiplicado por
        1000 y el valor NODATA donde no hay dato, y el crs_transform de 6
        valores [res, 0, oeste, 0, -res, norte] que define la grilla.
    """
    oeste, _, _, norte = bbox
    transform = [res, 0, oeste, 0, -res, norte]
    lado_m = res * M_POR_GRADO
    max_pix = int(math.ceil((lado_m / RES_NATIVA_M) ** 2 * 1.5))
    if max_pix > 65536:
        raise ValueError("--res demasiado grueso para reduceResolution; usar un valor menor")
    agregada = (img.reduceResolution(reducer=ee.Reducer.mean(), maxPixels=max(max_pix, 64))
                .reproject(crs="EPSG:4326", crsTransform=transform))
    return agregada.multiply(1000).round().toInt16().unmask(NODATA), transform


def descargar(img, transform, ancho, alto, ruta):
    """Descarga una imagen de Earth Engine como GeoTIFF multibanda.

    Con crs y crs_transform, Earth Engine ignora "region" y recorta con
    "dimensions". El transform va en orden de filas: [res, 0, oeste, 0, -res, norte].

    Parámetros
    ----------
    img : ee.Image
        Imagen a descargar (salida de agregar_a_grilla).
    transform : list de 6 float
        Transformación afín de la grilla en EPSG:4326.
    ancho : int
        Número de columnas de la imagen descargada.
    alto : int
        Número de filas de la imagen descargada.
    ruta : str
        Ruta del archivo .tif de salida.

    Devuelve
    --------
    float
        Tamaño del archivo descargado, en megabytes.
        Lanza RuntimeError si el servicio responde con un error.
    """
    url = img.getDownloadURL({
        "crs": "EPSG:4326",
        "crs_transform": transform,
        "dimensions": [ancho, alto],
        "format": "GEO_TIFF",
    })
    r = requests.get(url, timeout=900)
    if r.status_code != 200 or "json" in r.headers.get("Content-Type", ""):
        raise RuntimeError(f"Fallo la descarga ({r.status_code}): {r.text[:500]}")
    with open(ruta, "wb") as f:
        f.write(r.content)
    return len(r.content) / 1e6


def comprobar_extension(ruta, bbox):
    """Verifica que la extensión del GeoTIFF coincida con la caja solicitada.

    Si rasterio no está instalado, no hace nada.

    Parámetros
    ----------
    ruta : str
        Ruta del archivo .tif descargado.
    bbox : tuple de 4 float
        (oeste, sur, este, norte) esperados, en grados.

    Devuelve
    --------
    None
        Imprime una confirmación o un aviso si las extensiones difieren en más
        de 1e-6 grados.
    """
    try:
        import rasterio
    except ImportError:
        return
    with rasterio.open(ruta) as src:
        b = src.bounds
        esperado = (bbox[0], bbox[1], bbox[2], bbox[3])
        real = (b.left, b.bottom, b.right, b.top)
        if any(abs(x - y) > 1e-6 for x, y in zip(real, esperado)):
            print(f"AVISO: extensión del tif {real} distinta de la caja {esperado}")
        else:
            print(f"Extensión verificada: {real}, {src.count} bandas")


def teselas(bbox, tesela):
    """Divide una caja en teselas cuadradas, de oeste a este y de norte a sur.

    Parámetros
    ----------
    bbox : tuple de 4 float
        (oeste, sur, este, norte) en grados.
    tesela : float
        Lado de cada tesela en grados. Las teselas del borde este y sur pueden
        ser más pequeñas si la caja no es múltiplo de este valor.

    Devuelve
    --------
    list de tuple (int, int, tuple)
        Una entrada por tesela: (fila, columna, bbox_tesela), donde fila y
        columna cuentan desde la esquina noroeste (empiezan en 0) y bbox_tesela
        es (oeste, sur, este, norte) en grados.
    """
    oeste, sur, este, norte = bbox
    n_col = math.ceil((este - oeste) / tesela - 1e-9)
    n_fil = math.ceil((norte - sur) / tesela - 1e-9)
    lista = []
    for f in range(n_fil):
        for c in range(n_col):
            o = oeste + c * tesela
            n = norte - f * tesela
            caja = (round(o, 9), round(max(n - tesela, sur), 9),
                    round(min(o + tesela, este), 9), round(n, 9))
            lista.append((f, c, caja))
    return lista


def ensamblar(rutas, bbox, res, year, ruta_salida):
    """Une las teselas descargadas en un único GeoTIFF multibanda.

    Requiere rasterio. Cada tesela se coloca según su esquina noroeste dentro
    de la grilla de la caja completa.

    Parámetros
    ----------
    rutas : list de str
        Rutas de los GeoTIFF de las teselas (12 bandas, enteros de 16 bits).
    bbox : tuple de 4 float
        (oeste, sur, este, norte) de la caja completa, en grados.
    res : float
        Tamaño de píxel en grados, igual al usado para descargar las teselas.
    year : int
        Año procesado; se usa para nombrar las bandas (dNBR_AAAAMM).
    ruta_salida : str
        Ruta del GeoTIFF final.

    Devuelve
    --------
    None
        Escribe `ruta_salida` con 12 bandas, CRS EPSG:4326, nodata declarado
        (-32768) y compresión deflate.
    """
    import numpy as np
    import rasterio
    from rasterio.transform import from_origin

    oeste, sur, este, norte = bbox
    ancho = round((este - oeste) / res)
    alto = round((norte - sur) / res)
    datos = np.full((12, alto, ancho), NODATA, dtype="int16")
    for ruta in rutas:
        with rasterio.open(ruta) as src:
            col = round((src.bounds.left - oeste) / res)
            fila = round((norte - src.bounds.top) / res)
            datos[:, fila:fila + src.height, col:col + src.width] = src.read()
    with rasterio.open(
            ruta_salida, "w", driver="GTiff", height=alto, width=ancho, count=12,
            dtype="int16", crs="EPSG:4326", transform=from_origin(oeste, norte, res, res),
            nodata=NODATA, compress="deflate") as dst:
        dst.write(datos)
        for i in range(12):
            dst.set_band_description(i + 1, f"dNBR_{year}{i + 1:02d}")


def procesar(year, nombre, bbox, res=0.001, meses_previos=3, nubes_max=70.0,
             salida="salidas_dnbr", project=None, tesela=0.25):
    """Calcula el dNBR mensual de un año y lo descarga como GeoTIFF.

    La caja se divide en teselas que se descargan una por una y se ensamblan al
    final. Si el proceso se interrumpe, las teselas ya descargadas se reutilizan
    al volver a llamar a la función con los mismos argumentos.

    Parámetros
    ----------
    year : int
        Año a procesar (por ejemplo, 2020).
    nombre : str
        Nombre corto de la región; se usa en el nombre del archivo.
    bbox : tuple de 4 float
        (oeste, sur, este, norte) en grados, con bordes múltiplos de 0.25.
    res : float, opcional
        Tamaño de píxel de salida en grados. 0.001 equivale a unos 110 m.
    meses_previos : int, opcional
        Número de meses anteriores que forman la línea base (3 por defecto).
    nubes_max : float, opcional
        Porcentaje máximo de nubes por imagen en el prefiltro (70 por defecto).
    salida : str, opcional
        Carpeta donde se guarda el archivo; se crea si no existe.
    project : str, opcional
        Proyecto de Google Cloud habilitado para Earth Engine.
    tesela : float, opcional
        Lado de cada tesela de descarga en grados (0.25 por defecto). Debe ser
        múltiplo de `res`. Reducirlo si Earth Engine responde "User memory
        limit exceeded".

    Devuelve
    --------
    str
        Ruta del GeoTIFF generado: <salida>/dnbr_<nombre>_<year>.tif, con 12
        bandas (enero a diciembre), enteros de 16 bits con el dNBR multiplicado
        por 1000 y -32768 donde no hay dato.
        Lanza ValueError si la caja, la resolución o la tesela no son válidas, o
        si falta rasterio y hay más de una tesela. Lanza RuntimeError si una
        descarga falla.
    """
    validar_caja(bbox, res)
    pix_por_tesela = tesela / res
    if abs(pix_por_tesela - round(pix_por_tesela)) > 1e-9:
        raise ValueError(f"tesela={tesela} debe ser múltiplo de res={res}")

    lista = teselas(bbox, tesela)
    if len(lista) > 1:
        try:
            import rasterio  # noqa: F401
        except ImportError:
            raise ValueError("Para ensamblar teselas se requiere rasterio "
                             "(pip install rasterio) o usar tesela >= la caja.")
    for _, _, caja in lista:
        n_bytes, _, _ = estimar_bytes(caja, res)
        if n_bytes > LIMITE_BYTES:
            raise ValueError("Una tesela supera el límite de descarga: "
                             "aumentar res o reducir tesela.")

    ee.Initialize(project=project)
    os.makedirs(salida, exist_ok=True)
    carpeta = os.path.join(salida, f"teselas_{nombre}_{year}")
    os.makedirs(carpeta, exist_ok=True)

    print(f"{nombre} {year}: {len(lista)} teselas de {tesela} grados, "
          f"~{res * M_POR_GRADO:.0f} m por píxel")
    rutas = []
    for k, (f, c, caja) in enumerate(lista, start=1):
        ruta_t = os.path.join(carpeta, f"f{f}_c{c}.tif")
        rutas.append(ruta_t)
        if os.path.exists(ruta_t) and os.path.getsize(ruta_t) > 0:
            print(f"  [{k}/{len(lista)}] f{f} c{c}: ya descargada")
            continue
        o, s, e, n = caja
        region = ee.Geometry.Rectangle([o, s, e, n], "EPSG:4326", False)
        img = construir_dnbr(year, region, meses_previos, nubes_max)
        img, transform = agregar_a_grilla(img, caja, res)
        _, ancho_t, alto_t = estimar_bytes(caja, res)
        mb = descargar(img, transform, ancho_t, alto_t, ruta_t)
        print(f"  [{k}/{len(lista)}] f{f} c{c}: {mb:.1f} MB")

    ruta = os.path.join(salida, f"dnbr_{nombre}_{year}.tif")
    ensamblar(rutas, bbox, res, year, ruta)
    shutil.rmtree(carpeta)
    print(f"Guardado {ruta}")
    comprobar_extension(ruta, bbox)
    return ruta


def main():
    """Punto de entrada de la línea de comandos.

    Parámetros
    ----------
    Ninguno. Lee los argumentos con parsear_argumentos().

    Devuelve
    --------
    None
        Ejecuta procesar() y termina el programa con un mensaje si hay un ValueError.
    """
    args = parsear_argumentos()
    try:
        procesar(args.year, args.nombre, tuple(args.bbox), res=args.res,
                 meses_previos=args.meses_previos, nubes_max=args.nubes_max,
                 salida=args.salida, project=args.project, tesela=args.tesela)
    except ValueError as error:
        sys.exit(str(error))


if __name__ == "__main__":
    main()
