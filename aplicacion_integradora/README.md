# dNBR mensual de Sentinel-2 y fracción quemada en la grilla de 0.25°

Flujo de trabajo reutilizable para estimar área quemada a partir de imágenes Sentinel-2 y llevarla a la misma grilla que el modelo de área quemada Sparky v1.0 (ECMWF), con lectura opcional de MapBiomas Fuego como referencia. Material del curso Geocomputación (1GEO20), PUCP, 2026-II.

## Qué hace

1. **Descarga y procesa** con Google Earth Engine el dNBR mensual de un año sobre una caja geográfica, y lo guarda como GeoTIFF de 12 bandas.
2. **Lee y visualiza** los GeoTIFF: paneles mensuales y animaciones (GIF).
3. **Estima la fracción quemada** observada por celda de 0.25° y mes, y la lleva a la misma grilla que el modelo y que MapBiomas Fuego.

## Contenido

| Archivo | Descripción |
|---|---|
| `dnbr_mensual.py` | Calcula y descarga el dNBR mensual de un año con Earth Engine. Se usa como módulo o desde la línea de comandos. |
| `dnbr_vis.py` | Lectura de los GeoTIFF (`leer_dnbr`), panel de 12 mapas (`panel_mensual`) y animación GIF (`animar`). |
| `dnbr_analisis.py` | Anomalía espacial, clasificación quemado / no quemado, primer mes de detección, agregación a 0.25°, lectura del modelo y de MapBiomas, y alineación de grillas. |
| `integrador_dnbr_modelo.ipynb` | Cuaderno que recorre el flujo completo con una región de ejemplo. |
| `requirements.txt` | Librerías necesarias. |

Todas las funciones tienen docstring con la descripción de cada parámetro y de lo que devuelven (`help(nombre_de_funcion)`).

## Estructura de carpetas esperada

```
.
├── dnbr_mensual.py
├── dnbr_vis.py
├── dnbr_analisis.py
├── integrador_dnbr_modelo.ipynb
├── salidas_dnbr/        GeoTIFF generados: dnbr_<nombre>_<año>.tif
├── datos_modelo/        NetCDF del modelo: BAF_pred_<año>*.nc
└── mapbiomas/           (opcional) GeoTIFF anuales de MapBiomas Fuego
```

Las carpetas de datos se crean o se llenan al ejecutar el flujo, y las rutas se cambian en la primera celda del cuaderno.

## Requisitos

- Python 3 reciente (desarrollado y probado con 3.13).
- Librerías: `numpy`, `xarray`, `matplotlib`, `rasterio`, `netCDF4`, `requests` y `earthengine-api`.
- Una cuenta de Google con un proyecto de Google Cloud habilitado para Earth Engine.

```bash
pip install -r requirements.txt
```

Para ejecutar solo el cuaderno de análisis (sin descargar datos nuevos) no hace falta `earthengine-api`.

## Uso

### 1. Procesar el dNBR de un año (Earth Engine)

Primero se autentica Earth Engine, una vez por sesión:

```python
import ee
ee.Authenticate(auth_mode='notebook')
```

Luego se llama a la función desde Python:

```python
import dnbr_mensual as dn

ruta = dn.procesar(2020, 'ucayali', (-74.5, -9.5, -73.5, -8.5),
                   project='MI-PROYECTO-EE')
```

Los argumentos principales son el año, un nombre corto de la región, la caja `(oeste, sur, este, norte)` en grados y el proyecto de Earth Engine. Los bordes de la caja deben ser múltiplos de 0.25 para que coincidan con la grilla del modelo. El resto de parámetros (resolución, meses de línea base, nubosidad máxima, tamaño de tesela, carpeta de salida) se describen en la docstring de `procesar`.

El mismo proceso se puede ejecutar desde la terminal:

```bash
python dnbr_mensual.py --year 2020 --nombre ucayali \
    --bbox -74.5 -9.5 -73.5 -8.5 --project MI-PROYECTO-EE
```

Para repetirlo en varios años o regiones basta un bucle:

```python
import os

cajas = {"ucayali": (-74.5, -9.5, -73.5, -8.5)}
for nombre, bbox in cajas.items():
    for year in (2019, 2020, 2021):
        if os.path.exists(f"salidas_dnbr/dnbr_{nombre}_{year}.tif"):
            continue
        dn.procesar(year, nombre, bbox, project="MI-PROYECTO-EE")
```

### 2. Visualizar

```python
import dnbr_vis as dv

da = dv.leer_dnbr("salidas_dnbr/dnbr_ucayali_2020.tif")
dv.panel_mensual(da)
dv.animar(da, "ucayali_2020.gif")
```

### 3. Analizar con el cuaderno

Se abre `integrador_dnbr_modelo.ipynb` y se ajusta la primera celda (región, años, rutas y parámetros). El cuaderno sigue estos pasos.

1. Lectura y exploración del dNBR mensual.
2. Separación de la señal estacional común y elección del umbral.
3. Clasificación quemado / no quemado y primer mes de detección.
4. Fracción quemada observada en la grilla de 0.25°.
5. Lectura y exploración del modelo.
6. Alineación de la observación y el modelo, y guardado de los datos en un NetCDF.

### En Google Colab

1. Subir `dnbr_vis.py`, `dnbr_analisis.py` y, si se va a procesar con Earth Engine, `dnbr_mensual.py` a la carpeta de trabajo (`/content`), junto con los datos.
2. Comprobar las librerías al inicio:

   ```python
   try:
       import rasterio, netCDF4
   except ImportError:
       %pip install -q rasterio netCDF4
   ```

3. Los archivos de Colab desaparecen al terminar la sesión, por lo que conviene copiarlos a Drive.
4. Si se modifica un módulo después de importarlo, hay que reiniciar el entorno.

## Cómo se calcula

**dNBR mensual** (criterio heurístico, no corresponde a un método publicado):
- Sentinel-2 L2A armonizado (`COPERNICUS/S2_SR_HARMONIZED`) con máscara de nubes, sombras, cirros, nieve y agua a partir de la banda SCL.
- NBR = (B8A − B12) / (B8A + B12), con mediana mensual.
- dNBR del mes = mediana del NBR mensual de los tres meses anteriores menos el NBR del mes.

**Salida de `dnbr_mensual.py`:**
- Un GeoTIFF por región y año, de 12 bandas (enero a diciembre) en EPSG:4326, nombradas `dNBR_AAAAMM`.
- Enteros de 16 bits con el valor multiplicado por 1000 (el dNBR se recupera dividiendo por 1000).
- Sin dato igual a −32768, declarado como nodata.
- La caja se descarga en teselas, porque Earth Engine limita la memoria de las descargas directas. Las teselas ya descargadas se conservan y se reutilizan si el proceso se interrumpe.

**Fracción quemada observada** (`dnbr_analisis.py`):
- La anomalía resta a cada mes la mediana espacial de la caja, para quitar la señal estacional común.
- Un umbral sobre la anomalía clasifica cada píxel como quemado o no. El umbral es un parámetro y debe calibrarse para cada región.
- Cada píxel se cuenta solo en el primer mes en que aparece marcado.
- La fracción quemada de cada celda de 0.25° es la proporción de píxeles marcados entre los píxeles con dato.

## Datos de entrada

**Modelo Sparky v1.0** (variable `BAF_pred`, fracción quemada mensual a 0.25°). Se descarga del Climate Data Store de ECMWF: <https://xds.ecmwf.int/datasets/projections-fire-fuel-burned-area>. Los archivos se esperan con nombres de la forma `BAF_pred_<año>*.nc`. Es un producto modelado, no una observación. El atributo `product_type` de cada archivo indica si corresponde al periodo histórico o a proyecciones.

**MapBiomas Fuego, Monthly Burned Area** (opcional). Mapa anual de una sola banda, con el mes de la quema (1 a 12) en cada píxel y 0 donde no se detectó quema. El 0 no distingue entre ausencia de quema y ausencia de dato.

## Limitaciones

- El dNBR usa una línea base de tres meses que no es estacional, por lo que el paisaje completo sube y baja con las estaciones. La anomalía espacial corrige esa señal bajo el supuesto de que menos de la mitad de la caja está quemada en un mes.
- El resultado depende del umbral y de la ventana de detección, que son parámetros provisionales.
- Cada píxel procesado representa unos 110 m, por lo que las cicatrices pequeñas quedan diluidas.
- Las nubes y sombras que escapan a la máscara pueden producir falsos positivos, y los meses con pocas imágenes pueden dejar quemas sin detectar.

## Referencias

- McNorton, J. R., Comyn-Platt, E., Di Giuseppe, F. et al. *Modelled Global Burned Area from 1950 to 2100 under Shared Socioeconomic Pathways*, 2026. Preprint EGUsphere-2026-5374 (en revisión).
- Gorelick, N., Hancher, M., Dixon, M., Ilyushchenko, S., Thau, D. y Moore, R. (2017). Google Earth Engine: Planetary-scale geospatial analysis for everyone. *Remote Sensing of Environment*, 202, 18-27.
- MapBiomas Fuego: <https://mapbiomas.org>.
- Imágenes Sentinel-2: programa Copernicus, ESA.

## Autoría

Laura Paccini, Sección de Geografía y Medio Ambiente, PUCP.
