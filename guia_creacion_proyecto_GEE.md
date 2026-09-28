# Guía paso a paso: Creación y configuración de un proyecto gratuito en Google Earth Engine (Uso educativo / no comercial)

Esta guía detalla el proceso para habilitar el acceso a **Google Earth Engine (GEE)** mediante un proyecto gratuito de Google Cloud.

Se requiere:
- Una cuenta personal o institucional de Google (`@gmail.com` o correo universitario).
- Navegador web (se recomienda desactivar temporalmente bloqueadores de ventanas emergentes o de anuncios).



##  Paso 1: Crear el proyecto en Google Cloud Console

1. Ingresa a la consola de Google Cloud: **[console.cloud.google.com](https://console.cloud.google.com/)** e inicia sesión con tu cuenta de Google.
2. En la barra superior, haz clic en el selector desplegable de proyectos (junto al logo de *Google Cloud*).
3. En la ventana emergente, haz clic en el botón **Proyecto nuevo** (*New Project*) situado en la esquina superior derecha.
4. Completa los campos:
   - **Nombre del proyecto:** `nombre-geocomputacion` *(los nombres solo aceptan letras, números, espacios y guiones)*.
   - **ID del proyecto:** Haz clic en **Editar** si deseas personalizarlo (por ejemplo, `nombre-geocomputacion-edu`), o anota el ID generado automáticamente (ej. `nombre-geocomputacion-435102`).  
     > **Importante:** Este **Project ID** es el identificador único que usarás en tus scripts de Python, no el nombre general.
   - **Organización / Ubicación:** Dejarlo en `Sin organización` o indicar la cuenta institucional administrada.
5. Haz clic en **Crear** y espera unos segundos hasta que la notificación confirme que el proyecto fue creado.



## 3. Paso 2: Registrar y vincular el proyecto a Earth Engine

Para usar el proyecto con la API de Earth Engine sin costo:

1. Ve a la página de registro de Earth Engine: **[code.earthengine.google.com/register](https://code.earthengine.google.com/register)**.
2. Selecciona la opción de uso:
   - Marca **Noncommercial / Academic / Research** (Uso no comercial / Académico).
3. Selecciona el tipo de perfil:
   - Elige **Student** (Estudiante) o **Education** (Educación).
4. Elige cómo vincular el proyecto de Cloud:
   - Marca **Use an existing Google Cloud project** (*Usar un proyecto de Google Cloud existente*).
   - En la lista desplegable, selecciona el proyecto creado en el paso anterior (`nombre-geocomputacion`) o introduce su **Project ID**.
5. Lee y acepta los términos del servicio de Earth Engine.
6. Haz clic en **Confirm and continue** / **Submit**.

> **Nota:** La habilitación suele ser inmediata. Una vez completado, serás redirigido al Code Editor web de Earth Engine (`code.earthengine.google.com`).



## 4. Paso 3: Inicialización en Python (Jupyter / Google Colab)

Una vez registrado, ya puedes autenticarte y conectar tu entorno de trabajo en Python.

### En Google Colab

```python
# 1. Instalar la librería cliente (si no está instalada)
!pip install -q earthengine-api

import ee

# 2. Autenticación (se recomienda 'notebook')
ee.Authenticate(auth_mode='notebook')

# 3. Inicialización usando tu Project ID
PROJECT_ID = 'nombre-geocomputacion'  # Reemplaza con tu Project ID exacto
ee.Initialize(project=PROJECT_ID)

# 4. Prueba rápida de conexión
print("Conexión exitosa:", ee.String("Google Earth Engine activo").getInfo())