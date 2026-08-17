# Guía de instalación — Windows


## 1. Miniforge

Descarga `Miniforge3-Windows-x86_64.exe` desde:
`https://github.com/conda-forge/miniforge/releases/latest`

Durante la instalación:
- **Install for:** "Just Me" (no "All Users").
- Ruta de instalación por defecto: aceptar.
- **No marcar** "Add Miniforge3 to my PATH environment variable" (el propio instalador lo desaconseja).

## 2. Git for Windows

Descarga desde `https://git-scm.com/download/win`.

Durante la instalación:
- En "Adjusting your PATH environment", deja la opción recomendada por defecto: **"Git from the command line and also from 3rd-party software"**.
- El resto de opciones: valores por defecto.

## 3. Configurar Git Bash como terminal por defecto en VSCode

1. Paleta de comandos (`Ctrl + Shift + P`).
2. Busca y selecciona `Terminal: Select Default Profile` → elige **Git Bash**.
3. Cierra cualquier terminal abierta y abre una nueva (`` Ctrl + ` ``): debe abrir Git Bash, no PowerShell.

Verifica:
```bash
git --version
```

## 4. Conectar conda con Git Bash

Desde la terminal de Git Bash en VSCode:

```bash
"$USERPROFILE/miniforge3/Scripts/conda.exe" init bash
```

Cierra VSCode por completo y ábrelo de nuevo. Verifica:

```bash
conda --version
```

## 5. Configurar identidad de Git (una sola vez)

```bash
git config --global user.name "Tu Nombre"
git config --global user.email "tu-correo@ejemplo.com"
```

Usa el mismo correo asociado a tu cuenta de GitHub.

## 6. VSCode: extensiones

Desde el panel de Extensiones (`Ctrl + Shift + X`), instala:
- **Python** (`ms-python.python`)
- **Jupyter** (`ms-toolsai.jupyter`)

## 7. Clonar el repositorio del curso

1. En GitHub.com, entra al repositorio `1GEO20`, clic en **Code**, pestaña **HTTPS** (no SSH), copia la URL.
2. En VSCode: paleta de comandos (`Ctrl + Shift + P`) → `Git: Clone` → pega la URL HTTPS.
3. Elige dónde guardar la carpeta localmente.

## 8. Crear el entorno del curso

Desde Git Bash, dentro de la carpeta del repositorio ya clonado:

```bash
conda env create -f environment.yml
conda activate 1geo20
```

Para actualizarlo más adelante (por ejemplo, si se agrega una librería nueva), sin rehacer todo el entorno:

```bash
conda env update -f environment.yml --prune
```

<!-- ## 9. Primer commit y push

Desde el panel de Source Control (`Ctrl + Shift + G`) en VSCode: marca el archivo con el ícono "+" (stage), escribe un mensaje de commit, clic en el ✓ (Commit), luego en "Push" o "Sync Changes".

En este primer push se abre el navegador pidiendo iniciar sesión en GitHub. Autoriza el acceso.-->)

## Verificación final

```bash
conda --version
git --version
python -c "import numpy, pandas, xarray; print('entorno OK')"
```
