"""
Funciones de apoyo para la clasificación de cobertura de suelo (semana 11).

Se usan desde el notebook clasificacion_cobertura_demo.ipynb y se reutilizan en la práctica.
Cada función tiene un docstring: se puede consultar con help(uc.nombre_de_la_funcion).
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import pandas as pd
from sklearn.metrics import ConfusionMatrixDisplay, accuracy_score, precision_recall_fscore_support

# Colores por defecto de las clases de la demo (las clases nuevas usan una paleta de respaldo)
COLORES_CLASES = {
    "Bosque": "#1f8d49",
    "Bosque inundable": "#026975",
    "Pasto, agropecuario y herbazal": "#e6d27a",
    "Suelo sin vegetación": "#e08e5a",
    "Urbano": "#d4271e",
    "Minería": "#7a1b3a",
    "Agua": "#2a52e0",
}


def agrupar_clases(mapa, agrupacion):
    """
    Agrupa los códigos de un mapa de cobertura en clases nuevas.

    Parámetros
    ----------
    mapa : array 2D con los códigos originales (por ejemplo, de MapBiomas).
    agrupacion : dict {nombre de la clase: lista de códigos originales}.

    Devuelve
    --------
    grupos : array 2D de enteros. 1, 2, ... según el orden del diccionario, y 0 donde el
        código no pertenece a ninguna clase.
    nombres : lista con los nombres de las clases (la clase k tiene el nombre nombres[k - 1]).
    """
    grupos = np.zeros(mapa.shape, dtype="int16")
    for k, codigos in enumerate(agrupacion.values(), start=1):
        grupos[np.isin(mapa, codigos)] = k
    return grupos, list(agrupacion)


def muestrear_por_clase(y, n_por_clase, min_pixeles, semilla=0, nombres=None):
    """
    Elige hasta n_por_clase píxeles de cada clase. Las clases con menos de min_pixeles
    píxeles en total se descartan.

    Parámetros
    ----------
    y : array 1D con la clase de cada píxel.
    n_por_clase : máximo de píxeles por clase (si hay menos, se toman todos).
    min_pixeles : mínimo de píxeles que debe tener una clase en la imagen para conservarla.
    semilla : semilla del generador aleatorio.
    nombres : lista de nombres (opcional), solo para informar qué clases se descartan.

    Devuelve
    --------
    idx : índices (ordenados) de los píxeles elegidos.
    clases : lista de las clases conservadas.
    """
    rng = np.random.default_rng(semilla)
    valores, cuentas = np.unique(y, return_counts=True)
    clases = [int(v) for v, n in zip(valores, cuentas) if n >= min_pixeles]
    for v, n in zip(valores, cuentas):
        if n < min_pixeles:
            etiqueta = nombres[int(v) - 1] if nombres else int(v)
            print(f"Clase descartada por tener pocos píxeles ({n}): {etiqueta}")
    elegidos = []
    for clase in clases:
        posiciones = np.flatnonzero(y == clase)
        n = min(n_por_clase, len(posiciones))
        elegidos.append(rng.choice(posiciones, size=n, replace=False))
    return np.sort(np.concatenate(elegidos)), clases


def dividir_por_bloques(filas, columnas, tam_bloque, prop_prueba, semilla=0, y=None, min_por_clase=30):
    """
    Divide los píxeles en entrenamiento y prueba por zonas. La imagen se parte en bloques
    cuadrados de tam_bloque píxeles y cada bloque va completo a entrenamiento o a prueba.

    Parámetros
    ----------
    filas, columnas : posición (fila y columna) de cada píxel de la muestra.
    tam_bloque : lado del bloque, en píxeles.
    prop_prueba : proporción de bloques que se destinan a prueba.
    semilla : semilla del generador aleatorio.
    y : clase de cada píxel (opcional). Si se entrega, se vuelve a sortear (hasta 100 veces)
        hasta que cada clase tenga al menos min_por_clase píxeles en entrenamiento y en prueba.

    Devuelve
    --------
    es_prueba : array de booleanos, True para los píxeles de prueba.
    """
    rng = np.random.default_rng(semilla)
    id_bloque = (filas // tam_bloque) * 100000 + (columnas // tam_bloque)
    bloques = np.unique(id_bloque)
    n_prueba = max(1, int(round(prop_prueba * len(bloques))))
    for _ in range(100):
        es_prueba = np.isin(id_bloque, rng.choice(bloques, size=n_prueba, replace=False))
        if y is None:
            return es_prueba
        if all((y[es_prueba] == c).sum() >= min_por_clase and (y[~es_prueba] == c).sum() >= min_por_clase
               for c in np.unique(y)):
            return es_prueba
    print("Aviso: alguna clase quedó con muy pocos píxeles en entrenamiento o en prueba; "
          "prueba con otro tamaño de bloque o con otra semilla")
    return es_prueba


def evaluar(modelo, X_prueba, y_prueba, clases, nombres, ax=None, titulo=""):
    """
    Evalúa un modelo con los datos de prueba: imprime la exactitud global y el reporte de
    precisión y exhaustividad (las exactitudes del usuario y del productor) por clase, y dibuja la matriz de confusión si se entrega un eje.
    Las filas de la matriz son la clase de referencia y las columnas la clase predicha.

    Parámetros
    ----------
    modelo : modelo de scikit-learn ya entrenado.
    X_prueba, y_prueba : datos de prueba.
    clases : lista de clases (enteros) en el orden en que se muestran.
    nombres : nombres de esas clases, en el mismo orden.
    ax : eje de matplotlib para la matriz de confusión (opcional).
    titulo : título de la matriz.

    Devuelve
    --------
    La exactitud global en prueba.
    """
    prediccion = modelo.predict(X_prueba)
    exactitud = accuracy_score(y_prueba, prediccion)
    if ax is not None:
        ConfusionMatrixDisplay.from_predictions(
            y_prueba, prediccion, labels=clases, display_labels=nombres,
            xticks_rotation=45, ax=ax, colorbar=False)
        ax.set_xlabel("Clase predicha")
        ax.set_ylabel("Clase de referencia")
        ax.set_title(titulo)
    usuario, productor, _, n_prueba = precision_recall_fscore_support(
        y_prueba, prediccion, labels=clases, zero_division=0)
    tabla = pd.DataFrame({"Exactitud del usuario": usuario,
                          "Exactitud del productor": productor,
                          "Píxeles de prueba": n_prueba}, index=nombres).round(3)
    print(f"Exactitud global en prueba: {exactitud:.3f}\n")
    print(tabla)
    return exactitud


def predecir_imagen(modelo, pila):
    """
    Clasifica todos los píxeles de una imagen con un modelo ya entrenado.

    Parámetros
    ----------
    modelo : modelo de scikit-learn ya entrenado.
    pila : array (bandas, alto, ancho). Los píxeles con NaN en alguna banda no se clasifican.

    Devuelve
    --------
    Array 2D (alto, ancho) con la clase predicha, y NaN donde no hay datos.
    """
    n_bandas, alto, ancho = pila.shape
    X = pila.reshape(n_bandas, -1).T
    valido = ~np.isnan(X).any(axis=1)
    salida = np.full(X.shape[0], np.nan, dtype="float32")
    salida[valido] = modelo.predict(X[valido])
    return salida.reshape(alto, ancho)


# ----- Apoyo para las figuras -----
def colores_para(nombres):
    """Lista de colores, uno por clase, según COLORES_CLASES (o una paleta de respaldo)."""
    respaldo = iter(plt.get_cmap("tab10").colors)
    return [COLORES_CLASES.get(nombre, next(respaldo)) for nombre in nombres]


def leyenda_clases(ax, nombres, colores, **kwargs):
    """Agrega una leyenda con un cuadro de color por clase."""
    ax.legend(handles=[Patch(color=c, label=n) for n, c in zip(nombres, colores)], **kwargs)


def mapa_particion(filas, columnas, es_prueba, forma, ax, titulo, tam_bloque=None):
    """
    Dibuja dónde caen los píxeles de entrenamiento y de prueba dentro de la imagen.
    forma = (alto, ancho) de la imagen. Si se entrega tam_bloque, se marca la cuadrícula de bloques.
    """
    ax.scatter(columnas[~es_prueba], filas[~es_prueba], s=1, c="tab:blue", label="Entrenamiento")
    ax.scatter(columnas[es_prueba], filas[es_prueba], s=1, c="tab:orange", label="Prueba")
    ax.set_xlim(0, forma[1])
    ax.set_ylim(forma[0], 0)
    ax.set_aspect("equal")
    if tam_bloque:
        ax.set_xticks(np.arange(0, forma[1], tam_bloque))
        ax.set_yticks(np.arange(0, forma[0], tam_bloque))
        ax.grid(True, linewidth=0.5, color="gray")
        ax.tick_params(labelbottom=False, labelleft=False)
    ax.set_title(titulo)
    ax.legend(markerscale=8, loc="lower right")


def guardar_figura(fig, nombre, activo=True, carpeta="figuras"):
    """Guarda la figura como PNG en la carpeta indicada (solo si activo es True)."""
    if activo:
        Path(carpeta).mkdir(exist_ok=True)
        fig.savefig(Path(carpeta) / f"{nombre}.png", dpi=200, bbox_inches="tight")
