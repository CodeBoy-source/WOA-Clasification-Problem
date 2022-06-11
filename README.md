# ** _Brian Sena Simons 3ºA - A2_ **
# Práctica 3 - MetaHeurística

### - Utilize el "cmake CMakeLists.txt && make" para compilar;
### - Los algoritmos pueden ser ejecutados utilizando los scripts en ./scripts {Ej. runAll.sh}
        {EMPTY}
### - Los resultados se guardan en ./results;
## Ejecución individual:
### Algoritmo WOA:
        {EMPTY}

## Descripción breve del Problema
La idea es comparar distintos tipos de algoritmos para clasificar datos pertenecientes
a una base de datos públicas que nos provee el profesorado. Partiremos primero
de la implementación del típico algoritmo de clasficación K-NN dónde K representa
el número de vecinos a mirar y la idea conceptual es buscar los K vecinos más
cercanos para realizar una predicción sobre que clase pertenece el objeto a predecir.

Una vez implementado el algoritmo 1NN intentaremos mejorar el porcentaje de aciertos
utilizando técnicas de ponderación de características mediante un vector de pesos.
El grueso de la práctica está en el cálculo de esos pesos. En esta parte de la
práctica compararemos los algoritmos empleados anteriormente (1NN, Greedy,
búsqueda local) con un algoritmo bio-inspirado basado en las ballenas jorobadas.

El comportamiento del algoritmo es como una variación del algoritmo basado en
partículas, dónde el movimiento de las partículas en este caso es una descripción
de la caza que realizan las ballenas jorobadas. Simulan un tipo de espiral o círculo
alrededor de su preza, en este caso el óptimo de la función, en el cual se nota
un factor explorativo explotativo que demarca la gran eficiencia de este algoritmo
respecto a la optimización de función por su compromiso equilibrado entre ambas
características.

La idea es implementar y describir el algoritmo, analizar los resultados obtenidos
y planificar/implementar una mejora del algoritmo basandose en lo estudiado
en las clases de teoría sobre heurísticas y sus caracterísicas.
