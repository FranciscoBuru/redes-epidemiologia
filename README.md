# Instalación de Julia

`https://docs.junolab.org/stable/man/installation/`

An la liga de arriba está todo para instalar y usar con Atom.

# Introducción

Modelos SIS y SIR para metapoblaciones representadas por redes.

# Generado de grafos
Usamos la librería `Graphs` para generar las topologías a evaluar. Creamos
3 métodos que facilitan la obtención de los grafos generados. Cada
uno de los métodos regresa el grafo que indica su nombre.

Se pueden generar muchos más grafos, toda la documentación está disponible en
el siguiente [Link](https://juliagraphs.org/Graphs.jl/stable/generators/#All-Generators)

# Datos

Hay grafos de muchas cosas disponibles en el siguiente [Link](https://snap.stanford.edu/data/index.html)


# Importante:

El método para generar grafos small world, de la libreria (watts_strogatz) no es
determinístico con la seed puesta. Ejemplo:
 ```julia
 grafo1 = watts_strogatz(5, 4, 0.3, seed=2)
 grafo2 = watts_strogatz(5, 4, 0.3, seed=2)
 print(grafo1==grafo2)
 ```
 ```
 Output = true
 ```
 y
 ```julia
 grafo1 = watts_strogatz(10, 4, 0.3, seed=2)
 grafo2 = watts_strogatz(10, 4, 0.3, seed=2)
 print(grafo1==grafo2)
 ```
 ```
 Output = false
 ```
