# Instalación de Julia

`https://docs.junolab.org/stable/man/installation/`

An la liga de arriba está todo para instalar y usar con Atom.

# Introducción

Modelos SIS y SIR para metapoblaciones representadas por redes.

# Correr modelos.

Se necesitan una serie de paquetes de julia para poder ejecutar el sistema de modelado. A continuación están los comandos para añadirlos a julia.

`Pkg.add("Random")`
`Pkg.add("StatsBase")`
`Pkg.add("GraphPlot")`
`Pkg.add("Plots")`
`Pkg.add("DelimitedFiles")`
`Pkg.add("Compose")`
`Pkg.add("Cairo")`
`Pkg.add("Graphs")`

## Modelos sintéticos

Para solamente correr los modelos ejecutar los archivos `run-[topologia].jl`. Si lo que se quiere son las gráficas de promedio de infectados alternando los parametros `r` y `beta` correro `rung-[topologia].jl` CUIDADO: Las gráficas nuevas sobreescriben a las viejas en la carpeta `plots`.


En el modelo SIR todos los nodos se inicializan con 0.0 en la probabilidad de
estár recuperado.

# Generado de grafos
Usamos la librería `Graphs` para generar las topologías a evaluar. Creamos
3 métodos que facilitan la obtención de los grafos generados. Cada
uno de los métodos regresa el grafo que indica su nombre.

Se pueden generar muchos más grafos, toda la documentación está disponible en
el siguiente [Link](https://juliagraphs.org/Graphs.jl/stable/generators/#All-Generators)

# Datos

Hay grafos de muchas cosas disponibles en el siguiente [Link](https://snap.stanford.edu/data/index.html)


# Importante:

## Esto ya está corregido.

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

 # Compilar e instalar módulo de python local

 Para general los puntos de HitandRun usamos el paquete de python mhar. Pero estamos en julia entonces debemos correr python y meter el output a julia. Para correr python en julia debemos instalar el paquete PyCall. Ya con Pycall podemos llamar a nuestro módulo.

 Para instalar nuestro módulo globalmente primero debemos instalar las dos dependencias que usa:
```
pip install mhar
```
```
pip install numpy
```

 Y luego hay que ir al directorio general del código de python y correr:

 ```
python hitAndRun.py install
 ```
 o
 ```
pip3 install -e hitAndRun
 ```

 con eso podremos importar nuestro nuevo paquete. Con lo anterior podemos traer datos desde python pero cuidado, tarda un rato.

 Nota: El paquete está hardoceado para generar los puntos en el hipercubo de 12 dimensiones con restricciones >0 y <1 en todas las dimensiones.
