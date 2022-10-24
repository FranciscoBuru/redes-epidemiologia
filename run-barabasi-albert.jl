include("modelado.jl")
using Plots

NODOS = 500          #Cantidad de nodos
mu = 0.08           #Proba de recuperacion
beta = 0.23         #Proba de contacto con vecino
r = 0.7             #Proba de infección por concatco con vecino
eta = 0.0           #Proba de no contagio por interaccion con vecinos
w=0.02              #Proba de pasar de recuperado a suceptible (solo sir)
epocas = 20000      #Iteraciones
n0 = 4              #Nodos del grafo incial.
k_bar = 4           #Vecinos de cada nodo.
proba = 0.1         #Proba de empezar infectado.
seed = 6

x,y = SIS_Barabasi(NODOS, mu, beta, r, eta, w, epocas, proba, n0, k_bar, seed)
# y tiene: (promedio sus, promedio inf, suma de promedios (siempre 1))
#SIR_Barabasi(NODOS, mu, beta, r, eta, w, epocas, proba, n0, k_bar, seed)

plt = gplot(x)
draw(PNG("./plots/BA-Graph-50.png", 16cm, 16cm), gplot(x))
