include("modelado.jl")


NODOS = 15          #Cantidad de nodos
mu = 0.08           #Proba de recuperacion
beta = 0.23         #Proba de contacto con vecino
r = 0.7             #Proba de infección por concatco con vecino
eta = 0.0           #Proba de no contagio por interaccion con vecinos
w=0.02              #Proba de pasar de recuperado a suceptible (solo sir)
epocas = 30         #Iteraciones
seed = 2
#barabasi Opcionales
n0 = 4              #Nodos del grafo incial.
k_bar = 3           #Vecinos de cada nodo.
proba = 0.1         #Proba de empezar infectado.

#regular Opcionales
k_reg = 4           #Grado de los vertices.

#small-world
k_sw = 4            #Grado de los vertices.
beta_sw = 0.4       #Probabilidad de aristas.


SIS_smallWorld(NODOS, mu, beta, r, eta, w, epocas, proba, k_sw, beta_sw, seed)

SIR_smallWorld(NODOS, mu, beta, r, eta, w, epocas, proba, k_sw, beta_sw, seed)

SIS_Barabasi(NODOS, mu, beta, r, eta, w, epocas, proba, n0, k_bar, seed)

SIR_Barabasi(NODOS, mu, beta, r, eta, w, epocas, proba, n0, k_bar, seed)

SIR_regular(NODOS, mu, beta, r, eta, w, epocas, proba, k_reg, seed)

SIS_regular(NODOS, mu, beta, r, eta, w, epocas, proba, k_reg, seed)


#=
#Pruebas, muestra no determinisitico small world, con 5 nodos da true pero
#con 20 nodos da false.

grafo1 = watts_strogatz(5, 4, 0.3, seed=2)
grafo2 = watts_strogatz(5, 4, 0.3, seed=2)
print(grafo1==grafo2)
=#
