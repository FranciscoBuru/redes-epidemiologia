include("modelado.jl")
#using Plots

NODOS = 50          #Cantidad de nodos
mu = 0.08           #Proba de recuperacion
beta = 0.23         #Proba de contacto con vecino
r = 0.7             #Proba de infección por concatco con vecino
eta = 0.0           #Proba de no contagio por interaccion con vecinos
w=0.02              #Proba de pasar de recuperado a suceptible (solo sir)
epocas = 20000         #Iteraciones
k_sw = 4            #Grado de los vertices.
beta_sw = 0.3       #Probabilidad de aristas.
seed = 6

grafo1=SIS_smallWorld(NODOS, mu, beta, r, eta, w, epocas, proba, k_sw, beta_sw, seed)

#seed=9
#grafo2=SIS_smallWorld(NODOS, mu, beta, r, eta, w, epocas, proba, k_sw, beta_sw, seed)
#print(grafo1==grafo2)

#SIR_smallWorld(NODOS, mu, beta, r, eta, w, epocas, proba, k_sw, beta_sw, seed)
