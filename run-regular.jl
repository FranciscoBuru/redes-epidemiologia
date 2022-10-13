include("modelado.jl")
#using Plots

NODOS = 50          #Cantidad de nodos
mu = 0.08           #Proba de recuperacion
beta = 0.23         #Proba de contacto con vecino
r = 0.7             #Proba de infección por concatco con vecino
eta = 0.0           #Proba de no contagio por interaccion con vecinos
w=0.02              #Proba de pasar de recuperado a suceptible (solo sir)
epocas = 20000         #Iteraciones
k_reg = 4           #Grado de los vertices.
proba = 0.1         #Proba de empezar infectado.
seed = 6

g1=SIS_regular(NODOS, mu, beta, r, eta, w, epocas, proba, k_reg, seed)

#seed = 7
#g2=SIR_regular(NODOS, mu, beta, r, eta, w, epocas, proba, k_reg, seed)
#print(g1==g2)
#SIR_regular(NODOS, mu, beta, r, eta, w, epocas, proba, k_reg, seed)
