include("modelado.jl")
#using Plots

NODOS = 10000          #Cantidad de nodos
mu = 0.08           #Proba de recuperacion
beta = 0.23         #Proba de contacto con vecino
r = 0.7             #Proba de infección por concatco con vecino
eta = 0.0           #Proba de no contagio por interaccion con vecinos
w=0.02              #Proba de pasar de recuperado a suceptible (solo sir)
epocas = 20000      #Iteraciones
k_reg = 4           #Grado de los vertices.
proba = 0.9         #Proba de empezar infectado.
seed = 6

x,y=SIS_regular(NODOS, mu, beta, r, eta, w, epocas, proba, k_reg, seed)

#seed = 7
#g2=SIR_regular(NODOS, mu, beta, r, eta, w, epocas, proba, k_reg, seed)
#print(g1==g2)
#SIR_regular(NODOS, mu, beta, r, eta, w, epocas, proba, k_reg, seed)
maximum(broadcast(abs, adjacency_spectrum(x)))

plt = gplot(x)
draw(PNG("./plots/REG-Graph-50.png", 16cm, 16cm), gplot(x))
