include("modelado.jl")
using Plots
using Compose
using Cairo

NODOS = 50          #Cantidad de nodos
mu = 0.08           #Proba de recuperacion
beta = 0.23         #Proba de contacto con vecino
r = 0.7             #Proba de infección por concatco con vecino
eta = 0.0           #Proba de no contagio por interaccion con vecinos
w=0.02              #Proba de pasar de recuperado a suceptible (solo sir)
epocas = 500        #Iteraciones
k_sw = 4            #Grado de los vertices.
beta_sw = 0.3       #Probabilidad de aristas.
seed = 6
mallado = 50


#Guardamos grafo a tratar.
x, y = SIS_smallWorld(NODOS, mu, beta, r, eta, w, epocas, proba, k_sw, beta_sw, seed)
plt = gplot(x)
draw(PNG("./plots/WS-Graph.png", 16cm, 16cm), gplot(x))

#Variamos r
arrer = ones(mallado)
for i in 1:mallado
    x, y = SIS_smallWorld(NODOS, mu, beta, i/mallado, eta, w, epocas, proba, k_sw, beta_sw, seed)
    #Saco proba infectado
    arrer[i]=y[2]
end
plt = plot((1:mallado)./mallado, arrer, dpi=300, title="r-parameter variation in SIS-WattsStrogratz", lw=3, xlabel="r", ylabel="Infection probability", leg=false, xlims=(0,1), ylims=(0,1))
savefig(plt, "./plots/WS-rVar")

#Variamos beta
arreb = ones(mallado)
for i in 1:mallado
    x, y = SIS_smallWorld(NODOS, mu, i/mallado, r, eta, w, epocas, proba, k_sw, beta_sw, seed)
    #Saco proba infectado
    arreb[i]=y[2]
end
plt = plot((1:mallado)./mallado, arreb, dpi=300, title="beta-parameter variation in SIS-WattsStrogratz", lw=3, xlabel="beta", ylabel="Infection probability", leg=false, xlims=(0,1), ylims=(0,1))
savefig(plt, "./plots/WS-bVar")

# Variamos r y beta
arrerb = ones(mallado*mallado)
xx = ones(mallado*mallado)
yy = ones(mallado*mallado)
for i in 1:mallado
    for ii in 1:mallado
        x, y = SIS_smallWorld(NODOS, mu, i/mallado, ii/mallado, eta, w, epocas, proba, k_sw, beta_sw, seed)
        arrerb[(i-1)*mallado + ii]=y[2]
        xx[(i-1)*mallado + ii] = i/mallado
        yy[(i-1)*mallado + ii] = ii/mallado
    end
end
display(gplot(x))
plt = scatter(xx,yy,arrerb,dpi=300, title="beta-r-variation SIS-WattsStrogratz", xlabel="beta", ylabel="r", zlabel="Infection probability", leg=false, xlims=(0,1), ylims=(0,1), zlims=(0,1))
savefig(plt, "./plots/WS-Scatter")
