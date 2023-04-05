include("modelado.jl")
using Plots
using Compose
using Cairo

NODOS = 500          #Cantidad de nodos
mu = 0.08           #Proba de recuperacion
beta = 0.23         #Proba de contacto con vecino
r = 0.7             #Proba de infección por concatco con vecino
eta = 0.0           #Proba de no contagio por interaccion con vecinos
w=0.02              #Proba de pasar de recuperado a suceptible (solo sir)
epocas = 500        #Iteraciones
n0 = 4              #Nodos del grafo incial.
k_bar = 4           #Vecinos de cada nodo.
proba = 0.1         #Proba de empezar infectado.
seed = 6
mallado = 50

#Guardamos grafo a tratar.
x, y = SIS_Barabasi(NODOS, mu, beta, r, eta, w, epocas, proba, n0, k_bar, seed)
#plt = gplot(x)
#draw(PNG("./plots/BA-Graph.png", 16cm, 16cm), gplot(x))
maximum(broadcast(abs, adjacency_spectrum(x)))
#WW = adjacency_matrix(x[, T=Int; dir=:out])
#eigvals(adjacency_matrix(WW; permute::Bool=true, scale::Bool=true, sortby)

#Variamos r
arrer = ones(mallado)
for i in 1:mallado
    x, y = SIS_Barabasi(NODOS, mu, beta, i/mallado, eta, w, epocas, proba, n0, k_bar, seed)
    #Saco proba infectado
    arrer[i]=y[2]
end
plt = plot((1:mallado)./mallado, arrer, dpi=300, lw=3, xlabel="r", ylabel="Infection probability", leg=false, xlims=(0,1), ylims=(0,1))
savefig(plt, "./plots/BA-rVar")

#Variamos beta
arreb = ones(mallado)
for i in 1:mallado
    x, y = SIS_Barabasi(NODOS, mu, i/mallado, r, eta, w, epocas, proba, n0, k_bar, seed)
    #Saco proba infectado
    arreb[i]=y[2]
end
plt = plot((1:mallado)./mallado, arreb, dpi=300, lw=3, xlabel="beta", ylabel="Infection probability", leg=false, xlims=(0,1), ylims=(0,1))
savefig(plt, "./plots/BA-bVar")

# Variamos r y beta
mtz = zeros((mallado,mallado))
for i in 1:mallado
    for ii in 1:mallado
        x, y = SIS_Barabasi(NODOS, mu, i/mallado, ii/mallado, eta, w, epocas, proba, n0, k_bar, seed)
        mtz[i,ii] = y[2]
    end
end
#display(gplot(x))
gr()
plt = heatmap(1/mallado:1/mallado:1, 1/mallado:1/mallado:1, mtz,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="r", ylabel="beta", colorbar_title="Infection Probability")


#plt = scatter(xx,yy,arrerb,dpi=300, xlabel="beta", ylabel="r", zlabel="Infection probability", leg=false, xlims=(0,1), ylims=(0,1), zlims=(0,1))
savefig(plt, "./plots/BA-Heatmap")
