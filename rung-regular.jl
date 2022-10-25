include("modelado.jl")
using Plots
using Compose
using Cairo

NODOS = 50          #Cantidad de nodos
mu = 0.08           #Proba de recuperacion
beta = 0.23         #Proba de contacto con vecino
r = 0.7            #Proba de infección por concatco con vecino
eta = 0.0           #Proba de no contagio por interaccion con vecinos
w=0.02              #Proba de pasar de recuperado a suceptible (solo sir)
epocas = 500        #Iteraciones
k_reg = 4           #Grado de los vertices.
proba = 0.1         #Proba de empezar infectado.
seed = 6
mallado = 50

#Guardamos grafo a tratar.
#x, y = SIS_regular(NODOS, mu, beta, r, eta, w, epocas, proba, k_reg, seed)
#plt = gplot(x)
#draw(PNG("./plots/R-Graph.png", 16cm, 16cm), gplot(x))


#Variamos r
arrer = ones(mallado)
for i in 1:mallado
    x, y = SIS_regular(NODOS, mu, beta, i/mallado, eta, w, epocas, proba, k_reg, seed)
    #Saco proba infectado
    arrer[i]=y[2]
end
plt = plot((1:mallado)./mallado, arrer, dpi=300, lw=3, xlabel="r", ylabel="Infection probability", leg=false, xlims=(0,1), ylims=(0,1))
savefig(plt, "./plots/R-rVar")
#Variamos beta
arreb = ones(mallado)
for i in 1:mallado
    x, y = SIS_regular(NODOS, mu, i/mallado, r, eta, w, epocas, proba, k_reg, seed)
    #Saco proba infectado
    arreb[i]=y[2]
end
plt = plot((1:mallado)./mallado, arreb, dpi=300, lw=3, xlabel="beta", ylabel="Infection probability", leg=false, xlims=(0,1), ylims=(0,1))
savefig(plt, "./plots/R-bVar")

# Variamos r y beta
mtz = zeros((mallado,mallado))
for i in 1:mallado #beta
    for ii in 1:mallado #r
        x, y = SIS_regular(NODOS, mu, i/mallado, ii/mallado, eta, w, epocas, proba, k_reg, seed)
        mtz[i,ii] = y[2]
    end
end
gr()
plt = heatmap(1/mallado:1/mallado:1, 1/mallado:1/mallado:1, mtz,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="r", ylabel="beta", colorbar_title="Infection Probability")

savefig(plt, "./plots/R-Heatmap")
