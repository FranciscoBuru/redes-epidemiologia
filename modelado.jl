using Random
using StatsBase
using GraphPlot
include("constructor-listas.jl")
include("topologias.jl")
#using .topologias

NODOS = 1000
mu = 0.08
beta = 0.23
r = 0.7
eta = 0.0
w=0.02
epocas = 5
seed = 2
#barabasi Opcionales
n0 = 4
k_bar = 3
proba = 0.1

#regular Opcionales
k_reg = 4

#small-world
k_sw = 4
beta_sw = 0.4

#No determinístico
function SIS_smallWorld(Nodos::Int64, mu::Float64, beta::Float64, r::Float64, eta::Float64, w::Float64, epocas::Int64, proba::Float64=0.1, k::Int64=2, beta_sw::Float64=0.4, seed::Int64=1)
    grafo = smallWorld(Nodos, k, beta_sw)
    #grafo = watts_strogatz(Nodos, k, beta_sw, seed=seed)
    nodos = construyeListaNodos(Nodos, mu, beta, r, eta, w, "SIS", proba, seed)
    #print(nodos)
    modelo = SIS(grafo,nodos)
    for i in 1:epocas
        println(probas(modelo))
        modelo = etapa(modelo)
    end
    return grafo
end

#No determinístico
function SIR_smallWorld(Nodos::Int64, mu::Float64, beta::Float64, r::Float64, eta::Float64, w::Float64, epocas::Int64, proba::Float64=0.1, k::Int64=2, beta_sw::Float64=0.4, seed::Int64=1)
    grafo = smallWorld(Nodos, k, beta_sw, seed)
    #grafo = watts_strogatz(Nodos, k, beta_sw, seed=seed)
    nodos = construyeListaNodos(Nodos, mu, beta, r, eta, w, "SIR", proba, seed)
    #print(nodos)
    modelo = SIR(grafo,nodos)
    for i in 1:epocas
        println(probas(modelo))
        modelo = etapa(modelo)
    end
    return grafo
end

function SIS_regular(Nodos::Int64, mu::Float64, beta::Float64, r::Float64, eta::Float64, w::Float64, epocas::Int64, proba::Float64=0.1, k::Int64=2, seed::Int64=1)
    grafo = regular(Nodos, k, seed)
    #grafo = random_regular_graph(Nodos, k, seed=seed)
    nodos = construyeListaNodos(Nodos, mu, beta, r, eta, w, "SIS", proba, seed)
    #print(nodos)
    modelo = SIS(grafo,nodos)
    for i in 1:epocas
        println(probas(modelo))
        modelo = etapa(modelo)
    end
    return grafo
end

function SIR_regular(Nodos::Int64, mu::Float64, beta::Float64, r::Float64, eta::Float64, w::Float64, epocas::Int64, proba::Float64=0.1, k::Int64=2, seed::Int64=1)
    grafo = regular(Nodos, k, seed)
    #grafo = random_regular_graph(Nodos, k, seed=seed)
    nodos = construyeListaNodos(Nodos, mu, beta, r, eta, w, "SIR", proba, seed)
    #print(nodos)
    modelo = SIR(grafo,nodos)
    for i in 1:epocas
        println(probas(modelo))
        modelo = etapa(modelo)
    end
    return grafo
end


function SIS_Barabasi(Nodos::Int64, mu::Float64, beta::Float64, r::Float64, eta::Float64, w::Float64, epocas::Int64, proba::Float64=0.1, n0::Int64=2, k::Int64=2, seed::Int64=1)
    grafo = barabasi(Nodos, n0, k, seed)
    #grafo = barabasi_albert(Nodos, n0, k, seed = seed)
    nodos = construyeListaNodos(Nodos, mu, beta, r, eta, w, "SIS", proba, seed)
    #print(nodos)
    modelo = SIS(grafo,nodos)
    for i in 1:epocas
        println(probas(modelo))
        modelo = etapa(modelo)
    end
    return grafo
end

function SIR_Barabasi(Nodos::Int64, mu::Float64, beta::Float64, r::Float64, eta::Float64, w::Float64, epocas::Int64, proba::Float64=0.1, n0::Int64=2, k::Int64=2, seed::Int64=1)
    grafo = barabasi(Nodos, n0, k, seed)
    #grafo = barabasi_albert(Nodos, n0, k, seed = seed)
    nodos = construyeListaNodos(Nodos, mu, beta, r, eta, w, "SIR", proba, seed)
    #print(nodos)
    modelo = SIR(grafo,nodos)
    for i in 1:epocas
        println(probas(modelo))
        modelo = etapa(modelo)
    end
    return grafo
end

gplot(SIS_smallWorld(NODOS, mu, beta, r, eta, w, epocas, proba, k_sw, beta_sw, seed))

gplot(SIR_smallWorld(NODOS, mu, beta, r, eta, w, epocas, proba, k_sw, beta_sw, seed))

gplot(SIS_Barabasi(NODOS, mu, beta, r, eta, w, epocas, proba, n0, k_bar, seed))

gplot(SIR_Barabasi(NODOS, mu, beta, r, eta, w, epocas, proba, n0, k_bar, seed))

gplot(SIR_regular(NODOS, mu, beta, r, eta, w, epocas, proba, k_reg, seed))

gplot(SIS_regular(NODOS, mu, beta, r, eta, w, epocas, proba, k_reg, seed))

#=
Pruebas
gplot(barabasi(20, 1, 1, 1))

grafo = smallWorld(5, k_sw, beta_sw, 1)

grafo1 = watts_strogatz(5, 4, 0.3, seed=2)
grafo2 = watts_strogatz(5, 4, 0.3, seed=2)
print(grafo1==grafo2)
#grafo = smallWorld(20, k_sw, beta_sw, 1)

#grafo = watts_strogatz(Nodos, k, beta_sw, seed=seed)
nodos = construyeListaNodos(20, mu, beta, r, eta, "SIS", proba, 1)
#print(nodos)
modelo = SIS(grafo,nodos)
for i in 1:epocas
    println(probas(modelo))
    modelo = etapa(modelo)
end
=#
