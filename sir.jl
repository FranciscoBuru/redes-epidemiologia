using Random
using StatsBase
using Graphs
include("modelo.jl")

mutable struct NodoSIR
    s::Float64      #Proba de estar suceptible
    i::Float64      #Proba de estar infectaod
    rec::Float64    #Proba de estar recuperado
    mu::Float64     #Proba de recuperación
    beta::Float64   #Proba por contacto con vecino
    r::Float64      #Proba de contacto con vecino
    eta::Float64    #Proba de no contagio por interacción con todos los vecinos
    w::Float64      #Proba de pasar de reccuperado a suceptible
    #auxiliares
    iaux::Float64   #Auxilliar de proba de infectado
    saux::Float64   #Auxiliar de proba de suceptible
    raux::Float64   #Auxiliar de proba de infectado

    function NodoSIR(s::Int64, i::Float64, rec::Float64, mu::Float64, beta::Float64,r::Float64, eta::Float64, w::Float64 )
        return new(s,i,rec,mu,beta,r,eta,w,0,0,0)
    end

    function NodoSIR(s::Int64, i::Int64, rec::Float64, mu::Float64, beta::Float64,r::Float64, eta::Float64, w::Float64 )
        return new(s,i,rec,mu,beta,r,eta,w,0,0,0)
    end
end


mutable struct SIR<:Modelo
    grafo::SimpleGraph{Int64}
    nodos::Array{NodoSIR}

    function SIR(grafo::SimpleGraph, nodos::Array{NodoSIR})
        return new(grafo, nodos)
    end
end

function etapa(mod::SIR)
    A = adjacency_matrix(mod.grafo)
    N = length(mod.nodos)
    Threads.@threads for j in 1:N
        mod.nodos[j].eta = 1
        for jj in 1:N
            if(A[j,jj] == 1 && j != jj)
                mod.nodos[j].eta = mod.nodos[j].eta*(1-(mod.nodos[j].r*mod.nodos[j].beta*mod.nodos[jj].i))
            end
        end
    end
    Threads.@threads for j in 1:N
        mod.nodos[j].saux = mod.nodos[j].eta*mod.nodos[j].s + mod.nodos[j].w*mod.nodos[j].rec

        mod.nodos[j].iaux = (1-mod.nodos[j].eta)*mod.nodos[j].s + (1-mod.nodos[j].mu)*mod.nodos[j].i

        mod.nodos[j].raux = (1-mod.nodos[j].w)*mod.nodos[j].rec + mod.nodos[j].mu*mod.nodos[j].i
    end
    Threads.@threads for j in 1:N
        mod.nodos[j].s = mod.nodos[j].saux
        mod.nodos[j].i = mod.nodos[j].iaux
        mod.nodos[j].rec = mod.nodos[j].raux
    end
    return mod
    end

function probas(mod::SIR)
    n = length(mod.nodos)
    r=0
    s=0
    i=0
    for j in 1:n
        r += mod.nodos[j].rec
        s += mod.nodos[j].s
        i += mod.nodos[j].i
    end
    return s/n, i/n, r/n, (s+i+r)/n
end
