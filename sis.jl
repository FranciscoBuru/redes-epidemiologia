using Random
using StatsBase
using Graphs
include("modelo.jl")

mutable struct NodoSIS
    s::Float64      #Proba de estar suceptible
    i::Float64      #Proba de estar infectado
    mu::Float64     #Proba de recuperación
    beta::Float64   #Proba de infeccion por contacto con vecino
    r::Float64      #Proba de contacto con vecino
    eta::Float64    #Proba de no contagio por interacción con todos los vecinos
    #auxiliares
    iaux::Float64   #Auxilliar de proba de infectado
    saux::Float64   #Auxiliar de proba de suceptible

    function NodoSIS(s::Float64, i::Float64, mu::Float64, beta::Float64,r::Float64, eta::Float64)
        return new(s,i,mu,beta,r,eta,0,0)
    end

    function NodoSIS(s::Int64, i::Int64, mu::Float64, beta::Float64,r::Float64, eta::Float64)
        return new(s,i,mu,beta,r,eta,0,0)
    end
end

mutable struct SIS<:Modelo
    grafo::SimpleGraph{Int64}
    nodos::Array{NodoSIS}

    function SIS(grafo::SimpleGraph, nodos::Array{NodoSIS})
        return new(grafo, nodos)
    end
end


function etapa(mod::SIS)
    A = adjacency_matrix(mod.grafo)
    N = length(mod.nodos)
    Threads.@threads for j in 1:N
        mod.nodos[j].eta = 1
        for jj in 1:N
            if(A[j,jj] == 1 && j!=jj)
                mod.nodos[j].eta = mod.nodos[j].eta*(1-(mod.nodos[j].r*mod.nodos[j].beta*mod.nodos[jj].i))
            end
        end
    end
    Threads.@threads for j in 1:N
        mod.nodos[j].iaux = ((1-mod.nodos[j].mu)*mod.nodos[j].i) + ((1-mod.nodos[j].eta)*mod.nodos[j].s)
        mod.nodos[j].saux = (mod.nodos[j].mu*mod.nodos[j].i) + (mod.nodos[j].eta*mod.nodos[j].s)
    end
    Threads.@threads for j in 1:N
        mod.nodos[j].i = mod.nodos[j].iaux
        mod.nodos[j].s = mod.nodos[j].saux
    end
    return mod
end

function probas(mod::SIS)
    n = length(mod.nodos)
    s=0
    i=0
    for j in 1:n
        s += mod.nodos[j].s
        i += mod.nodos[j].i
    end
    return s/n, i/n, (s+i)/n
end
