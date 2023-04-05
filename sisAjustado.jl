using Random
using StatsBase
using Graphs
include("modelo.jl")

mutable struct NodoSISAjustado
    s::Float64      #Proba de estar suceptible
    i::Float64      #Proba de estar infectaod
    mu::Float64     #Proba de recuperación
    beta::Array{Float64,1}   #Proba por contacto con vecino
    r::Array{Float64,1}      #Proba de contacto con vecino
    eta::Float64    #Proba de no contagio por interacción con todos los vecinos
    precip::Array{Float64,1} #Precipitacion.
    #auxiliares
    iaux::Float64   #Auxilliar de proba de infectado
    saux::Float64   #Auxiliar de proba de suceptible

    function NodoSISAjustado(s::Float64, i::Float64, mu::Float64, beta::Array{Float64,1},r::Array{Float64,1}, eta::Float64, precip::Array{Float64,1})
        return new(s,i,mu,beta,r,eta,precip,0,0)
    end

    function NodoSISAjustado(s::Int64, i::Int64, mu::Float64, beta::Array{Float64,1},r::Array{Float64,1}, eta::Float64, precip::Array{Float64,1})
        return new(s,i,mu,beta,r,eta,precip, 0,0)
    end
end

mutable struct SISAjustado<:Modelo
    grafo::SimpleGraph{Int64}
    nodos::Array{NodoSISAjustado,1}

    function SISAjustado(grafo::SimpleGraph, nodos::Array{NodoSISAjustado, 1})
        return new(grafo, nodos)
    end
end


function etapa(mod::SISAjustado, semana::Integer)
    A = adjacency_matrix(mod.grafo)
    N = length(mod.nodos)
    sem = semana
    if (semana > 52 && semana<=104)
            sem = semana-52;
    end
    Threads.@threads for j in 1:N
        mod.nodos[j].eta = 1
        for jj in 1:N
            if(A[j,jj] == 1 || j==jj)
                # Lluvia y beta del año
                mod.nodos[j].eta = mod.nodos[j].eta*(1-(mod.nodos[j].precip[Int(floor(semana/4.333333)+1)]*mod.nodos[j].r[jj]*mod.nodos[j].beta[Int(floor(semana/52)+1)]*mod.nodos[jj].i))
                #Luvia del año
                #mod.nodos[j].eta = mod.nodos[j].eta*(1-(mod.nodos[j].precip[Int(floor(semana/4.333333)+1)]*mod.nodos[j].r[jj]*mod.nodos[j].beta[1]*mod.nodos[jj].i))
                #Beta del año
                #mod.nodos[j].eta = mod.nodos[j].eta*(1-(mod.nodos[j].precip[Int(floor(sem/4.333333)+1)]*mod.nodos[j].r[jj]*mod.nodos[j].beta[Int(floor(semana/52)+1)]*mod.nodos[jj].i))
                #Solo los pasados
                #mod.nodos[j].eta = mod.nodos[j].eta*(1-(mod.nodos[j].precip[Int(floor(sem/4.333333)+1)]*mod.nodos[j].r[jj]*mod.nodos[j].beta[1]*mod.nodos[jj].i))

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

function probas(mod::SISAjustado)
    n = length(mod.nodos)
    s=0
    i=0
    for j in 1:n
        s += mod.nodos[j].s
        i += mod.nodos[j].i
    end
    return s/n, i/n, (s+i)/n
end

function valoresActuales(mod::SISAjustado)
    return [mod.nodos[1].i, mod.nodos[2].i, mod.nodos[3].i, mod.nodos[4].i]
end
