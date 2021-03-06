using Random
using StatsBase
using GraphPlot
include("constructor-listas.jl")
include("topologias.jl")


function modelo_mexico(puntos::Integer)

end

function corre_modelo(modelo::Modelo, datos::Array{Float64,2}, epocas::Integer)
    arre = zeros(Float64, epocas,4)
    difs = zeros(Float64, epocas)
    for i in 1:epocas
        modelo = etapa(modelo, i)
        arre[i,:] = valoresActuales(modelo)
        difs[i] = abs(sum(datos[i,:] .- valoresActuales(modelo)))
    end
    #println(probas(modelo))
    return sum(difs), arre
end

function corre_modelo(modelo::Modelo)
    for i in 1:epocas
        modelo = etapa(modelo)
    end
    println(probas(modelo))
    return modelo
end

function grafo_mexico(A::Array{Float64,2}, epocas::Int64, data::Array{Float64,2}, mu::Float64, beta::Array{Float64,2}, eta::Float64, precip::Array{Float64,2})
    grafo = mexico()
    #display(gplot(grafo))
    nodos = construyeListaNodos(A, mu, beta, eta, precip)
    modelo = SISAjustado(grafo,nodos)
    diff, arre = corre_modelo(modelo, data, epocas)
    return diff, arre
end

#No determinístico
function SIS_smallWorld(Nodos::Int64, mu::Float64, beta::Float64, r::Float64, eta::Float64, w::Float64, epocas::Int64, proba::Float64=0.1, k::Int64=2, beta_sw::Float64=0.4, seed::Int64=1)
    grafo = smallWorld(Nodos, k, beta_sw)
    display(gplot(grafo))
    nodos = construyeListaNodos(Nodos, mu, beta, r, eta, w, "SIS", proba, seed)
    modelo = SIS(grafo,nodos)
    modelo = corre_modelo(modelo)
    return grafo
end

#No determinístico
function SIR_smallWorld(Nodos::Int64, mu::Float64, beta::Float64, r::Float64, eta::Float64, w::Float64, epocas::Int64, proba::Float64=0.1, k::Int64=2, beta_sw::Float64=0.4, seed::Int64=1)
    grafo = smallWorld(Nodos, k, beta_sw, seed)
    display(gplot(grafo))
    nodos = construyeListaNodos(Nodos, mu, beta, r, eta, w, "SIR", proba, seed)
    modelo = SIR(grafo,nodos)
    modelo = corre_modelo(modelo)
    return grafo
end

function SIS_regular(Nodos::Int64, mu::Float64, beta::Float64, r::Float64, eta::Float64, w::Float64, epocas::Int64, proba::Float64=0.1, k::Int64=2, seed::Int64=1)
    grafo = regular(Nodos, k, seed)
    display(gplot(grafo))
    nodos = construyeListaNodos(Nodos, mu, beta, r, eta, w, "SIS", proba, seed)
    modelo = SIS(grafo,nodos)
    modelo = corre_modelo(modelo)
    return grafo
end

function SIR_regular(Nodos::Int64, mu::Float64, beta::Float64, r::Float64, eta::Float64, w::Float64, epocas::Int64, proba::Float64=0.1, k::Int64=2, seed::Int64=1)
    grafo = regular(Nodos, k, seed)
    display(gplot(grafo))
    nodos = construyeListaNodos(Nodos, mu, beta, r, eta, w, "SIR", proba, seed)
    modelo = SIR(grafo,nodos)
    modelo = corre_modelo(modelo)
    return grafo
end


function SIS_Barabasi(Nodos::Int64, mu::Float64, beta::Float64, r::Float64, eta::Float64, w::Float64, epocas::Int64, proba::Float64=0.1, n0::Int64=2, k::Int64=2, seed::Int64=1)
    grafo = barabasi(Nodos, n0, k, seed)
    display(gplot(grafo))
    nodos = construyeListaNodos(Nodos, mu, beta, r, eta, w, "SIS", proba, seed)
    modelo = SIS(grafo,nodos)
    modelo = corre_modelo(modelo)
    return grafo
end

function SIR_Barabasi(Nodos::Int64, mu::Float64, beta::Float64, r::Float64, eta::Float64, w::Float64, epocas::Int64, proba::Float64=0.1, n0::Int64=2, k::Int64=2, seed::Int64=1)
    grafo = barabasi(Nodos, n0, k, seed)
    display(gplot(grafo))
    nodos = construyeListaNodos(Nodos, mu, beta, r, eta, w, "SIR", proba, seed)
    modelo = SIR(grafo,nodos)
    modelo = corre_modelo(modelo)
    return grafo
end
