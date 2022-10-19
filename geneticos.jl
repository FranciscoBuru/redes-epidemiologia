using LinearAlgebra
using Random
using StatsBase
using Distributions
using Plots

abstract type Genoma end


function mutar!(x::Array{Float64}, Pm::Float64, rng::MersenneTwister)
    n=length(x)
    indiceMutar=rand(rng, 1:12)
    x[indiceMutar] = rand(rng,Uniform(0,1))
end


mutable struct Poblacion{T<:Genoma}
    pob::Array{T} #Poblacion en si
    rng::MersenneTwister #Seed
    funCalif::Function #Funcion calificadora, debe ser positiva entera y a maximizar
    keepbest::Bool #Numero de genomas que se guardan de la poblacion pasada
    random::Int64 #Numero de genomas que se inicializan al azar cada vuelta
    Pm::Float64 #Probabilidad de mutacion
    #Constructor entero
    function Poblacion{T}(n::Int64, funCalif::Function; keepbest::Bool=true, random::Int64=2, Pm::Float64=0.05, seed::Int64=1234) where T<:Genoma
        rng=MersenneTwister(seed)
        pob=Array{T}(undef, n)
        for i in 1:n
            pob[i]=T(rng, funCalif)
        end
        new{T}(pob, rng, funCalif, keepbest, random, Pm)
    end
    #Constructor default
    function Poblacion{T}(pob::Array{T}, rng::MersenneTwister, funCalif::Function, keepbest::Bool, random::Int64, Pm::Float64) where T<:Genoma
        new(pob, rng, funCalif, keepbest, random, Pm)
    end
end

function rouletteselector(pob::Array{T}, n::Int64, rng::MersenneTwister) where T<:Genoma
    f=x->x.calif; w=f.(pob)
    W=FrequencyWeights(-w.+(maximum(w)+1))
    sample(rng, pob, W, n; replace=false)
end

function reproduce(p::Poblacion{T}) where T<:Genoma
    n=length(p.pob)
    pob=Array{T}(undef, n)
    m=n-p.keepbest-p.random
    seeds=sample(p.rng, Vector(1:1000), ceil(Int,m/2))
    Threads.@threads for i in 1:ceil(Int,m/2)
        rng=MersenneTwister(seeds[i])
        s=rouletteselector(p.pob, 2, rng)
        c=crossover(s[1], s[2], rng, p.funCalif, p.Pm)
        #(p1::OX, p2::OX, rng::MersenneTwister, funCalif::Function, Pm::Float64)
        pob[2*i-1]=c[1]
        if 2*i<=m
            pob[2*i]=c[2]
        end
    end
    if p.random>0
        for i in 1:p.random
            pob[m+i]=T(p.rng, p.funCalif)
        end
    end
    if p.keepbest
        f=x->x.calif
        pob[n]=p.pob[findfirst(isequal(minimum(f,p.pob)),f.(p.pob))]
    end
    Poblacion{T}(pob, p.rng, p.funCalif, p.keepbest, p.random, p.Pm)
end


function getbest(pob::Poblacion)
    f=x->x.calif;
    return pob.pob[findfirst(isequal(minimum(f,pob.pob)),f.(pob.pob))]
end

function algoritmoGenetico(funCalif::Function, tipo, pobsize, generations; Pm=0.01, intStart=0, intEnd=0, binlen=10, seed=1234, keepbest=true,random=2)
    if tipo=="WAC"
        pob=Poblacion{WAC}(pobsize, funCalif)#; keepbest=keepbest, random=random, Pm=Pm, seed=seed)
        #WAC(rng::MersenneTwister, funCalif::Function)
    end

    for i in 2:generations
        println(getbest(pob).calif)
        pob=reproduce(pob)
        #println(i)
    end
    res=getbest(pob)
    return res.calif, res.genoma
end




function algoritmoGeneticoTiempo(funCalif::Function, tipo, pobsize, tiempo_sec; Pm=0.01, best=718, intStart=0, intEnd=0, binlen=10, seed=1234, keepbest=true,random=2)
    if tipo=="PMX"
        pob=Poblacion{PMX}(pobsize, intStart, intEnd, funCalif)#; keepbest=keepbest, random=random, Pm=Pm, seed=seed)
        #Poblacion{T}(n::Int64, starter::Int64, ending::Int64, funCalif::Function; keepbest::Bool=true, random::Int64=2, Pm::Float64=0.05, seed::Int64=1234)
    elseif tipo=="OX"
        pob=Poblacion{OX}(pobsize, intStart, intEnd, funCalif; keepbest=keepbest, random=random, Pm=Pm, seed=seed)
    elseif tipo=="CX"
        # println("Dentro")
        pob=Poblacion{CX}(pobsize, intStart, intEnd, funCalif; keepbest=keepbest, random=random, Pm=Pm, seed=seed)
    end
    ini =time_ns()
    dif=0
    i=1
    while dif<tiempo_sec
        #println(i,getbest(pob))
        pob=reproduce(pob)
        #println(i)
        dif=(time_ns()-ini)/1e9
        i+=1
        if getbest(pob).calif<best
            println(getbest(pob).calif)
            best = getbest(pob).calif
        end
    end
    res=getbest(pob)
    return res.calif, res.genoma, i
end


# Ini WAC
mutable struct WAC<:Genoma
    genoma::Array{Float64}
    calif::Float64
    function WAC(rng::MersenneTwister, funCalif::Function)
        genoma=vec(rand(rng, Uniform(0,1),12,1))
        calif=funCalif(genoma)
        new(genoma, calif)
    end
    function WAC(genoma::Array{Float64}, funCalif::Function)
        new(genoma, funCalif(genoma))
    end
end

function crossover(p1::WAC, p2::WAC, rng::MersenneTwister, funCalif::Function, Pm::Float64)
    n=length(p1.genoma)
    a = rand(rng, Uniform(0,1))
    # Primer Hijo, sobreescribimos en y
    x1 = a.*p1.genoma + (1-a).*p2.genoma
    x2 = a.*p2.genoma + (1-a).*p1.genoma

    # Segundo Hijo


    #Mutamos
    mutar!(x1, Pm, rng); mutar!(x2, Pm, rng)
    return WAC(x2,funCalif), WAC(x1,funCalif)
end
# Fin WAC
