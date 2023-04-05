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

function reproduceH2(p::Poblacion{T}) where T<:Genoma
    f=x->x.calif
    indexAlpha = findfirst(isequal(minimum(f,p.pob)),f.(p.pob))
    alpha = p.pob[indexAlpha]
    n=length(p.pob)
    pob=Array{T}(undef, n)
    m=n-p.keepbest-p.random
    seeds=sample(p.rng, Vector(1:1000), ceil(Int,m/2))
    Threads.@threads for i in 1:ceil(Int,m/2)
        rng=MersenneTwister(seeds[i])
        c=crossover(alpha, p.pob[rand(rng,1:10)], rng, p.funCalif, p.Pm)
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


function reproduceH(p::Poblacion{T}) where T<:Genoma
    f=x->x.calif
    indexAlpha = findfirst(isequal(minimum(f,p.pob)),f.(p.pob))
    alpha = p.pob[indexAlpha]
    n=length(p.pob)
    pob=Array{T}(undef, n*2 - 1)
    seeds=sample(p.rng, Vector(1:1000), ceil(Int,n))
    Threads.@threads for i in 1:(indexAlpha-1)
        rng=MersenneTwister(seeds[i])
        c=crossover(alpha, p.pob[i], rng, p.funCalif, p.Pm)
        #(p1::OX, p2::OX, rng::MersenneTwister, funCalif::Function, Pm::Float64)
        pob[2*(i)-1]=c[1]
        pob[2*(i)]=c[2]
    end
    Threads.@threads for i in (indexAlpha+1):ceil(Int,n)
        rng=MersenneTwister(seeds[i])
        c=crossover(alpha, p.pob[i], rng, p.funCalif, p.Pm)
        #(p1::OX, p2::OX, rng::MersenneTwister, funCalif::Function, Pm::Float64)
        pob[2*(i-1)-1]=c[1]
        pob[2*(i-1)]=c[2]
    end
    pob[n*2 - 1] = alpha
    #Ordena de mejor a peor. sort!(v, by = x -> x[1]);
    sort!(pob, by = x -> x.calif)

    if p.random>0
        for i in 1:p.random
            pob[10-(i-1)]=T(p.rng, p.funCalif)
        end
    end
    Poblacion{T}(pob[1:10], p.rng, p.funCalif, p.keepbest, p.random, p.Pm)
end


function getbest(pob::Poblacion)
    f=x->x.calif;
    return pob.pob[findfirst(isequal(minimum(f,pob.pob)),f.(pob.pob))]
end

function algoritmoGenetico(funCalif::Function, tipo, pobsize, generations; Pm=0.01, intStart=0, intEnd=0, binlen=10, seed=1234, keepbest=true,random=2)
    if tipo=="WAC"
        pob=Poblacion{WAC}(pobsize, funCalif, random=seed)#; keepbest=keepbest, random=random, Pm=Pm, seed=seed)
        #WAC(rng::MersenneTwister, funCalif::Function)
    elseif tipo=="SAC"
        pob=Poblacion{SAC}(pobsize, funCalif, random=seed)#; keepbest=keepbest, random=random, Pm=Pm, seed=seed)
    end

    for i in 2:generations
        println(i)
        pob=reproduce(pob)
        #pob=reproduceH(pob)  #SDC
        #pob=reproduceH2(pob) #SDA
        println(getbest(pob).calif)
    end
    res=getbest(pob)
    return res.calif, res.genoma
end




function algoritmoGeneticoTiempo(funCalif::Function, tipo, pobsize, tiempo_sec; Pm=0.01, best=718, intStart=0, intEnd=0, binlen=10, seed=1234, keepbest=true,random=2)
    if tipo=="WAC"
        pob=Poblacion{WAC}(pobsize, funCalif, random=seed)#; keepbest=keepbest, random=random, Pm=Pm, seed=seed)
        #WAC(rng::MersenneTwister, funCalif::Function)
    elseif tipo=="SAC"
        pob=Poblacion{SAC}(pobsize, funCalif, random=seed)#; keepbest=keepbest, random=random, Pm=Pm, seed=seed)
    end
    ini =time_ns()
    dif=0
    i=1
    while dif<tiempo_sec
        #println(i,getbest(pob))
        #pob=reproduce(pob)
        pob=reproduceH(pob)  #SDC
        #pob=reproduceH2(pob) #SDA **
        println(getbest(pob).calif)
        dif=(time_ns()-ini)/1e9
        i+=1
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


# Ini SAC
mutable struct SAC<:Genoma
    genoma::Array{Float64}
    calif::Float64
    function SAC(rng::MersenneTwister, funCalif::Function)
        genoma=vec(rand(rng, Uniform(0,1),12,1))
        calif=funCalif(genoma)
        new(genoma, calif)
    end
    function SAC(genoma::Array{Float64}, funCalif::Function)
        new(genoma, funCalif(genoma))
    end
end

function crossover(p1::SAC, p2::SAC, rng::MersenneTwister, funCalif::Function, Pm::Float64)
    n=length(p1.genoma)
    a = rand(rng, 1:12)
    # Primer Hijo, sobreescribimos en y
    x1 = vcat(p1.genoma[1:a-1],p2.genoma[a:12])
    x2 = vcat(p2.genoma[1:a-1],p1.genoma[a:12])

    #Mutamos
    mutar!(x1, Pm, rng); mutar!(x2, Pm, rng)
    return SAC(x2,funCalif), SAC(x1,funCalif)
end
# Fin SAC
