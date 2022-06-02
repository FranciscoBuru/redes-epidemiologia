include("sis.jl")
include("sir.jl")
include("sisAjustado.jl")


function construyeListaNodos(N::Int64, mu::Float64, beta::Float64, r::Float64, eta::Float64, w::Float64, tipo::String="SIS", proba::Float64=0.1, seed::Int64=1)
    rng = MersenneTwister(seed)
    seeds=sample(rng, Vector(1:1000), N)
    if tipo == "SIS"
        arre = Array{NodoSIS}(undef, N)
        Threads.@threads for j in 1:N
            in_rng = MersenneTwister(seeds[j])
            infectado = 0
            if rand(in_rng,Float64) < proba
                infectado = 1
            end
            arre[j] = NodoSIS(1-infectado, infectado, mu, beta, r, eta)
        end
    elseif tipo =="SIR"
        arre = Array{NodoSIR}(undef, N)
        Threads.@threads for j in 1:N
            in_rng = MersenneTwister(seeds[j])
            infectado = 0
            if rand(in_rng,Float64) < proba
                infectado = 1
            end
            arre[j] = NodoSIR(1-infectado, infectado, 0.0, mu, beta, r, eta, w)
        end
    end
    return arre
end

function construyeListaNodos(A::Array{Float64,2}, mu::Float64, beta::Array{Float64,1}, eta::Float64, precip::Array{Float64,2})
    arre = Array{NodoSISAjustado}(undef, 4)
    for j in 1:4
        if j == 1
            arre[j] = NodoSISAjustado(1, 0, mu, beta[j], A[j,:], eta, precip[j,:])
        else
            arre[j] = NodoSISAjustado(1-0.00000042, 0.00000042, mu, beta[j], A[j,:], eta, precip[j,:])
        end
    end
    return arre
end
