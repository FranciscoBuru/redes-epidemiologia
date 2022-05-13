include("sis.jl")
include("sir.jl")


function construyeListaNodos(N::Int64, mu::Float64, beta::Float64, r::Float64, eta::Float64, tipo::String="SIS", proba::Float64=0.1, seed::Int64=1)
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
            arre[j] = NodoSIR(1-infectado, infectado, 0.0, mu, beta, r, eta)
        end
    end
    return arre
end
