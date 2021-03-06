using Graphs
#https://juliagraphs.org/Graphs.jl/stable/generators/#All-Generators

#https://juliagraphs.org/Graphs.jl/stable/generators/#Graphs.SimpleGraphs.watts_strogatz-Tuple{Integer,%20Integer,%20Real}
function smallWorld(N::Int64=20, K::Int64=4, beta::Float64=0.4, seed::Int64=1)
    return watts_strogatz(N, K, beta, seed=seed)
end

#https://juliagraphs.org/Graphs.jl/stable/generators/#Graphs.SimpleGraphs.random_regular_graph-Tuple{Integer,%20Integer}
function regular(N::Int64=20, k::Int64=3, seed::Int64=1)
    return random_regular_graph(N, k, seed=seed)
end

#https://juliagraphs.org/Graphs.jl/stable/generators/#Graphs.SimpleGraphs.barabasi_albert-Tuple{Integer,%20Integer,%20Integer}
function barabasi(N::Int64=20, n0::Int64=2, k::Int64=2, seed::Int64=1)
    return barabasi_albert(N, n0, k, seed = seed)
end

function mexico()
    return SimpleGraph([0 0 1 1; 0 0 1 0; 1 1 0 1; 1 0 1 0])
end
