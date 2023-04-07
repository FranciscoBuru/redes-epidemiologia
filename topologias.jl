using Graphs

#https://juliagraphs.org/Graphs.jl/dev/core_functions/simplegraphs_generators/

#https://juliagraphs.org/Graphs.jl/dev/core_functions/simplegraphs_generators/#Graphs.SimpleGraphs.watts_strogatz-Tuple{Integer,%20Integer,%20Real}
function smallWorld(N::Int64=20, K::Int64=4, beta::Float64=0.4, seed::Int64=1)
    print(seed)
    return watts_strogatz(N, K, beta, seed=seed)
end

#https://juliagraphs.org/Graphs.jl/dev/core_functions/simplegraphs_generators/#Graphs.SimpleGraphs.random_regular_graph-Tuple{Integer,%20Integer}
function regular(N::Int64=20, k::Int64=3, seed::Int64=1)
    return random_regular_graph(N, k, seed=seed)
end

#https://juliagraphs.org/Graphs.jl/dev/core_functions/simplegraphs_generators/#Graphs.SimpleGraphs.barabasi_albert-Tuple{Integer,%20Integer,%20Integer}
function barabasi(N::Int64=20, n0::Int64=2, k::Int64=2, seed::Int64=1)
    return barabasi_albert(N, n0, k, seed = seed)
end

function mexico()
    return SimpleGraph([0 0 1 1; 0 0 1 0; 1 1 0 1; 1 0 1 0])
end

using Random
using StatsBase
using GraphPlot
using Compose
using Cairo

x = regular(6, 4)
display(gplot(x))
draw(PNG("./plots/R-chica.png", 16cm, 16cm), gplot(x))
x = smallWorld(6,4,0.4)
display(gplot(x))
draw(PNG("./plots/SW-chica.png", 16cm, 16cm), gplot(x))
x = barabasi(8,2,2)
display(gplot(x))
draw(PNG("./plots/BAR-chica.png", 16cm, 16cm), gplot(x))
