include("geneticos.jl")
include("modelado.jl")
using Plots
using DelimitedFiles


# Beta primer anio.
beta = [0.231857 0.427434 0.010000 0.558327 0.221011 0.294170;
        0.145923 0.396093 0.010000 0.397411 0.010000 0.334287;
        0.221471 0.310545 0.420411 0.570422 0.236681 0.190461;
        0.576013 0.596109 0.568935 0.010000 0.169991 0.356684] #5

# Datos de casos
data = readdlm("./data.txt")

#Precipitacion primer anio
precip = [0.24940130 0.13479302 0.15121451 0.18918919 0.64762231 0.80465275 0.76565173 0.55319877 1.00000000 0.83099555 0.27129661 0.14471433 0.03838534;
          0.15494466 0.00071403 0.01285255 0.00963941 0.30024991 0.86826134 0.82506248 1.00000000 0.92466976 0.53052481 0.02392003 0.00464120 0.03721910;
          0.12779813 0.08424908 0.09686610 0.15873016 0.51037851 1.00000000 0.79894180 0.78591779 0.91982092 0.64794465 0.12779813 0.08506309 0.06378927;
          0.23616089 0.13898845 0.28793309 0.29669454 0.59418558 1.00000000 0.61967344 0.78733572 0.88410992 0.79291119 0.28474711 0.22182397 0.14065282]

epocas = 52    # Numero de semanas a ver. 52 es de 1 anio.
mu = 0.1428     # Probabilidad de recuperacion.
eta=0.0

function costoRuta(M, epocas, data, mu, beta, eta, precip)
    A = [M[1] 0 M[2] M[3]; 0 M[4] M[5] 0; M[6] M[7] M[8] M[9]; M[10] 0 M[11] M[12]]
    A = convert(Array{Float64,2}, A)
    costo, kk, kk2  = grafo_mexico(A, epocas, data, mu, beta, eta, precip,0)
    return costo
end


calif=x->costoRuta(x, epocas, data, mu, beta, eta, precip)

genes = 10000
println(genes)
ini = time_ns()
califFinal, genomaFinal = algoritmoGenetico(calif, "WAC", 10, genes; intStart=1, intEnd=100, random=2, seed=2, keepbest=true)
print((time_ns() - ini )/ 1e9)
print(genomaFinal)


##Grafica en grid.
A = [genomaFinal[1] 0 genomaFinal[2] genomaFinal[3]; 0 genomaFinal[4] genomaFinal[5] 0; genomaFinal[6] genomaFinal[7] genomaFinal[8] genomaFinal[9]; genomaFinal[10] 0 genomaFinal[11] genomaFinal[12]]
A = convert(Array{Float64,2}, A)
sjahdfg, mtz, arreDifs = grafo_mexico(A, epocas, data, mu, beta, eta, precip, 0)


x = 1:(epocas-1); y = mtz[1:(epocas-1)] # These are the plotting data
plot(x, y, legend = false)
yy = data[1:(epocas-1)]
p1=plot!(1:(epocas-1),yy, title = "Chiapas", legend = false)

#Col2
x = 1:(epocas-1); y = mtz[(epocas+1):(2*epocas-1)] # These are the plotting data
plot(x, y, legend = false)
yy = data[floor(Int, length(data)/4 +1):floor(Int, length(data)/4 +1)+(epocas-2)]
p2=plot!(x, yy,title = "Guerrero", legend = false)

#Col3
x = 1:(epocas-1); y = mtz[(2*epocas+1):(3*epocas-1)] # These are the plotting data
plot(x, y, legend = false)
yy = data[2*floor(Int, length(data)/4)+1:2*floor(Int, length(data)/4 +1)+(epocas-3)]
p3=plot!(x, yy,title = "Oaxaca", legend = false)

#Col4
x = 1:(epocas-1); y = mtz[(3*epocas+1):(4*epocas-1)] # These are the plotting data
plot(x, y, legend = false)
yy = data[3*floor(Int, length(data)/4)+1:3*floor(Int, length(data)/4 +1)+(epocas-4)]
p4=plot!(x, yy,title = "Veracruz", legend = false)

plot(p1, p2, p3, p4, layout = (2, 2), legend = false)
