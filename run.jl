include("modelado.jl")
include("hitAndRun.jl")
using DelimitedFiles
using Plots

NODOS = 15          #Cantidad de nodos
mu = 0.08           #Proba de recuperacion
beta = 0.23         #Proba de contacto con vecino
r = 0.7             #Proba de infección por concatco con vecino
eta = 0.0           #Proba de no contagio por interaccion con vecinos
w=0.02              #Proba de pasar de recuperado a suceptible (solo sir)
epocas = 53        #Iteraciones
seed = 2
#barabasi Opcionales
n0 = 4              #Nodos del grafo incial.
k_bar = 3           #Vecinos de cada nodo.
proba = 0.1         #Proba de empezar infectado.

#regular Opcionales
k_reg = 4           #Grado de los vertices.

#small-world
k_sw = 4            #Grado de los vertices.
beta_sw = 0.4       #Probabilidad de aristas.

beta = [0.231857 0.427434 0.010000 0.558327 0.221011 0.294170;
        0.145923 0.396093 0.010000 0.397411 0.010000 0.334287;
        0.221471 0.310545 0.420411 0.570422 0.236681 0.190461;
        0.576013 0.596109 0.568935 0.010000 0.169991 0.356684]

data = readdlm("./data.txt")

#Precipitacion primer mes
precip = [0.24940130 0.13479302 0.15121451 0.18918919 0.64762231 0.80465275 0.76565173 0.55319877 1.00000000 0.83099555 0.27129661 0.14471433 0.03838534;
          0.15494466 0.00071403 0.01285255 0.00963941 0.30024991 0.86826134 0.82506248 1.00000000 0.92466976 0.53052481 0.02392003 0.00464120 0.03721910;
          0.12779813 0.08424908 0.09686610 0.15873016 0.51037851 1.00000000 0.79894180 0.78591779 0.91982092 0.64794465 0.12779813 0.08506309 0.06378927;
          0.23616089 0.13898845 0.28793309 0.29669454 0.59418558 1.00000000 0.61967344 0.78733572 0.88410992 0.79291119 0.28474711 0.22182397 0.14065282]

epocas = 300
mu = 0.1428
eta=0.0

#Precipitacion tiempo completo
precip = [0.24940130 0.13479302 0.15121451 0.18918919 0.64762231 0.80465275 0.76565173 0.55319877 1.00000000 0.83099555 0.27129661 0.14471433 0.03838534 0.02575532 0.09039128 0.10252600 0.36428925 1.00000000 0.74244676 0.71471025 0.84720158 0.79692917 0.13150074 0.11614661 0.24354244 0.11240420 0.11268805 0.12063582 0.66193585 1.00000000 0.82685211 0.62248084 0.80982118 0.80158955 0.26142492 0.23417542 0.30572082 0.10389016 0.08604119 0.09816934 0.19931350 0.43112128 0.45148741 0.84988558       0.67597254 1.00000000 0.12151030 0.05469108 0.18650712 0.12585104 0.06663916 0.18465030 0.36723747 0.76191459 0.67113679 0.58345368 1.00000000 0.44914380 0.05178461 0.08231896 0.15584866 0.04547032 0.11003124 0.11801458 0.42519958 1.00000000 0.54251996 0.82575495 0.88753905 0.63658452 0.33772995 0.14960083;
                    0.15494466 0.00071403 0.01285255 0.00963941 0.30024991 0.86826134 0.82506248 1.00000000 0.92466976 0.53052481 0.02392003 0.00464120 0.03721910 0.00245787 0.01123596 0.00316011 0.05582865 0.56776685 0.99473315 1.00000000 0.83110955 0.34515449 0.01404494 0.00456461 0.00417537 0.00104384 0.00487126 0.02157272 0.23347251 0.86430063 0.88413361 0.88239388 1.00000000 0.64752958 0.26826722 0.00347947 0.01248581 0.01986379 0.00681044 0.08030647 0.13110102 0.47814983 0.63308740 1.00000000 0.54710556 0.34733258 0.02099886 0.00539160 0.00000000 0.00000000 0.00067249 0.04572966 0.06086079 0.63584398 0.89509079 0.74445192 1.00000000 0.42030935 0.00437122 0.00000000 0.00891819 0.01008143 0.00193874 0.00853044 0.25746413 0.76812718 0.57696782 0.69329197 1.00000000 0.46956185 0.02946879 0.00542846;
                    0.12779813 0.08424908 0.09686610 0.15873016 0.51037851 1.00000000 0.79894180 0.78591779 0.91982092 0.64794465 0.12779813 0.08506309 0.06378927 0.04760394 0.03713107 0.06696287 0.30371311 0.66772453 0.88638527 1.00000000 0.73087909 0.64677880 0.17200889 0.08251349 0.14789272 0.06104725 0.02886335 0.04699872 0.27075351 0.66130268 1.00000000 0.74201788 0.74201788 0.44163474 0.31545338 0.13997446 0.07111001 0.10555005 0.03666997 0.10579782 0.12413280 0.37388503 0.73909812 1.00000000 0.65535183 0.50346878 0.13206145 0.05921705 0.13169488 0.03895201 0.02990958 0.10201716 0.21956875 0.83167169 0.83329469 0.49200093 1.00000000 0.37305820 0.08277301 0.04544401 0.06737120 0.06032585 0.01012770 0.07617790 0.30162924 0.93306913 0.69044474 0.73139586 1.00000000 0.57815940 0.08410392 0.10656099;
                    0.23616089 0.13898845 0.28793309 0.29669454 0.59418558 1.00000000 0.61967344 0.78733572 0.88410992 0.79291119 0.28474711 0.22182397 0.14065282 0.13115727 0.07685460 0.03204748 0.21810089 0.87596439 0.75934718 1.00000000 0.72878338 0.95133531 0.36587537 0.16973294 0.26446551 0.12880026 0.09316770 0.10264792 0.30369402 0.68649886 1.00000000 0.70415168 0.80385747 0.70153645 0.54593004 0.40503432 0.17140083 0.21485798 0.05462491 0.11556203 0.11337703 0.26074290 0.35324108 1.00000000 0.79825200 0.65986890 0.26074290 0.13522700 0.14970666 0.07809023 0.05098119 0.13959134 0.18409873 0.64313170 0.61723650 0.32551082 1.00000000 0.49949423 0.13675905 0.06595185 0.13127188 0.13681447 0.03208868 0.09655776 0.19778296 0.40227538 0.43232205 0.54929988 1.00000000 0.58838973 0.28646441 0.20157526]



#x = 1:54; y = mtz[1:54] # These are the plotting data
#plot(x, y, legend = false)
#yy = data[1:51]
#plot(1:51,yy, legend = false)

#writedlm("200k", M)


puntos = 500000
arre = ones(puntos)
M = hit_and_run(puntos)
writedlm("./puntos/500k", M)
#M = readdlm("./puntos/5k")
for q in 1:puntos
    A = [M[q,1] 0 M[q,2] M[q,3]; 0 M[q,4] M[q,5] 0; M[q,6] M[q,7] M[q,8] M[q,9]; M[q,10] 0 M[q,11] M[q,12]]
    A = convert(Array{Float64,2}, A)
    arre[q], kk, kk2  = grafo_mexico(A, epocas, data, mu, beta, eta, precip,0)
end

findmin(arre)
min = findmin(arre)[2]
min = 40562
print(M[min,:])
#writedlm("./minimos/min1Mmax", M[min,:])


nuevos = M[min,:]

A = [nuevos[1] 0 nuevos[2] nuevos[3]; 0 nuevos[4] nuevos[5] 0; nuevos[6] nuevos[7] nuevos[8] nuevos[9]; nuevos[10] 0 nuevos[11] nuevos[12]]
A = convert(Array{Float64,2}, A)
#A = [0.31 0.00 0.13 0.93; 0.0 0.99 0.99 0.0; 0.60 0.99 0.20 0.90; 0.38 0.0 0.2 0.92]
diff, mtz, arreDifs = grafo_mexico(A, epocas, data, mu, beta, eta, precip,0)
diff, errInd, arreDifs = grafo_mexico(A, epocas, data, mu, beta, eta, precip, 1)


## Seccion para error por estado dado un punto minimo
x = 1:(epocas-1); y = errInd[1:(epocas-1)] # These are the plotting data
p1=plot(x, y, title="Chiapas", legend = false)

#Col2
x = 1:(epocas-1); y = errInd[(epocas+1):(2*epocas-1)] # These are the plotting data
p2=plot(x, y,title = "Guerrero", legend = false)

#Col3
x = 1:(epocas-1); y = errInd[(2*epocas+1):(3*epocas-1)] # These are the plotting data
p3=plot(x, y,title = "Oaxaca", legend = false)

#Col4
x = 1:(epocas-1); y = errInd[(3*epocas+1):(4*epocas-1)] # These are the plotting data
p4=plot(x, y, title = "Veracruz", legend = false)

plot(p1, p2, p3, p4, layout = (2, 2), legend = false)

## Seccion para graficar el comportamiento del error total dado un punto minimo.
x = 1:(epocas);
plot(x, arreDifs, legend = false)



## Seccion para graficar los 4 estados en un grid
#Grafica col 1
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


## Esto es para correr la primer parte del proyecto, ver como se comportan los
#  modelos en distintas topologias.


SIS_smallWorld(NODOS, mu, beta, r, eta, w, epocas, proba, k_sw, beta_sw, seed)

SIR_smallWorld(NODOS, mu, beta, r, eta, w, epocas, proba, k_sw, beta_sw, seed)

SIS_Barabasi(NODOS, mu, beta, r, eta, w, epocas, proba, n0, k_bar, seed)

SIR_Barabasi(NODOS, mu, beta, r, eta, w, epocas, proba, n0, k_bar, seed)

SIR_regular(NODOS, mu, beta, r, eta, w, epocas, proba, k_reg, seed)

SIS_regular(NODOS, mu, beta, r, eta, w, epocas, proba, k_reg, seed)


#=
#Pruebas, muestra no determinisitico small world, con 5 nodos da true pero
#con 20 nodos da false.

grafo1 = watts_strogatz(5, 4, 0.3, seed=2)
grafo2 = watts_strogatz(5, 4, 0.3, seed=2)
print(grafo1==grafo2)
=#
