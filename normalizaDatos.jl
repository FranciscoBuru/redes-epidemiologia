using DelimitedFiles

semanas  = 52
estados = 4
data = readdlm("./data.txt")

data2 = zeros(Float64,(semanas,estados))
inicio=105
fin=156
offset=104
for r in 1:estados
    for rr in 1:semanas
        data2[rr,r] = (data[rr+offset,r]-minimum(data[inicio:fin,r]))/(maximum(data[inicio:fin,r])-minimum(data[inicio:fin,r]))
    end
end
data2 = replace!(data2, Inf=>0.0)
data2 = replace!(data2, NaN=>0.0)
writedlm("./datosNorm3", data2)
