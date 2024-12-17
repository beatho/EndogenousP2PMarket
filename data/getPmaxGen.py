import os
# ID;name;country;origin;latitude;longitude;status;primaryfuel;secondaryfuel;capacity;lincost;cyclecost;minuptime;mindowntime;minonlinecapacity
indice =0
print("begining")
with open("PowerMaxGen.txt", "w") as out:
    with open("generator_info.csv") as file:
       for line in file:
            if indice>0:
               print(indice)
               ID = line.split(";")[0]
               capacity = line.split(";")[9]
               print(ID, capacity, file=out)
            indice = indice+1