import os
# ID;name;country;origin;latitude;longitude;status;primaryfuel;secondaryfuel;capacity;lincost;cyclecost;minuptime;mindowntime;minonlinecapacity
indice =0
with open("genCarac.txt", "w") as out :
    with open("generator_info.csv") as file:
        for line in file :
            if indice>0:
                gene = line.split(",") # il arrive qu'une virgule apparaisse dans le nom, ce qui décalle le split !
                if(len(gene)==15):
                    offset=0
                elif (len(gene)==16):
                    offset=1
                else:
                    print("probleme de split ligne ", indice)
                ID = gene[0]
                Bus = gene[3+offset]
                lincost = gene[10+offset] 
                capacity = gene[9+offset]
                print(ID, lincost, capacity, Bus, file=out)
            indice = indice+1