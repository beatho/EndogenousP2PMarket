import os
debut = 0
with open("costGen.txt", "w") as out :
    with open("generator_info.csv") as file:
        for line in file :
            gene = line.split(";")
            id = gene[0]
            lincost = gene[10] 
            if(debut<2):
                debut = debut+1
            else:
                print(id, lincost,file=out)
                   