import os
indice =0
Pays =[]
with open("network_nodes.csv") as file:
    for line in file :
        if indice>0:
            bus = line.split(",")
            Pays.append(bus[2]);
        indice=indice+1

with open("BusAgent.txt", "w") as out :
    with open("load_signal.csv") as file:
        for line in file : # seule la première ligne est interressante
            nums = line.split(",")[1:]
            agent =0
            for num in nums:
                print(agent, int(num), Pays[agent], file=out)
                agent = agent+1
            break
                   