import os
# On va d'abord trouver les noeuds qui sont dans un pays
indice =0
Pays = [] # pays de chaque bus
Pay = "ITA"
Pay2 = "Italy"
Noeud = [] # liste des noueds/bus dans le bon pays
Agents =[]
PosAgent=[]

#ID,name,country,voltage,latitude,longitude
with open("network_nodes.csv") as file:
    for line in file :
        if indice>0:
            bus = line.split(",")
            Pays.append(bus[2]);
        indice=indice+1


# on fait la correspondance entre les noeuds et les agents
with open("BusAgent"+ Pay2 +".txt", "w") as out :
    with open("load_signal.csv") as file:
        iter = 0
        for line in file : # seule la première ligne est interressante
            nums = line.split(",")[1:] # liste de tous les bus
            agent = 0
            for num in nums:
                if(Pays[iter]==Pay):
                    print(agent, int(num), Pays[iter], file=out)
                    Agents.append(agent);
                    agent = agent + 1
                    Noeud.append(int(num))
                    PosAgent.append(iter)
                iter = iter + 1
            break

# on regarde les liens qui sont exclusivement en Pay
# fromNode;toNode;X;Y;numLines;limit;length;
for node in Noeud:
    print(node)
print("fin noeud")

indice =0
nb0 = 0
nbLine = 0
with open("Network" + Pay2 +".txt", "w") as out:
    with open("network_edges.csv") as file:
       for line in file:
            if indice>0:
                fromNode = line.split(",")[0]
                toNode = line.split(",")[1]
                Y = (float(line.split(",")[3]) * float(line.split(",")[4]))
                limit = line.split(",")[5]
                if(int(fromNode) in Noeud):
                    if(int(toNode) in Noeud):
                        nbLine = nbLine+1 
                        if(float(limit)==0):
                            nb0 = nb0+1
                        #print(fromNode, toNode, Y, limit)
                        print(fromNode, toNode, Y, limit, file=out)
            indice = indice+1   
print(nbLine, nb0)
# On regarde les générateurs Français dont le noeud le plus proche est aussi en France
# ID;name;country;origin;latitude;longitude;status;primaryfuel;secondaryfuel;capacity;lincost;cyclecost;minuptime;mindowntime;minonlinecapacity
indice =0
with open("genCarac"+ Pay2 +".txt", "w") as out :
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
                country = gene[2+offset]
                Bus = gene[3+offset]
                lincost = gene[10+offset] 
                capacity = gene[9+offset]
                if(int(Bus) in Noeud):
                    if(country==Pay or country==Pay2):
                        print(ID, lincost, capacity, Bus, file=out)
                    else:
                        print(Bus, country) # generateur pas dans le pays mais le noeud le plus proche est dans le pays, possible ? oui

            indice = indice+1

# il faut générer l'ensemble des charges en fonction du temps

indice = 0
loadTotal =[Noeud];
dayMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
month = 0
nomFile =""
print(loadTotal)
with open("load_signal.csv") as file:
    for line in file:
        if indice>0:
            load=[]
            charge = line.split(",")
             
            date = charge[0].split(" ")
            heure = date[1]
            mois = date[0].split("-")[0:2]
            NewnomFile = "load"+ Pay2 +"/Month/"+ mois[0] + "-" + mois[1] + ".txt"
            if(nomFile !=NewnomFile):
                if(nomFile !=""):
                    with open(nomFile, "w") as out:
                        n = len(loadTotal[0])
                        m = len(loadTotal)
                        for i in range(0,n):
                            for j in range(0,m):
                                print(float(loadTotal[j][i]), end=' ', file=out)
                            print("", file=out)
                        loadTotal =[Noeud];
                nomFile = NewnomFile
            for node in PosAgent:
                load.append(charge[node+1])#+1 car le premier c'est la date
            loadTotal.append(load);
        else:
            AllNode = line.split(",")[1:]
            indice=1
    with open(NewnomFile, "w") as out:
        n = len(loadTotal[0])
        m = len(loadTotal)
        #for i in range(0,n):
            #for j in range(0,m):
                #print(loadTotal[j][i], end=' ', file=out)
            #print("", file=out)
        #loadTotal =[Noeud];
#print(loadTotal)