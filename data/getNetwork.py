import os
# fromNode;toNode;X;Y;numLines;limit;length;
indice =0
nb0 = 0
nbLine = 0
print("begining")
with open("Network.txt", "w") as out:
    with open("network_edges.csv") as file:
       for line in file:
            if indice>0:
               fromNode = line.split(",")[0]
               toNode = line.split(",")[1]
               Y = (float(line.split(",")[3]) * float(line.split(",")[4]))
               limit = line.split(",")[5]
               nbLine = nbLine+1
               if(float(limit)==0):
                  nb0 = nb0+1
               else:
                  print(fromNode, toNode, Y, limit)
                  #print(fromNode, toNode, Y, limit, file=out)
            indice = indice+1   

print(nbLine, nb0)