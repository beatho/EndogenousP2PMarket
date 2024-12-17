import os
with open("PowerMaxCountry.csv", "w") as out:
    for name in os.listdir("./"):
        if ("inter.txt" in name):
            with open(name) as file:
                for line in file:
                    print(name.split("-inter.txt")[0], line.split("]")[0].split(",")[-1].split(" ")[-1])
                    print(line.split("]")[0].split(",")[-1].split(" ")[-1], file=out)
                    #print(name.split("-inter.txt")[0], line.split("]")[0].split(",")[-1].split(" ")[-1], file=out)
                    break

with open("CoefPoly.csv","w") as out:
    for name in os.listdir("./"):
        if ("-poly.txt" in name):
            with open(name) as file:
                for line in file:
                   print(name.split("-poly.txt")[0], line.split(",")[0].split("[")[-1],line.split(",")[1])
                   print(line.split(",")[0].split("[")[-1],line.split(",")[1], file=out)
                   #print(name.split("-poly.txt")[0], line.split(",")[0].split("[")[-1],line.split(",")[1], file=out)
                   break 