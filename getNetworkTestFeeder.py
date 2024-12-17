from math import pi
import os
## Input :
# Lines.csv : name; fromNode;toNode;Phase, length, unit, line Code
# LineCodes.csv : Code, nphase, R1, X1, R0, X0, C1, C0, Units
# LinesMaxCurrent.csv : code, Ima
# BusCoords.csv : bus, x, y -> just to have the number of bus
# Loads.csv: name, numphase (1), bus, phase (A,B ou C), kV, Modele, Connection, kW, PF, NameShape
# LoadShapes.csv : Name; npts, minterval, File, useactual -> not used, useless because Shape_N -> Load_profile_N.csv

## Output :
# BranchTestFeeder.txt : from, to, Ys real, Ys Im, Yp, tau, thetha , limit, zs real, zs Imag
# CaseTestFeeder.txt : Sbase, Vbase, Zbase, nAgent, nBus, nLine, V0, theta0
# AgentTestFeeder.txt : bus, Pmult, facteur puissance
# AgentConsumptionTestFeeder.txt : Mat[i][j] = P_i à t=j


indice =0
nb0 = 0
nbLine = 0
dictLine = {};
LineCode =[];
Imax = [];
V0 = 230; # 230V
S0 = 1; # 1kW
Z0 = V0*V0/(S0*1000);
f = 50; # 50Hz
w0 = 2 * pi * f;

print("begining")

with open("European_LV_Test_Feeder_v2\European_LV_CSV\BusCoords.csv") as file:
   for line in file:
      indice = indice + 1;
nBus = indice - 2;
indice = 0;
with open("European_LV_Test_Feeder_v2\European_LV_CSV\LinesMaxCurrent.csv") as file:
   for line in file:
      if indice>1:
         lineCode = line.split(",")[0];
         try:
            I = float(line.split(",")[1])
         except:
            print("------------")
            break;
         else:
            dictLine[lineCode] = [];
            dictLine[lineCode].append(I);
            LineCode.append(lineCode);
            Imax.append(I);
            #print(lineCode,I);      
      indice = indice+1;  
indice =0;
NtypeLine = len(LineCode);
with open("European_LV_Test_Feeder_v2\European_LV_CSV\LineCodes.csv") as file:
   for line in file: # Code, nphase, R1, X1, R0, X0, C1, C0, Units   
      if indice>1:
         lineCode = line.split(",")[0];
         try:
            R = float(line.split(",")[2])
            X = float(line.split(",")[3])
            C = float(line.split(",")[6])
            Units = line.split(",")[8]
            #print(Units, " km")
            if Units == "km" or Units == "km\n" :
               #print("changement d'echelle")
               R = R/1000;
               X = X/1000;
               C = C/1000;
         except:
            if indice<NtypeLine:
               print("********ERROR*********");
            break;
         else:
            dictLine[lineCode].append(R);
            dictLine[lineCode].append(X);
            dictLine[lineCode].append(C);
            #print(lineCode,R ,X, C);      
      indice = indice+1;  

#print(NtypeLine);
#print(dictLine);
indice = 0;
offset = 1;
with open("BranchTestFeeder.txt", "w") as out: #from, to, Ys real, Ys Im, Yp, tau, theta , limit, zs real, zs Imag
   with open("European_LV_Test_Feeder_v2\European_LV_CSV\Lines.csv") as file:
      for line in file: # name; fromNode;toNode;Phase, length, unit, line Code 
         if indice>1:
            nbLine = nbLine +1;
            BusFrom = int(line.split(",")[1]) - offset;
            BusTo = int(line.split(",")[2]) - offset;
            if BusTo != (indice - 1):
               print("ERROR : wrong order for the line \n")
            tau = 0;
            theta = 0;
            length = float(line.split(",")[4])
            unit = line.split(",")[5]
            if(unit != "m"):
               print("ERROR : wrong unit for the line \n")
            lineCode = line.split(",")[6][0:-1];
            limit = 230 * dictLine[lineCode][0];
            R = dictLine[lineCode][1] * length / Z0;
            X = dictLine[lineCode][2] * length / Z0;
            C = dictLine[lineCode][3] * length / Z0;
            Zs = complex(R,X);
            Ys = 1/Zs;
            Yp = 0;
            Yp = 2*C;

            YsRe = Ys.real;
            YsIm = Ys.imag;
            ZsRe = Zs.real;
            ZsIm = Zs.imag;
            print(BusFrom, BusTo, YsRe, YsIm, Yp, tau, theta, limit, ZsRe, ZsIm, file=out);
         indice = indice + 1;

indice = 0;
with open("AgentTestFeeder.txt", "w") as out:
   with open("European_LV_Test_Feeder_v2\European_LV_CSV\Loads.csv") as file:
      for line in file:
         if indice>2:
            bus = int(line.split(",")[2]) - offset;
            Pmult = line.split(",")[7];
            PF = line.split(",")[8];
            print(bus, Pmult, PF, file=out);
         indice = indice + 1;
nAgent = indice - 3;

nMaxShape = 100;
nMinute = 24 * 60;
PowerMat = [[0]*nMaxShape for i in range(nMinute)] # matrice de b lignes et m colonnes
for i in range(1,nMaxShape+1):
   name = "European_LV_Test_Feeder_v2\European_LV_CSV\Load Profiles\Load_profile_" + str(i) + ".csv";
   indice = 0;
   with open(name) as file:
      for line in file:
         if indice>0:
            PowerMat[indice-1][i-1] = float(line.split(",")[1]);
         indice = indice + 1;

with open("AgentConsumptionTestFeeder.txt", "w") as out:
    for i in range(0,nMinute):
      for j in range(0,nMaxShape):
         print(float(PowerMat[i][j]), end=' ', file=out)
      print("", file=out)
      # CaseTestFeeder.txt : Sbase, Vbase, Zbase, nAgent, nBus, nLine, V0, theta0

with open("CaseTestFeeder.txt", "w") as out:
   print(S0, V0, Z0, nAgent, nBus, nbLine, 1, 0, file=out);