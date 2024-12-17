#pragma once


#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "Agent.h"
#include "MatrixCPU.h"
#include "MatrixGPU.cuh"
#include "MatrixCPUD.h"
#include <math.h>
#include <iostream>
#include <vector>

#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}


enum class Qlinks {FULLY, ASP, ONELINK};

class StudyCaseAgent
{
 
    int _nAgent = 0;
    int _nPro = 0;
    int _nGen = 0;
    int _nGenNFle = 0;
    int _nGenFle = 0;
    int _nCons = 0;
    int temporalStep = 0;
    bool _AC = false;

    Qlinks qlink = Qlinks::FULLY;

    std::string _name = "None";
    
    float _timeInit;
    Agent* _agents = nullptr;
    MatrixCPU _BETA;
    MatrixCPU _connect;
    MatrixCPU _a;
    MatrixCPU _b;
    MatrixCPU _Ub;
    MatrixCPU _Lb;
    MatrixCPU _Pmin;
    MatrixCPU _Pmax;
    MatrixCPU _nVoisin;
    MatrixCPU _Pobj; // -b/a for consumer,  a standard value for producer
    MatrixCPUD _PobjD; // -b/a for consumer,  a standard value for producer
    MatrixCPU _PobjTemp;
    MatrixCPU _factor; // P = factor * Pobj
    MatrixCPU _PF; // P = PF*S; Q = sqrt(1 - PF) S
    void setMatFromFile(const std::string& path, const std::string& date, MatrixCPU* Pgen, MatrixCPU* P0, MatrixCPU* costGen);
    void setGenFromFile(const std::string& path, MatrixCPU* Pgen, MatrixCPU* costGen, MatrixCPU* BusGen);
    void initCaseFromPobj();

    float rand1() const; 
    int randab(int a, int b) const; // return a random number between a and b
    float randabfloat(float a, float b) const;
    void initMat();
    void initMatAC();
    
    MatrixCPU GenBus; // when from files

public:
    int getNFileline(std::string nameFile);
    
    float _Sbase = 1;


    StudyCaseAgent();
    
    StudyCaseAgent(int nAgent, float P, float dP, float a, float da, float b, float db, float ppropCons = 0.4, float propPro = 0.125);
    StudyCaseAgent(int nAgent, float P, float dP, float Q, float dQ, float a, float da, float aQ, float daQ, float b, float db, float Gamma, float dGamma, float propCons, float propGenNFle, float propPro); // AC
    StudyCaseAgent(int nAgent, float P, float dP, float a, float da, float b, float db, float Gamma, float dGamma, float propCons, float propGenNFle, float propPro); // DC
    StudyCaseAgent(int nAgent, float P0, float dP, float b, float db, float propCons = 0.4);


    StudyCaseAgent(const StudyCaseAgent& s);
    StudyCaseAgent(std::string fileName);
    StudyCaseAgent& operator= (const StudyCaseAgent& s);

    // generator
    void UpdateP0(MatrixCPU* P0);
    void genConnec(MatrixCPU* connec);
    
    void genBetaUniforme(float beta);
    void genBetaDistance(float s, MatrixCPU* Distance);
    void genBetaZone(MatrixCPU* Mats, MatrixCPU* Distance, MatrixCPU zones);
    void genBetaRandom(float beta, float dbeta);

    void genAgents(int nAgent, float propCons, float Pconso, float dPconso,  float bProd, float dbProd, float Pprod, float dPpord, float gamma, float dGamma); // gamma = -1 distance ?
    void genAgentsAC(int nAgent, float propCons, float propGenNFle, float Pconso, float dPconso, float bProd, float dbProd, float dQconso, float Pprod, float dPpord, float gamma, float dGamma);
    void genAgentsFullRandom(int nAgent, float aMin, float aMax, float P0Min, float P0Max, float gammaMin, float gammaMax, float propConsoMin, float propConsoMax, float borneMin, float borneMax);
    
    void nextStepPobj();

    //Set case
    MatrixCPU Set29node(bool AC=false);
    MatrixCPU Set4Agent();


    MatrixCPU Set3BusOld(bool AC=false);
    MatrixCPU Set3Bus(bool AC = false);
    MatrixCPU Set2node(bool AC = false);
    MatrixCPU Set4nodeBis(bool AC = false);
    void SetEuropeP0(const std::string& path, MatrixCPU* P0);
    void SetStudyCaseAgent(std::string path, std::string name, MatrixCPU* P0);
 
    MatrixCPU SetACFromFile(std::string name, std::string path = "data/ACGrid/");
    MatrixCPU SetACFromFileSimplify(std::string name, std::string path = "data/ACGrid/");
    MatrixCPU SetEuropeTestFeeder(std::string path = "data/ACGrid/", int beggining = 0);

    // getter 
    MatrixCPU getBeta() const;
    MatrixCPU getC() const;
    MatrixCPU geta() const;
    MatrixCPU getb() const;
    MatrixCPU getUb() const;
    MatrixCPU getLb() const;
    MatrixCPU getPmin() const;
    MatrixCPU getPmax() const;
    MatrixCPU getNvoi() const;
    
    MatrixCPU getPobj();
    MatrixCPUD getPobjD();

    MatrixCPU getGenBus() const;
   
    float getTimeInit() const;
    int getNagent() const;
    int getNCons() const;
  
    MatrixCPU getVoisin(int agent) const;
    Agent getAgent(int agent) const;
    std::string getName() const;



    // modification
    void removeLink(int i, int j);
    Agent removeAgent(int agent);
    void restoreAgent(Agent& agent, bool all = false);
    void addLink(int i, int j);
    
    void saveCSV(const std::string& fileName, bool all = true);
    void display(int type=0);
  
    ~StudyCaseAgent();

};
