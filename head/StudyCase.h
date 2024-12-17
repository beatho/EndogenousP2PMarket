#pragma once


#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include "MatrixCPU.h"
#include "MatrixGPU.cuh"
#include "MatrixCPUD.h"
#include <math.h>
#include <iostream>
#include <vector>
#include "StudyCaseAgent.h"
#include "StudyCaseDCGrid.h"
#include "StudyCaseACGrid.h"

// rajouter securité pour que les Sbases soit les mêmes, et que les réseaux DC utilise un marché DC (et inversement)
// refaire les simulations sans le alreadyDefined pour refaire les matrices
// generation europe il faut remettre le caclul de power Sensi dans le linkage 
// genBetaDistanceZone


#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}



class StudyCase
{
    
    StudyCaseAgent SCAg; // Agents
    StudyCaseDCGrid* SCGrid = nullptr; // AC or DC-Grid
    //StudyCaseACGrid SCACG; // AC-Grid   
    
    float _Sbase = 1;
    bool DC = true;
    bool toReduce = false;
    int _nAgent = 0;
    int _nBus = 1;
    int _nLine = 0;
   
    std::string _name = "None";
    
    // melange grid/agent
    float _timeInit;
    MatrixCPU _SensiPower; // G = A*I
    MatrixCPU _CoresBusAgent; // I
    MatrixCPU _SensiPowerReduce; // Gred
    MatrixCPU _Distance; // Dij = sum |Pl|

    // Permettre de calculer W0
    MatrixCPU _CoresBusAgentLin;
    MatrixCPU _CoresAgentBusLin;
    MatrixCPU _CoresAgentBusLinBegin;
    MatrixCPU _nAgentByBus;

    

    float rand1() const; 
    int randab(int a, int b) const; // return a random number between a and b
    void genCoresBusAgent(bool all = false);
    void initMat();
    void createGrid(bool _DC);

public:

  
    ~StudyCase();
    StudyCase();
    StudyCase(int nAgent, float P, float dP, float a, float da, float b, float db, float propCons = 0.4, float propPro = 0.125);
    StudyCase(int nAgent, float P, float dP, float Q, float dQ, float a, float da, float aQ, float daQ, float b, float db, float Gamma, float dGamma, float propCons, float propGenNFle, float propPro);
    StudyCase(int nAgent, float P, float dP, float a, float da, float b, float db, float Gamma, float dGamma, float propCons, float propGenNFle, float propPro);
    StudyCase(int nAgent, float P0, float dP, float b, float db, float propCons = 0.4);


    StudyCase(const StudyCase& s);
    StudyCase(std::string fileName);
    StudyCase& operator= (const StudyCase& s);

    // generator
    void UpdateP0(MatrixCPU* P0);
        
    void genBetaUniforme(float beta);
    void genBetaDistance(float s); //agent + grid
    void genBetaDistanceByZone(MatrixCPU* s);
    void genGrid(int _nBus, int _nMajorLine, int _minorLine, float ReacMajor, float DeltaReacMajor, float ReacMinor, float DeltaReacMinor, float LlimitMajor, float dLlimitMajor, float LlimitMinor, float dLlimitMinor);
    void genGridBT(int _nBus, int Nbranch, int Ndeep, float length, float dlength);
    void genGridBTSpecial(int nBus, int Nbranch, int Ndeep, float length, float dlength, RadialType type);
    void genGridHTB(int nBus, int nLine, int dnLine, float length, float dlength);

    void genGridFromFile(std::string path, bool alreadyDefine=true);

    void genAgents(int nAgent, float propCons, float Pconso, float dPconso,  float bProd, float dbProd, float Pprod, float dPpord, float gamma, float dGamma); // gamma = -1 distance ?
    void genAgentsAC(int nAgent, float propCons, float propGenNFle, float Pconso, float dPconso, float bProd, float dbProd, float dQconso, float Pprod, float dPprod, float Gamma, float dGamma);
    void genAgentsFullRandom(int nAgent, float aMin, float aMax, float P0Min, float P0Max, float gammaMin, float gammaMax, float propConsoMin, float propConsoMax, float borneMin, float borneMax);
    
    void genLinkGridAgent();
    void computeSensiPower();
    void genLineLimit(int nLine, float limit, float dlimit );
    void genDCGridFromAC();

    //Set case
    void Set29node();
    void Set39Bus(std::string path = "data/", bool alreadyDefine=false);
    void Set3Bus(std::string path = "data/");
    void Set4Agent();
    void Set2node();
    void Set4nodeBis(std::string path);
    void Set2nodeConstraint(float lim = 0.8);
    void SetEuropeP0(const std::string& path, MatrixCPU* P0, bool alreadyDefine=0);
    void SetEuropeP0WithoutConstraint(const std::string& path, MatrixCPU* P0);
    void SetStudyCase(std::string path, std::string name, MatrixCPU* P0, bool alreadyDefine=0);
   
    void SetAC39Bus(std::string path = "data/", bool alreadyDefine = false);
    void SetAC3Bus(std::string path = "data/");
    void SetAC2node();
    void SetACFromFile(std::string name, std::string path = "data/ACGrid/");
    void SetACFromFileSimplify(std::string name, std::string path = "data/ACGrid/");
    void SetEuropeTestFeeder(std::string path = "data/ACGrid/", int typeOfAgentGen = 0, int beggining = 0);
   
    // set Parametters
    void setDistance(bool alreadyDefine = false, std::string path = "data/distance.txt");
    void setLineLimitMin(float lineMin);
    void setLineLimitRelaxation(float eps);
    void setLineLimit(int line, float limit);
    void setReduce(bool toReduce);
    
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

    MatrixCPU getPowerSensi() const;
    MatrixCPU getLineLimit() const;
    MatrixCPU getCurrentLimit() const;

    MatrixCPU getCoresLineBus(bool force=false) const;
    MatrixCPU getCoresBusAgent() const;
    MatrixCPU getCoresBusAgentLin() const;
    MatrixCPU getCoresAgentBusLin() const;
    MatrixCPU getCoresAgentBusLinBegin() const;
    MatrixCPU getNagentByBus() const;
    MatrixCPU getLineSuceptance() const;
    MatrixCPU getLineReactance() const;
    MatrixCPUD getLineSuceptanceD() const;
    MatrixCPUD getLineReactanceD() const;
    MatrixCPU getUpperBound() const;
    MatrixCPU getLowerBound() const;
    MatrixCPU getPobj();
    MatrixCPUD getPobjD();
    MatrixCPUD getSolPF() const;
    double getV0() const;
    double gettheta0() const;
    MatrixCPU getZsRe() const;
    MatrixCPU getZsImag() const;
    MatrixCPU getYd() const;
    MatrixCPU getGlin() const;
    MatrixCPU getBlin() const;
    MatrixCPU getGlin2() const;
    MatrixCPU getBlin2() const;
    MatrixCPU getVoltageInit() const;
    MatrixCPUD getVoltageInitD() const;
    MatrixCPUD getGlinD() const;
    MatrixCPUD getBlinD() const;
    MatrixCPU getNLines() const;
    MatrixCPU getNLinesBegin() const;
    MatrixCPU getCoresBusLin() const;
    MatrixCPU getCoresVoiLin() const;

    
    float getTimeInit() const;
    int getNagent() const;
    int getNCons() const;
    bool isAC() const;
    bool isRadial() const;
    bool isCurrentLimit() const;
    int getNLine(bool force=false) const;
    int getNBus() const;
    float getSbase() const;
    int getLastBus() const;
    MatrixCPU getVoisin(int agent) const;
    Agent getAgent(int agent) const;
    std::string getName() const;



    // modification
    void removeLink(int i, int j);
    Agent removeAgent(int agent);
    void restoreAgent(Agent& agent, bool all = false);
    void addLink(int i, int j);
    void nextStepPobj();
   
    
    void saveCSV(const std::string& fileName, bool all = true);
    void display(int type=0);
    void displayLineCores(MatrixCPU* g, bool all = true);


};
