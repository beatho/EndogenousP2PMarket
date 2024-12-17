#pragma once

#include "StudyCase.h"
#include "Simparam.h"




/* interface*/
#include "MethodP2P.cuh"
#include "Method.h"
#include "MethodOPF.h"

/* Market */
#include "PAC.h"


/* Market endogene DC*/
#include "ADMMConst.h"
#include "ADMMConst1.h"
#include "OSQPConst.h"
#include "PACConst.h"

#include "ADMMGPUConst1.cuh"
#include "ADMMGPUConst1T.cuh"
#include "ADMMGPUConst2.cuh"
#include "ADMMGPUConst3.cuh"
#include "ADMMGPUConst4.cuh"
#include "ADMMGPUConst5.cuh"

#include "ADMMGPUConstCons.cuh"
#include "ADMMGPUConstCons2.cuh"
#include "ADMMGPUConstCons3.cuh"


/* OPF */
#include "OPFADMM.h"
#include "OPFADMMGPU.cuh"
#include "OPFADMM2.h"
#include "OPFADMMGPU2.cuh"
#include "OPFPDIPM.h"



/* Other*/
#include "ADMMACConst1.h"
#include "DCOPFOSQP.h"



#include <iostream>
#include <string>

#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}

class System
{
public:

	System();
	System(float rho, int iterMaxGlobal, int iterMaxLocal, float epsGlobal, float epsLocal, std::string nameMethode, int nAgent, float P = 0, float dP = 0, float a = 0, float da = 0, float b = 0, float db = 0);
	~System();
	const std::string sADMMConst = "ADMMConst";
	const std::string sADMMConst1 = "ADMMConst1";
	const std::string sADMMGPUConst1 = "ADMMGPUConst1";
	const std::string sADMMGPUConst1T = "ADMMGPUConst1T";
	const std::string sADMMGPUConst2 = "ADMMGPUConst2";
	const std::string sADMMGPUConst3 = "ADMMGPUConst3";
	const std::string sADMMGPUConst4 = "ADMMGPUConst4";
	const std::string sADMMGPUConst5 = "ADMMGPUConst5";
	const std::string sADMMGPUConstCons = "ADMMGPUConstCons";
	const std::string sADMMGPUConstCons2 = "ADMMGPUConstCons2";
	const std::string sADMMGPUConstCons3 = "ADMMGPUConstCons3";
	const std::string sDCOPFOSQP = "DCOPFOSQP";
	const std::string sADMMACConst1 = "ADMMACConst1";
	const std::string sPAC = "PAC";
	const std::string sPACConst = "PACConst";
	const std::string sOSQPConst = "OSQPConst";
	const std::string sOPFADMM = "OPFADMM";
	const std::string sOPFADMMGPU = "OPFADMMGPU";
	const std::string sOPFPDIPM = "OPFPDIPM";
	Simparam solve();
	void solveIntervalle(std::string path, MatrixCPU* interval, int nCons, int nGen); 
	void solveIntervalle(std::string path, std::string name, MatrixCPU* interval);
	void solveIntervalle(std::string path, int begin, int end, int chosenAgentGen);
	void UpdateP0();
	void resetMethod(); // permet de forcer l'initialisation, même si ce n'est pas la première itération
	void resetParam();
	void removeLink(int i, int j);
	void addLink(int i, int j);
	Agent removeAgent(int agent);
	void restoreAgent(Agent& agent, bool all = false);
	void setBestRho(float rhoMax = 0, bool rhoVar = 0, float rhoTest = 0);


	// mettre tous les set permettant de modifier les paramètres...
	void setStudyCase(const StudyCase& cas);
	void setSimparam(const Simparam& param);
	void setMethod(std::string nameMethode);
	void setMethod(Method* method);
	void setRho(float rho);
	void setRho1(float rho1);
	void setRhoL(float rho);
	void setIter(int iterG, int iterL);
	void setStep(int stepG, int stepL);
	void setEpsG(float epsG);
	void setEpsGC(float epsgC);
	void setEpsL(float epsL);
	void setTrade(MatrixCPU* trade);
	void setLineLimitMin(float lineMin);
	void setWarmStart(bool warmstart = true);
	void setConstraintRelaxation(float factor = 1);

	MatrixCPU getRes() const;
	MatrixCPU getTrade() const;
	MatrixCPU getTemps() const;
	MatrixCPU getIter() const;
	MatrixCPU getConv() const;
	MatrixCPU getFc() const;
	MatrixCPU getResR() const;
	MatrixCPU getResS() const;
	MatrixCPU getResX() const;
	MatrixCPU getPn() const;
	int getNbSimu(MatrixCPU* interval) const;
	int getNagent() const;
	
	int getNTrade() const;

	std::string generateDate(int year, int month, int day, int hour);
	void generateP0(MatrixCPU* P0, std::string path, std::string month);
	std::string generateMonth(int year, int month);
	void display(int type=0);
	void displayTime(std::string fileName = "SimulationFB.csv") const;
	void displayTradesAgent();

private:
	int dayMonth[12] = { 31, 28, 31, 30, 31, 30 , 31, 31, 30, 31, 30, 31 };
	int m[12] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 };
	bool useOPF = false;
	StudyCase _case;
	Simparam _simparam;
	Simparam* _result = nullptr;
	Method* _methode = nullptr;
	MatrixCPU _temps;
	MatrixCPU _iter;
	MatrixCPU _conv;
	MatrixCPU _fc;
	MatrixCPU _ResR;
	MatrixCPU _ResS;
	MatrixCPU _ResX;
	
	int getNFileline(std::string nameFile);


	
};




