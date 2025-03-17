#pragma once

#include "Simparam.h"
#include "StudyCase.h"


/* Market */
#include "ADMMMarket.h"
#include "ADMMMarketOpenMP.h"


/* OPF */
#include "MethodOPF.h"
#include "OPFADMM.h"
#include "OPFADMM2.h"

/* PF */
#include "CPUPF.h"
#include "CPUPFdist.h"
#include "CPUPFdistPQ.h"
#include "CPUPFGS.h"


/* Other*/
#include "ADMMACConst1.h"
#ifdef OSQP
	#include "DCOPFOSQP.h"
#endif


#include <iostream>
#include <string>

#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}

class GPUPF; // inclu dans le .cuh !!!


class System
{
public:

	System();
	System(float rho, int iterMaxGlobal, int iterMaxLocal, float epsGlobal, float epsLocal, std::string nameMethode, int nAgent, float P = 0, float dP = 0, float a = 0, float da = 0, float b = 0, float db = 0);
	~System();
	
	const std::string sNR 	     = "NR";
	const std::string sGS        = "GS";
	const std::string sDistPQ    = "DistPQ";
	const std::string sNRGPU 	 = "NRGPU";
	const std::string sGSGPU     = "GSGPU";
	const std::string sDistPQGPU = "DistPQGPU";

	const std::string sADMMMarket 	 = "ADMM";
	const std::string sADMMMarketMP  = "ADMMMP";
	const std::string sADMMMarketGPU = "ADMMGPU";

	const std::string sADMMConst = "ADMMConst";
	const std::string sADMMConst1 = "DCEndoMarket"; //ADMMConst1
	const std::string sADMMGPUConst1 = "ADMMGPUConst1";
	const std::string sADMMGPUConst1T = "ADMMGPUConst1T";
	const std::string sADMMGPUConst2 = "ADMMGPUConst2";
	const std::string sADMMGPUConst3 = "ADMMGPUConst3";
	const std::string sADMMGPUConst4 = "DCEndoMarketGPU"; //ADMMGPUConst4
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
	Simparam solvePF(); // on pourra rajouter des paramètres
	ResultInterface* solve(ResultInterface* res, ParamInterface* param, StudyCaseInterface* caseInter, bool AC);
	ResultInterface* solvePF(ResultInterface* res, ParamInterface* param, StudyCaseInterface* caseInter);
	void solveIntervalle(std::string path, MatrixCPU* interval, int nCons, int nGen); 
	void solveIntervalle(std::string path, std::string name, MatrixCPU* interval);
	void solveIntervalle(std::string path, int begin, int end, int chosenAgentGen);
	void UpdateP0();
	void resetMethod(); // permet de forcer l'initialisation, meme si ce n'est pas la premiere iteration
	void resetParam();
	void removeLink(int i, int j);
	void addLink(int i, int j);
	Agent removeAgent(int agent);
	void restoreAgent(Agent& agent, bool all = false);
	void setBestRho(float rhoMax = 0, bool rhoVar = 0, float rhoTest = 0);


	// mettre tous les set permettant de modifier les param�tres...
	void setStudyCase(const StudyCase& cas);
	void setStudyCase(std::string caseName);
	void setSimparam(const Simparam& param);
	void setMethod(std::string nameMethode);
	void setMethodPF(std::string nameMethode, bool isDouble);
	void setMethod(Method* method);
	void setRho(float rho);
	void setRho1(float rho1);
	void setRhoL(float rho);
	void setIter(int iterG, int iterL);
	void setItIntern(int iter);
	void setStep(int stepG, int stepL);
	void setStep(int stepG, int stepL, int stepIntern);
	void setEpsG(float epsG);
	void setEpsGC(float epsgC);
	void setEpsIntern(float eps);
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
	bool usePFGPU = false;
	bool useDoublePF = false;
	StudyCase _case;
	Simparam _simparam;
	Simparam* _result = nullptr;
	Method* _methode = nullptr;
	CPUPF* _methodePF = nullptr;
	GPUPF* _methodePFGPU = nullptr;
	MatrixCPU _temps;
	MatrixCPU _iter;
	MatrixCPU _conv;
	MatrixCPU _fc;
	MatrixCPU _ResR;
	MatrixCPU _ResS;
	MatrixCPU _ResX;
	
	int getNFileline(std::string nameFile);


	
};




