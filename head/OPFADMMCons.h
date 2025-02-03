#pragma once
#include "MethodOPF.h"
#include <iostream>
#include <string>
#include <chrono>
#include "Utilities.h"



class OPFADMMCons : public MethodOPF // OPFADMM2 mais on consid�re que l'on peut avoir un fonction cout sur les pertes, il y a donc un bus fictif
{
public:
	OPFADMMCons();
	OPFADMMCons(float rho);
	virtual ~OPFADMMCons();
	void setParam(float rho);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	virtual void solveConsensus(float eps, MatrixCPU* PSO);
	virtual void initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO);
	virtual void updateConsensus(MatrixCPU* Pmarket);
	
	
	
	bool chekcase();
	std::string NAME ="OPFADMMCons";
	void updateGlobalProb();
	
	virtual void updateX();
	void updateXWOCurrent();
	virtual float getPLoss();
	virtual float getQLoss();

	void updateMu();
	void updateChat();
	void CommunicationX();
	float updateRes(int indice);
	float updateResRhoFixe(int indice);
	virtual int feasiblePoint();

	void display();
private:
	float DFSP(int j); // recursive fonction to initialize
	float DFSQ(int j); // P and Q flows

	// ne change pas avec P0
	int _nBus = 0;
	int _nBusWLoss = 0;
	int _nLine = 0;
	int _nAgent = 0;
	float _rho = 0;
	float _rhoInv = 0;
	double coefPoly2[2];
	double root2[4];
	double root3[4];
	double root4[4];
	double root5[4];
	double root6[4];
	double coefPoly3[3];
	
	int _iterGlobal = 0;
	int _iterG = 0;
	int _stepG = 0;
	float _Ploss = 0;
	float _Qloss = 0;
	clock_t timeOPF = 0;
	bool isCurrentLimited = false;
	LossType losstype = LossType::CURRENT;

	// parameter agent and iteration dependant (but not here for now)


	MatrixCPU tempN1; // Matrix temporaire pour aider les calculs
	MatrixCPU tempNN; // plut�t que de re-allouer de la m�moire � chaque utilisation
	MatrixCPU* tempM1 = nullptr; //
	MatrixCPU* tempM = nullptr; //
	

	MatrixCPU Cost1;
	MatrixCPU Cost2;
	MatrixCPU _CoresBusAgent;
	MatrixCPU _CoresAgentBus;
	MatrixCPU _CoresAgentBusBegin;

	MatrixCPU PosAgent;
	MatrixCPU _nAgentByBus;

	MatrixCPU ZsRe;
	MatrixCPU ZsIm;
	MatrixCPU ZsNorm;
	MatrixCPU VoltageLimit; // (vmin^2, vmax^2) * sqrt(Nchild + 1 / 2)
	MatrixCPU VoltageLimitReal; // vmin, vmax
	MatrixCPU FluxLimit;
	
	MatrixCPU* X = nullptr; // (Pi,Qi,vi,li,pi,qi, pn, vai,pji,qji,lji)
	MatrixCPU* Ypre = nullptr;
	MatrixCPU* Y = nullptr; // (Pi,Qi,vi,li,pi,qi, pn, vai,pji,qji,lji)

	MatrixCPU* Mu = nullptr;

	MatrixCPU* Chat = nullptr;

	MatrixCPU* Hinv = nullptr;
	MatrixCPU* Q = nullptr;
	MatrixCPU* A = nullptr;

	MatrixCPU nChild;
	MatrixCPU* Childs = nullptr;
	MatrixCPU PosChild; 
	MatrixCPU Ancestor;
	MatrixCPU CoresLineBus;

	MatrixCPU sizeOPFADMMCons;

	MatrixCPU resF;

	// Local resolution
	MatrixCPU tempN2; // size : (_nAgent*2, 1)
	MatrixCPU tempB2; // size : (_nBus  *2, 1)
	MatrixCPU CoresSoloBusAgent;
	MatrixCPU Pn;
	MatrixCPU Pmin;
	MatrixCPU Pmax;
	MatrixCPU Pbmax;
	MatrixCPU Pbmin;
	MatrixCPU Pb;



	// connsensus
	MatrixCPU etaSO; // consensus with market
	float _rhoSO = 0;



};


	




