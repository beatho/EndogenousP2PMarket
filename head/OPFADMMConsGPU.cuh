#pragma once
#include "MethodOPFGPU.cuh"
#include <iostream>
#include <string>
#include <chrono>
#include "Utilities.cuh"
#include "Utilities.h"

#include "kernelFunction.cuh"


class OPFADMMConsGPU : public MethodOPFGPU // OPFADMM2 mais on consid�re que l'on peut avoir un fonction cout sur les pertes, il y a donc un bus fictif
{
public:
	OPFADMMConsGPU();
	OPFADMMConsGPU(float rho);
	virtual ~OPFADMMConsGPU();
	void setParam(float rho);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	
	virtual void initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO);
	virtual void updateConsensus(MatrixGPU* Pmarket);
	virtual void solveConsensus(float eps, MatrixGPU* PSO);
	virtual void updateConsensus(MatrixCPU* Pmarket);
	virtual void solveConsensus(float eps, MatrixCPU* PSO);
	
	bool chekcase();
	std::string NAME ="OPFADMMConsGPU";
	void updateGlobalProb();
	
	virtual void updateX();
	void updateXWOCurrent();
	void updateXWOCurrentOnCPU();
	void updateXWOCurrentOnCPUBis();
	void updateXWOCurrentOnCPUBis(bool first);

	void computeLoss();

	virtual float getPLoss();
	virtual float getQLoss();

	void updateMu();
	void updateChat();
	void CommunicationX();
	float updateRes(int indice);
	float updateResRhoFixe(int indice);
	virtual int feasiblePoint();

	void ComputePFromAgentToBus();
	virtual MatrixCPU getPb();
	virtual MatrixCPU getPhi();
	virtual MatrixCPU getE();


	void display();
private:
	// ne change pas avec P0
	int _numBlocksB;
	int _numBlocksN;
	int _numBlocksH;

	int _nBus = 0;
	int _nBusWLoss = 0;
	int _nLine = 0;
	int _nAgent = 0;
	int _nAgentOn0 = 0;

	int _sizeOPFADMMConsTotal = 0;
	int _debutloss = 0;
	int _sizeOPFADMMConsMax = 0;
	int _sizeChat = 0;

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


	MatrixGPU tempN1; // Matrix temporaire pour aider les calculs
	MatrixGPU tempNN; // plut�t que de re-allouer de la m�moire � chaque utilisation
	MatrixGPU tempM1; //
	MatrixGPU tempM; //
	
	// Reste sur CPU
	MatrixCPU _nAgentByBusCPU;
	MatrixCPU resF;
	MatrixCPU nChildCPU;
	MatrixCPU CoresLineBusCPU;
	MatrixCPU nVoisinCPU;
	MatrixCPU ZsNorm;
	MatrixCPU VoltageLimitCPU;
	MatrixCPU indiceBusBeginCPU; // size : nBus
	MatrixCPU CoresChatBeginCPU;


	MatrixGPU Cost1;
	MatrixGPU Cost2;
	MatrixGPU _CoresBusAgent;
	MatrixGPU _CoresAgentBus;
	MatrixGPU _CoresAgentBusBegin;

	MatrixGPU PosAgent;
	MatrixGPU _nAgentByBus;

	MatrixGPU ZsRe;
	MatrixGPU ZsIm;
	
	MatrixGPU VoltageLimit; // (vmin^2, vmax^2) * sqrt(Nchild + 1 / 2)
	MatrixGPU VoltageLimitReal; // vmin, vmax
	MatrixGPU FluxLimit;
	
	MatrixGPU X; // (Pi, Qi, li, vi, pn..., qn..., vai, Pci ..., Qci... , lci...) !!!!!
	MatrixGPU Ypre;
	MatrixGPU Y;// (Pi, Qi, li, vi, pn..., qn..., vai, Pci ..., Qci... , lci...) !!!!!
	MatrixGPU Mu;

	MatrixGPU Chat;

	MatrixGPU Hinv;
	MatrixGPU Q;

	MatrixGPU nChild;
	MatrixGPU Childs;
	MatrixGPU PosChild;
	MatrixGPU Ancestor;
	MatrixGPU CoresLineBus;

	MatrixGPU sizeOPFADMMConsGPU;
	MatrixGPU sizeOPFADMMConsGPUBig;

	// Local resolution
	MatrixGPU tempN2; // size : (_nAgent*2, 1)
	MatrixGPU tempB2; // size : (_nBus  *2, 1)
	MatrixGPU CoresSoloBusAgent;
	MatrixGPU Pn;
	MatrixGPU Pmin;
	MatrixGPU Pmax;
	MatrixGPU Pbmax;
	MatrixGPU Pbmin;
	MatrixGPU Pb;

	// special GPU
	MatrixGPU _CoresChatBegin;
	MatrixGPU tempL;
	MatrixGPU _indiceBusBegin; // size : nBus
	MatrixGPU _indiceBusBeginBig; // size : sizeOPFTotal
	MatrixGPU _indiceChildBegin;
	MatrixGPU PSOGPU; // dans le cas o� on onterface avec un march� sur CPU



	// consensus
	MatrixGPU etaSO; // consensus with market
	float _rhoSO = 0;



};


	

template <unsigned int _blockSizeSmall>
__global__ void updateChatGPU4(float* Chat, float* Y, float* MU, float* nChild, float* Ancestor, float* posChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, float* CoresChatBegin, float* indiceAgentBegin, float* CoresAgentBus, float* nAgentByBus, float _rho, int losstype, int nBus);


__global__ void updateConsensusGPU(float* Cost2, float* etaSO, float* Pn, float* Pmarket, float _rhoSO, int nAgent);