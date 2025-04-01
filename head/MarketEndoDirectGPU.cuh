#pragma once
#include <device_launch_parameters.h>
#include "MethodP2PGPU.cuh"
#include <iostream>
#include <string>
#include <chrono>
#include "Utilities.cuh"
#include "Utilities.h"

#include "ADMMMarket.h"
#include "ADMMMarketGPU.cuh"
#include "kernelFunction.cuh"


class MarketEndoDirectGPU : public MethodP2PGPU // OPFADMMCons mais on met des trades sur le bus fictif
{
public:

	
	MarketEndoDirectGPU();
	MarketEndoDirectGPU(float rho);
	virtual ~MarketEndoDirectGPU();
	void setParam(float rho);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	
	bool chekcase();
	std::string NAME ="MarketEndoDirectGPU";
	void updateGlobalProb();
	
	virtual void updateX();
	void updateXWOCurrent();
	void updateXWOCurrentCPU();

	virtual float getPLoss();
	virtual float getQLoss();

	void computeLoss();

	void updateMu();
	void updateChat();
	void CommunicationX();
	virtual float updateResEndo(int indice);
	float updateResRhoFixe(int indice);
	virtual int feasiblePoint();
	void ComputePFromAgentToBus();

	//Market
	void initMarket(const Simparam& sim, const StudyCase& cas);
	
	void updateBp2();
	
	void updateLocalProb();
	
	void updatePMarket(); // suite de updateX
	
	bool initWithMarketClear = true; //false
	void display();
	virtual MatrixCPU getPb();
	virtual MatrixCPU getPhi();
	virtual MatrixCPU getE();
private:

	// ne change pas avec P0
	
	int _blockSizeSmall = 32;
	int _numBlocksB;
	int _numBlocksH;

	int _nBus = 0;
	int _nBusWLoss = 0;
	
	int _nAgentOn0 = 0;
	int _sizeEndoMarketTotal = 0;
	int _debutloss = 0;
	int _sizeEndoMarketMax = 0;
	int _sizeChat = 0;
	int _numLineByBlockY = 1;
	
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
	float _epsL = 0;
	int _nIterL = 0;

	
	float _Ploss = 0;
	float _Qloss = 0;
	clock_t timeMarketEndo = 0;

	// market
	int _nAgentTrue = 0; // _nAgent = _nAgentTrue + (isAc)*_nAgent
	int _nTradeP = 0;
	int _nTradeQ = 0;

	
	LossType losstype = LossType::CURRENT;

	// Reste sur CPU
	
	MatrixCPU _nAgentByBusCPU;
	MatrixCPU resF;
	MatrixCPU nChildCPU;
	MatrixCPU CoresLineBusCPU;
	
	MatrixCPU ZsNorm;
	MatrixCPU indiceBusBeginCPU;
	MatrixCPU CoresChatBeginCPU;
	MatrixCPU VoltageLimitCPU;

	MatrixGPU tempB2; // size : (_nBus  *2, 1)


	MatrixGPU Ap3;
	MatrixGPU Ap123;

	MatrixGPU Bt2;
	MatrixGPU Bp1;
	MatrixGPU Bp2;
	

	// Reseau
	MatrixGPU tempM1; //
	MatrixGPU tempM; //
	MatrixGPU Pbmax;
	MatrixGPU Pbmin;
	MatrixGPU Pb;
	MatrixGPU CoresSoloBusAgent;

	MatrixGPU _CoresBusAgent;
	

	MatrixGPU PosAgent;
	MatrixGPU _nAgentByBus;


	
	MatrixGPU VoltageLimit; // (vmin^2, vmax^2) * sqrt(Nchild + 1 / 2)
	MatrixGPU VoltageLimitReal; // vmin, vmax
	//MatrixGPU FluxLimit;
	
	MatrixGPU X; // (Pi, Qi, li, vi, pn..., qn..., vai, Pci ..., Qci... , lci...) !!!!!
	MatrixGPU Ypre;
	MatrixGPU Y;// (Pi, Qi, li, vi, pn..., qn..., vai, Pci ..., Qci... , lci...) !!!!!
	MatrixGPU Mu;

	MatrixGPU Chat;

	MatrixGPU Hinv;
	MatrixGPU Q;

	MatrixGPU ZsRe;
	MatrixGPU ZsIm;

	
	
	MatrixGPU Childs;
	MatrixGPU PosChild; 
	MatrixGPU Ancestor;
	MatrixGPU CoresLineBus;

	MatrixGPU sizeMarketEndoDirectGPU;
	MatrixGPU sizeMarketEndoDirectGPUBig;



	// special GPU
	MatrixGPU _CoresAgentBus;
	MatrixGPU _CoresAgentBusBegin;
	MatrixGPU _CoresChatBegin;
	MatrixGPU tempL;
	MatrixGPU _indiceBusBegin; // size : nBus
	MatrixGPU _indiceBusBeginBig; // size : sizeOPFTotal
	MatrixGPU _indiceChildBegin;
	MatrixGPU nChild;

};

template <unsigned int _blockSizeSmall>
__global__ void updateChatGPU3(float* Chat, float* Y, float* MU, float* nChild, float* Ancestor, float* posChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, float* CoresChatBegin, float* nAgentByBus, float _rho, int LossType, int nBus);
	

__global__ void updateBp2GPU(float* Bp2, float* Y, float* MU, float* indiceBusBegin, float* indiceAgentBegin, float* CoresAgentBus, float* nAgentByBus, int losstype, float rho, int nBus, int nAgent);


