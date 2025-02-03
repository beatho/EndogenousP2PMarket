#pragma once
#include "MethodOPFGPU.cuh"
#include <iostream>
#include <string>
#include <chrono>
#include "Utilities.cuh"
#include "Utilities.h"

#include "kernelFunction.cuh"


class OPFADMMGPU2 : public MethodOPFGPU
{
public:
	OPFADMMGPU2();
	OPFADMMGPU2(float rho);
	virtual ~OPFADMMGPU2();
	void setParam(float rho);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);

	virtual void solveConsensus(float eps, MatrixCPU* PSO);
	virtual void initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO);
	virtual void updateConsensus(MatrixCPU* Pmarket);

	virtual void solveConsensus(float eps, MatrixGPU* PSO);
	virtual void updateConsensus(MatrixGPU* Pmarket);


	bool chekcase();
	std::string NAME ="OPFADMMGPU2";
	void updateGlobalProb();

	void updateMu();

	virtual float getPLoss();
	virtual float getQLoss();

	void ComputePFromAgentToBus();

	void updateChat();
	void CommunicationX();
	float updateRes(int indice);

	virtual int feasiblePoint();

	virtual void display();
private:
	// ne change pas avec P0
	int _numBlocksB = 1;
	int _numBlocksN = 1;
	int _numBlocksM = 1;

	int _nBus = 0;
	int _nLine = 0;
	int _nAgent = 0;
	int _nAgentOn0 = 0;
	int _sizeOPFTotal = 0;
	int _sizeChat = 0;
	int _sizeOPFMax = 0;
	float _rho = 0;
	float _rhoInv = 0;
	MatrixGPU root; // peut �tre juste defini en local c'est suffisant
	MatrixGPU coefPoly2;
	MatrixGPU coefPoly3;
	float _rhol = 0;
	int _iterLocal = 0;
	int _iterGlobal = 0;
	int _iterG = 0;
	int _stepG = 0;
	
	bool consensus = false;
	clock_t timeOPF = 0;
	
	MatrixGPU _apt1;

	// parameter agent and iteration dependant (but not here for now)


	MatrixGPU tempN1; // Matrix temporaire pour aider les calculs
	MatrixGPU tempNN; // plut�t que de re-allouer de la m�moire � chaque utilisation
	MatrixGPU tempM1; //
	MatrixGPU tempM; //

	

	MatrixGPU Cost1;
	MatrixGPU Cost2;
	MatrixGPU _CoresBusAgent;
	MatrixGPU _nAgentByBus;
	MatrixCPU _nAgentByBusCPU;

	MatrixGPU ZsRe;
	MatrixGPU ZsIm;
	MatrixGPU ZsNorm;
	MatrixGPU VoltageLimit; // (vmin^2, vmax^2) * sqrt(Nchild + 1 / 2)
	MatrixGPU VoltageLimitReal; // vmin, vmax
	
	MatrixGPU X; // (Pi, Qi, li, vi, pn..., qn..., vai, Pci ..., Qci... , lci...) !!!!!
	MatrixGPU Ypre;
	MatrixGPU Y;// (Pi, Qi, li, vi, pn..., qn..., vai, Pci ..., Qci... , lci...) !!!!!
	MatrixGPU Mu;

	MatrixGPU Chat;

	MatrixGPU Hinv;
	MatrixGPU Q;

	MatrixCPU nChildCPU;
	MatrixGPU Childs;
	MatrixGPU PosChild; 
	MatrixGPU Ancestor;
	MatrixGPU CoresLineBus;

	MatrixGPU sizeOPFADMMGPU2;
	MatrixGPU sizeOPFADMMGPU2Big; // size : sizeOPFTotal

	MatrixCPU resF;

	// Local resolution
	MatrixGPU tempN2; // size : (_nAgent*2, 1)
	MatrixGPU tempB2; // size : (_nBus  *2, 1)
	MatrixGPU CoresSoloBusAgent;
	MatrixGPU Pn;
	MatrixGPU Pmin;
	MatrixGPU Pmax;
	MatrixGPU Pbmin;
	MatrixGPU Pbmax;
	MatrixGPU Pb;


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
__global__ void updateChatGPU2(float* Chat, float* Y, float* MU, float* nChild, float* Ancestor, float* posChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, float* CoresChatBegin, float* nAgentByBus, float _rho, int nBus);
