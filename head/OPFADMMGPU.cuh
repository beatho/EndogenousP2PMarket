#pragma once
#include <device_launch_parameters.h>
#include "MethodOPFGPU.cuh"
#include <iostream>
#include <string>
#include <chrono>
#include "Utilities.cuh"
#include "Utilities.h"

#include "kernelFunction.cuh"


class OPFADMMGPU : public MethodOPFGPU
{
public:
	OPFADMMGPU();
	OPFADMMGPU(float rho);
	virtual ~OPFADMMGPU();
	void setParam(float rho);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);

	virtual void solveConsensus(float eps, MatrixCPU* PSO);
	virtual void solveConsensus(float eps, MatrixGPU* PSO);
	virtual void initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO);
	virtual void updateConsensus(MatrixCPU* Pmarket);
	virtual void updateConsensus(MatrixGPU* Pmarket);

	bool chekcase();
	std::string NAME ="OPFADMMGPU";
	void updateGlobalProb();
	void updateLocalProb(float epsL, int nIterL);
	void updateXWOCurrentOnCPU();
	void updateMu();

	void ComputePFromAgentToBus();
	virtual MatrixCPU getPb();
	virtual MatrixCPU getPhi();
	virtual MatrixCPU getE();


	void updateChat();
	void CommunicationX();
	float updateRes(int indice);
	virtual int feasiblePoint();

	void display();
private:
	// ne change pas avec P0
	
	int _numBlocksB = 1;
	int _numBlocksN = 1;
	int _numBlocksM = 1;

	int _nBus = 0;
	int _nLine = 0;
	
	int _sizeOPFTotal = 0;
	int _sizeOPFMax = 0;
	float _rho = 0;
	float _rhoInv = 0;
	double coefPoly2[2];
	double root2[4];
	double root3[4];
	double root4[4];
	double root5[4];
	double root6[4];
	double coefPoly3[3];

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

	MatrixGPU ZsRe;
	MatrixGPU ZsIm;
	MatrixGPU ZsNorm;
	MatrixGPU VoltageLimit; // (vmin^2, vmax^2) * sqrt(Nchild + 1 / 2)
	MatrixGPU VoltageLimitReal; // vmin, vmax
	
	MatrixGPU X; // (Pi, Qi, li, vi, pi, qi, vai, Pci ..., Qci... , lci...) !!!!!
	MatrixGPU Ypre;
	MatrixGPU Y; // (Pi, Qi, li, vi, pi, qi, vai, Pci ..., Qci... , lci...) !!!!!
	//MatrixGPU YTrans; // (Pi,Qi,vi,li,pi,qi,pai,qai,lai,vij) !!!
	MatrixGPU Mu;

	MatrixGPU Chat;

	MatrixGPU Hinv;
	MatrixGPU Q;

	MatrixGPU Childs;
	MatrixGPU PosChild; 
	MatrixGPU Ancestor;
	

	MatrixGPU sizeOPFADMMGPU;
	MatrixGPU sizeOPFADMMGPUBig; // size : sizeOPFTotal


	MatrixCPU nChildCPU;
	MatrixCPU resF;
	MatrixCPU VoltageLimitCPU;
	MatrixCPU indiceBusBeginCPU; // size : nBus
	MatrixCPU nAgentByBusCPU;
	MatrixCPU CoresLineBus;

	// Local resolution
	MatrixGPU tempN2; // size : (_nAgent*2, 1)
	MatrixGPU tempB2; // size : (_nBus  *2, 1)
	MatrixGPU CoresSoloBusAgent;
	
	MatrixGPU Pmin;
	MatrixGPU Pmax;
	MatrixGPU PnTmin;
	MatrixGPU PnTmax;
	MatrixGPU PnMoy;
	MatrixGPU PnPre;
	MatrixGPU MuL;
	MatrixGPU PnTilde;
	MatrixGPU Bp1;
	MatrixGPU Bpt1;
	MatrixGPU Bpt2;
	MatrixGPU Apt1;
	MatrixGPU Apt2;
	

	// special GPU
	MatrixGPU _CoresAgentBus;
	MatrixGPU _CoresAgentBusBegin;
	MatrixGPU tempL;
	MatrixGPU _indiceBusBegin; // size : nBus
	MatrixGPU _indiceBusBeginBig; // size : sizeOPFTotal
	MatrixGPU _indiceChildBegin;
	MatrixGPU nChild;


};

//__global__ void setAncestorChild(float* coresLineBus, float* mat2, int N);





template <unsigned int blockSize>
__global__ void updateChatGPU(float* Chat, float* Y, float* MU, float* nChild, float* Ancestor, float* posChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin,
	float _rho, int nBus);

__global__ void updateBpt2(float* Bpt2, float* Chat, float* nAgentByBus, int nBus);

template <unsigned int blockSize>
__global__ void updatePnPGPUSharedResidual(float* Pn, float* PnPre, float* PnMoy, float* PnTilde, float* MUL, float* nAgentByBus, float _rhol, float* Ap2, float* Cp, float* Pmin,
	float* Pmax, float* Apt1, float* Apt2, float* Bpt2, float* CoresSoloBusAgent, float* CoresBusAgent, float* CoresBusAgentBegin, float eps, int nIterLMax, int nAgent, int nBus);


