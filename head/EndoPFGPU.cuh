#pragma once

#include <iostream>
#include <string>
#include <chrono>

// PF
#include "GPUPF.cuh"
#include "GPUPFdistPQ.cuh"
#include "ADMMMarketGPU.cuh"

class EndoPFGPU : public MethodP2P
{
public:
	EndoPFGPU();
	EndoPFGPU(float rho);
	virtual ~EndoPFGPU();
	void setParam(float rho);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="EndoPFGPU";
	void updateGlobalProbGPU();
	void updateLocalProbGPU(float epsL, int nIterL);
	
	void updateCp2();
	void updateSensi();
	float updateResBis(MatrixCPU* res, MatrixGPU* Tlocal, int iter, MatrixGPU* tempNN);
	
	
	void display();
	
	bool initWithMarketClear = true;
private:
	// ne change pas avec P0
	float _mu = 40;
	float _tau = 2;

	int _blockSize = 256;
	int _numBlocksN = 0;
	int _numBlocksM = 0;
	int _numBlocksL = 0;
	int _numBlocksBL = 0;

	float _rho = 0;
	float _rhol = 0;
	int _nAgent = 0;
	int _nAgentTrue = 0; // _nAgent = _nAgentTrue + (isAc)*_nAgent
	int _nTrade = 0;
	int _nTradeP = 0;
	int _nTradeQ = 0;
	float _rhog = 0;
	float _at1 = 0;
	float _at2 = 0;
	int _iterGlobal = 0;
	int _iterG = 0;
	int _stepG = 0;
	clock_t timeEndoPF = 0;

	MatrixGPU tempNN; // Matrix temporaire pour aider les calculs
	MatrixGPU tempN1; // plutôt que de re-allouer de la mémoire à chaque utilisation


	MatrixGPU Tlocal;
	MatrixGPU P; // moyenne des trades
	MatrixGPU Pn; // somme des trades
	MatrixGPU Pnpre;
	MatrixGPU dP;

	MatrixGPU a;
	MatrixGPU Ap1;
	MatrixGPU Ap2;
	MatrixGPU Ap12;
	MatrixGPU Bt1;
	MatrixGPU Bt2;
	MatrixGPU Bp1;
	MatrixGPU Ct;
	MatrixGPU matUb;

	MatrixGPU nVoisin;
	MatrixGPU tradeLin;
	MatrixGPU Tmoy;
	MatrixGPU LAMBDALin;
	MatrixGPU Tlocal_pre;
	MatrixGPU MU;

	MatrixGPU CoresMatLin;
	MatrixGPU CoresLinAgent;
	MatrixGPU CoresAgentLin;
	MatrixGPU CoresLinVoisin;
	MatrixGPU CoresLinTrans;

	// Matrices for the result
	MatrixCPU LAMBDA;
	MatrixCPU trade;
	MatrixCPU resF;
	MatrixCPU resX;
	MatrixCPU nVoisinCPU;



	// change avec P0
	MatrixGPU b;
	MatrixGPU matLb;
	MatrixGPU Cp;
	MatrixGPU Pmin;
	MatrixGPU Pmax;

	// Pour le réseau
	int _nLine;
	int _nBus;
	int _nVarPF; // _nLine + 2 * _nBus
	float _rho1;
	bool isRadial;
		
	MatrixGPU delta1; // \delta_1^ { k + 1 } = \delta_1^ { k } - G(P) + \overline{ l } - z_1
	MatrixGPU delta2; // \delta_2^ { k + 1 } = \delta_2^ { k } + G(P) + \overline{ l } - z_2
	MatrixGPU Z1; //z_1^ { k + 1 } = max(0, -Y + \overline{ Y } + \delta_1) 
	MatrixGPU Z2; //z_2^ { k + 1 } = max(0,  Y + \overline{ Y } + \delta_2)
	MatrixGPU Y;
	MatrixGPU Ypre;
	MatrixGPU Cp1;
	MatrixGPU Cp2;

	GPUPF* pf = nullptr;

	MatrixGPU tempL1;
	MatrixGPU G;
	MatrixGPU Ylimit;
	MatrixGPU dY;
	MatrixGPU SensiBis; // (z_1^ k - z_2 ^ k + \delta_2 ^ k - \delta_1 ^ k + 2 Y ^ k)
	MatrixGPU YOffset; // |Y - Yoffset| < Ylimit
};





__global__ void initLimits(float* Ylimit, float* Yoffset, float* limitsLb, float* limitsUb, int nVarPF);

__global__ void updateZGPU(float* Z1, float* Z2, float* Ylimit, float* delta1, float* delta2, float* Y, int nVarPF);
__global__ void updateDeltaGPU(float* delta1, float* delta2, float* Z1, float* Z2, float* Ylimit, float* Y, int nVarPF);
__global__ void updateZDeltaGPU(float* Z1, float* Z2, float* Ylimit, float* delta1, float* delta2, float* Y, int nVarPF);

template <unsigned int blockSize>
__global__ void updateCp2GPU(float* Cp2, float* SensiBis, float* G, float* nVoisin, float rho1, int nVarPF);
__global__ void updateSensiBis(float* sensiBis, float* Y, float* Z1, float* Z2, float* delta1, float* delta2, int nVarPF);
__global__ void updateSensiGPU(float* G, float* dY, float* dP, int nVar);
	




