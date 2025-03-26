#pragma once
#include <device_launch_parameters.h>
#include "MethodP2P.cuh"
#include "kernelFunction.cuh"
#include <iostream>
#include <string>
#include <chrono>


class ADMMMarketGPU : public MethodP2P
{
public:
	ADMMMarketGPU();
	ADMMMarketGPU(float rho);
	virtual ~ADMMMarketGPU();
	void setParam(float rho);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void solveWithMinPower(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="ADMMMarketGPU";
	void updateGlobalProb();
	void updateLocalProbGPU(float epsL, int nIterL);
	float updateRes(MatrixCPU* res, int iter, MatrixGPU* tempNN);

	void display();
private:
	// ne change pas avec P0
	float _mu = 50;
	float _tau = 2;

	float _rho = 0;
	float _rhol = 0;
	int _blockSize = 512;
	int _numBlocksN = 0;
	int _numBlocksM = 0;
	
	int _nAgent = 0;
	int _nAgentTrue = 0; // _nAgent = _nAgentTrue + (isAc)*_nAgent
	int _nTrade = 0;
	int _nTradeP = 0;
	int _nTradeQ = 0;
	float _rhog = 0;
	float _at1 = 0;
	float _at2 = 0;
	bool isAC = false;
	int _iterGlobal = 0;
	int _iterG = 0;
	int _stepG = 0;
	clock_t tMarket;

	MatrixGPU tempNN; // Matrix temporaire pour aider les calculs
	MatrixGPU tempN1; // plut�t que de re-allouer de la m�moire � chaque utilisation


	MatrixGPU Tlocal;
	MatrixGPU P; // moyenne des trades
	MatrixGPU Pn; // somme des trades

	MatrixGPU a;
	MatrixGPU Ap2;
	MatrixGPU Ap1;
	MatrixGPU Ap12;
	MatrixGPU Bt1;
	//MatrixGPU Bt2;
	//MatrixGPU Bp1;
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

	// Matrices kept on CPU
	MatrixCPU nVoisinCPU;
	MatrixCPU LAMBDA;
	MatrixCPU trade;
	MatrixCPU resF;
	MatrixCPU resX;



	// change avec P0
	MatrixGPU b;
	MatrixGPU matLb;
	MatrixGPU Cp;
	MatrixGPU Pmin;
	MatrixGPU Pmax;

};

__global__ void setMinPowerForSolve(float* Pmax, float* Pmin, int nCons);




