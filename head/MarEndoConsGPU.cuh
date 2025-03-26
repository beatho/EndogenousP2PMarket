#pragma once
#include <device_launch_parameters.h>
#include "MethodP2P.cuh"
#include <iostream>
#include <string>
#include <chrono>

#include "OPFADMMConsGPU.cuh"
#include "OPFADMMCons.h"
#include "ADMMMarketGPU.cuh"


class MarEndoConsGPU : public MethodP2P
{
public:
	MarEndoConsGPU();
	MarEndoConsGPU(float rho);
	virtual ~MarEndoConsGPU();
	void setParam(float rho);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="MarEndoConsGPU";
	void updateGlobalProb();
	void updateLocalProbGPU(float epsL2, int iterL);
	 
	void updateEtaSO();
	void updateBp3();
	  
	float updateResBis(int iter );
	 
	 
	void display();
	bool initWithMarketClear = true;
	bool OPFonCPU = false; // false true
private:
	// ne change pas avec P0
	float _mu = 40;
	float _mu1 = 40;
	float _mul = 40;
	float _tau = 2;
	float _delta = 1;
	float _epsLim = 0.1;
	int _stepIntern = 0;

	int _blockSize = 512;
	int _numBlocksN = 0;
	int _numBlocksM = 0;


	float _rho = 0;
	float _rhol = 0;
	int _nAgent = 0;
	int _nAgentTrue = 0;
	int _nTrade = 0;
	int _nTradeP = 0;
	int _nTradeQ = 0;
	float _rhog = 0;
	float _at1 = 0;
	float _at2 = 0;
	float _resG = 0;
	int _iterGlobal = 0;
	int _stepG = 0;
	int _iterG = 0;
	int _stepL = 0;
	clock_t timeMarketEndo;

	MatrixGPU tempNN; // Matrix temporaire pour aider les calculs
	MatrixGPU tempN1; // plut�t que de re-allouer de la m�moire � chaque utilisation


	MatrixGPU Tlocal;
	MatrixGPU P; // moyenne des trades
	MatrixGPU Pn; // somme des trades

	MatrixGPU a;
	MatrixGPU Ap2;
	MatrixGPU Ap1;
	MatrixGPU Ap123;
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
	//MatrixCPU resX;
	MatrixCPU nVoisinCPU;
	MatrixCPU PSOCPU;
	MatrixCPU PnCPU;



	// change avec P0
	MatrixGPU b;
	MatrixGPU matLb;
	MatrixGPU Cp;
	MatrixGPU Pmin;
	MatrixGPU Pmax;

	// Pour le r�seau

	float _rhoSO = 0;
	MatrixGPU Ap3;
	bool _radial = false;
	
	MethodOPFGPU* OPF = nullptr;
	MethodOPF* OPFCPU = nullptr;
	MatrixGPU PSO;
	MatrixGPU etaSO;
	MatrixGPU Bp3;
	Simparam paramOPF;
};


	




