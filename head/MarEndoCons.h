#pragma once
#include "MethodP2P.cuh"
#include <iostream>
#include <string>
#include <chrono>

#include "OPFADMMCons.h"
#include "OPFPDIPM.h"


class MarEndoCons : public MethodP2P
{
public:
	MarEndoCons();
	MarEndoCons(float rho);
	virtual ~MarEndoCons();
	void setParam(float rho);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="MarEndoCons";
	void updateGlobalProb();
	void updateLocalProb();
	void updateLambda();
	void updateEtaSO();
	void updateBp3();
	void updateBt1();
	void updateBt2();
	void updateBp1();
	void updateTl();
	float calcRes();
	float updateResBis(MatrixCPU* res, int iter, MatrixCPU* tempNN);
	void updateP();
	void updateMU();
	void display();
	bool initWithMarketClear = true;
private:
	// ne change pas avec P0
	float _mu = 40;
	float _mu1 = 40;
	float _mul = 40;
	float _tau = 2;
	float _delta = 0.01;
	float _epsLim = 0.1;

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
	int _stepIntern = 0;

	clock_t timeMarketEndo;

	MatrixCPU tempNN; // Matrix temporaire pour aider les calculs
	MatrixCPU tempN1; // plutôt que de re-allouer de la mémoire à chaque utilisation


	MatrixCPU Tlocal;
	MatrixCPU P; // moyenne des trades
	MatrixCPU Pn; // somme des trades

	MatrixCPU a;
	MatrixCPU Ap2;
	MatrixCPU Ap1;
	MatrixCPU Ap123;
	MatrixCPU Bt1;
	MatrixCPU Bt2;
	MatrixCPU Bp1;
	MatrixCPU Ct;
	MatrixCPU matUb;

	MatrixCPU nVoisin;
	MatrixCPU tradeLin;
	MatrixCPU Tmoy;
	MatrixCPU LAMBDALin;
	MatrixCPU Tlocal_pre;
	MatrixCPU MU;

	MatrixCPU CoresMatLin;
	MatrixCPU CoresLinAgent;
	MatrixCPU CoresAgentLin;
	MatrixCPU CoresLinVoisin;
	MatrixCPU CoresLinTrans;

	// Matrices for the result
	MatrixCPU LAMBDA;
	MatrixCPU trade;
	MatrixCPU resF;
	//MatrixCPU resX;



	// change avec P0
	MatrixCPU b;
	MatrixCPU matLb;
	MatrixCPU Cp;
	MatrixCPU Pmin;
	MatrixCPU Pmax;

	// Pour le réseau

	float _rhoSO = 0;
	MatrixCPU Ap3;
	bool _radial = false;
	
	MethodOPF* OPF = nullptr;
	MatrixCPU PSO;
	MatrixCPU etaSO;
	MatrixCPU Bp3;
	Simparam paramOPF;
};


	




