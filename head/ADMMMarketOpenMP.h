#pragma once
#include "MethodP2P.h"
#include <iostream>
#include <string>
#include <chrono>

#ifdef _OPENMP
	#include "omp.h"
#endif


class ADMMMarketOpenMP : public MethodP2P
{
public:
	ADMMMarketOpenMP();
	ADMMMarketOpenMP(float rho);
	virtual ~ADMMMarketOpenMP();
	void setParam(float rho);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void solveWithMinPower(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	
	std::string NAME ="ADMMMarketOpenMP";
	
	void display();
private:
	// ne change pas avec P0
	float _mu = 50;
	float _tau = 2;

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
	bool isAC = false;
	int _iterGlobal = 0;
	int _iterG = 0;
	int _stepG = 0;
	clock_t tMarket;

	MatrixCPU tempNN; // Matrix temporaire pour aider les calculs
	MatrixCPU tempN1; // plut�t que de re-allouer de la m�moire � chaque utilisation


	MatrixCPU Tlocal;
	MatrixCPU P; // moyenne des trades
	MatrixCPU Pn; // somme des trades

	MatrixCPU a;
	MatrixCPU Ap2;
	MatrixCPU Ap1;
	MatrixCPU Ap12;
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
	MatrixCPU resX;



	// change avec P0
	MatrixCPU b;
	MatrixCPU matLb;
	MatrixCPU Cp;
	MatrixCPU Pmin;
	MatrixCPU Pmax;

};


	




