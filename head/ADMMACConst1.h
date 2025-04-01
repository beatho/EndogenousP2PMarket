#pragma once
#include "MethodP2P.h"
#include "CPUPF.h"
#include <iostream>
#include <string>
#include <chrono>


class ADMMACConst1 : public MethodP2P
{
public:
	ADMMACConst1();
	ADMMACConst1(float rho);
	virtual ~ADMMACConst1();
	void setParam(float rho);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="ADMMConst";
	void updateGlobalProb();
	void updateLocalProb();
	void updateLambda();
	void updateBt1();
	void updateBt2();
	void updateBp1();
	void updateTl();
	float calcRes();
	float updateResBis(MatrixCPU* res, int iter, MatrixCPU* tempNN);
	void updateP();
	void updateMU();
	void updateQ();
	void display();
private:
	// ne change pas avec P0
	float _mu = 40;
	float _mu1 = 40;
	float _mul = 40;
	float _tau = 2;

	float _rho = 0;
	float _rhol = 0;
	int _nAgent = 0;
	int _nTrade = 0;
	float _rhog = 0;
	float _at1 = 0;
	float _at2 = 0;

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

	// Pour le r�seau

	CPUPF PF;

	int _nLine = 0;
	int _nBus = 0;
	int _nConstraint = 0;
	float _rho1 = 0;
	int _nAgent2 = 0; // 2*_nAgent
	int _nTradeQ = 0; // N(N-1)
	int _nTradeP = 0;
	int _nVoisinRef = 0;
	float Ploss = 0;
	float Qloss = 0;
	int iterGlobal = 0;

	MatrixCPU PQ;


	MatrixCPU Kappa1;
	MatrixCPU Kappa2;
	MatrixCPU Kappa1_pre;
	MatrixCPU Kappa2_pre;
	MatrixCPU Cp1;
	MatrixCPU Cp2;
	MatrixCPU Ap3;
	MatrixCPU Ap123;
	MatrixCPU Qpart;
	MatrixCPU Qtot;
	MatrixCPU rhoMn;
	MatrixCPU rhoMn2;
	MatrixCPU UpperBound;
	MatrixCPU LowerBound;
	MatrixCPU DiffBound;
	MatrixCPU* Constraint;

	MatrixCPU tempL1;
	MatrixCPU tempC1;
	MatrixCPU tempCN;

	MatrixCPU* G;
	MatrixCPU G2; // G*G
	
};


	




