#pragma once
#include "MethodP2P.cuh"
#include <iostream>
#include <string>
#include <chrono>


class ADMMConst1 : public MethodP2P
{
public:
	ADMMConst1();
	ADMMConst1(float rho);
	virtual ~ADMMConst1();
	void setParam(float rho);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="ADMMConst1";
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
	MatrixCPU tempN1; // plutôt que de re-allouer de la mémoire à chaque utilisation


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

	// Pour le réseau
	int _nLine;
	int _nBus;
	float _rho1;

	MatrixCPU Kappa1;
	MatrixCPU Kappa2;
	MatrixCPU Kappa1_pre;
	MatrixCPU Kappa2_pre;
	MatrixCPU Cp1;
	MatrixCPU Cp2;
	MatrixCPU Qpart;
	MatrixCPU Qtot;
	MatrixCPU alpha;

	MatrixCPU tempL1;
	MatrixCPU G;
	MatrixCPU GTrans;
	MatrixCPU G2; // G^t.*G^t
	MatrixCPU lLimit;
};


	




