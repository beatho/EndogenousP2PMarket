#pragma once

#include <iostream>
#include <string>
#include <chrono>

// PF
#include "CPUPF.h"
#include "CPUPFdistPQ.h"
#include "ADMMMarket.h"

class EndoPF : public MethodP2P
{
public:
	EndoPF();
	EndoPF(float rho);
	virtual ~EndoPF();
	void setParam(float rho);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="EndoPF";
	void updateGlobalProb();
	void updateLocalProb();
	void updateLambda();
	void updateBt1();
	void updateBt2();
	void updateBp1();
	void updateTl();
	void updateSensi();
	void updateZ();
	void updateDelta();
	float calcRes();
	float updateResBis(MatrixCPU* res, int iter, MatrixCPU* tempNN);
	void updateP();
	void updateMU();
	virtual void updateCp2();
	void display();
	
	bool initWithMarketClear = true;
private:
	// ne change pas avec P0
	float _mu = 40;
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
	int _iterGlobal = 0;
	int _iterG = 0;
	int _stepG = 0;

	MatrixCPU tempNN; // Matrix temporaire pour aider les calculs
	MatrixCPU tempN1; // plutôt que de re-allouer de la mémoire à chaque utilisation


	MatrixCPU Tlocal;
	MatrixCPU P; // moyenne des trades
	MatrixCPU Pn; // somme des trades
	MatrixCPU Pnpre;
	MatrixCPU dP;

	MatrixCPU a;
	MatrixCPU Ap1;
	MatrixCPU Ap2;
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
	int _nVarPF; // _nLine + 2 * _nBus
	float _rho1;
	bool isRadial;
		
	MatrixCPU delta1; // \delta_1^ { k + 1 } = \delta_1^ { k } - G(P) + \overline{ l } - z_1
	MatrixCPU delta2; // \delta_2^ { k + 1 } = \delta_2^ { k } + G(P) + \overline{ l } - z_2
	MatrixCPU Z1; //z_1^ { k + 1 } = max(0, -Y + \overline{ Y } + \delta_1) 
	MatrixCPU Z2; //z_2^ { k + 1 } = max(0,  Y + \overline{ Y } + \delta_2)
	MatrixCPU Y;
	MatrixCPU Ypre;
	MatrixCPU Cp1;
	MatrixCPU Cp2;

	CPUPF* pf = nullptr;

	MatrixCPU tempL1;
	MatrixCPU G;
	MatrixCPU Ylimit;
	MatrixCPU dY;
	MatrixCPU SensiBis; // (z_1^ k - z_2 ^ k + \delta_2 ^ k - \delta_1 ^ k + 2 Y ^ k)
	MatrixCPU YOffset; // |Y - Yoffset| < Ylimit
};


	




