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
	void updateBt1();
	void updateBt2();
	void updateBp1();
	void updateTl();
	void updateSensi();
	void updateZ();
	void updateDelta();
	virtual float updateResEndo(int iter);
	void updateP();
	void updateMU();
	virtual void updateCp2();
	void display();
	
	bool initWithMarketClear = true;
private:
	// ne change pas avec P0
	
	int _nAgentTrue = 0; // _nAgent = _nAgentTrue + (isAc)*_nAgent
	int _nTradeP = 0;
	int _nTradeQ = 0;
	

	int _iterGlobal = 0;
	int _iterG = 0;
	int _stepG = 0;

	
	MatrixCPU Pnpre;
	MatrixCPU dP;

	// Pour le r�seau
	
	int _nVarPF; // _nLine + 2 * _nBus

	bool isRadial;
		
	MatrixCPU delta1; // \delta_1^ { k + 1 } = \delta_1^ { k } - G(P) + \overline{ l } - z_1
	MatrixCPU delta2; // \delta_2^ { k + 1 } = \delta_2^ { k } + G(P) + \overline{ l } - z_2
	MatrixCPU Z1; //z_1^ { k + 1 } = max(0, -Y + \overline{ Y } + \delta_1) 
	MatrixCPU Z2; //z_2^ { k + 1 } = max(0,  Y + \overline{ Y } + \delta_2)
	MatrixCPU Y;
	MatrixCPU Ypre;
	
	CPUPF* pf = nullptr;


	MatrixCPU Ylimit;
	MatrixCPU dY;
	MatrixCPU SensiBis; // (z_1^ k - z_2 ^ k + \delta_2 ^ k - \delta_1 ^ k + 2 Y ^ k)
	MatrixCPU YOffset; // |Y - Yoffset| < Ylimit
};


	




