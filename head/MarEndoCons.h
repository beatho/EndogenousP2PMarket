#pragma once
#include "MethodP2P.h"
#include <iostream>
#include <string>
#include <chrono>

#include "OPFADMMCons.h"
#include "ADMMMarket.h"



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
	float _delta = 0.01f;
		
	// Pour le r�seau

	float _rhoSO = 0;
	bool _radial = false;
	
	MethodOPF* OPF = nullptr;
	MatrixCPU PSO;
	MatrixCPU etaSO;
	MatrixCPU Bp3;
	Simparam paramOPF;
};


	




