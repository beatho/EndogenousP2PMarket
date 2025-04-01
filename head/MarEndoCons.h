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
	float _delta = 0.01;
	float _epsLim = 0.1;


	int _nAgentTrue = 0;
	int _nTradeP = 0;
	int _nTradeQ = 0;


	float _resG = 0;
	int _iterGlobal = 0;
	int _stepG = 0;
	int _iterG = 0;
	int _stepL = 0;
	int _stepIntern = 0;

	clock_t timeMarketEndo;

	
	MatrixCPU Ap123;


	// Pour le r�seau

	float _rhoSO = 0;
	MatrixCPU Ap3;
	bool _radial = false;
	
	MethodOPF* OPF = nullptr;
	MatrixCPU PSO;
	MatrixCPU etaSO;
	MatrixCPU Bp3;
	Simparam paramOPF;
};


	




