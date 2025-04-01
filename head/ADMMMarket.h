#pragma once
#include "MethodP2P.h"
#include <iostream>
#include <string>
#include <chrono>

#ifdef _OPENMP
	#include "omp.h"
#endif

class ADMMMarket : public MethodP2P
{
public:
	ADMMMarket();
	ADMMMarket(float rho);
	virtual ~ADMMMarket();
	void setParam(float rho);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void solveWithMinPower(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="ADMMMarket";
	void updateGlobalProb();
	void updateLocalProb();
	void updateBt1();
	void updateBt2();
	void updateBp1();
	void updateTl();
	
	void updateP();
	void updateMU();
	virtual void display();
private:
	// ne change pas avec P0
	int _nAgentTrue = 0; // _nAgent = _nAgentTrue + (isAc)*_nAgent

	int _nTradeP = 0;
	int _nTradeQ = 0;
	
	bool isAC = false;
	
	int _iterGlobal = 0;
	int _iterG = 0;
	int _stepG = 0;
	clock_t tMarket;


};


	




