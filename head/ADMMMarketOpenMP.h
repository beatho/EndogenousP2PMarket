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
	virtual void setParam(float rho);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void solveWithMinPower(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	
	std::string NAME ="ADMMMarketOpenMP";
private:
	clock_t tMarket;
};


	




