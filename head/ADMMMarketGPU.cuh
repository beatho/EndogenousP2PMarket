#pragma once
#include <device_launch_parameters.h>
#include "MethodP2PGPU.cuh"
#include "kernelFunction.cuh"
#include <iostream>
#include <string>
#include <chrono>


class ADMMMarketGPU : public MethodP2PGPU
{
public:
	ADMMMarketGPU();
	ADMMMarketGPU(float rho);
	virtual ~ADMMMarketGPU();
	virtual void setParam(float rho);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void solveWithMinPower(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="ADMMMarketGPU";
	void updateGlobalProb();
	void updateLocalProbGPU(float epsL, int nIterL);
	
	virtual void display();
private:
	// ne change pas avec P0
	
	int _iterGlobal = 0;
	clock_t tMarket;


};

__global__ void setMinPowerForSolve(float* Pmax, float* Pmin, int nCons);




