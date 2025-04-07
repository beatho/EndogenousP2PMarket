#pragma once
#include <device_launch_parameters.h>
#include "MethodP2PGPU.cuh"
#include <iostream>
#include <string>
#include <chrono>

#include "OPFADMMConsGPU.cuh"
#include "OPFADMMCons.h"
#include "ADMMMarketGPU.cuh"


class MarEndoConsGPU : public MethodP2PGPU
{
public:
	MarEndoConsGPU();
	MarEndoConsGPU(float rho);
	virtual ~MarEndoConsGPU();
	void setParam(float rho);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="MarEndoConsGPU";
	void updateGlobalProb();
	void updateLocalProbGPU(float epsL2, int iterL);
	 
	void updateEtaSO();
	void updateBp3();
	  
	float updateResBis(int iter );
	 
	 
	void display();
	bool initWithMarketClear = true;
	bool OPFonCPU = false; // false true
private:
	// ne change pas avec P0
	
	float _delta = 1;
	float _resG = 0;
	

	// Matrices kept on CPU
	MatrixCPU PSOCPU;
	MatrixCPU PnCPU;

	// Pour le r�seau

	float _rhoSO = 0;
	bool _radial = false;
	
	MethodOPFGPU* OPF = nullptr;
	MethodOPF* OPFCPU = nullptr;
	MatrixGPU PSO;
	MatrixGPU etaSO;
	MatrixGPU Bp3;
	Simparam paramOPF;
};


	




