#pragma once
#include <device_launch_parameters.h>
#include "MethodP2PGPU.cuh"
#include "MatrixGPU.cuh"
#include "MatrixCPU.h"
#include "kernelFunction.cuh"

#include <iostream>
#include <string>
#include <cuda_runtime.h>
#include <chrono>


class ADMMGPUConst4 : public MethodP2PGPU
{
public:
	ADMMGPUConst4();
	ADMMGPUConst4(float rho);
	virtual ~ADMMGPUConst4();
	void setParam(float rho);
	void setTau(float tau);

	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	void updateGlobalProbGPU();
	void updateLocalProbGPU(float epsL, int nIterL);
	void init(const Simparam& sim, const StudyCase& cas);
	void updateP0(const StudyCase& cas);
	std::string NAME ="ADMMGPUConst4";
	
private:
	MatrixGPU Ap2a; // stocke la premi�re partie de Ap2a qui ne change par avec rho1
	MatrixGPU Ap2b; // Mn^2 * sum(G2)

};







