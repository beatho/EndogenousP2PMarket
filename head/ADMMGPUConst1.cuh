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


class ADMMGPUConst1 : public MethodP2PGPU
{
public:
	ADMMGPUConst1();
	ADMMGPUConst1(float rho);
	virtual ~ADMMGPUConst1();
	void setParam(float rho);
	void setTau(float tau);

	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	void updateGlobalProbGPU();
	void updateLocalProbGPU(MatrixGPU* Tlocal, MatrixGPU* P);
	void init(const Simparam& sim, const StudyCase& cas);
	void updateP0(const StudyCase& cas);
	std::string NAME ="ADMMGPUConst1";

private:

	MatrixGPU Ap1Copy;
};


