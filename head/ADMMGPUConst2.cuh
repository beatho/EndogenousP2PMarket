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


class ADMMGPUConst2 : public MethodP2PGPU
{
public:
	ADMMGPUConst2();
	ADMMGPUConst2(float rho);
	virtual ~ADMMGPUConst2();
	virtual void setParam(float rho);
	void setTau(float tau);

	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	void updateGlobalProbGPU();
	void updateLocalProbGPU(MatrixGPU* Tlocal, MatrixGPU* P);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="ADMMGPUConst2";
		

private:
	

};


