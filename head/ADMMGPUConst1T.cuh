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


class ADMMGPUConst1T : public MethodP2PGPU
{
public:
	ADMMGPUConst1T();
	ADMMGPUConst1T(float rho);
	virtual ~ADMMGPUConst1T();
	virtual void setParam(float rho);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	void updateGlobalProbGPU();
	void updateLocalProbGPU(MatrixGPU* Tlocal, MatrixGPU* P);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="ADMMGPUConst1T";


private:
	bool Ap1Changed = false;
	MatrixGPU Ap1Copy;

};














