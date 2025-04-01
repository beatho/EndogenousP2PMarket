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
	void setParam(float rho);
	void setTau(float tau);

	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	void updateGlobalProbGPU();
	void updateLocalProbGPU(MatrixGPU* Tlocal, MatrixGPU* P);
	void init(const Simparam& sim, const StudyCase& cas);
	void updateP0(const StudyCase& cas);
	std::string NAME ="ADMMGPUConst1T";


private:
	bool Ap1Changed = false;
	MatrixGPU Ap2a; // stocke la premi�re partie de Ap2a qui ne change par avec rho1
	MatrixGPU Ap1Copy;

};














