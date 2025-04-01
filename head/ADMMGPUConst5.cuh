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


class ADMMGPUConst5 : public MethodP2PGPU
{
public:
	ADMMGPUConst5();
	ADMMGPUConst5(float rho);
	virtual ~ADMMGPUConst5();
	void setParam(float rho);
	void setTau(float tau);

	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	void updateGlobalProbGPU();
	void updateLocalProbGPU(float epsL, int nIterL);
	void init(const Simparam& sim, const StudyCase& cas);
	void updateP0(const StudyCase& cas);
	std::string NAME ="ADMMGPUConst5";
		
private:
	
	MatrixGPU Ap2a; // a *Mn^2
	MatrixGPU Ap2b; // Mn^2 * sum(G2)
	
};

__device__ float warpReduceMax5(volatile float* r);

template <unsigned int blockSize>
__global__ void maxMonoBlock5(float* g_idata, float* g_odata, unsigned int n, unsigned int id);



