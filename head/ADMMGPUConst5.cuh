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
	virtual void setParam(float rho);
	void setTau(float tau);

	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	void updateGlobalProbGPU();
	void updateLocalProbGPU(float epsL, int nIterL);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="ADMMGPUConst5";
		
private:
		
};

__device__ float warpReduceMax5(volatile float* r);

template <unsigned int blockSize>
__global__ void maxMonoBlock5(float* g_idata, float* g_odata, unsigned int n, unsigned int id);



