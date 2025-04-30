#pragma once

#include <device_launch_parameters.h>
#include <iostream>
#include <string>
#include <chrono>

// PF
#include "GPUPF.cuh"
#include "GPUPFdistPQ.cuh"
#include "ADMMMarketGPU.cuh"

class EndoPFGPU : public MethodP2PGPU
{
public:
	EndoPFGPU();
	EndoPFGPU(float rho);
	virtual ~EndoPFGPU();
	void setParam(float rho);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="EndoPFGPU";
	void updateGlobalProbGPU();
	void updateLocalProbGPU(float epsL, int nIterL);
	
	void updateCp2();
	void updateSensi();
	virtual float updateResEndo(int iter);
	
	
	void display();
	
	bool initWithMarketClear = true;
private:
	// ne change pas avec P0
	int _numBlocksBL = 0;
	
	MatrixGPU Pnpre;
	MatrixGPU dP;
	
	// Pour le reseau
	
	int _nVarPF; // _nLine + 2 * _nBus
	bool isRadial;
		
	MatrixGPU delta1; // \delta_1^ { k + 1 } = \delta_1^ { k } - G(P) + \overline{ l } - z_1
	MatrixGPU delta2; // \delta_2^ { k + 1 } = \delta_2^ { k } + G(P) + \overline{ l } - z_2
	MatrixGPU Z1; //z_1^ { k + 1 } = max(0, -Y + \overline{ Y } + \delta_1) 
	MatrixGPU Z2; //z_2^ { k + 1 } = max(0,  Y + \overline{ Y } + \delta_2)
	MatrixGPU Y;
	MatrixGPU Ypre;
	


	GPUPF* pf = nullptr;

	MatrixGPU Ylimit;
	MatrixGPU dY;
	MatrixGPU SensiBis; // (z_1^ k - z_2 ^ k + \delta_2 ^ k - \delta_1 ^ k + 2 Y ^ k)
	MatrixGPU YOffset; // |Y - Yoffset| < Ylimit
};





__global__ void initLimits(float* Ylimit, float* Yoffset, float* limitsLb, float* limitsUb, int nVarPF);

__global__ void updateZGPU(float* Z1, float* Z2, float* Ylimit, float* delta1, float* delta2, float* Y, int nVarPF);
__global__ void updateDeltaGPU(float* delta1, float* delta2, float* Z1, float* Z2, float* Ylimit, float* Y, int nVarPF);
__global__ void updateZDeltaGPU(float* Z1, float* Z2, float* Ylimit, float* delta1, float* delta2, float* Y, int nVarPF);

template <unsigned int blockSize>
__global__ void updateCp2GPU(float* Cp2, float* SensiBis, float* G, float* nVoisin, float rho1, int nVarPF);
__global__ void updateSensiBis(float* sensiBis, float* Y, float* Z1, float* Z2, float* delta1, float* delta2, int nVarPF);
__global__ void updateSensiGPU(float* G, float* dY, float* dP, int nVar);
	




