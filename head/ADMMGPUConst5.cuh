#pragma once
#include <device_launch_parameters.h>
#include "MethodP2P.cuh"
#include "MatrixGPU.cuh"
#include "MatrixCPU.h"
#include "kernelFunction.cuh"

#include <iostream>
#include <string>
#include <cuda_runtime.h>
#include <chrono>


class ADMMGPUConst5 : public MethodP2P
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
	
	virtual float updateRes(MatrixCPU* res, MatrixGPU* Tlocal, int iter, MatrixGPU* tempNN);
	
	
	void display();

private:
	// ne change pas avec P0
	float _mu = 1000;
	float _mu1 = 1000;
	float _muInit = 1000;
	float _tau = 1.1;
	float _rho = 0;
	int _blockSize = 512;
	int _numBlocksN = 0;
	int _numBlocksM = 0;
	int _numBlocksL = 0;
	int _numBlocksNL = 0;
	int _nAgent = 0;
	int _nTrade = 0;
	float _rhog = 0;
	float _rhol = 0;
	float _at1 = 0;
	float _at2 = 0;

	MatrixGPU tempNN; // Matrix temporaire pour aider les calculs
	MatrixGPU tempN1; // plut�t que de re-allouer de la m�moire � chaque utilisation
	MatrixGPU tempL1;
	MatrixGPU tempL2; // il faut deux vecteurs de taille L
	
	MatrixGPU Tlocal;
	MatrixGPU P; // moyenne des trades
	MatrixGPU Pn; // somme des trades

	MatrixGPU a;
	MatrixGPU Ap2;
	MatrixGPU Ap2a; // a *Mn^2
	MatrixGPU Ap2b; // Mn^2 * sum(G2)
	MatrixGPU Ap1;
	MatrixGPU Ap12;
	MatrixGPU Bt1;
	MatrixGPU Ct;
	MatrixGPU matUb;
	
	MatrixGPU nVoisin;
	MatrixGPU tradeLin;
	MatrixGPU Tmoy;
	MatrixGPU LAMBDALin;
	MatrixGPU Tlocal_pre;
	MatrixGPU MU;

	MatrixGPU CoresMatLin;
	MatrixGPU CoresLinAgent;
	MatrixGPU CoresAgentLin;
	MatrixGPU CoresLinVoisin;
	MatrixGPU CoresLinTrans;

	// Matrices kept on CPU
	MatrixCPU nVoisinCPU;
	MatrixCPU LAMBDA;
	MatrixCPU trade;
	MatrixCPU resF;
	MatrixCPU resX;

	// change avec P0
	MatrixGPU b;
	MatrixGPU matLb;
	MatrixGPU Cp;
	MatrixGPU Pmin;
	MatrixGPU Pmax;
	
	// Pour le r�seau
	int _nLine;
	int _nBus;
	float _rho1;

	MatrixGPU Kappa1;
	MatrixGPU Kappa2;
	MatrixGPU Kappa1_pre;
	MatrixGPU Kappa2_pre;
	MatrixGPU Cp1;
	MatrixGPU Cp2;
	MatrixGPU Qpart;
	MatrixGPU Qtot;
	MatrixGPU alpha;

	MatrixGPU G;
	MatrixGPU GTrans;
	MatrixGPU G2; // G^t.*G^t
	MatrixGPU lLimit;


};

__device__ float warpReduceMax5(volatile float* r);

template <unsigned int blockSize>
__global__ void maxMonoBlock5(float* g_idata, float* g_odata, unsigned int n, unsigned int id);



