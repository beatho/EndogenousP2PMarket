#pragma once
#include <device_launch_parameters.h>
#ifdef OSQP
#include "MethodP2P.h"
#include "MatrixGPU.cuh"
#include "MatrixCPU.h"
#include "kernelFunction.cuh"

//#include "C:\Program Files\OSQP\osqp\include\osqp.h"

	#include <osqp.h>

#include <iostream>
#include <string>
#include <cuda_runtime.h>
#include <chrono>

/// <summary>
/// But : utiliser une autre m�thode pour g�rer les contraintes : avoir un SO qui r�soud un DC-OPF pour respecter les contraintes
/// Change : Init car pas m�me probl�me local, et Fb3a et b (bref la partie r�seau)
/// Ici la resolution se fait avec un OSQP (CPU)
/// </summary>
class ADMMGPUConstCons : public MethodP2P
{
public:
	ADMMGPUConstCons();
	ADMMGPUConstCons(float rho);
	virtual ~ADMMGPUConstCons();
	void setParam(float rho);
	void setTau(float tau);

	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	void updateGlobalProbGPU();
	void updateLocalProbGPU(float epsL, int nIterL);
	void init(const Simparam& sim, const StudyCase& cas);
	void updateP0(const StudyCase& cas);
	void solveOPF();
	void updateQ();
	std::string NAME ="ADMMGPUConstCons";
	
	virtual float updateRes(MatrixCPU* res, MatrixGPU* Tlocal, int iter, MatrixGPU* tempNN);
	
	
	void display();

private:
	// ne change pas avec P0
	float _rho = 0;
	float _tau = 1.1;
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
	MatrixGPU Ap2; // a * Mn^2
	MatrixGPU Ap1; // rhol * Mn
	MatrixGPU Ap123; // Ap1 + Ap2 + Ap3
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
	
	// Pour le consensus
	MatrixGPU etaP;
	MatrixGPU Ap3; // rho1*Mn^2
	MatrixGPU Bp3; // 1/Mn * (Pso + P)/2 - eta/rho1
	MatrixGPU Pso;

		

	// pour osqp

	OSQPSettings* settings;
	OSQPWorkspace* work;
#ifndef OSQPGPU
	OSQPData* data;
#endif
	
	MatrixCPU Aosqp; // l <A*P<u
	MatrixCPU Hosqp;
	MatrixCPU qTosqp;
	MatrixCPU PsoCPU;
	c_float* Q = nullptr;
	c_float* xResult = nullptr;

	// Pour le r�seau
	int _nLine;
	int _nBus;
	float _rho1;
	MatrixGPU G;
	MatrixGPU GTrans;
	MatrixGPU lLimit;

};

/*__device__ float warpReduceMax5(volatile float* r);

template <unsigned int blockSize>
__global__ void maxMonoBlock5(float* g_idata, float* g_odata, unsigned int n, unsigned int id);*/


#endif

