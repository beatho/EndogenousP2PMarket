#pragma once
#include "MethodP2P.cuh"
#include "MatrixGPU.cuh"
#include "MatrixCPU.h"
#include "kernelFunction.cuh"

#include <iostream>
#include <string>
#include <cuda_runtime.h>
#include <chrono>


class ADMMGPUConst1 : public MethodP2P
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
	
	virtual float updateRes(MatrixCPU* res, MatrixGPU* Tlocal, int iter, MatrixGPU* tempNN);
	
	float calcRes(MatrixGPU* Tlocal, MatrixGPU* P);
	
	void display();

private:
	// ne change pas avec P0
	float _mu = 40.0f;
	float _mu1 = 40.0f;
	float _mul = 40.0f;
	float _tau = 2.0f;
	
	float _rho = 0.0f;
	float _rhol = 0.0f;
	int _blockSize = 512;
	int _numBlocksN = 0;
	int _numBlocksM = 0;
	int _numBlocksL = 0;
	int _numBlocksNL = 0;
	int _nAgent = 0;
	int _nTrade = 0;
	float _rhog = 0.0f;
	float _at1 = 0.0f;
	float _at2 = 0.0f;

	MatrixGPU tempNN; // Matrix temporaire pour aider les calculs
	MatrixGPU tempN1; // plut�t que de re-allouer de la m�moire � chaque utilisation
	MatrixGPU tempL1; 
	MatrixGPU tempL2; 
	

	MatrixGPU Tlocal;
	MatrixGPU P; // moyenne des trades
	MatrixGPU Pn; // somme des trades

	MatrixGPU a;
	MatrixGPU Ap2;
	MatrixGPU Ap2a; // stocke la premi�re partie de Ap2a qui ne change par avec rho1
	MatrixGPU Ap1;
	MatrixGPU Ap1Copy;
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
	int _nLine = 0;
	int _nBus = 0;
	float _rho1 = 0.0f;

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


