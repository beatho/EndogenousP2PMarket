#pragma once

#include "Utilities.cuh"
#include "Utilities.h"
#include <cuda_runtime.h>
#include "MatrixGPU.cuh"
#include "MethodP2P.cuh"
#include "StudyCase.h"
#include "MatrixCPU.h"

#include "Simparam.h"
#include "Method.h"


#include <iostream>
#include <string>

class MethodP2P : public Method
{
public:
	MethodP2P();
	virtual ~MethodP2P();
	virtual void setParam(float param) = 0;
	//virtual void setTau(float tau) = 0;

	void updateLAMBDA(MatrixCPU* LAMBDA, MatrixCPU* trade, float rho);
	void updateLAMBDA(MatrixGPU* LAMBDA, MatrixGPU* trade, float rho, MatrixGPU* tempNN);
	void updateKappa(MatrixCPU* Kappa1, MatrixCPU* Kappa2, MatrixCPU* L, MatrixCPU* Qtot);
	void updateKappa(MatrixGPU* Kappa1, MatrixGPU* Kappa2, MatrixGPU* L, MatrixGPU* Qtot);
	void updateCp2(MatrixCPU* Cp2, float rho1, MatrixCPU* Kappa1, MatrixCPU* Kappa2, MatrixCPU* G, MatrixCPU* temp1L, MatrixCPU* Qpart, MatrixCPU* nVoisin, int nLine, int nAgent);
	virtual float updateRes(MatrixCPU* res, MatrixGPU* Tlocal, MatrixGPU* trade, int iter, MatrixGPU* tempNN);
	
	float updateRes(MatrixCPU* res, MatrixCPU* Tlocal, MatrixCPU* trade, int iter);
	float updateRes(MatrixCPU* res, MatrixCPU* Tlocal, MatrixCPU* trade, int iter, MatrixCPU* Kappa1, MatrixCPU* Kappa2, MatrixCPU* Kappa1_pre, MatrixCPU* Kappa2_pre);
	
	virtual void updatePn(MatrixCPU* Pn, MatrixCPU* Tmoy, MatrixCPU* nVoisin);
	void updatePn(MatrixGPU* Pn, MatrixGPU* Tmoy, MatrixGPU* nVoisin);
	void updatePn(MatrixCPU* Pn, MatrixCPU* trade);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas)=0;
	virtual void solveWithMinPower(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas) = 0;
	virtual void init(const Simparam& sim, const StudyCase& cas) = 0;
	virtual void display() = 0;
	
	void resetId();

	

};





// Regroupement des appels kernels pour ne pas devoir d�finir plusieurs fois le m�me avec des noms differents
__global__ void updateLAMBDAGPU(float* LAMBDALin, float* tradeLin, float rho, float* CoresLinTrans, int const N);
__global__ void updateBt1GPU(float* Bt1, float* tradeLin, float rho, float* LAMBDA, float* CoresLinTrans, int const N);
__global__ void updateLAMBDABt1GPU(float* Bt1, float* LAMBDALin, float* tradeLin, float rho, float* CoresLinTrans, int const N);



__global__ void updateDiffGPU(float* tempNN, float* Tlocal, float* CoresLinTrans, int const N);
__global__ void updateResKappa(float* result, float* Kappa1, float* Kappa2, float* Kappapre1, float* Kappapre2, float ratio, int const L);
__global__ void selectResidual(float* res, unsigned int id1, unsigned int id2, unsigned int id3, float* output);
__global__ void updateKappaGPU(float* Kappa1, float* Kappa2, float* Llimit, float* Qtot, int nLine);
__global__ void diffKappa(float* tempL1, float* Kappa1, float* Kappa2, const int nLine);
__global__ void updateCpOld(float* Cp, float* Cp1, float* Cp2, float* tempN1, float* nVoisin, const float rho1, const int nAgent);
__global__ void updateCp(float* Cp, float* Cp1, float* Cp2, const int nAgent);

__global__ void updateQpart(float* Qpart, float* alpha, const int nAgent);
__global__ void updateQpartTrans(float* Qpart, float* alpha, const int N, const int nLine);

__global__ void updateQtot(float* Qtot, float* Qpart, float* alpha, const int nLine, const int nAgent);
__global__ void updateQtotTrans(float* Qtot, float* Qpart, float* alpha, const int nLine);


__global__ void updatePnGPU(float* Pn, float* Tmoy, float* nVoisin, const int nAgent);

__global__ void updateAlpha(float* alpha, float* G, float* Pn, const int nLine, const int nAgent);
__global__ void updateAlphaTrans(float* alpha, float* GTrans, float* Pn, const int nLine, const int nAgent);

__global__ void updateResX(float* res, float* Kappa1, float* Kappa2, float* KappaPre1, float* KappaPre2, const int nLine);




// Methode pour le dual ascent


/*
__global__ void updateQt(float* qt, float* Pso, float* Pn, float* etaSO, float rho1, int N);
__global__ void updateUAiq(float* UAiq, float* u, float* Aiq, int N, int size);
__global__ void updateRu(float* Ru, float* U, float* g, float epsi, int N, int L2);
__global__ void updateV(float* v, float* pas, float* alpha, int offset);*/


__global__ void updatePI(float* PI, float* c, float mu, float valMin, int L);
__global__ void updatePso(float* Pso, float* pas, float* alpha, int N);
__global__ void updateU(float* U, float* pas, float* alpha, int N, int L);

 
__global__ void updateEtaPBp3(float* Bp3, float* etaP, float* nVoisin, float* Pso, float* Pn, float rho, const int nAgent);
