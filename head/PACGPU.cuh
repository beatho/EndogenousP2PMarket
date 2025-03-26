#pragma once
#include <device_launch_parameters.h>
#include "MethodP2P.cuh"
#include "kernelFunction.cuh"
#include <iostream>
#include <string>
#include <chrono>


class PACGPU : public MethodP2P // g�re P et Q comme si tous les Q �taient des agents diff�rents
	//et donc P et Q sont 2 probl�mes compl�tement distincts
{
public:
	PACGPU();
	PACGPU(float rho);
	virtual ~PACGPU();
	void setParam(float rho);
	void setGamma(float gamma); 
	void setGammahat(float gammahat);
	void setInitCoef(float alpha, float phi, float theta);
	void setBestRhoGamma(float lambdaMax, float lambdaMin, const StudyCase& cas);
	void setBestRhoGammaHeuristic(const StudyCase& cas);
	void updateCoef();
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	virtual void setBestParam(const StudyCase& cas);
	std::string NAME ="PACGPU";
	void updateGlobalProb();
	void updateLocalProb();

	void updateXhat();
	void updateMu();
	void updateNu();
	void updateQ();

	float updateRes(int indice);

	void display();
private:
	bool augmente = true;
	// ne change pas avec P0
	float _gamma = 0.1;
	float _gammahat = 0.1;
	bool isAC = false;
	int _blockSize = 512;
	int _numBlocksN = 0;
		
	int _nAgent = 0;
	int _nAgentTrue = 0;
	int _nTrade = 0;
	int _nTradeP = 0;
	int _nTradeQ = 0;
	float _rhog = 0;
	float _rho = 1.5;
	float _rhoInv = 0;
	int _sizePACGPU = 0;
	int _sizePACX = 0;
	int _sizePACMu = 0;
	int _sizePACNu = 0;
	int _sizeHinv = 0;
	
	// parameter agent and iteration dependant (but not here for now)
	float _alpha = 0; // 0.6
	float _phi = 0; // 0.2
	float _theta = 0; // 0.2
	

	MatrixGPU tempN1; // Matrix temporaire pour aider les calculs
	MatrixGPU tempNN; // plut�t que de re-allouer de la m�moire � chaque utilisation
	MatrixGPU tempM1; //
	MatrixGPU tempM; //
	

	MatrixGPU Cost1;
	MatrixGPU Cost2;
	MatrixGPU BETA;
	MatrixGPU BETALin;

	MatrixGPU X;
	MatrixGPU Xpre;
	MatrixGPU Xhat;
	MatrixGPU Mu;
	MatrixGPU Muhat;
	MatrixGPU Nu;
	MatrixGPU Nuhat;


	MatrixGPU Hinv;
	MatrixCPU H;
	MatrixGPU _sizeQ;

	MatrixGPU nVoisin;
	MatrixGPU CoresMatLin;
	MatrixGPU CoresLinAgent;
	MatrixGPU CoresAgentLin;
	MatrixGPU CoresLinVoisin;
	MatrixGPU CoresLinTrans;
	MatrixGPU CoresLinTransLocal;
	MatrixGPU CoresIndiceNu;
	MatrixGPU CoresAgentLinBig;

	// Matrices for the result
	MatrixCPU Pn;
	MatrixCPU trade;
	MatrixCPU resF;
	MatrixCPU resX;
	MatrixCPU nVoisinCPU;



	// change avec P0
	
	MatrixGPU matLb;
	MatrixGPU Q;
	MatrixGPU Qinit;

	MatrixGPU matUb;

};


	
__global__ void updateNuAug(float* Nu, float* Nuhat, float* Xhat, float* nVoisin, float* CoresAgentLin, float* CoresLinTrans, float* CoresindiceNu, float rho, float gamma, float theta);
__global__ void updateNuGPU(float* Nu, float* Nuhat, float* Xhat, float* nVoisin, float* CoresAgentLin, float* CoresLinTrans, float* CoresindiceNu, float rho, float gamma, float gammahat);
__global__ void updateQAug(float* Q, float* Qinit, float* Xhat, float* Muhat, float* Nuhat, float* nVoisin, float* CoresAgentLin, float* CoresLinTrans, float* CoresIndiceNu, float rhoInv, float rho, float gamma, bool augmente);

template <unsigned int blockSize>
__global__ void updateLocalProblPAC(float* X, float* Q, float* Hinv, float* matLb, float* matUb, float* CoresAgentLinBig, float* sizeQ, int sizeOPFmax);


__global__ void calculFcPAC(float* tempN, float* tempNN, float* a, float* b, float* betaLin, float* X, float* CoresAgentLin, float* nVoisin);
__global__ void updateP0PAC(float* matLb, float* matUb, float* Q, float* Qinit, float* Pmin, float* Pmax, float* Cost2, float* Lb, float* Ub, float* CoresAgentLin, float* nVoisin);