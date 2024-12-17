#pragma once
#include "MethodP2P.cuh"
#include "MatrixGPU.cuh"
#include "MatrixCPU.h"
#include "kernelFunction.cuh"


#include <iostream>
#include <string>
#include <cuda_runtime.h>
#include <chrono>

/// <summary>
/// But : utiliser une autre méthode pour gérer les contraintes : avoir un SO qui résoud un DC-OPF pour respecter les contraintes
/// Change : Init car pas même problème local, et Fb3a et b (bref la partie réseau)
/// Utilise la methode du dual ascent GPU
/// </summary>
class ADMMGPUConstCons3 : public MethodP2P
{
public:
	ADMMGPUConstCons3();
	ADMMGPUConstCons3(float rho);
	virtual ~ADMMGPUConstCons3();
	void setParam(float rho);
	void setTau(float tau);

	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	void updateGlobalProbGPU(float epsOPF);
	void updateLocalProbGPU(float epsL, int nIterL);
	void init(const Simparam& sim, const StudyCase& cas);
	void updateP0(const StudyCase& cas);
	void solveOPF(float epsOPF);
	void findalpha();
	std::string NAME ="ADMMGPUConstCons3";
	
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
	MatrixGPU tempN1; // plutôt que de re-allouer de la mémoire à chaque utilisation
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


	// Pour le réseau
	int _nLine;
	int _nBus;
	float _rho1;
	MatrixGPU G;
	MatrixGPU GTrans;
	MatrixGPU lLimit;

	// dual ascent version eva
	float _epsi;
	int L2; // 2*_nLine;
	int _Msize; // nAgent+2*nLine +1;
	int _Asize; // L2*nAgent
	float* alpha = nullptr; // a mettre sur GPU !
	float err1; // 

	float mu; // relaxation, fonction barriere

	MatrixGPU tempL21; // taille 2*L+1 vecteur colonne.

	MatrixGPU H;
	MatrixGPU q; // 0.5x^THx + q^T*x
	MatrixGPU diffPso;

	MatrixGPU c; // contrainte Ax+b>0 ou = 0 pour egalité
	MatrixGPU Ai;
	MatrixGPU bi;
	MatrixGPU Apas; // Ai* pas

	MatrixGPU M; // M = (H -Atrans ZA W)
	MatrixGPU pas;// M*pas = R
	MatrixGPU Minv; // Minv = M^-1 

	MatrixGPU ZA; 
	MatrixGPU Z;// Z = diag(Zvect)
	MatrixGPU Zvect; 
	MatrixGPU W; // W = diag(Wvect)
	MatrixGPU Wvect; // c inegalites, mu egalite
	MatrixGPU Atrans;
	
	
	
	MatrixGPU R; // R = (Rx1+Rx2 Ru)
	MatrixGPU Rx1; // Hx+q
	MatrixGPU Rx2; // Ai*U
	
	MatrixGPU Ru; // Ru = W*(U-PI)
	
	MatrixGPU U;
	MatrixGPU PI;

};





