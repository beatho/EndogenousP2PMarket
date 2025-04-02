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

/// <summary>
/// But : utiliser une autre m�thode pour g�rer les contraintes : avoir un SO qui r�soud un DC-OPF pour respecter les contraintes
/// Change : Init car pas m�me probl�me local, et Fb3a et b (bref la partie r�seau)
/// Ici on utilise la m�thode du dual ascent (CPU)
/// </summary>
class ADMMGPUConstCons2 : public MethodP2PGPU
{
public:
	ADMMGPUConstCons2();
	ADMMGPUConstCons2(float rho);
	virtual ~ADMMGPUConstCons2();
	virtual void setParam(float rho);
	void setTau(float tau);

	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	void updateGlobalProbGPU();
	void updateLocalProbGPU(float epsL, int nIterL);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	void solveOPF();
	void findalpha();
	std::string NAME ="ADMMGPUConstCons2";
	
	virtual float updateResEndo(int iter);
	


private:

	MatrixGPU Ap123; // Ap1 + Ap2 + Ap3
	
		
	// Pour le consensus
	MatrixGPU etaP;
	MatrixGPU Ap3; // rho1*Mn^2
	MatrixGPU Bp3; // 1/Mn * (Pso + P)/2 - eta/rho1
	MatrixGPU Pso;




	// dual ascent version eva
	float _epsi;
	int L2; // 2*_nLine;
	int _Msize; // nAgent+2*nLine +1;
	int _Asize; // L2*nAgent
	float alpha; // a mettre sur GPU !
	float err1; // a mettre sur GPU !
	float v; //  Ploss = 0 a mettre sur GPU !
	float mu; // relaxation, fonction barriere

	MatrixGPU H;
	MatrixGPU q; // 0.5x^THx + q^T*x
	MatrixGPU diffPso;

	MatrixGPU c; // contrainte Ax+b>0 ou = 0 pour egalit�
	MatrixGPU Ai;
	MatrixGPU bi;

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





