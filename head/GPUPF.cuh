#pragma once

#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}


#include <stdio.h>
#include <stdlib.h>
#include "MatrixGPU.cuh"
#include "MatrixGPUD.cuh"
#include "StudyCase.h"
#include <chrono>
#include "kernelFunction.cuh"
#include "Utilities.cuh"
#include "Utilities.h"

class GPUPF
{
	
public:
	GPUPF();
	~GPUPF();
	virtual void init(const StudyCase& cas, MatrixGPU* Pn);
	virtual void init(const StudyCase& cas, MatrixGPU* Pn, MatrixGPUD* PnD, bool useDouble);
	virtual void solve();
	virtual void updatePQ(MatrixGPU* PQ);
	virtual void calculW0(MatrixGPU* PQ);
	virtual void calculW0D(MatrixGPUD* PQD);
	virtual void calculW0Bis(MatrixGPU* PQ);
	virtual void calculW0DBis(MatrixGPUD* PQD);
	virtual void calcW(bool end=false);
	virtual void calcJac();
	virtual int calcVoltage();
	
	
	
	virtual void calcPhi();
	virtual void calcJacPhiE();
	virtual void calcE();
	MatrixGPU* calcG();
	virtual MatrixGPU  getY();


	virtual void setE(MatrixGPU* Enew);
	virtual void setE(MatrixGPUD* Enew);
	
	void setW(MatrixGPU* Wnew);

	float getPloss();
	float getQloss();
	float getRes();
	int getIter();
	float getP0();
	float getQ0();
	float getTime();
	int getConv();
	MatrixCPU getE();
	MatrixCPU getW();

	virtual void display(bool all);

	void saveTimeBlock(std::string fileName);

protected:

	bool _useDouble = false;
	int _blockSize = 32; // petit car chaque noeud a peu de voisin, � voir...
	int numBlock = 0;
	int Nagent = 0;
	int Nbus = 0;
	int Nline = 0;
	int Nconstraint = 0;
	int B2 = 0; // 2*Nbus =0
	int N2 = 0; // 2*NAgent =0
	int BL2 = 0; // Nbus+2*Nline
	double V0 = 1; // voltage at the ref bus
	double theta0 = 0; // voltage angle at the ref bus
	int iterM = 10;
	double epsPF = 0.0005;
	double err = 0;
	int iter = 0;
	clock_t time =0;
	std::string _name = "Newton GPU";
	int status;


	MatrixGPU I_aug; // agent to bus
	MatrixGPU I; // agent to bus

	MatrixGPU W0;
	MatrixGPU W;
	MatrixGPU dW;
	MatrixGPU E;
	
	MatrixGPU dE;
	MatrixGPU Jac;
	MatrixGPU JacInv;
	MatrixGPU Bgrid;
	MatrixGPU Ggrid;

	MatrixGPU G;
	MatrixGPU Phi;
	MatrixGPU Y; // [E Phi]

	MatrixGPU tempLN2; // JacPhiE * Jacinv * Iaug
	MatrixGPU tempB2N2; // Jacinv * Iaug
	MatrixGPU JacPhiE; // L * B2
	MatrixGPU CoresLineBus;
	MatrixGPU Ggrid2Bgrid2; // Ggrid^2+Bgrid^2

	MatrixGPU _Blin;
	MatrixGPU _Glin;
	MatrixGPU _CoresVoiLin;
	MatrixGPU _CoresBusLin;
	MatrixGPU _nLines;
	
	MatrixGPU _Pintermediate;
	MatrixGPU _Qintermediate;

	MatrixGPU A;
	MatrixGPU P;
	
	MatrixGPU CoresAgentBus;
	MatrixGPU CoresAgentBusBegin;
	MatrixGPU NagentByBus;

	// si double
	MatrixGPUD ED;
	MatrixGPUD dED;

	MatrixGPUD _BlinD;
	MatrixGPUD _GlinD;
	// les autres sont des entiers, pas besoin que cela soit des float ?

	MatrixGPUD JacD;
	MatrixGPUD JacInvD;

	MatrixGPUD W0D;
	MatrixGPUD WD;
	MatrixGPUD dWD;

	MatrixGPUD _PintermediateD;
	MatrixGPUD _QintermediateD;


	MatrixGPUD AD;
	MatrixGPUD PD;

	// Pour interactionmarche endog�ne
	MatrixGPU _Blin2; // que les lignes
	MatrixGPU _Glin2;
	MatrixGPU CoresLineBusGPU; // doit prendre la transpos�e 


	// instrumentation
	MatrixCPU timePerBlock = MatrixCPU(7, 1); // Fb0, Fb1abc, Fb2, Fb3abc, Fb4, Fb5, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	MatrixCPU occurencePerBlock = MatrixCPU(7, 1);; //nb de fois utilis� pendant la simu

	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;
};


template <unsigned int blockSize>
__global__ void calcW0GPU(float* W0, float* PQ, float* Cores, int N, int B);

template <unsigned int blockSize>
__global__ void calcW0GPUD(double* W0D, double* PQD, float* Cores, int N, int B);

template <unsigned int blockSize>
__global__ void calcW0GPUBis(float* W0, float* PQ, float* Cores, float* nAgentByBus, float* beginBus, int N, int B);

template <unsigned int blockSize>
__global__ void calcW0GPUDBis(double* W0D, double* PQD, float* Cores, float* nAgentByBus, float* beginBus, int N, int B);

__global__ void initED(double* ED, double theta0, double V0, int B);

__global__ void initE(float* E, float theta0, float V0, int B);

template <unsigned int blockSize>
__global__ void calcWGPU(float* W, float* Pinter, float* Qinter, float* CoresBusLin, float* nLines, int B);
template <unsigned int blockSize>
__global__ void calcWGPU(float* W, float* W0, float* Pinter, float* Qinter, float* CoresBusLin, float* nLines, int B);

__global__ void calcWinter(float* Pinter, float* Qinter, float* E, float* Glin, float* Blin, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B);

__global__ void calcJacGPU(float* Jac, float* W, float* E, float* Glin, float* Blin, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B);

template <unsigned int blockSize>
__global__ void calcWGPUD(double* W, double* W0, double* Pinter, double* Qinter, float* CoresBusLin, float* nLines, int B);

template <unsigned int blockSize>
__global__ void calcWGPUD(double* W, double* Pinter, double* Qinter, float* CoresBusLin, float* nLines, int B);

__global__ void calcWinterD(double* Pinter, double* Qinter, double* E, double* Glin, double* Blin, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B);

__global__ void calcJacGPUD(double* Jac, double* W, double* E, double* Glin, double* Blin, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B);