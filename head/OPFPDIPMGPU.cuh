#pragma once
#include "Method.h"
#include <iostream>
#include <string>
#include <chrono>
#include "Utilities.cuh"
#include "CPUPF.h"
#include "CPUPFdist.h"
#include "ADMMMarket.h"


class OPFPDIPMGPU : public Method
{
public:
	OPFPDIPMGPU();
	OPFPDIPMGPU(float rho);
	virtual ~OPFPDIPMGPU();
	void setParam(float rho);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="OPFPDIPMGPU";
	
	
	float getPLoss();
	float getQLoss();

	void updatePerturbedFactor();
	void correctionEquation();
	void updateStep();
	void updateVariable();
	void setCenteringParameter();

	void updateJ();
	void updateH();
	void setSystem();
	void computeConstraint(); // L. & PSI(.,0)
	void updatePSILlu();
	
	
	float updateRes(int indice);

	void display();
private:
	void initWithMarket(const Simparam& sim, const StudyCase& cas);
	void initWithPF(const StudyCase& cas);

	// ne change pas avec P0
	int _blockSize = 512;
	int _blockSizeSmall = 64;
	int _numBlocksB = 0;
	int _numBlocksN = 0;
	int _numBlocksM = 0;
	int _nBus = 0;
	int _nLine = 0;
	int _nAgent = 0;
	int _sizeEq = 0;   // R
	int _sizeInEq = 0; // S
	int _sizeVar = 0;  // M
	float* _sigma;
	float _mu = 0;
	float _V0 = 0;
	
	int _iterGlobal = 0;
	int _iterG = 0;
	int _stepG = 0;
	float* _stepGPU = nullptr;
	float Cgap = 0;
	float epsMin = 0.0000001;


	bool _initWithMarket = true;
	bool _initWithPF = true;
	float valMin = 0.01;
	float valMax = 10;

	clock_t timeOPF = 0;
	
	
	MatrixGPU X; // (pn, qn, ei, fi)
	MatrixGPU Xpre;
	MatrixGPU Y; // PF : real for all bus, imag for all bus
	MatrixGPU u;
	MatrixGPU l;
	MatrixGPU z;
	MatrixGPU w;

	MatrixGPU dX; // (pn, qn, ei, fi)
	MatrixGPU dY; // PF : real for all bus, imag for all bus
	MatrixGPU du;
	MatrixGPU dl;
	MatrixGPU dz;
	MatrixGPU dw;
	MatrixGPU dXY;

	MatrixGPU PSI;
	MatrixGPU MatSys; // Matsys*dXY = PSI
	MatrixGPU d; // U^-1 * W - L^-1 * Z
	MatrixGPU VectSys;

	MatrixGPU XHess;
	MatrixGPU hHess; 
	MatrixGPU gHess;
	MatrixGPU Hess;

	MatrixGPU XJac;
	MatrixGPU hJac;
	MatrixGPU gJac;

	MatrixGPU Lx;
	MatrixGPU Ly;
	MatrixGPU Lz;
	MatrixGPU Lw;
	MatrixGPU Ll;
	MatrixGPU Lu;

	//MatrixGPU tempN1; // Matrix temporaire pour aider les calculs
	//MatrixGPU tempNN; // plutôt que de re-allouer de la mémoire à chaque utilisation
	//MatrixGPU tempM1; //
	MatrixGPU tempM1; //
	MatrixGPU tempS1;
	MatrixGPU tempS1bis;
	MatrixGPU tempN2;
	//MatrixGPU tempR1;

	
	MatrixGPU _CoresBusAgent;
	MatrixGPU _nAgentByBus;


	MatrixGPU BgridLin;
	MatrixGPU GgridLin;
	MatrixGPU B2; // only for the l lines
	MatrixGPU G2;// only for the l lines
	MatrixGPU CoresVoiLin;
	MatrixGPU CoresBusLin;
	MatrixGPU nLines;
	MatrixGPU nLinesBegin;
	
	MatrixGPU UpperBound; //
	MatrixGPU LowerBound; //
	MatrixGPU Pmin;
	MatrixGPU Pmax;
	MatrixGPU Pn;

	MatrixGPU Cost1;
	MatrixGPU Cost2;
	
	MatrixGPU PnBus;
	MatrixGPU CoresLineBus;

	// for the system
	MatrixGPU AD;
	MatrixGPU PD;


	// Matrix that stay on CPU

	MatrixCPU resF;
	MatrixCPU _nAgentByBusCPU;

	
	MatrixGPU E; // voltage : [angle, norme]
	Simparam paramMarketInit;

};

__global__ void initBoundGPU(float* lowerBound, float* upperBound, float* lowerBoundBusLine, float* upperBoundBusLine, float* Pmin, float* Pmax, int nAgent, int nBus, int nLine);

__global__ void initFlatVoltageGPU(float* voltage, int nBus);

__global__ void initForAgent(float* X, float* XHess, float* l, float* u, float* Pn, float* Cost1, float* lowerBound, float* upperBound, float valMin, float valMax, int nAgent, int sizeVar);
__global__ void initForBus(float* X, float* l, float* u, float* E, float* lowerBound, float* upperBound, float valMin, float valMax, int nBus, int nAgent);
__global__ void initForLine(float* l, float* u, float* E, float* CoresLineBus, float* lowerBound, float* upperBound, float* B2, float* G2, float valMin, float valMax, int nBus, int nAgent, int nLine);


template <unsigned int blockSize>
__global__ void updateStepGPU(float* step, float* l, float* u, float* z, float* w, float* dl, float* du, float* dz, float* dw, int sizeInEq);
__global__ void updatedXX(float* X, float* dX, float* Xpre, float* stepGPU, int sizeVar);
__global__ void updatedYY(float* Y, float* dY, float* stepGPU, int sizeEq);

__global__ void updateSlack(float* l, float* u, float* w, float* z, float* dl, float* du, float* dw, float* dz, float* stepGPU, int sizeInEq);


__global__ void updateJacAgent(float* XJac, float* gJac, float* hJac, float* X, float* Cost1, float* Cost2, int nAgent);

template <unsigned int blockSize>
__global__ void updateCenteringParameter(float* _sigma, float* l, float* u, float* w, float* z, float* dl, float* du, float* dw, float* dz, float* stepGPU, float Cgap, int sizeInEq);

