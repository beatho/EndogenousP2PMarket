#pragma once
#include "MethodOPF.h"
#include <iostream>
#include <string>
#include <chrono>
#include "Utilities.cuh"
#include "Utilities.h"

#include "CPUPF.h"
#include "CPUPFdist.h"
#include "CPUPFdistPQ.h"
#include "ADMMMarket.h"


class OPFPDIPM : public MethodOPF
{
public:
	OPFPDIPM();
	OPFPDIPM(float rho);
	virtual ~OPFPDIPM();
	void setParam(float rho);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);

	virtual void solveConsensus(float eps, MatrixCPU* PSO);
	virtual void initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO);
	virtual void updateConsensus(MatrixCPU* Pmarket);

	std::string NAME ="OPFPDIPM";
	
	
	virtual float getPLoss();
	virtual float getQLoss();

	void updatePerturbedFactor();
	void correctionEquation();
	void updateStep();
	void updateVariable();
	bool setCenteringParameter();

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
	int _nBus = 0;
	int _nLine = 0;
	int _nAgent = 0;
	int _sizeEq = 0;   // R
	int _sizeInEq = 0; // S
	int _sizeVar = 0;  // M
	float _sigma = 0;
	float _mu = 0;
	float _V0 = 0;
	
	int _iterGlobal = 0;
	int _iterBest = 0;
	int _iterG = 0;
	int _stepG = 0;
	float _stepD = 0;
	float _stepP = 0;
	float Cgap = 0;
	float CgapBest = 0;
	float epsMin = 0.0000001;


	bool _initWithMarket = false;
	bool _initWithPF = false;
	bool _initFlatPower = false;
	bool _saveInterSol = true;
	float valMin = 0.01;
	float valMax = 10;

	clock_t timeOPF = 0;
	
	
	MatrixCPU X; // (pn, qn, ei, fi)
	MatrixCPU Xbest; // (pn, qn, ei, fi)
	MatrixCPU Xpre;
	MatrixCPU Y; // PF : real for all bus, imag for all bus
	MatrixCPU u;
	MatrixCPU l;
	MatrixCPU z;
	MatrixCPU w;

	MatrixCPU dX; // (pn, qn, ei, fi)
	MatrixCPU dY; // PF : real for all bus, imag for all bus
	MatrixCPU du;
	MatrixCPU dl;
	MatrixCPU dz;
	MatrixCPU dw;
	MatrixCPU dXY;

	MatrixCPU PSI;
	MatrixCPU MatSys; // Matsys*dXY = PSI
	MatrixCPU d; // U^-1 * W - L^-1 * Z
	MatrixCPU VectSys;

	MatrixCPU XHess;
	MatrixCPU hHess; 
	MatrixCPU gHess;
	MatrixCPU Hess;

	MatrixCPU XJac;
	MatrixCPU hJac;
	MatrixCPU gJac;

	MatrixCPU Lx;
	MatrixCPU Ly;
	MatrixCPU Lz;
	MatrixCPU Lw;
	MatrixCPU Ll;
	MatrixCPU Lu;

	//MatrixCPU tempN1; // Matrix temporaire pour aider les calculs
	//MatrixCPU tempNN; // plut�t que de re-allouer de la m�moire � chaque utilisation
	//MatrixCPU tempM1; //
	MatrixCPU tempM1; //
	MatrixCPU tempS1;
	MatrixCPU tempS1bis;
	MatrixCPU tempN2;
	//MatrixCPU tempR1;

	
	MatrixCPU _CoresBusAgent;
	MatrixCPU _nAgentByBus;


	MatrixCPU BgridLin;
	MatrixCPU GgridLin;
	MatrixCPU CoresVoiLin;
	MatrixCPU CoresBusLin;
	MatrixCPU nLines;
	
	MatrixCPU UpperBound; //
	MatrixCPU LowerBound; //
	MatrixCPU Pmin;
	MatrixCPU Pmax;
	MatrixCPU Pn;

	MatrixCPU Cost1;
	MatrixCPU Cost2;
	
	MatrixCPU PnBus;
	MatrixCPU CoresLineBus;

	MatrixCPU resF;

	
	MatrixCPU E; // voltage : [angle, norme]
	Simparam paramMarketInit;


	MatrixCPU etaSO; // consensus with market
	float _rhoSO = 0;

};


	




