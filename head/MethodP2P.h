#pragma once


#include "Utilities.h"

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
	virtual void init(const Simparam& sim, const StudyCase& cas) = 0;
	virtual void updateP0(const StudyCase& cas);
	//virtual void setTau(float tau) = 0;

	virtual void updateLAMBDA();
	void updateLambda();
	void updateKappa();
	void updateCp2();

	float updateResMat(int iter);
	float updateRes(int iter);
	float updateResEndo(int iter);
	float calcRes();

	void initLinForm(const StudyCase& cas);
	void initSize(const StudyCase& cas);
	void initSimParam(const Simparam& sim);
	void initDCEndoGrid(const StudyCase& cas); // init with G transposed
	void initCaseParam(const Simparam& sim, const StudyCase& cas);
	void initDCEndoMarket();
	void initP2PMarket();

	void setResult(Simparam* result, bool casAC);
	// compute : fc = 0.5 * cost1 * pn^2 + cost2 * pn + BETA * trade
	float calcFc();

	virtual void updatePn();
	virtual void solveWithMinPower(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void display();
protected:
	// ne change pas avec P0
	clock_t tMarket;

	float _mu = 50;
	float _mu1 = 50;
	float _mul = 50;
	float _tau = 2;

	float _rho = 0;
	float _rhol = 0;
	float _rhog = 0;

	int _iterGlobal = 0;
	int _iterG = 0;
	int _iterL = 0;
	int _iterIntern = 0;

	int _stepG = 1;
	int _stepL = 1;
	int _stepIntern = 1;

	float _epsL = 0;
	float _epsG = 0;
	float _epsX = 0;
	float _epsIntern = 0;
	
	bool isAC = false;
	int _nAgent = 0;
	int _nAgentTrue = 0;

	int _nTrade = 0;
	int _nTradeP = 0;
	int _nTradeQ = 0;
	
	float _at1 = 0;
	float _at2 = 0;

	MatrixCPU tempNN; // Matrix temporaire pour aider les calculs
	MatrixCPU tempN1; // plut�t que de re-allouer de la m�moire � chaque utilisation
	

	MatrixCPU Tlocal;
	MatrixCPU P; // moyenne des trades
	MatrixCPU Pn; // somme des trades

	MatrixCPU a;
	MatrixCPU Ap1;
	MatrixCPU Ap2;
	MatrixCPU Ap2a;
	MatrixCPU Ap2b;
	MatrixCPU Ap3;
	MatrixCPU Ap12;
	MatrixCPU Ap123;


	MatrixCPU Bt1;
	MatrixCPU Bt2;
	MatrixCPU Bp1;
	MatrixCPU Ct;
	MatrixCPU matUb;

	MatrixCPU nVoisin;
	MatrixCPU tradeLin;
	MatrixCPU Tmoy;
	MatrixCPU LAMBDALin;
	MatrixCPU Tlocal_pre;
	MatrixCPU MU;

	MatrixCPU CoresMatLin;
	MatrixCPU CoresLinAgent;
	MatrixCPU CoresAgentLin;
	MatrixCPU CoresLinVoisin;
	MatrixCPU CoresLinTrans;

	// Matrices for the result
	MatrixCPU LAMBDA;
	MatrixCPU trade;
	MatrixCPU resF;
	MatrixCPU resX;



	// change avec P0
	MatrixCPU b;
	MatrixCPU matLb;
	MatrixCPU Cp;
	MatrixCPU Pmin;
	MatrixCPU Pmax;

	// Pour le reseau
	int _nLine;
	int _nBus;
	float _rho1;

	MatrixCPU tempL1;
	MatrixCPU tempL2;

	MatrixCPU Kappa1;
	MatrixCPU Kappa2;
	MatrixCPU Kappa1_pre;
	MatrixCPU Kappa2_pre;
	MatrixCPU Cp1;
	MatrixCPU Cp2;
	MatrixCPU Qpart;
	MatrixCPU Qtot;
	MatrixCPU alpha;


	MatrixCPU G;
	MatrixCPU GTrans;
	MatrixCPU G2; // G^t.*G^t
	MatrixCPU lLimit;

};



