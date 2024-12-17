
#pragma once
#include "MethodP2P.cuh"
#include <iostream>
#include <string>

#include <osqp.h>
#include <chrono>
#ifdef OSQPGPU
	//#include <csc_utils.h>
	#include <csc_type.h>
#endif // OSQPGPU


class OSQP : public MethodP2P
{ // method ou P fait parti de l'inconnue x
public:
	OSQP();
	OSQP(float alpha);
	virtual ~OSQP();
	void setParam(float alpha);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	void updateTrade(c_float* xResult, int agent);
	void updateQ(c_float* Q, int agent);
	void updateLAMBDA();
	float calcFc(MatrixCPU* cost1, MatrixCPU* cost2, MatrixCPU* trade, MatrixCPU* Pn, MatrixCPU* BETA, MatrixCPU* tempN1, MatrixCPU* tempNN);
	float updateResBis(MatrixCPU* res, int iter, MatrixCPU* tempNN);

	std::string NAME = "OSQP";
	void display();
private:
	float _alpha = 1.6;
	int _nAgent = 0;
	float _rho;
	int _nAgentTrue = 0; // _nAgent = _nAgentTrue + (isAc)*_nAgent
	int _nTrade = 0;
	int _nTradeP = 0;
	int _nTradeQ = 0;
	bool isAC = false;

	MatrixCPU nVoisin;
	MatrixCPU LAMBDA;
	MatrixCPU trade;
	MatrixCPU a;
	MatrixCPU b;
	MatrixCPU GAMMA;
	MatrixCPU Tlocal;
	
	// Faire la correspondance entre agent + num trade et voisin
	MatrixCPU CoresAgentLin;
	MatrixCPU CoresLinVoisin;


	//MatrixCPU* OMEGA;

	MatrixCPU Ub;

	c_float** Q = nullptr;
	c_float** u = nullptr;
	c_float** l = nullptr;
	c_float** xResult = nullptr;

	c_float** Pdata = nullptr;
	c_int** Pidx = nullptr;
	c_int** Pptr = nullptr;

	c_float** Adata = nullptr;
	c_int** Aidx = nullptr;
	c_int** Aptr = nullptr;

	OSQPSettings* settings = nullptr;
#ifdef OSQPGPU
	OSQPSolver** work = nullptr;
	csc** Pcsv = nullptr;
	csc** Acsv = nullptr;
	csc** Pu = nullptr;
#else
	OSQPWorkspace** work = nullptr;
	OSQPData** data = nullptr;
#endif // OSQPGPU
};









