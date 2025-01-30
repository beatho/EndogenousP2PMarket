
#pragma once
#include "Method.h"
#include <iostream>
#include <string>

#ifdef OSQP
	#include <osqp.h>



#include <chrono>
#ifdef OSQPGPU
	#include <csc_utils.h>
	#include <cs.h>
#endif // OSQPGPU


class OSQPCentralized2 : public Method
{ // methode centralis� ou P ne fait pas partie des inconnus, car sinon on a plein de ligne nulle dans H
	// Donc q = bn + gamma	
	// Et H = a_n par bloc 

public:
	OSQPCentralized2();
	OSQPCentralized2(float alpha);
	virtual ~OSQPCentralized2();
	void setParam(float alpha);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	void updatePn(MatrixCPU* Pn, MatrixCPU* trade);
	virtual float calcFc(MatrixCPU* cost1, MatrixCPU* cost2, MatrixCPU* trade, MatrixCPU* Pn, MatrixCPU* BETA, MatrixCPU* tempN1, MatrixCPU* tempNN);
	void updateTrade();
	std::string NAME = "OSQPCentralized2";
	void display();
private:
	float _alpha;
	int _nAgent;
	int _nAgentTrue = 0; // _nAgent = _nAgentTrue + (isAc)*_nAgent
	
	int _nTrade = 0;
	int _nTradeP = 0;
	int _nTradeQ = 0;

	bool isAC = false;

	int nConstraint;

	MatrixCPU nVoisin;
	MatrixCPU trade;
	MatrixCPU a;
	MatrixCPU b;
	MatrixCPU GAMMA;
	
	MatrixGPU CoresMatLin;
	MatrixGPU CoresLinAgent;
	MatrixGPU CoresAgentLin;
	MatrixGPU CoresLinVoisin;
	MatrixGPU CoresLinTrans;
	MatrixCPU GAMMALin;
	//MatrixCPU* OMEGA;

	MatrixCPU Ub;
	
	c_float* Q = nullptr;
	c_float* u = nullptr;
	c_float* l = nullptr;
	c_float* xResult = nullptr;
	c_float* Pdata = nullptr;
	c_int* Pidx = nullptr;
	c_int* Pptr = nullptr;
	c_float* Adata = nullptr;
	c_int* Aidx = nullptr;
	c_int* Aptr = nullptr;

	OSQPSettings* settings = nullptr;
#ifdef OSQPGPU
	OSQPSolver* work;
	csc* Pcsv = new csc();
	csc* Acsv = new csc();
	csc* Pu = new csc();
#else
	OSQPWorkspace* work = nullptr;
	OSQPData* data = nullptr;
#endif // OSQPGPU
};


#endif






