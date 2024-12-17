
#pragma once
#include "MethodP2P.cuh"
#include <iostream>
#include <string>

#include <osqp.h>
#include <chrono>
#ifdef OSQPGPU
	#include <csc_utils.h>
#endif // OSQPGPU


class OSQPConst : public MethodP2P
{ // method ou P fait parti de l'inconnue x
public:
	OSQPConst();
	OSQPConst(float alpha);
	virtual ~OSQPConst();
	void setParam(float alpha);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	virtual void updatePn();
	virtual void updateLAMBDA();
	float updateResBis(MatrixCPU* res, int iter);
	void updatePhi();
	void updateTrade(c_float* xResult, int agent, MatrixCPU* omega);

	void updateQ();
	std::string NAME = "OSQPConst";
	void display();
private:
	float _alpha = 1.6;
	int _nAgent = 0;
	float _rho;

	MatrixCPU nVoisin;
	MatrixCPU LAMBDA;
	MatrixCPU trade;
	MatrixCPU a;
	MatrixCPU b;
	MatrixCPU Beta;
	MatrixCPU Tlocal;
	
	//MatrixCPU* OMEGA;

	MatrixCPU Ub;

	int** IdVoisin = nullptr;
	MatrixCPU PosVoisin; //Pos(n,m)=i; pos(m,n)=j with n the i^th neigbourh of m and m the j^th neigbourh of n

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

	// Pour le réseau
	int _nLine;
	int _nBus;
	float _rho1;

	MatrixCPU Kappa1;
	MatrixCPU Kappa2;
	MatrixCPU Kappa1_pre;
	MatrixCPU Kappa2_pre;
	MatrixCPU Phipart;
	MatrixCPU Phi;
	MatrixCPU Ggrid;
	MatrixCPU Ggrid2; //Ggrid*Ggrid
	
	MatrixCPU tempL1;
	MatrixCPU lLimit;
	MatrixCPU Pn;
};









