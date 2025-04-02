#pragma once
#include "MethodP2P.h"
#include "CPUPF.h"
#include <iostream>
#include <string>
#include <chrono>


class ADMMACConst1 : public MethodP2P
{
public:
	ADMMACConst1();
	ADMMACConst1(float rho);
	virtual ~ADMMACConst1();
	void setParam(float rho);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="ADMMConst";
	void updateGlobalProb();
	void updateLocalProb();
	void updateLambda();
	void updateBt1();
	void updateBt2();
	void updateBp1();
	void updateTl();
	float calcRes();
	float updateResBis(MatrixCPU* res, int iter, MatrixCPU* tempNN);
	void updateP();
	void updateMU();
	void updateQ();
	void display();
private:

	// Pour le r�seau

	CPUPF PF;

	int _nConstraint = 0;
	

	int _nVoisinRef = 0;
	float Ploss = 0;
	float Qloss = 0;
	int iterGlobal = 0;

	MatrixCPU PQ;



	MatrixCPU Ap3;
	MatrixCPU Ap123;
	MatrixCPU rhoMn;
	MatrixCPU rhoMn2;
	MatrixCPU UpperBound;
	MatrixCPU LowerBound;
	MatrixCPU DiffBound;
	MatrixCPU* Constraint;


	MatrixCPU tempC1;
	MatrixCPU tempCN;
	
};


	




