#pragma once
#include "MethodOPF.h"
#include <iostream>
#include <string>
#include <chrono>
#include "Utilities.h"


class OPFADMM : public MethodOPF
{
public:
	OPFADMM();
	OPFADMM(float rho);
	
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	virtual void initLocalProb(const StudyCase& cas);
	virtual void solveConsensus(float eps, MatrixCPU* PSO){};
	virtual void initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO){};
	virtual void updateConsensus(MatrixCPU* Pmarket){};
	
	std::string NAME ="OPFADMM";
	
	void updateLocalProb();
	void updateXPQ();

	void updateBp1();
	void updatepl();
	void updatePTilde();
	float calcRes();
	void updateChat();
	void CommunicationX();
	

	virtual float updateRes(int indice);
	
private:

	const int indPc = 7;
	const int indQc = 8;
	const int indlc = 9;
	// ne change pas avec P0
	MatrixCPU _apt1;
	
	// Local resolution
	MatrixCPU PnTmin;
	MatrixCPU PnTmax;
	MatrixCPU PnMoy;
	MatrixCPU PnPre;
	MatrixCPU MuL;
	MatrixCPU PnTilde;
	MatrixCPU Ap12; // a_n + rhol
	MatrixCPU Bp1;
	MatrixCPU Bpt1;
	MatrixCPU Bpt2;
	MatrixCPU Apt1;
	MatrixCPU Apt2;
	MatrixCPU Apt12;

};


	




