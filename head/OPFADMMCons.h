#pragma once
#include "MethodOPF.h"
#include <iostream>
#include <string>
#include <chrono>
#include "Utilities.h"



class OPFADMMCons : public MethodOPF // OPFADMM2 mais on consid�re que l'on peut avoir un fonction cout sur les pertes, il y a donc un bus fictif
{
public:
	OPFADMMCons();
	OPFADMMCons(float rho);
	
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	virtual void solveConsensus(float eps, MatrixCPU* PSO);
	virtual void initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO);
	virtual void updateConsensus(MatrixCPU* Pmarket);
	
	std::string NAME ="OPFADMMCons";
	
	virtual void updateX();
	void updateXPQ();
	
	void updateChat();
	void CommunicationX();

	
private:
	// connsensus
	int _nBusWLoss = 0;
	MatrixCPU etaSO; // consensus with market
	float _rhoSO = 0;

	const int indPc = 5;
	const int indQc = 6;
	const int indlc = 7;

};


	




