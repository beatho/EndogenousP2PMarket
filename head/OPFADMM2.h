#pragma once
#include "MethodOPF.h"
#include <iostream>
#include <string>
#include <chrono>
#include "Utilities.h"


class OPFADMM2 : public MethodOPF // on remonte les puissances des agents sur les variables duales et pi n'est plus une variable
{
public:
	OPFADMM2();
	OPFADMM2(float rho);
	
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	virtual void solveConsensus(float eps, MatrixCPU* PSO){};
	virtual void initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO){};
	virtual void updateConsensus(MatrixCPU* Pmarket){};

	std::string NAME ="OPFADMM2";
	
	void updateXPQ();

	void updateChat();
	void CommunicationX();
		
private:
	const int indPc = 5;
	const int indQc = 6;
	const int indlc = 7;

};


	




