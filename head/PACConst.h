#pragma once
#include "MethodP2P.h"
#include <iostream>
#include <string>
#include <chrono>


class PACConst : public MethodP2P
{
public:
	PACConst();
	PACConst(float rho);
	virtual ~PACConst();
	void setParam(float rho);
	void setGamma(float gamma); 
	void setGammahat(float gammahat);
	void setInitCoef(float alpha, float phi, float theta);
	void setBestRhoGamma(float lambdaMax, float lambdaMin, const StudyCase& cas);
	void updateCoef();
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="PACConst";
	void updateGlobalProb();
	void updateLocalProb();

	void updateXhat();
	void updateMu();
	void updateNu();
	void updateQ();
	void updateDelta();

	float updateRes(int indice);

	void display();
private:
	bool augmente = true;
	// ne change pas avec P0
	float _gamma = 0.01f;
	float _gammahat = 0.01f;
		
	int _nAgent = 0;
	int _nTrade = 0;
	int _nLine = 0;
	float _rhog = 0;
	float _rho = 1.5f;
	float _rhoInv = 0;
	int _sizePACConst = 0;
	
	// parameter agent and iteration dependant (but not here for now)
	float _alpha = 0.6f;
	float _phi = 0.2f;
	float _theta = 0.2f;
	

	MatrixCPU* tempM1 = nullptr; //
	MatrixCPU* tempM = nullptr;; //
	MatrixCPU tempL;

	MatrixCPU BETA;
	MatrixCPU Phi;


	MatrixCPU* X = nullptr;;
	MatrixCPU* Xpre = nullptr;;
	MatrixCPU* Xhat = nullptr;;
	MatrixCPU* Mu = nullptr;;
	MatrixCPU* Muhat = nullptr;;
	MatrixCPU* Nu = nullptr;;
	MatrixCPU* Nuhat = nullptr;;
	MatrixCPU delta;
	MatrixCPU deltahat;

	MatrixCPU* Hinv = nullptr;;
	MatrixCPU* H = nullptr;;

	MatrixCPU nVoisin;
	MatrixCPU CoresMatLin;
	MatrixCPU CoresLinAgent;
	MatrixCPU CoresAgentLin;
	MatrixCPU CoresLinVoisin;
	MatrixCPU CoresLinTrans;
	MatrixCPU CoresLinTransLocal;

	// Matrices for the result
	MatrixCPU Pn;
	

	// change avec P0
	
	MatrixCPU* matLb = nullptr;;
	MatrixCPU* Q = nullptr;;
	MatrixCPU* Qinit = nullptr;;

	MatrixCPU* matUb = nullptr;;

};


	




