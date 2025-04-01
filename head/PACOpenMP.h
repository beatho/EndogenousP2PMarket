#pragma once
#include "MethodP2P.h"
#include <iostream>
#include <string>
#include <chrono>

#ifdef _OPENMP
	#include "omp.h"
#endif


class PACOpenMP : public MethodP2P // g�re P et Q comme si tous les Q �taient des agents diff�rents
	//et donc P et Q sont 2 probl�mes compl�tement distincts
{
public:
	PACOpenMP();
	PACOpenMP(float rho);
	virtual ~PACOpenMP();
	void setParam(float rho);
	void setGamma(float gamma); 
	void setGammahat(float gammahat);
	void setInitCoef(float alpha, float phi, float theta);
	void setBestRhoGamma(float lambdaMax, float lambdaMin, const StudyCase& cas);
	void setBestRhoGammaHeuristic(const StudyCase& cas);
	void updateCoef();
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	virtual void setBestParam(const StudyCase& cas);
	std::string NAME ="PACOpenMP";
	void updateGlobalProb();
	void updateLocalProb();

	void updateXhat();
	void updateMu();
	void updateNu();
	void updateQ();

	float updateRes(int indice);

	void display();
private:
	bool augmente = true;
	
	// ne change pas avec P0
	float _gamma = 1; // 0.1
	float _gammahat = 0.1;
	bool isAC = false;
		
	int _nAgent = 0;
	int _nAgentTrue = 0;
	int _nTrade = 0;
	int _nTradeP = 0;
	int _nTradeQ = 0;
	float _rhog = 0;
	float _rho = 0;
	float _rhoInv = 0;
	int _sizePACOpenMP = 0;
	
	// parameter agent and iteration dependant (but not here for now)
	float _alpha = 0; // 0.6
	float _phi = 0; // 0.2
	float _theta = 0; // 0.2
	

	MatrixCPU tempN1; // Matrix temporaire pour aider les calculs
	MatrixCPU tempNN; // plut�t que de re-allouer de la m�moire � chaque utilisation
	MatrixCPU* tempM1 = nullptr; //
	MatrixCPU* tempM = nullptr; //
	

	MatrixCPU Cost1;
	MatrixCPU Cost2;
	MatrixCPU BETA;

	MatrixCPU* X = nullptr;
	MatrixCPU* Xpre = nullptr;
	MatrixCPU* Xhat = nullptr;
	MatrixCPU* Mu = nullptr;
	MatrixCPU* Muhat = nullptr;
	MatrixCPU* Nu = nullptr;
	MatrixCPU* Nuhat = nullptr;


	MatrixCPU* Hinv = nullptr;
	MatrixCPU* H = nullptr;

	MatrixCPU nVoisin;
	MatrixCPU CoresMatLin;
	MatrixCPU CoresLinAgent;
	MatrixCPU CoresAgentLin;
	MatrixCPU CoresLinVoisin;
	MatrixCPU CoresLinTrans;
	MatrixCPU CoresLinTransLocal;

	// Matrices for the result
	MatrixCPU Pn;
	MatrixCPU trade;
	MatrixCPU resF;
	MatrixCPU resX;



	// change avec P0
	
	MatrixCPU* matLb = nullptr;
	MatrixCPU* Q = nullptr;
	MatrixCPU* Qinit = nullptr;

	MatrixCPU* matUb = nullptr;

};


	




