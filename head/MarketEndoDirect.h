#pragma once
#include "Method.h"
#include <iostream>
#include <string>
#include <chrono>
#include "Utilities.h"

#include "ADMMMarket.h"



class MarketEndoDirect : public Method // OPFADMMCons mais on met des trades sur le bus fictif
{
public:

	
	MarketEndoDirect();
	MarketEndoDirect(float rho);
	virtual ~MarketEndoDirect();
	void setParam(float rho);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	
	
	
	
	bool chekcase();
	std::string NAME ="MarketEndoDirect";
	void updateGlobalProb();
	
	virtual void updateX();
	void updateXWOCurrent();

	virtual float getPLoss();
	virtual float getQLoss();

	void updateMu();
	void updateChat();
	void CommunicationX();
	float updateRes(int indice);
	float updateResRhoFixe(int indice);
	virtual int feasiblePoint();

	//Market
	void initMarket(const Simparam& sim, const StudyCase& cas);
	void updateBt1();
	void updateBt2();
	void updateBp1();
	void updateBp2();
	void updateTl();
	void updateP();
	void updateMU();
	void updateLocalProb();
	void updateLambda();
	void updatePMarket(); // suite de updateX
	float calcRes();
	bool initWithMarketClear = true;


	void display();
	void computePb();
	virtual MatrixCPU getPb();
	virtual MatrixCPU getPhi();
	virtual MatrixCPU getE();

private:
	float DFSP(int j); // recursive fonction to initialize
	float DFSQ(int j); // P and Q flows

	// ne change pas avec P0
	int _nBus = 0;
	int _nBusWLoss = 0;
	int _nLine = 0;
	int _nAgent = 0;
	
	float _rho = 0;
	float _rhoInv = 0;
	
	double coefPoly2[2];
	double root2[4];
	double root3[4];
	double root4[4];
	double root5[4];
	double root6[4];
	double coefPoly3[3];
	
	int _iterGlobal = 0;
	int _iterG = 0;
	int _stepG = 0;
	float _mu = 40;
	float _tau = 2;
	float _Ploss = 0;
	float _Qloss = 0;
	clock_t timeMarketEndo = 0;

	// market
	int _nAgentTrue = 0; // _nAgent = _nAgentTrue + (isAc)*_nAgent
	int _nTrade = 0;
	int _nTradeP = 0;
	int _nTradeQ = 0;
	float _rhol = 0;
	float _at1 = 0;
	float _at2 = 0;
	float _epsL = 0;
	int _stepL = 0;
	int _iterL = 0;
	
	LossType losstype = LossType::CURRENT;

	
	MatrixCPU LAMBDA;
	MatrixCPU trade;

	MatrixCPU Tlocal;
	MatrixCPU Tlocal_pre;
	MatrixCPU TradeLin;
	MatrixCPU P; // moyenne des trades
	MatrixCPU Pn; // somme des trades
	MatrixCPU tempN2; // size : (_nAgent*2, 1)
	MatrixCPU tempB2; // size : (_nBus  *2, 1)

	MatrixCPU Pmin;
	MatrixCPU Pmax;
	MatrixCPU Ap1;
	MatrixCPU Ap2;
	MatrixCPU Ap3;
	MatrixCPU Ap123;
	MatrixCPU Bt1;
	MatrixCPU Bt2;
	MatrixCPU Bp1;
	MatrixCPU Bp2;
	MatrixCPU Ct;
	MatrixCPU Cp;
	MatrixCPU matUb;
	MatrixCPU matLb;

	MatrixCPU nVoisin;
	MatrixCPU Tmoy;
	MatrixCPU LAMBDALin;
	
	MatrixCPU MU;

	//corespondance
	MatrixCPU CoresMatLin;
	MatrixCPU CoresLinAgent;
	MatrixCPU CoresAgentLin;
	MatrixCPU CoresLinVoisin;
	MatrixCPU CoresLinTrans;


	MatrixCPU tempN1; // Matrix temporaire pour aider les calculs
	MatrixCPU tempNN; // plut�t que de re-allouer de la m�moire � chaque utilisation


	// Reseau
	MatrixCPU* tempM1 = nullptr; //
	MatrixCPU* tempM = nullptr; //
	MatrixCPU Pbmax;
	MatrixCPU Pbmin;
	MatrixCPU Pb;

	MatrixCPU Cost1;
	MatrixCPU Cost2;
	MatrixCPU _CoresBusAgent;
	MatrixCPU _CoresAgentBus;
	MatrixCPU _CoresAgentBusBegin;

	MatrixCPU PosAgent;
	MatrixCPU _nAgentByBus;

	MatrixCPU ZsRe;
	MatrixCPU ZsIm;
	MatrixCPU ZsNorm;
	MatrixCPU VoltageLimit; // (vmin^2, vmax^2) * sqrt(Nchild + 1 / 2)
	MatrixCPU VoltageLimitReal; // vmin, vmax
	MatrixCPU FluxLimit;
	
	MatrixCPU* X = nullptr; // (Pi,Qi,vi,li,pi,qi, pn, vai,pji,qji,lji)
	MatrixCPU* Ypre = nullptr;
	MatrixCPU* Y = nullptr; // (Pi,Qi,vi,li,pi,qi, pn, vai,pji,qji,lji)

	MatrixCPU* Mu = nullptr;

	MatrixCPU* Chat = nullptr;

	MatrixCPU* Hinv = nullptr;
	MatrixCPU* Q = nullptr;
	MatrixCPU* A = nullptr;

	MatrixCPU nChild;
	MatrixCPU* Childs = nullptr;
	MatrixCPU PosChild; 
	MatrixCPU Ancestor;
	MatrixCPU CoresLineBus;

	MatrixCPU sizeMarketEndoDirect;

	MatrixCPU resF;


};


	




