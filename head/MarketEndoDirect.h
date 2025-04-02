#pragma once
#include "MethodP2P.h"
#include <iostream>
#include <string>
#include <chrono>
#include "Utilities.h"

#include "ADMMMarket.h"



class MarketEndoDirect : public MethodP2P // OPFADMMCons mais on met des trades sur le bus fictif
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
	int _nBusWLoss = 0;
	float _rhoInv = 0;
	
	double coefPoly2[2];
	double root2[4];
	double root3[4];
	double root4[4];
	double root5[4];
	double root6[4];
	double coefPoly3[3];
	
	
	float _Ploss = 0;
	float _Qloss = 0;
	clock_t timeMarketEndo = 0;


	
	// problem : TradeLin tempN2 cost1 cost2
	LossType losstype = LossType::CURRENT;

	MatrixCPU tempB2; // size : (_nBus  *2, 1)

	
	MatrixCPU Ap3;
	MatrixCPU Ap123;
	MatrixCPU Bp2;
	MatrixCPU Cp;


	// Reseau
	MatrixCPU* tempM1 = nullptr; //
	MatrixCPU* tempM = nullptr; //
	MatrixCPU Pbmax;
	MatrixCPU Pbmin;
	MatrixCPU Pb;

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

};


	




