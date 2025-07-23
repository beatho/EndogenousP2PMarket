#pragma once

#include "Method.h"


class MethodOPF : public Method
{
public:
	MethodOPF();
	MethodOPF(float rho);
	virtual ~MethodOPF();

	void initSimParam(const Simparam& sim);
	void initSize(const StudyCase& cas);
	void initGrid(const StudyCase& cas);
	void initAgent(const StudyCase& cas);
	void allocateTab();
	
	void setParam(float rho);
	virtual float getPLoss();
	virtual float getQLoss();
	virtual float getPLoss2(); // compute with the other method
	virtual float getQLoss2(); // idem
	double calcFc(MatrixCPUD* cost1, MatrixCPUD* cost2, MatrixCPUD* Pn, MatrixCPUD* tempN1);
	float calcFc();
	virtual int feasiblePoint();
	void setLagrange(bool lagrange);
	float getMurhoVar();
	float getTaurhoVar();
	bool chekcase();

	void updateGlobalProb();
	void updateMu();
	void updateXWOCurrent(); // (Pi, Qi, li, vi)
	

	float updateRes(int indice);
	float updateResRhoFixe(int indice);

	virtual void computePb();

	virtual void solveConsensus(float eps, MatrixCPU* PSO) = 0;
	virtual void initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO) = 0;
	virtual void updateConsensus(MatrixCPU* Pmarket) = 0;
	virtual MatrixCPU getPb();
	virtual MatrixCPU getPhi();
	virtual MatrixCPU getE();

	float DFSP(int j);
	float DFSQ(int j);

	virtual void display();
protected:
	LossType losstype = LossType::POWER;
	bool Lagrange = false; // true  false
	bool isCurrentLimited = false;

	// (Pi,Qi,li, vi, vai, pi, qi, pji, qji,lji) (j'ai modifé)
	const int indPi  = 0;
	const int indQi  = 1;
	const int indli  = 2;
	const int indvi  = 3;
	const int indvai = 4;
	const int indpi  = 5;
	const int indqi  = 6;
	const int indChatpi  = 4;
	const int indChatqi  = 5;
	

	float _mu = 100;
	float _tau = 2;
	float _rho = 0;
	float _rhol = 0;
	float _rhoInv = 0;

	int _iterGlobal = 0;
	int _iterLocal = 0;

	int _iterG  = 0;
	int _stepG  = 0;
	int _iterL  = 0;
	int _stepL  = 0;
	float _epsG = 0.0f;
	float _epsL = 0.0f;

	int _nAgent = 0;
	int _nBus = 0;
	int _nLine = 0;
	int _sizeProbGlob = 0;

	float _Ploss = 0;
	float _Qloss = 0;
	clock_t timeOPF = 0;

	
	double coefPoly2[2];
	double root2[4];
	double root3[4];
	double root4[4];
	double root5[4];
	double root6[4];
	double coefPoly3[3];
	
	

	MatrixCPU* X = nullptr; 
	MatrixCPU* Ypre = nullptr;
	MatrixCPU* Y = nullptr; // (Pi,Qi,vi,li,pi,qi, pn, vai,pji,qji,lji)

	MatrixCPU* Mu = nullptr;

	MatrixCPU* Chat = nullptr;

	MatrixCPU* Hinv = nullptr;
	MatrixCPU* Q = nullptr;
	MatrixCPU* A = nullptr;

	MatrixCPU Pn;
	
	MatrixCPU tempN2; // size : (_nAgent*2, 1)
	MatrixCPU tempN1; // Matrix temporaire pour aider les calculs
	MatrixCPU* tempM1 = nullptr; //
	MatrixCPU* tempM = nullptr; //
	

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
	
	

	MatrixCPU nChild;
	MatrixCPU* Childs = nullptr;
	MatrixCPU PosChild; 
	MatrixCPU Ancestor;
	MatrixCPU CoresLineBus;

	MatrixCPU sizeOPFADMM;

	MatrixCPU resF;

	// Local resolution
	
	MatrixCPU tempB2; // size : (_nBus  *2, 1)
	MatrixCPU CoresSoloBusAgent;
	
	MatrixCPU Pmin;
	MatrixCPU Pmax;
	MatrixCPU Pbmax;
	MatrixCPU Pbmin;
	MatrixCPU Pb;

};

