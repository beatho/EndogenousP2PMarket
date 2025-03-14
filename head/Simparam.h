#pragma once
#include "MatrixCPU.h"
#include "ParamInterface.h"
#include <iostream>
#include <string>


class Simparam
{
public:
	static constexpr float RHOG = 1.5f;
	static const int ITERMAXGLOBAL = 100;
	static const int ITERMAXLOCAL = 10000;
	static constexpr float EPSGLOBAL = 0.0001f; // 1e-4
	static constexpr float EPSLOCAL = 0.00001f; // 1e-5 (pas plus car float !!!!!)
	static constexpr float EPSGLOBALCONST = 0.01f; // 1e-2
	static const int STEPL = 5;
	static const int STEPG = 5;
	float _lineLimitMin = 0;
	float offsetConstraint = 0;
	bool _AC = false;
	bool _warmstart = true;
	Simparam();
	Simparam(const Simparam& sim);
	Simparam(int nAgent);
	Simparam(int nAgent, int nLine, bool AC = false);
	Simparam(int nAgent, int nLine, int nLineConstraint, bool AC = false);
	Simparam(float rho, int iterMaxGlobal, int iterMaxLocal, float epsGlobal, float epsLocal, int _nAgent);

	void reset(int oldN=0, int oldL=0, bool AC = false);

	Simparam& operator= (const Simparam& sim);
	void setFromInterface(ParamInterface* param, bool AC);
	void convertToResultInterface(ResultInterface* res);
	~Simparam();
	
	// getter

	float getRho() const;
	float getRho1() const;

	int getIterG() const;
	int getIter() const;
	int getIterLTot() const;
	int getIterL() const;
	int getIterIntern() const;

	float getEpsG() const;
	float getEpsGC() const;
	float getEpsL() const;
	float getEpsIntern() const;
	int getNAgent() const;
	int getNLine() const;

	int getStepL() const;
	int getStepG() const;
	int getStepIntern() const;


	MatrixCPU getTrade() const;
	MatrixCPU getTradeSym() const;
	MatrixCPU getLambda() const;
	MatrixCPU getRes() const;
	MatrixCPU getPn() const;
	MatrixCPU getDelta1() const;
	MatrixCPU getDelta2() const;
	MatrixCPU getPb() const;
	MatrixCPU getPhi() const;
	MatrixCPU getE() const;

	float getTime() const;
	MatrixCPU getMU() const;
	float getFc() const;
	float getFcSym() const;
	
	// setter
	void setItG(int iter);
	void setItL(int iter);
	void setItIntern(int iter);

	void setLAMBDA(MatrixCPU* m);
	void setTrade(MatrixCPU* m);
	void setTradeSym(MatrixCPU* m);
	void setResF(MatrixCPU* m);
	void setMU(MatrixCPU* m);
	void setPn(MatrixCPU* m);
	void setDelta1(MatrixCPU* delta1);
	void setDelta2(MatrixCPU* delta2);
	void setPb(MatrixCPU* Pb);
	void setPhi(MatrixCPU* Phi);
	void setE(MatrixCPU* E);
	void setIter(int c);
	void setTime(float f);
	void setTimeBloc(MatrixCPU* time, MatrixCPU* occurrence);
	void setFc(float f);
	void setFcSym(float f);
	void setNagent(int n);
	void setNLine(int l);
	void setNAgentLine(int n, int l, bool AC = false);
	void setItLTot(int iterLocalTotal);

	void setRho(float rho);
	void setRho1(float rho1);
	void setStep(int stepG, int stepL);
	void setStep(int stepG, int stepL, int stepIntern);
	void setEpsG(float epsG);
	void setEpsGC(float epsGConst);
	void setEpsL(float epsL);
	void setEpsIntern(float eps);

	// other
	void saveCSV(const std::string& filename, int all=1);
	void display(int type=0) const;
	void displayTime(std::string fileName = "SimulationFB.csv") const;

private:
	/* Param */
	float _rho;
	float _rho1;
	int _iterMaxGlobal;
	int _iterMaxLocal;
	int _iterIntern;
	int _stepG;
	int _stepL;
	int _stepIntern = 1;
	float _epsGlobal;
	float _epsGlobalConst;
	float _epsLocal;
	float _epsIntern;
	int _nAgent;
	int _nLine;
	int _nLineConst = 0;
	int _nBus;
	/* Result Market*/
	MatrixCPU _LAMBDA;
	MatrixCPU _trade;
	MatrixCPU _tradeSym;
	MatrixCPU _MU;
	float _fcSym;
	/*Result DC-Endo*/
	MatrixCPU _delta1;
	MatrixCPU _delta2;

	/*Result OPF and PF*/
	MatrixCPU _Pb;
	MatrixCPU _Phi;
	MatrixCPU _E;

	/* Result All */
	int _iter;
	int _iterLocalTotal;
	float _time;
	MatrixCPU _Pn;
	MatrixCPU _resF;
	
	float _fc;
	
	MatrixCPU timePerBlock; // Fb0, Fb1abc, Fb2, Fb3, Fb4, Fb5, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	MatrixCPU occurencePerBlock; //nb de fois utilis� pendant la simu

};




