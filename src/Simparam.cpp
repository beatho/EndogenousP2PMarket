


#include "../head/Simparam.h"
Simparam::Simparam()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "simparam constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_rho = RHOG;
	_iterMaxGlobal = ITERMAXGLOBAL;
	_iterMaxLocal = ITERMAXLOCAL;
	
	_epsGlobal = EPSGLOBAL;
	_epsGlobalConst = EPSGLOBALCONST;
	_epsLocal = EPSLOCAL;
	
	_stepL = STEPL;
	_stepG = STEPG;
	_nAgent = 0;
	_nLine = 0;
	_iter = 0;
	_time = 0;
	_fc = 0;
	_fcSym = 0;
	_iterLocalTotal = 0;
	_rho1 = RHOG;


	_iterIntern = ITERMAXLOCAL;
	_epsIntern = EPSLOCAL;
	_stepIntern = STEPL;


	_resF = MatrixCPU(3, (_iterMaxGlobal / STEPG) + 1);
}

Simparam::Simparam(const Simparam& sim)
{
	_rho = sim._rho;
	_iterMaxGlobal = sim._iterMaxGlobal;
	_iterMaxLocal = sim._iterMaxLocal;
	_epsGlobal = sim._epsGlobal;
	_epsGlobalConst = sim._epsGlobalConst;
	_epsLocal = sim._epsLocal;
	_nAgent = sim._nAgent;
	_nLine = sim._nLine;
	_stepL = sim._stepL;
	_stepG = sim._stepG;
	_iter = sim._iter;
	_time = sim._time;
	_fc = sim._fc;
	_fcSym = sim._fcSym;
	_iterLocalTotal = sim._iterLocalTotal;
	_rho1 = sim._rho1;
	_lineLimitMin = sim._lineLimitMin;
	_warmstart = sim._warmstart;
	_AC = sim._AC;

	_iterIntern = sim._iterIntern;
	_epsIntern = sim._epsIntern;
	_stepIntern = sim._stepIntern;


	_LAMBDA = MatrixCPU(sim._LAMBDA);
	_trade = MatrixCPU(sim._trade);
	_tradeSym = MatrixCPU(sim._tradeSym);
	_Pn = MatrixCPU(sim._Pn);
	_resF = MatrixCPU(sim._resF);
	_MU = MatrixCPU(sim._MU);
	_delta1 = MatrixCPU(sim._delta1);
	_delta2 = MatrixCPU(sim._delta2);

}



Simparam::Simparam(int nAgent) {
#ifdef DEBUG_CONSTRUCTOR
	std::cout << "simparam constructor 1" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	
	_nAgent = nAgent;
	_rho = RHOG;
	_iterMaxGlobal = ITERMAXGLOBAL;
	_iterMaxLocal = ITERMAXLOCAL;
	_epsGlobal = EPSGLOBAL;
	_epsGlobalConst = EPSGLOBALCONST;
	_epsLocal =  EPSLOCAL;
	_stepL = STEPL;
	_stepG = STEPG;
	_rho1 = RHOG;
	_nLine = 0;


	_iterIntern = ITERMAXLOCAL;
	_epsIntern = EPSLOCAL;
	_stepIntern = STEPL;

	_LAMBDA = MatrixCPU(nAgent, nAgent);
	_trade = MatrixCPU(nAgent, nAgent);
	_tradeSym = MatrixCPU(nAgent, nAgent);
	_Pn = MatrixCPU(nAgent, 1);
	_resF = MatrixCPU(3, (_iterMaxGlobal/STEPG)+1);
	_MU = MatrixCPU(nAgent, 1);
	_delta1 = MatrixCPU(_nLine, 1);
	_delta2 = MatrixCPU(_nLine, 1);
	_iter = 0;
	_time = 0;
	_fc = 0;
	_fcSym = 0;
	_iterLocalTotal = 0;
}

Simparam::Simparam(int nAgent, int nLine, bool AC)
{
	_nAgent = nAgent;
	_rho = RHOG;
	_iterMaxGlobal = ITERMAXGLOBAL;
	_iterMaxLocal = ITERMAXLOCAL;
	_epsGlobal = EPSGLOBAL;
	_epsGlobalConst = EPSGLOBALCONST;
	_epsLocal = EPSLOCAL;
	_stepL = STEPL;
	_stepG = STEPG;
	_rho1 = RHOG;
	_nLine = nLine;
	_AC = AC;

	_iterIntern = ITERMAXLOCAL;
	_epsIntern = EPSLOCAL;
	_stepIntern = STEPL;
	
	if (AC) {
		_Pn = MatrixCPU(2 * nAgent, 1);
		_MU = MatrixCPU(2 * nAgent, 1);
		_LAMBDA = MatrixCPU(2 * nAgent, nAgent);
		_trade = MatrixCPU(2 * nAgent, nAgent);
		_tradeSym = MatrixCPU(2 * nAgent, nAgent);
	}
	else {
		_Pn = MatrixCPU(nAgent, 1);
		_MU = MatrixCPU(nAgent, 1);
		_LAMBDA = MatrixCPU(nAgent, nAgent);
		_trade = MatrixCPU(nAgent, nAgent);
		_tradeSym = MatrixCPU(nAgent, nAgent);
	}
	
	_resF = MatrixCPU(3, (_iterMaxGlobal / STEPG) + 1);
	
	_delta1 = MatrixCPU(_nLine, 1);
	_delta2 = MatrixCPU(_nLine, 1);
	_iter = 0;
	_time = 0;
	_fc = 0;
	_fcSym = 0;
	_iterLocalTotal = 0;
}

Simparam::Simparam(int nAgent, int nLine, int nLineConstraint, bool AC)
{
	_nAgent = nAgent;
	_rho = RHOG;
	_iterMaxGlobal = ITERMAXGLOBAL;
	_iterMaxLocal = ITERMAXLOCAL;
	_epsGlobal = EPSGLOBAL;
	_epsGlobalConst = EPSGLOBALCONST;
	_epsLocal = EPSLOCAL;
	_stepL = STEPL;
	_stepG = STEPG;
	_rho1 = RHOG;
	_nLine = nLine;

	_iterIntern = ITERMAXLOCAL;
	_epsIntern = EPSLOCAL;
	_stepIntern = STEPL;

	_AC = AC;

	if (AC) {
		_Pn = MatrixCPU(2 * nAgent, 1);
		_MU = MatrixCPU(2 * nAgent, 1);
		_LAMBDA = MatrixCPU(2 * nAgent, nAgent);
		_trade = MatrixCPU(2 * nAgent, nAgent);
		_tradeSym = MatrixCPU(2 * nAgent, nAgent);
	}
	else {
		_Pn = MatrixCPU(nAgent, 1);
		_MU = MatrixCPU(nAgent, 1);
		_LAMBDA = MatrixCPU(nAgent, nAgent);
		_trade = MatrixCPU(nAgent, nAgent);
		_tradeSym = MatrixCPU(nAgent, nAgent);
	}

	_resF = MatrixCPU(3, (_iterMaxGlobal / STEPG) + 1);

	_delta1 = MatrixCPU(nLineConstraint, 1);
	_delta2 = MatrixCPU(nLineConstraint, 1);
	_iter = 0;
	_time = 0;
	_fc = 0;
	_fcSym = 0;
	_iterLocalTotal = 0;
}


Simparam::Simparam(float rho, int iterMaxGlobal, int iterMaxLocal, float epsGlobal, float epsLocal, int nAgent)
{
	#if DEBUG_CONSTRUCTOR
		std::cout << "constructeur simparam 2" << std::endl;
	#endif // DEBUG_CONSTRUCTOR
	_rho = rho;
	_iterMaxGlobal = iterMaxGlobal;
	_iterMaxLocal = iterMaxLocal;
	_nAgent = nAgent;
	_epsGlobal = epsGlobal;
	_epsGlobalConst = epsGlobal;
	_epsLocal =  epsLocal;

	_iterIntern = iterMaxLocal;
	_epsIntern = epsLocal;
	_stepIntern = STEPL;
	_nLine = 0;
	
	_stepL = STEPL;
	_stepG = STEPG;
	_LAMBDA = MatrixCPU(nAgent, nAgent);
	_trade = MatrixCPU(nAgent, nAgent);
	_tradeSym = MatrixCPU(nAgent, nAgent);
	_Pn = MatrixCPU(nAgent, 1);
	_resF = MatrixCPU(3, (_iterMaxGlobal / STEPG)+1);
	_MU = MatrixCPU(nAgent, 1);
	_delta1 = MatrixCPU(_nLine, 1);
	_delta2 = MatrixCPU(_nLine, 1);
	_iter = 0;
	_time = 0;
	_fc = 0;
	_fcSym = 0;
	_iterLocalTotal = 0;
}

void Simparam::reset(int oldN, int oldL, bool AC)
{

	if (_nAgent != oldN) {
		//std::cout << "reset Trade, Lambda, Pn, Mu" << std::endl;
		int n = _nAgent;
		if (AC) {
			n *= 2;
		}
		_trade = MatrixCPU(n, n);
		_tradeSym = MatrixCPU(n, n);
		_LAMBDA = MatrixCPU(n, n);
		_Pn = MatrixCPU(n, 1);
		_MU = MatrixCPU(n, 1);
	}
	if (_nLine != oldL) {
		//std::cout << "reset delta" << std::endl;
		_delta1 = MatrixCPU(_nLine, 1);
		_delta2 = MatrixCPU(_nLine, 1);
	}
	timePerBlock = MatrixCPU(1, 11, 0); 
	occurencePerBlock = MatrixCPU(1, 11, 0); 
	
}

Simparam& Simparam::operator=(const Simparam& sim)
{
#if DEBUG_CONSTRUCTOR
	std::cout << "operateur =" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_rho = sim._rho;
	_iterMaxGlobal = sim._iterMaxGlobal;
	_iterMaxLocal = sim._iterMaxLocal;
	_epsGlobal = sim._epsGlobal;
	_epsGlobalConst = sim._epsGlobalConst;
	_epsLocal = sim._epsLocal;
	_nAgent = sim._nAgent;
	_nLine = sim._nLine;
	_AC = sim._AC;


	_iterIntern = sim._iterIntern;
	_epsIntern = sim._epsIntern;
	_stepIntern = sim._stepIntern;

    _iter = sim._iter;
	_iterLocalTotal = sim._iterLocalTotal;
	_time = sim._time;
	_fc = sim._fc;
	_fcSym = sim._fcSym;

	_lineLimitMin = sim._lineLimitMin;
	_warmstart = sim._warmstart;

	_LAMBDA = MatrixCPU(sim._LAMBDA);
	_trade = MatrixCPU(sim._trade);
	_tradeSym = MatrixCPU(sim._tradeSym);
	_Pn = MatrixCPU(sim._Pn);
	_resF = MatrixCPU(sim._resF);
	_MU = MatrixCPU(sim._MU);
	_delta1 = MatrixCPU(sim._delta1);
	_delta2 = MatrixCPU(sim._delta2);
	


	return *this;
}

Simparam::~Simparam()
{
}

float Simparam::getRho() const
{
	return _rho;
}

float Simparam::getRho1() const
{
	return _rho1;
}

int Simparam::getIterG() const
{
	return _iterMaxGlobal;
}

int Simparam::getIter() const
{
	return _iter;
}

int Simparam::getIterLTot() const
{
	return _iterLocalTotal;
}

int Simparam::getIterL() const
{
	return _iterMaxLocal;
}

int Simparam::getIterIntern() const
{
	return _iterIntern;
}

float Simparam::getEpsG() const
{
	return _epsGlobal;
}

float Simparam::getEpsGC() const
{
	return _epsGlobalConst;
}

float Simparam::getEpsL() const
{
	return _epsLocal;
}

float Simparam::getEpsIntern() const
{
	return _epsIntern;
}

int Simparam::getNAgent() const
{
	return _nAgent;
}

int Simparam::getNLine() const
{
	return _nLine;
}

int Simparam::getStepL() const
{
	return _stepL;
}
int Simparam::getStepIntern() const
{
	return _stepIntern;
}
int Simparam::getStepG() const
{
	return _stepG;
}


MatrixCPU Simparam::getTrade() const {
	return _trade;
}
MatrixCPU Simparam::getTradeSym() const
{
	return _tradeSym;
}
MatrixCPU Simparam::getLambda() const
{
	return _LAMBDA;
}
MatrixCPU Simparam::getRes() const {
	return _resF;
}

MatrixCPU Simparam::getPn() const
{
	return _Pn;
}

MatrixCPU Simparam::getDelta1()
{
	return _delta1;
}

MatrixCPU Simparam::getDelta2()
{
	return _delta2;
}

float Simparam::getTime() const
{
	return _time;
}

MatrixCPU Simparam::getMU() const
{
	return _MU;
}

float Simparam::getFc() const
{
	return _fc;
}

float Simparam::getFcSym() const
{
	return _fcSym;
}

void Simparam::setItG(int iter) {
	_iterMaxGlobal = iter;
	_resF = MatrixCPU(3, (_iterMaxGlobal/_stepG)+1);
}

void Simparam::setItL(int iter) {
	_iterMaxLocal = iter;
}

void Simparam::setItIntern(int iter)
{
	_iterIntern = iter;
}


void Simparam::setLAMBDA(MatrixCPU* m)
{
	_LAMBDA = *m;
}

void Simparam::setTrade(MatrixCPU* m)
{
	_trade = *m;
}

void Simparam::setTradeSym(MatrixCPU* m)
{
	_tradeSym = *m;
}


void Simparam::setResF(MatrixCPU* m)
{
	_resF = *m;
}

void Simparam::setMU(MatrixCPU* m)
{
	_MU = *m;
}
void Simparam::setPn(MatrixCPU* m)
{
	_Pn = *m;
}

void Simparam::setDelta1(MatrixCPU* delta1)
{
	_delta1 = *delta1;
}

void Simparam::setDelta2(MatrixCPU* delta2)
{
	_delta2 = *delta2;
}

void Simparam::setIter(int c)
{
	_iter = c;
}

void Simparam::setTime(float f)
{
	_time = f;
}

void Simparam::setTimeBloc(MatrixCPU* time, MatrixCPU* occurrence)
{
	timePerBlock = *(time);
	occurencePerBlock = *(occurrence);
}

void Simparam::setFc(float f)
{
	_fc = f;
}

void Simparam::setFcSym(float f)
{
	_fcSym = f;
}

void Simparam::setNagent(int n)
{
	int oldN = _nAgent;
	_nAgent = n;
	reset(oldN,_nLine);
}

void Simparam::setNLine(int l)
{
	int oldL = _nLine;
	_nLine = l;
	reset(_nAgent,oldL);
}

void Simparam::setNAgentLine(int n, int l, bool AC)
{
	//std::cout << "setNAgentLine " << n << " " << l << " " << _nAgent << " " << _nLine << std::endl;
	int oldN = _nAgent;
	int oldL = _nLine;
	_nAgent = n;
	_nLine = l;
	reset(oldN, oldL, AC);
}

void Simparam::setItLTot(int iterLocalTotal)
{
	_iterLocalTotal = iterLocalTotal;
}

void Simparam::setRho(float rho) {
	_rho = rho;
}

void Simparam::setRho1(float rho1)
{
	_rho1 = rho1;
}

void Simparam::setStep(int stepG, int stepL)
{
	_stepG = stepG;
	_stepL = stepL;
	_resF = MatrixCPU(3, (_iterMaxGlobal / _stepG) + 1);
}

void Simparam::setStep(int stepG, int stepL, int stepIntern)
{
	_stepG = stepG;
	_stepL = stepL;
	_stepIntern = stepIntern;
	_resF = MatrixCPU(3, (_iterMaxGlobal / _stepG) + 1);
}

void Simparam::setEpsG(float epsG)
{
	_epsGlobal = epsG;
}

void Simparam::setEpsGC(float epsGConst)
{
	_epsGlobalConst = epsGConst;
}

void Simparam::setEpsL(float epsL)
{
	_epsLocal = epsL;
}

void Simparam::setEpsIntern(float eps)
{
	_epsIntern = eps;
}



void Simparam::saveCSV(const std::string& fileName, int all)
{
	/*
	float _rho;
	int _iterMaxGlobal;
	int _iterMaxLocal;
	float _epsGlobal;
	float _epsLocal;
	int _nAgent;
	MatrixCPU _LAMBDA;
	MatrixCPU _trade;
	MatrixCPU _Pn;
	MatrixCPU _resF;
	MatrixCPU MU;
	int _iter;
	float _time;
	float _fc;*/
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	

	int Nligne = 1 + _nAgent + _nAgent + 1 + 3 + 1;
	if (_MU.getNLin() > 0) {
		Nligne = Nligne + _nAgent;
	}

	if (all) {
		MatrixCPU param(1, 7);
		param.set(0, 0, _rho);
		param.set(0, 1, _iterMaxGlobal);
		param.set(0, 2, _iterMaxLocal);
		param.set(0, 3, _epsGlobal);
		param.set(0, 4, _epsLocal);
		param.set(0, 5, _nAgent);
		param.set(0, 6, Nligne);
		param.saveCSV(fileName, mode);
	}
	
	_LAMBDA.saveCSV(fileName, mode);
	_trade.saveCSV(fileName, mode);
	if (all) {

		MatrixCPU temp(1, _nAgent);
		temp.addTrans(&_Pn);
		temp.saveCSV(fileName, mode);
		_resF.saveCSV(fileName, mode);
		_MU.saveCSV(fileName, mode);

		MatrixCPU result(1, 3);
		result.set(0, 0, _iter);
		result.set(0, 1, _time);
		result.set(0, 2, _fc);
		result.saveCSV(fileName, mode);
	}
}

void Simparam::display(int type) const
{  // type =1 parameters, type = 2 all result, other just some results
	if (type==1) {
		std::cout << "Simulation's parameters : " << std::endl;
		std::cout << "Agents' count : " << _nAgent << std::endl;
		std::cout << "rho : " << _rho << std::endl;
		std::cout << "k_max : " << _iterMaxGlobal << std::endl;
		std::cout << "j_max : " << _iterMaxLocal << std::endl;
		std::cout << "eps_g : " << _epsGlobal << std::endl;
		std::cout << "eps_l : " << _epsLocal << std::endl;
		std::cout << "StepG/StepL : " << _stepG << " / " << _stepL << std::endl;

	}
	else if (type == 2) {
		std::cout << " Simulation result : " << std::endl;
		std::cout << " Agents' count : " << _nAgent << std::endl;
		std::cout << "f_c : " << _fc << std::endl;
		std::cout << "iter : " << _iter << std::endl;
		std::cout << "Pn : " << std::endl;
		_Pn.display();
		std::cout << "Residuals : " << _resF.get(0, (_iter - 1)/_stepG) << " " << _resF.get(1, (_iter - 1) / _stepG) << " " << _resF.get(2, (_iter - 1) / _stepG) << std::endl;
		std::cout << "computation time " << _time << std::endl;
		
		std::cout << "Trades : " << std::endl;
		_trade.display();
		std::cout << "LAMBDA : " << std::endl;
		_LAMBDA.display();
		std::cout << "delta1 : " << std::endl;
		_delta1.display();
		std::cout << "delta2 : " << std::endl;
		_delta2.display();
	}
	else {
		std::cout << " Simulation result : " << std::endl;
		std::cout << " Agents' count : " << _nAgent << std::endl;
		std::cout << "f_c : " << _fc << std::endl;
		std::cout << "iter : " << _iter << std::endl;
		std::cout << "Residuals : " << _resF.get(0, (_iter - 1) / _stepG) << " " << _resF.get(1, (_iter - 1) / _stepG) << " " << _resF.get(2, (_iter - 1) / _stepG) << std::endl;
		std::cout << "computation time : " << _time << std::endl;
	}
}

void Simparam::displayTime(std::string fileName) const
{
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	float factor = 1000000; // go from ns to ms fot the time

	
	/*if (occurencePerBlock.get(0, 0) != 0) {
		std::cout << "total resolution time :" << timePerBlock.sum() / (1000 * factor) << "s" << std::endl;
		std::cout << " Fb0 : " << timePerBlock.get(0, 0) / factor << "ms and occurence :" << occurencePerBlock.get(0, 0) << std::endl;
		if (occurencePerBlock.get(0, 2) != 0) {
			std::cout << " Fb1a : " << timePerBlock.get(0, 1) / factor << "ms and occurence :" << occurencePerBlock.get(0, 1) << std::endl;
			std::cout << " Fb1b : " << timePerBlock.get(0, 2) / factor << "ms and occurence :" << occurencePerBlock.get(0, 2) << std::endl;
			std::cout << " Fb1c : " << timePerBlock.get(0, 3) / factor << "ms and occurence :" << occurencePerBlock.get(0, 3) << std::endl;
		}
		else {
			std::cout << " Fb1 : " << timePerBlock.get(0, 1) / factor << "ms and occurence :" << occurencePerBlock.get(0, 1) << std::endl;
		}


		std::cout << " Fb2 : " << timePerBlock.get(0, 4) / factor << "ms and occurence :" << occurencePerBlock.get(0, 4) << std::endl;
		if (occurencePerBlock.get(0, 5) != 0) {
			std::cout << " Fb3a : " << timePerBlock.get(0, 5) / factor << "ms and occurence :" << occurencePerBlock.get(0, 5) << std::endl;
			std::cout << " Fb3b : " << timePerBlock.get(0, 6) / factor << "ms and occurence :" << occurencePerBlock.get(0, 6) << std::endl;
			std::cout << " Fb3c : " << timePerBlock.get(0, 7) / factor << "ms and occurence :" << occurencePerBlock.get(0, 7) << std::endl;
		}
		else {
			std::cout << " Fb3 : " << timePerBlock.get(0, 5) / factor << "ms and occurence :" << occurencePerBlock.get(0, 5) << std::endl;
		}

		std::cout << " Fb4 : " << timePerBlock.get(0, 8) / factor << "ms and occurence :" << occurencePerBlock.get(0, 8) << std::endl;
		std::cout << " Fb5 : " << timePerBlock.get(0, 9) / factor << "ms and occurence :" << occurencePerBlock.get(0, 9) << std::endl;
		if (occurencePerBlock.get(0, 10) != 0) {
			std::cout << " Fb0 update : " << timePerBlock.get(0, 10) / factor << "ms and occurence :" << occurencePerBlock.get(0, 10) << std::endl;
		}
	}
	else {
		std::cout << "pas de temps � afficher, ou alors il n'y a pas eut d'initialisation" << std::endl;
	}*/

	occurencePerBlock.saveCSV(fileName, mode);
	timePerBlock.saveCSV(fileName, mode);

}
