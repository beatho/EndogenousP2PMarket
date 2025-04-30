#pragma once
#include "../head/MethodP2P.h"
 


MethodP2P::MethodP2P() : Method()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "method constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1abcd, Fb2, Fb3, Fb5, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu
}

MethodP2P::~MethodP2P() 
{
}


void MethodP2P::updateP0(const StudyCase& cas){
	
	_id = _id + 1;
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	// Change : Power Limits, cost function


	MatrixCPU Lb(cas.getLb());
	MatrixCPU Ub(cas.getUb());
	MatrixCPU BETA(cas.getBeta());
	
	if (cas.isAC() && !isAC) {
		MatrixCPU aT = cas.geta();
		MatrixCPU bT = cas.getb();
		MatrixCPU PminT = cas.getPmin();
		MatrixCPU PmaxT = cas.getPmax();
		for (int n = 0; n < _nAgent; n++) {
			a.set(n, 0, aT.get(n, 0));
			b.set(n, 0, bT.get(n, 0));
			Pmin.set(n, 0, PminT.get(n, 0));
			Pmax.set(n, 0, PmaxT.get(n, 0));
		}
	}
	else if(!cas.isAC() && isAC){
		throw std::invalid_argument("updateP0 : Study Case is not AC, but this method require AC information");
	}
	else {
		a = (cas.geta());
		b = (cas.getb());
		Pmin = (cas.getPmin());
		Pmax = (cas.getPmax());
	}
	Cp1 = b;
	int indice = 0;

	// hypothese : ce sont les mêmes voisins !!!
	for (int idAgent = 0; idAgent < _nAgent; idAgent++) {
		int Nvoisinmax = (int) nVoisin.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			int idVoisin = (int) CoresLinVoisin.get(indice, 0);
			if(Lb.getNCol()==1){
				matLb.set(indice, 0, Lb.get(idAgent, 0));
				matUb.set(indice, 0, Ub.get(idAgent, 0));
			} else {
				matLb.set(indice, 0, Lb.get(idAgent, idVoisin));
				matUb.set(indice, 0, Ub.get(idAgent, idVoisin));
			}
			Ct.set(indice, 0, BETA.get(idAgent, idVoisin));
			indice = indice + 1;
		}
	}
	for (int idAgent = _nAgentTrue; idAgent < _nAgent; idAgent++) {
		for (int voisin = 0; voisin < (_nAgent - 1); voisin++) {
			int idVoisin = (int) CoresLinVoisin.get(indice, 0);
			if(Lb.getNCol()==1){
				matLb.set(indice, 0, Lb.get(idAgent, 0));
				matUb.set(indice, 0, Ub.get(idAgent, 0));
			} else {
				matLb.set(indice, 0, Lb.get(idAgent, idVoisin));
				matUb.set(indice, 0, Ub.get(idAgent, idVoisin));
			}
			Ct.set(indice, 0, BETA.get(idAgent, idVoisin));
			indice = indice + 1;
		}
	}

	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);
	Cp1.multiplyT(&nVoisin);

	Ap2a = a;
	Ap2.add(&Ap2a, &Ap2b);
	Ap2.multiplyT(&nVoisin);
	Ap2.multiplyT(&nVoisin);
	Ap12.add(&Ap1, &Ap2);
	Ap123.add(&Ap12, &Ap3);
	Cp.add(&Cp1, &Cp2);

#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION
	
}

void MethodP2P::updateLAMBDA()
{
	tempNN = trade;
	tempNN.addTrans(&trade);
	tempNN.multiply(_rho);
	tempNN.multiply(0.5);
	LAMBDA.add(&LAMBDA, &tempNN);
}
void MethodP2P::updateLambda()
{
	for (int t = 0; t < _nTrade; t++) {
		int k = CoresLinTrans.get(t, 0);
		float lamb = 0.5 * _rhog * (tradeLin.get(t, 0) + tradeLin.get(k, 0));
		LAMBDALin.set(t, 0, LAMBDALin.get(t, 0) + lamb);
	}
}



void MethodP2P::updateKappa()
{
	//
	Kappa1.projectNeg();
	Kappa1.add(&lLimit);
	Kappa1.subtract(&Qtot);
	
	Kappa2.projectNeg();
	Kappa2.add(&lLimit);
	Kappa2.add(&Qtot);
	//
}
void MethodP2P::updateCp2()
{
	tempL1.subtractAbs(&Kappa1, &Kappa2);
	//Cp2->multiplyTrans(G, tempL1, 0);

	float r = 0;
	for (int i = 0; i < _nAgent; ++i)
	{
		r = 0;
		for (int k = 0; k < _nLine; ++k)
		{
			r +=  G.get(k, i) * (tempL1.get(k, 0) + 2 * Qpart.get(k, i));
		}
		Cp2.set(i, 0, r);
	}

	Cp2.multiply(_rho1);
	Cp2.multiplyT(&nVoisin);
}
float MethodP2P::updateResMat(int iter)
{
	tempNN.subtract(&Tlocal, &trade);
	
	float resS = tempNN.max2();
	tempNN.set(&Tlocal);
	tempNN.addTrans(&Tlocal);
	float resR = tempNN.max2();
	
	if (iter > 0 && _tau > 1) {
		if (resR > _mu * resS) {
			_rhog = _tau * _rhog;
			_at1 = _rhog;
			//std::cout << iter << ", rho augmente :" << _rhog << std::endl;
		}
		else if (resS > _mu * resR) {// rho = rho / tau_inc;
			_rhog = _rhog / _tau;
			_at1 = _rhog;
			//std::cout << iter << ", rho diminue :" << _rhog << std::endl;
		}
	}

	resF.set(0, iter, resR);
	resF.set(1, iter, resS);
	

	return resR* (resR > resS) + resS * (resR <= resS);
}
float MethodP2P::updateRes(int iter)
{
	for (int t = 0; t < _nTrade; t++) {
		int k = CoresLinTrans.get(t, 0);
		tempNN.set(t, 0, tradeLin.get(t, 0) + tradeLin.get(k, 0));
	}
	float resR = tempNN.max2();

	
	tempNN.subtract(&Tlocal, &tradeLin);
	float resS = tempNN.max2();

	resF.set(0, iter, resR);
	resF.set(1, iter, resS);

	if (iter > 0 && _tau > 1) {
		if (resR > _mu * resS) {
			_rhog = _tau * _rhog;
			_at1 = _rhog;
			//std::cout << iter << ", rho augmente :" << _rhog << std::endl;
		}
		else if (resS > _mu * resR) {// rho = rho / tau_inc;
			_rhog = _rhog / _tau;
			_at1 = _rhog;
			//std::cout << iter << ", rho diminue :" << _rhog << std::endl;
		}
	}



	return resR * (resR > resS) + resS * (resR <= resS);

}
float MethodP2P::updateResEndo(int iter){
	for (int t = 0; t < _nTrade; t++) {
		int k = CoresLinTrans.get(t, 0);
		tempNN.set(t, 0, tradeLin.get(t, 0) + tradeLin.get(k, 0));
	}
	float resR = tempNN.max2();

	
	tempNN.subtract(&Tlocal, &tradeLin);
	float resS = tempNN.max2();

	// version de l'article
	tempL1.set(&Kappa1);
	tempL2.set(&Kappa2);
	Kappa1_pre.projectNeg();
	Kappa2_pre.projectNeg();
	tempL1.projectNeg();
	tempL2.projectNeg();
	tempL1.subtract(&Kappa1_pre);
	tempL2.subtract(&Kappa2_pre);
	tempL1.multiplyT(&tempL1);
	tempL2.multiplyT(&tempL2);
	tempL1.add(&tempL2);

	float resXf = _ratioEps * sqrt(tempL1.max2());
	resF.set(0, iter, resR);
	resF.set(1, iter, resS);
	resF.set(2, iter, resXf);
	return MYMAX(MYMAX(resXf, resS), resR);
}

float MethodP2P::calcRes()
{
	float d1 = Tlocal.max2(&Tlocal_pre);
	float d2 = P.max2(&Tmoy);

	return d1* (d1 > d2) + d2 * (d2 >= d1);
}

void MethodP2P::initSize(const StudyCase& cas){
	_nAgentTrue = cas.getNagent();
	_nAgent = _nAgentTrue + isAC * _nAgentTrue;
	if (cas.isAC() && !isAC) {
		MatrixCPU nVoisinT = cas.getNvoi();
		nVoisin = MatrixCPU(_nAgent, 1);
		for (int n = 0; n < _nAgent; n++) {
			nVoisin.set(n, 0, nVoisinT.get(n, 0));
		}
	}else if(!cas.isAC() && isAC){
		throw std::invalid_argument("initSize : Study Case is not AC, but this method require AC information");
	}
	else {
		nVoisin = cas.getNvoi();
	}
	if(isAC){
		_nLine = cas.getNLine(true);
	} else {
		_nLine = cas.getNLine();
	}
	
	_nBus = cas.getNBus();
	_nTrade = (int) nVoisin.sum();
	_nTradeP = 0;
	int testReduc = 0;
	if(!isAC){
		_nTradeP = _nTrade;
		_nTradeQ = 0;
	} else{
		#pragma omp parallel for reduction(+ : testReduc)
		for (int n = 0; n < _nAgentTrue; n++) {
			testReduc += (int) nVoisin.get(n, 0);
		}
		_nTradeP = testReduc;
		_nTradeQ = _nTrade - _nTradeP;
	}

}

void MethodP2P::initSimParam(const Simparam& sim){
	
	_rhog = sim.getRho();
	_rho1 = sim.getRho1();
	_rhol = _rho;
	if (_rho == 0) {
		_rhol = _rhog;
	}

	_iterG = sim.getIterG();
	_iterL = sim.getIterL();
	_iterIntern = sim.getIterIntern();

	_stepG = sim.getStepG();
	_stepL = sim.getStepL();
	_stepIntern = sim.getStepIntern();

	_epsG = sim.getEpsG();
	_epsX = sim.getEpsGC();
	_epsIntern = sim.getEpsIntern();
	_epsL = sim.getEpsL();
	_ratioEps = _epsG / _epsX;

	resF = MatrixCPU(3, (_iterG / _stepG) + 1);
	resX = MatrixCPU(4, (_iterG / _stepG) + 1);

	tempNN = MatrixCPU(_nTrade, 1, 0);
	tempN1 = MatrixCPU(_nAgent, 1, 0); // plut�t que de re-allouer de la m�moire � chaque utilisation
	tempL1 = MatrixCPU(_nLine, 1, 0);
	tempL2 = MatrixCPU(_nLine, 1, 0);
}

void MethodP2P::initLinForm(const StudyCase& cas){
	
	MatrixCPU BETA(cas.getBeta());
	MatrixCPU Ub(cas.getUb());
	MatrixCPU Lb(cas.getLb());
	
	//std::cout << "mise sous forme lin�aire" << std::endl;

	CoresMatLin = MatrixCPU(_nAgent, _nAgentTrue, -1);
	CoresAgentLin = MatrixCPU(_nAgent + 1, 1);
	CoresLinAgent = MatrixCPU(_nTrade, 1);
	CoresLinVoisin = MatrixCPU(_nTrade, 1);
	CoresLinTrans = MatrixCPU(_nTrade, 1);

	Tlocal_pre = MatrixCPU(_nTrade, 1);
	tradeLin = MatrixCPU(_nTrade, 1);
	LAMBDALin = MatrixCPU(_nTrade, 1);

	matLb = MatrixCPU(_nTrade, 1);
	matUb = MatrixCPU(_nTrade, 1);
	Ct = MatrixCPU(_nTrade, 1);
	
	int indice = 0;

	for (int idAgent = 0; idAgent < _nAgentTrue; idAgent++) { // P
		MatrixCPU omega(cas.getVoisin(idAgent));
		int Nvoisinmax = (int) nVoisin.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			int idVoisin = (int) omega.get(voisin, 0);
			if(Lb.getNCol()==1){
				matLb.set(indice, 0, Lb.get(idAgent, 0));
				matUb.set(indice, 0, Ub.get(idAgent, 0));
			} else {
				matLb.set(indice, 0, Lb.get(idAgent, idVoisin));
				matUb.set(indice, 0, Ub.get(idAgent, idVoisin));
			}
			Ct.set(indice, 0, BETA.get(idAgent, idVoisin));
			tradeLin.set(indice, 0, trade.get(idAgent, idVoisin));
			Tlocal_pre.set(indice, 0, trade.get(idAgent, idVoisin));
			LAMBDALin.set(indice, 0, LAMBDA.get(idAgent, idVoisin));
			CoresLinAgent.set(indice, 0, idAgent);
			CoresLinVoisin.set(indice, 0, idVoisin);
			CoresMatLin.set(idAgent, idVoisin, indice);
			indice = indice + 1;
		}
		CoresAgentLin.set(idAgent + 1, 0, indice);
	}
	//std::cout << " Q " << std::endl;
	for (int idAgent = _nAgentTrue; idAgent < _nAgent; idAgent++) { // Q
		for (int idVoisin = 0; idVoisin < _nAgentTrue; idVoisin++) {
			if (idVoisin != (idAgent - _nAgentTrue)) {
				if(Lb.getNCol()==1){
					matLb.set(indice, 0, Lb.get(idAgent, 0));
					matUb.set(indice, 0, Ub.get(idAgent, 0));
				} else {
					matLb.set(indice, 0, Lb.get(idAgent, idVoisin));
					matUb.set(indice, 0, Ub.get(idAgent, idVoisin));
				}
				tradeLin.set(indice, 0, trade.get(idAgent, idVoisin));
				Tlocal_pre.set(indice, 0, trade.get(idAgent, idVoisin));
				LAMBDALin.set(indice, 0, LAMBDA.get(idAgent, idVoisin));
				CoresLinAgent.set(indice, 0, idAgent);
				CoresLinVoisin.set(indice, 0, idVoisin + _nAgentTrue);
				CoresMatLin.set(idAgent, idVoisin, indice);
				indice = indice + 1;
			}
		}

		CoresAgentLin.set(idAgent + 1, 0, indice);
	}
	#pragma omp parallel for
	for (int lin = 0; lin < _nTrade; lin++) {
		int i = (int) CoresLinAgent.get(lin, 0);
		int j = (int) CoresLinVoisin.get(lin, 0);
		if (lin >= _nTradeP) {
			i -= _nAgentTrue;
		}

		int k = (int) CoresMatLin.get(j, i);
		CoresLinTrans.set(lin, 0, k);
	}

}

void MethodP2P::initDCEndoGrid(const StudyCase& cas){
	
	Kappa1 = MatrixCPU(_nLine, 1, 0);
	Kappa2 = MatrixCPU(_nLine, 1, 0);
	Kappa1_pre = MatrixCPU(_nLine, 1, 0);
	Kappa2_pre = MatrixCPU(_nLine, 1, 0);
	Qpart = MatrixCPU(_nLine, _nAgent, 0);
	Qtot = MatrixCPU(_nLine, 1, 0);
	alpha = MatrixCPU(_nLine, _nAgent, 0);

	G = cas.getPowerSensi();
	lLimit = cas.getLineLimit();

	GTrans = MatrixCPU(_nAgent, _nLine);
	GTrans.setTrans(&G);

	G2 = GTrans;
	G2.multiplyT(&GTrans);

}

void MethodP2P::initCaseParam(const Simparam& sim,const StudyCase& cas){

	Tlocal = MatrixCPU(_nTrade, 1, 0);
	P = MatrixCPU(_nAgent, 1, 0); // moyenne des trades
	LAMBDA = sim.getLambda();
	trade = sim.getTrade();

	// si cas AC, a, b ,Nvoisin, Pmin,Pma n'ont pas la bonne taille !!!
	if (cas.isAC() && !isAC) {
		MatrixCPU aT = cas.geta();
		MatrixCPU bT = cas.getb();
		MatrixCPU PminT = cas.getPmin();
		MatrixCPU PmaxT = cas.getPmax();
		MatrixCPU MUT = sim.getMU(); // facteur reduit i.e lambda_l/_rho
		MatrixCPU TmoyT = sim.getPn();
		a = MatrixCPU(_nAgent, 1);
		b = MatrixCPU(_nAgent, 1);
		Pmin = MatrixCPU(_nAgent, 1);
		Pmax = MatrixCPU(_nAgent, 1);
		MU = MatrixCPU(_nAgent, 1);
		Tmoy = MatrixCPU(_nAgent, 1);
		#pragma omp parallel for
		for (int n = 0; n < _nAgent; n++) {
			a.set(n, 0, aT.get(n, 0));
			b.set(n, 0, bT.get(n, 0));
			Pmin.set(n, 0, PminT.get(n, 0));
			Pmax.set(n, 0, PmaxT.get(n, 0));
			MU.set(n, 0, MUT.get(n, 0));
			Tmoy.set(n, 0, TmoyT.get(n, 0));
		}

	}
	else if(!cas.isAC() && isAC){
		throw std::invalid_argument("initCaseParam : Study Case is not AC, but this method require AC information");
	}
	else {
		a = cas.geta();
		b = cas.getb();
		Pmin = cas.getPmin();
		Pmax = cas.getPmax();
		MU = sim.getMU(); // facteur reduit i.e lambda_l/_rho
		Tmoy = sim.getPn();
	}
	Pn = Tmoy;
}

void MethodP2P::initDCEndoMarket(){
	initP2PMarket();

	Ap2a = a;
	Ap2b = MatrixCPU(_nAgent, 1);

	Cp2 = MatrixCPU(_nAgent, 1);
	Cp1 = b;

	Cp1.multiplyT(&nVoisin);
	

	Ap2b.sum(&G2);
	Ap2b.multiply(2 * _rho1);
	Ap2.add(&Ap2a, &Ap2b);

	Ap2.multiplyT(&nVoisin);
	Ap2.multiplyT(&nVoisin);
	Ap12.add(&Ap1, &Ap2);
	
	Cp = Cp1;
}
void MethodP2P::initP2PMarket(){
	_at1 = _rhog; 
	_at2 = _rhol;

	Ap1 = nVoisin;
	Ap2 = a;
	
	Ap12 = MatrixCPU(_nAgent, 1);

	Bt1 = MatrixCPU(_nTrade, 1);
	Bt2 = MatrixCPU(_nTrade, 1);
	Bp1 = MatrixCPU(_nAgent, 1);
	
	Cp = b;

	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);
	Tmoy.divideT(&nVoisin);

	Ap1.multiply(_rhol);
	Cp.multiplyT(&nVoisin);
	
	Ap2.multiplyT(&nVoisin);
	Ap2.multiplyT(&nVoisin);
	Ap12.add(&Ap1, &Ap2);

	/* not used by default but must exists for P0 */
	Ap2a   = MatrixCPU(_nAgent, 1, 0);
	Ap2b   = MatrixCPU(_nAgent, 1, 0);
	Ap3    = MatrixCPU(_nAgent, 1, 0); 
	Ap123  = MatrixCPU(_nAgent, 1, 0); 
	Cp2    = MatrixCPU(_nAgent, 1, 0);

}

void MethodP2P::setResult(Simparam* result, bool casAC){
	
	if(_nLine){
		Kappa1.projectNeg(); //delta1
		Kappa2.projectNeg(); // delta2
	}
	
	updatePn();
	
	#pragma omp parallel for
	for(int lin = 0; lin <_nTradeP; lin++){
		int idAgent  = (int) CoresLinAgent.get(lin, 0);
		int idVoisin = (int) CoresLinVoisin.get(lin,0); 
		trade.set(idAgent, idVoisin, tradeLin.get(lin, 0));
		LAMBDA.set(idAgent, idVoisin, LAMBDALin.get(lin, 0));
	}
	
	#pragma omp parallel for
	for(int lin = _nTradeP; lin < _nTrade; lin++){
		int idAgent  = (int) CoresLinAgent.get(lin, 0);
		int idVoisin = (int) CoresLinVoisin.get(lin,0) - _nAgentTrue; 
		trade.set(idAgent, idVoisin, tradeLin.get(lin, 0));
		LAMBDA.set(idAgent, idVoisin, LAMBDALin.get(lin, 0));
	}
	
	if(casAC == isAC){
		result->setLAMBDA(&LAMBDA);
		result->setTrade(&trade); 
		result->setMU(&MU);
		result->setPn(&Pn);
	} else if (!casAC && isAC){
		throw std::invalid_argument("setResult : method is AC and case is DC");
	} else{
		MatrixCPU tradeTot(2 * _nAgent, _nAgent);
		MatrixCPU LAMBDATot(2 * _nAgent, _nAgent);
		MatrixCPU PnTot(2 * _nAgent, 1);
		MatrixCPU MUTot(2 * _nAgent, 1);
		#pragma omp parallel for
		for (int n = 0; n < _nAgent; n++) {
			for (int m = 0; m < _nAgent; m++) {
				tradeTot.set(n, m, trade.get(n, m));
				LAMBDATot.set(n, m, LAMBDA.get(n, m));
			}
			PnTot.set(n, 0, Pn.get(n, 0));
			MUTot.set(n, 0, MU.get(n, 0));
		}
		result->setLAMBDA(&LAMBDATot);
		result->setTrade(&tradeTot);
		result->setMU(&MUTot);
		result->setPn(&PnTot);
	}
	
	float fc = calcFc();
	result->setResF(&resF);
	result->setIter(_iterGlobal);
	result->setFc(fc);

	if (_nLine) {
		result->setDelta1(&Kappa1);
		result->setDelta2(&Kappa2);
	}

	tMarket = clock() - tMarket;
	result->setTime((float)tMarket / CLOCKS_PER_SEC);
	
}

float MethodP2P::calcFc()
{

	float fc = 0;
	tempN1.set(&a);
	tempN1.multiply(0.5);
	tempN1.multiplyT(&Pn);
	tempN1.add(&b);
	tempN1.multiplyT(&Pn);
	
	fc = fc + tempN1.sum();
	tempNN.set(&tradeLin);
	tempNN.multiplyT(&Ct);
	

	fc = fc + tempNN.sum();

	return fc;

}

void MethodP2P::updatePn()
{
	Pn.set(&Tmoy);
	Pn.multiplyT(&nVoisin);
}



void MethodP2P::solveWithMinPower(Simparam* result, const Simparam& sim, const StudyCase& cas)
{
	std::cout << "solveWithMinPower : should not be called" << std::endl;
}

void MethodP2P::display(){
	if (_iterGlobal == 0) {
		std::cout << "algorithm not launch" << std::endl;
	}
	else if (_iterGlobal < _iterG) {
		std::cout << "method " << _name << " converged in " << _iterGlobal << " iterations." << std::endl;
		std::cout << "Converged in " << (float)tMarket / CLOCKS_PER_SEC << " seconds" << std::endl;

	}
	else {
		std::cout << "method " << _name << " not converged in " << _iterGlobal << " iterations." << std::endl;
		std::cout << "time taken " << (float) tMarket / CLOCKS_PER_SEC << " seconds" << std::endl;
	}
	resF.display();
	std::cout << "The power error of this state is (constraint) " << resF.get(0, _iterGlobal / _stepG - (_iterGlobal>0)) << " and convergence " << resF.get(1, _iterGlobal / _stepG - (_iterGlobal > 0)) << std::endl;
	std::cout << "===============================================================|" << std::endl;
	std::cout << "      System Summary                                           |" << std::endl;
	std::cout << "===============================================================|" << std::endl;
	std::cout << "Agent            " << _nAgentTrue << std::endl;
	std::cout << "Nombre d'echange " << _nTrade << std::endl;

	std::cout << std::endl << std::endl;

	std::cout << "==================================================================================================|" << std::endl;
	std::cout << "      Agent Data                                                                                  |" << std::endl;
	std::cout << "==================================================================================================|" << std::endl;
	std::cout << " Agent |  Cost   |  Cost   |        Power Injection          |           Power Injection          |" << std::endl;
	std::cout << "  #    |  a (pu) |  b (pu) |  P (pu) | Pmin (pu) | Pmax (pu) |  Q (pu)    | Qmin (pu) | Qmax (pu) |" << std::endl;
	std::cout << "-------|---------|---------|---------|-----------|-----------|------------|-----------|-----------|" << std::endl;

	for (int n = 0; n < _nAgentTrue; n++) {
		
		std::cout << std::setw(7) << n << "|" << std::setw(8) << a.get(n, 0) << " |" << std::setw(9)
			<< b.get(n, 0) << "|" << std::setw(9) << Pn.get(n, 0) << "|" << std::setw(11)
			<< Pmin.get(n, 0) * nVoisin.get(n,0) << "|" << std::setw(11) << Pmax.get(n, 0) * nVoisin.get(n, 0) << "|";
		if (isAC) {
			std::cout << std::setw(12) << Pn.get(n + _nAgentTrue, 0) << "|" << std::setw(11)
				<< Pmin.get(n + _nAgentTrue, 0) * nVoisin.get(n + _nAgentTrue, 0) << "|" << std::setw(11) << Pmax.get(n + _nAgentTrue, 0) * nVoisin.get(n + _nAgentTrue, 0) << "|" << std::endl;
		}
		else {
			std::cout << std::setw(10) << 0 << "|" << std::setw(11)
				<< 0 << "|" << std::setw(11) << 0 << "|" << std::endl;
		}
		
	}


	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "                      END PRINT                                                                         |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;


	std::cout << std::endl << std::endl;
}

