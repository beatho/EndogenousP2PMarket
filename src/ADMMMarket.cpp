#include "../head/ADMMMarket.h"
#define MAX(X, Y) X * (X >= Y) + Y * (Y > X)


ADMMMarket::ADMMMarket() : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " ADMMMarket Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
}


ADMMMarket::ADMMMarket(float rho) : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default ADMMMarket Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
}

ADMMMarket::~ADMMMarket()
{
}
void ADMMMarket::setParam(float rho)
{
	_rho = rho;
}

void ADMMMarket::setTau(float tau)
{
	if (tau < 1) {
		throw std::invalid_argument("tau must be greater than 1");
	}
	_tau = tau;
}



void ADMMMarket::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
{
#ifdef DEBUG_SOLVE
	cas.display();
	sim.display(1);
#endif // DEBUG_SOLVE
	
	tMarket =clock();
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;
#endif // INSTRUMENTATION
	
	// FB 0
	if (_id == 0) {
#ifdef INSTRUMENTATION
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		init(sim, cas);
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 0, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION
	}
	
	_rhog = sim.getRho();
	_at1 = _rhog;

	int iterL = sim.getIterL();
	int stepL = sim.getStepL();
	float epsL = sim.getEpsL();
	float epsG = sim.getEpsG();
	

	float fc = 0;

	int iterLocal = 0;
	_iterGlobal = 0;
	float resG = 2 * epsG;
	float resL = 2 * epsL;

	while ((_iterGlobal < _iterG) && (resG>epsG)) {
		/*Tlocal.display();
		P.saveCSV("testPCPU2.csv", 11, 1);
		Tlocal.saveCSV("testTCPU2.csv", 11, 1);
		Bt1.saveCSV("testBCPU2.csv", 11, 1);*/
		//P.display();
		//std::cout << "lambda" << std::endl;
		//LAMBDALin.display();
		resL = 2 * epsL;
		iterLocal = 0;
		while (iterLocal< iterL && resL>epsL) {
			updateLocalProb();
			// FB 2
			if (!(iterLocal % stepL)) {
#ifdef INSTRUMENTATION
				t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
				resL = calcRes();
#ifdef INSTRUMENTATION
				t2 = std::chrono::high_resolution_clock::now();
				timePerBlock.increment(0, 4, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
			}
			
			Tlocal_pre.swap(&Tlocal); 
			iterLocal++;
		}
		if (iterLocal == iterL) {
			//std::cout << _iterGlobal << " " << iterLocal << " " << resL << " " << resG << std::endl;
		}
		//std::cout << "*";
#ifdef INSTRUMENTATION
		occurencePerBlock.increment(0, 1, iterLocal);
		occurencePerBlock.increment(0, 2, iterLocal);
		occurencePerBlock.increment(0, 3, iterLocal);
		occurencePerBlock.increment(0, 4, iterLocal/stepL);
#endif // INSTRUMENTATION
		Tlocal_pre.swap(&Tlocal);
		tradeLin.swap(&Tlocal);
		
#ifdef INSTRUMENTATION
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		updateGlobalProb();
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 5, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		// FB 4
		if (!(_iterGlobal % _stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateResBis(&resF, (_iterGlobal / _stepG), &tempNN);
			//std::cout << _iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, _iterGlobal / _stepG) << " " << resF.get(1, _iterGlobal / _stepG) << std::endl;
#ifdef INSTRUMENTATION
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 6, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		}
		
		_iterGlobal++;
	}
	//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, (iterGlobal - 1) / stepG) << " " << resF.get(1, (iterGlobal - 1) / stepG) << " " << resF.get(2, (iterGlobal - 1) / stepG) << std::endl;
#ifdef INSTRUMENTATION	
	occurencePerBlock.increment(0, 5, _iterGlobal);
	occurencePerBlock.increment(0, 6, _iterGlobal / _stepG);
	

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	int indice = 0;
	for (int idAgent = 0; idAgent < _nAgentTrue; idAgent++) {
		MatrixCPU omega(cas.getVoisin(idAgent));
		int Nvoisinmax = nVoisin.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			int idVoisin = omega.get(voisin, 0);
			trade.set(idAgent, idVoisin, tradeLin.get(indice, 0));
			LAMBDA.set(idAgent, idVoisin, LAMBDALin.get(indice, 0));
			indice = indice + 1;
		}
	}
	for (int idAgent = _nAgentTrue; idAgent < _nAgent; idAgent++) {
		for (int idVoisin = 0; idVoisin < _nAgentTrue; idVoisin++) {
			if (idVoisin != (idAgent - _nAgentTrue)) {
				trade.set(idAgent, idVoisin, tradeLin.get(indice, 0));
				LAMBDA.set(idAgent, idVoisin, LAMBDALin.get(indice, 0));
				indice = indice + 1;
			}
			
		}
	}
	updatePn(&Pn, &P, &nVoisin);
	/*trade.display();
	std::cout << "Trade" << std::endl;
	tradeLin.display();
	std::cout << "power" << std::endl;
	Tmoy.display();
	P.display();
	Pn.display();*/

	//std::cout << "lambda" << std::endl;
	//LAMBDALin.display();

	fc = calcFc(&a, &b, &tradeLin, &Pn, &Ct, &tempN1, &tempNN);
	// FB 5
	
	
	result->setResF(&resF);
	
	result->setLAMBDA(&LAMBDA);
	
	result->setTrade(&trade); 
	
	result->setIter(_iterGlobal);
	
	result->setMU(&MU);
	
	result->setPn(&Pn);
	
	result->setFc(fc);

#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 7, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 7, 1);

	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
	tMarket = clock() - tMarket;

	result->setTime((float)tMarket / CLOCKS_PER_SEC);
}

void ADMMMarket::solveWithMinPower(Simparam* result, const Simparam& sim, const StudyCase& cas)
{
	init(sim, cas);
	_id = _id + 1;
	int nCons = cas.getNCons();
	for (int n = 0; n < nCons; n++) {
		Pmin.set(n, 0, Pmax.get(n, 0));
	}
	for (int n = _nAgentTrue; n < _nAgent; n++) {
		Pmin.set(n, 0, 0);
		Pmax.set(n, 0, 0);
	}
	
	solve(result, sim, cas);
}

void ADMMMarket::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	Pmin = cas.getPmin();
	Pmax = cas.getPmax();


	MatrixCPU Lb(cas.getLb());

	b = cas.getb();
	Cp = b;
	int indice = 0;

	for (int idAgent = 0; idAgent < _nAgentTrue; idAgent++) {
		int Nvoisinmax = nVoisin.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			matLb.set(indice, 0, Lb.get(idAgent, 0));
			indice = indice + 1;
		}
	}
	for (int idAgent = _nAgentTrue; idAgent < _nAgent; idAgent++) {
		for (int voisin = 0; voisin < (_nAgent - 1); voisin++) {
			matLb.set(indice, 0, Lb.get(idAgent, 0));
			indice = indice + 1;
		}
	}

	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);
	Cp.multiplyT(&nVoisin);
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION

}

void ADMMMarket::init(const Simparam& sim, const StudyCase& cas)
{
	// intitilisation des matrixs et variables 
	
	clock_t t = clock();

	isAC = cas.isAC();
	//std::cout << "init " << std::endl;
	_rhog = sim.getRho();
	
	_iterG = sim.getIterG();
	_stepG = sim.getStepG();
	float epsG = sim.getEpsG();
	
	_nAgentTrue = sim.getNAgent();
	_nAgent = _nAgentTrue + isAC * _nAgentTrue;

	_rhol = _rho; //*nAgent
	
				  //std::cout << "rho " << _rho << std::endl;
	if (_rho == 0) {
		_rhol = _rhog;
	}

	nVoisin = cas.getNvoi();
	//nVoisin.display();
	_nTrade = nVoisin.sum();
	_nTradeP = 0;
	if (isAC) {
		for (int n = 0; n < _nAgentTrue; n++) {
			_nTradeP += nVoisin.get(n, 0);
		}
		_nTradeQ = _nTrade - _nTradeP;
		if (_nTradeQ != (_nAgentTrue * (_nAgentTrue - 1))) {
			std::cout << "err ADMM : " << _nAgent << " " << _nAgentTrue << " " << _nTrade << " " << _nTradeP << " " << _nTradeQ << std::endl;
			throw std::invalid_argument("Agent must be fully conected for the Q echanges, WIP");
		}
	}
	else {
		_nTradeP = _nTrade;
	}

	//std::cout << _iterG << " " << _stepG << std::endl;
	//std::cout << isAC << " " <<  _nAgentTrue << " " << _nAgent << " " << _nTrade << " " << _nTradeP << " " << _nTradeQ << std::endl;
	


	_at1 = _rhog; 
	_at2 = _rhol;

	resF = MatrixCPU(3, (_iterG / _stepG) + 1);

	MatrixCPU BETA(cas.getBeta());
	MatrixCPU Ub(cas.getUb());
	MatrixCPU Lb(cas.getLb());
	LAMBDA = sim.getLambda();
	trade = sim.getTrade();

	//std::cout << "mise sous forme linéaire" << std::endl;
	
	
	CoresMatLin = MatrixCPU(_nAgent, _nAgentTrue, -1);
	CoresAgentLin = MatrixCPU( _nAgent + 1, 1);
	CoresLinAgent = MatrixCPU(_nTrade, 1);
	CoresLinVoisin = MatrixCPU(_nTrade, 1);
	CoresLinTrans = MatrixCPU(_nTrade, 1);

	Tlocal_pre = MatrixCPU(_nTrade, 1);
	tradeLin = MatrixCPU(_nTrade, 1);
	LAMBDALin = MatrixCPU(_nTrade, 1);

	matLb = MatrixCPU(_nTrade, 1);
	matUb = MatrixCPU(_nTrade, 1);
	Ct = MatrixCPU(_nTrade, 1);
	Bt2 = MatrixCPU(_nTrade, 1);

	int indice = 0;
	//std::cout << " P " << std::endl;
	for (int idAgent = 0; idAgent < _nAgentTrue; idAgent++) { // P
		MatrixCPU omega(cas.getVoisin(idAgent));
		int Nvoisinmax = nVoisin.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			int idVoisin = omega.get(voisin, 0);
			matLb.set(indice, 0, Lb.get(idAgent, 0));
			matUb.set(indice, 0, Ub.get(idAgent, 0));
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
				matLb.set(indice, 0, Lb.get(idAgent, 0));
				matUb.set(indice, 0, Ub.get(idAgent, 0));
				//Ct.set(indice, 0, BETA.get(idAgent, idVoisin));
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
	for (int lin = 0; lin < _nTrade; lin++) {
		int i = CoresLinAgent.get(lin, 0);
		int j = CoresLinVoisin.get(lin, 0);
		if (lin >= _nTradeP) {
			i -= _nAgentTrue;
		}
				
		int k = CoresMatLin.get(j, i);
		CoresLinTrans.set(lin, 0, k);
	}


	/*std::cout << "trade bound" << std::endl;
	matLb.display();
    matUb.display();*/
	//CoresLinTrans.saveCSV("Cores.csv", 11, 1);


	//std::cout << "autres donnée sur CPU" << std::endl;
	tempNN = MatrixCPU(_nTrade, 1, 0);
	tempN1 = MatrixCPU(_nAgent, 1, 0); // plutôt que de re-allouer de la mémoire à chaque utilisation
	//MatrixCPU temp1N(1, _nAgent, 0, 1);

	Tlocal = MatrixCPU(_nTrade, 1, 0);
	
	
	
	
	Pn = sim.getPn(); // somme des trades
	P = Pn; // moyenne des trades, ici c'est juste pour qu'il ait la même taille sans avoir besoin de se poser de question

	a = MatrixCPU(cas.geta());
	b = MatrixCPU(cas.getb());
	Ap2 = a;
	Ap1 = nVoisin;
	Ap12 = MatrixCPU(_nAgent, 1, 0);

	Bt1 = MatrixCPU(_nTrade, 1, 0);
	Cp = b;
	Bp1 = MatrixCPU(_nAgent, 1, 0);

	Pmin = cas.getPmin();
	Pmax = cas.getPmax();
	MU = sim.getMU(); // facteur reduit i.e lambda_l/_rho
	Tmoy = sim.getPn();


	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);

	/*std::cout << "Power bound" << std::endl;
	Pmin.display();
	Pmax.display();*/

	Ap1.multiply(_rhol);
	Cp.multiplyT(&nVoisin);
	Tmoy.divideT(&nVoisin);
	
	Ap2.multiplyT(&nVoisin);
	Ap2.multiplyT(&nVoisin);
	Ap12.add(&Ap1, &Ap2);


	updateGlobalProb();
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	

	//std::cout << _at1 << " " << _at2 << std::endl;
	//Ap1.display();
	//Ap2.display();
	//Cp.display();
	//std::cout << "************" << std::endl;


}

void ADMMMarket::updateGlobalProb() {
	updateLambda();
	updateBt1();
}

void ADMMMarket::updateLocalProb() {
	// FB 1a
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	updateBt2();
	updateTl();
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 1, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 1b
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION	

	for (int i = 0; i < _nAgent; i++) {
		int nVoisinLocal = nVoisin.get(i, 0);
		int beginLocal = CoresAgentLin.get(i, 0);
		int endLocal = beginLocal + nVoisinLocal;
		float m = 0;
		for (int j = beginLocal; j < endLocal; j++) {
			m += Tlocal.get(j, 0);
		}
		Tmoy.set(i, 0, m / nVoisinLocal);
	}

	
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 2, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	
	// FB 1c
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION	
	updateBp1();
	updateP();
	updateMU();
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 3, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
}

void ADMMMarket::updateLambda()
{
	for (int t = 0; t < _nTrade; t++) {
		int k = CoresLinTrans.get(t, 0);
		float lamb = 0.5 * _rhog * (tradeLin.get(t, 0) + tradeLin.get(k, 0));
		LAMBDALin.set(t, 0, LAMBDALin.get(t, 0) + lamb);
	}
}


void ADMMMarket::updateBt1()
{
	
	// subtractTrans
	for (int t = 0; t < _nTrade; t++) {
		int k = CoresLinTrans.get(t,0);
		Bt1.set(t, 0, tradeLin.get(t, 0) - tradeLin.get(k, 0));
	}
	Bt1.multiply(0.5*_rhog); 
	Bt1.subtract(&LAMBDALin);
	Bt1.divide(_rhog);
}

void ADMMMarket::updateBt2()
{
	for (int i = 0; i < _nAgent; i++) {
		int nVoisinLocal = nVoisin.get(i, 0);
		int beginLocal = CoresAgentLin.get(i,0);
		int endLocal = beginLocal + nVoisinLocal; 
		for (int j = beginLocal; j < endLocal; j++) {
			float m = Tlocal_pre.get(j,0) - Tmoy.get(i, 0) + P.get(i, 0) - MU.get(i, 0); 
			Bt2.set(j, 0, m);
		}
	}
}

void ADMMMarket::updateBp1()
{
	Bp1.add(&MU, &Tmoy);
}

void ADMMMarket::updateTl()
{

	float ada = _at1 / _at2; 
	float apa = _at1 + _at2;

	Tlocal.set(&Bt1);
	Tlocal.multiply(ada);
	Tlocal.add(&Bt2);
	Tlocal.multiply(_at2);

	Tlocal.subtract(&Ct);
	Tlocal.divide(apa);
	Tlocal.project(&matLb, &matUb);
}

float ADMMMarket::calcRes()
{
	float d1 = Tlocal.max2(&Tlocal_pre);
	float d2 = Tmoy.max2(&P);
	

	return d1 * (d1 > d2) + d2 * (d2 >= d1);
}

float ADMMMarket::updateResBis(MatrixCPU* res, int iter, MatrixCPU* tempNN)
{
	//std::cout << "tradeLin" << std::endl;
	//tradeLin.display();
	
	for (int t = 0; t < _nTrade; t++) {
		int k = CoresLinTrans.get(t, 0);
		tempNN->set(t, 0, tradeLin.get(t, 0) + tradeLin.get(k, 0));
	}
	//std::cout << "tempNN" << std::endl;
	//tempNN->display();
	float resR = tempNN->max2();


	float resS = Tlocal.max2(&tradeLin);
	//std::cout << iter << " " << resR << " " << resS << std::endl;
	if (iter > 0) {
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
	}/**/
	
	
	res->set(0, iter, resR);
	res->set(1, iter, resS);
	
	return MAX(resS, resR);
}


void ADMMMarket::updateP()
{
	P.multiplyT(&Ap1, &Bp1);
	P.subtract(&Cp);
	
	P.divideT(&Ap12);
	P.project(&Pmin, &Pmax);
}

void ADMMMarket::updateMU()
{
	MU.add(&Tmoy);
	MU.subtract(&P);
}


void ADMMMarket::display() {

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
	std::cout << "The power error of this state is (constraint) " << resF.get(0, _iterGlobal / _stepG - (_iterGlobal>1)) << " and convergence " << resF.get(1, _iterGlobal / _stepG - (_iterGlobal > 1)) << std::endl;
	std::cout << "===============================================================|" << std::endl;
	std::cout << "      System Summary                                           |" << std::endl;
	std::cout << "===============================================================|" << std::endl;
	std::cout << "Agent            " << _nAgentTrue << std::endl;
	std::cout << "Nombre d'échange " << _nTrade << std::endl;

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
