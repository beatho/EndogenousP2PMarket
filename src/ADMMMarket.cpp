#include "../head/ADMMMarket.h"
 


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
		timePerBlock.increment(0, 0, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION
	}
	
	_rhog = sim.getRho();
	_at1 = _rhog;

	
	float fc = 0;

	int iterLocal = 0;
	_iterGlobal = 0;
	float resG = 2 * _epsG;
	float resL = 2 * _epsL;

	while ((_iterGlobal < _iterG) && (resG>_epsG)) {
		/*Tlocal.display();
		P.saveCSV("testPCPU2.csv", 11, 1);
		Tlocal.saveCSV("testTCPU2.csv", 11, 1);
		Bt1.saveCSV("testBCPU2.csv", 11, 1);*/
		//P.display();
		//std::cout << "lambda" << std::endl;
		//LAMBDALin.display();
		resL = 2 * _epsL;
		iterLocal = 0;
		while (iterLocal< _iterL && resL>_epsL) {
			updateLocalProb();
			// FB 2
			if (!(iterLocal % _stepL)) {
#ifdef INSTRUMENTATION
				t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
				resL = calcRes();
#ifdef INSTRUMENTATION
				t2 = std::chrono::high_resolution_clock::now();
				timePerBlock.increment(0, 4, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
			}
			
			Tlocal_pre.swap(&Tlocal); 
			iterLocal++;
		}
		if (iterLocal == _iterL) {
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
		timePerBlock.increment(0, 5, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		// FB 4
		if (!(_iterGlobal % _stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			
			resG = updateRes((_iterGlobal / _stepG));
			
			//std::cout << _iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, _iterGlobal / _stepG) << " " << resF.get(1, _iterGlobal / _stepG) << std::endl;
#ifdef INSTRUMENTATION
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
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
	updatePn();
	/*trade.display();
	std::cout << "Trade" << std::endl;
	tradeLin.display();
	std::cout << "power" << std::endl;
	Tmoy.display();
	P.display();
	Pn.display();*/

	//std::cout << "lambda" << std::endl;
	//LAMBDALin.display();

	fc = calcFc();
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
	timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
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

void ADMMMarket::init(const Simparam& sim, const StudyCase& cas)
{
	// intitilisation des matrixs et variables 
	
	clock_t t = clock();
	isAC = cas.isAC();
	initSize(cas);
	initSimParam(sim);

	//std::cout << "mise sous forme lin�aire" << std::endl;
	initLinForm(sim, cas);
	
	//std::cout << "autres donnee sur CPU" << std::endl;
	initCaseParam(sim, cas);
	
	initP2PMarket();

	updateGlobalProb();
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	
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
	timePerBlock.increment(0, 1, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

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
	timePerBlock.increment(0, 2, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	
	// FB 1c
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION	
	updateBp1();
	updateP();
	updateMU();
	
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 3, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
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
