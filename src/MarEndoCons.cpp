#include "../head/MarEndoCons.h"



MarEndoCons::MarEndoCons() : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " MarEndoCons Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu
}


MarEndoCons::MarEndoCons(float rho) : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default MarEndoCons Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu
}

MarEndoCons::~MarEndoCons()
{
	DELETEB(OPF);
}


void MarEndoCons::setParam(float rho)
{
	_rho = rho;
}

void MarEndoCons::setTau(float tau)
{
	if (tau < 1) {
		throw std::invalid_argument("tau must be greater than 1");
	}
	_tau = tau;
}

void MarEndoCons::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
{
#ifdef DEBUG_SOLVE
	cas.display();
	sim.display(1);
#endif // DEBUG_SOLVE
	
	tMarket = clock();
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
		timePerBlock.set(0, 0, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		occurencePerBlock.set(0, 0, 1);
#endif // INSTRUMENTATION
	}
	_rhog = sim.getRho();
	_at1 = _rhog;
	
	float fc = 0;

	int iterLocal = 0;
	_resG = 2 * _epsG;
	float resL = 2 * _epsL;
	_iterGlobal = 0;
	while ((_iterGlobal < _iterG) && (_resG>_epsG) || (_iterGlobal <=_stepG)) { // || (_iterGlobal <= _stepG)
		/*std::cout << "---------------------------------" << std::endl;
		std::cout << " Pn " << std::endl;
		Pn.display();
		std::cout << " Pso " << std::endl;
		PSO.display();
		std::cout << " Bp3 " << std::endl;
		Bp3.display();*/

		/*std::cout << " Tlocal " << std::endl;
		Tlocal.display();
		std::cout << " Bt1 " << std::endl;
		Bt1.display();
		std::cout << " matlb " << std::endl;
		matLb.display();
		std::cout << " matUb " << std::endl;
		matUb.display();*/
		
		
		resL = 2 * _epsL;
		iterLocal = 0;

#ifdef INSTRUMENTATION
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		while (iterLocal< _iterL && resL>_epsL) {
			
			updateLocalProb();
			// FB 3
			if (!(iterLocal % _stepL)) {

				resL = calcRes();
			
			}
			Tlocal_pre.swap(&Tlocal); 
			iterLocal++;
		}

#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 1, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION

		//std::cout << _iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, _iterGlobal / _stepG) << " " << resF.get(1, _iterGlobal / _stepG) << std::endl;
		Tlocal_pre.swap(&Tlocal);
		tradeLin.swap(&Tlocal);
		
		updateGlobalProb();
		
		// FB 4
		if (!(_iterGlobal % _stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			_resG = updateResBis(&resF, (_iterGlobal / _stepG), &tempNN);
			//std::cout << _iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, _iterGlobal / _stepG)
				//<< " " << resF.get(1, _iterGlobal / _stepG) << " " << resF.get(2, _iterGlobal / _stepG) << std::endl;
#ifdef INSTRUMENTATION
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION		
		}

		_iterGlobal++;
	}
	//std::cout << _iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, (_iterGlobal - 1) / _stepG) << " " << resF.get(1, (_iterGlobal - 1) / _stepG) << " " << resF.get(2, (_iterGlobal - 1) / _stepG) << std::endl;
#ifdef INSTRUMENTATION
	occurencePerBlock.increment(0, 1, _iterGlobal);
	occurencePerBlock.increment(0, 3, _iterGlobal);
	occurencePerBlock.increment(0, 4, _iterGlobal);
	occurencePerBlock.increment(0, 5, _iterGlobal);
	occurencePerBlock.increment(0, 6, _iterGlobal / _stepG);

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	//trade.display();
	/*std::cout << "PSO, Pn" << std::endl;
	PSO.display();
	Pn.display();*/

	// FB 5
	MatrixCPU Pb(OPF->getPb());
	MatrixCPU Phi(OPF->getPhi());
	MatrixCPU E(OPF->getE());

	/*Pb.display();
	Phi.display();
	E.display();*/
	
	result->setE(&E);
	result->setPhi(&Phi);
	result->setPb(&Pb);

	setResult(result, cas.isAC());
	
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 7, 1);
	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
	
	//std::cout << "****" << std::endl;
	//OPF->display();
	//display();
}

void MarEndoCons::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	Pmin = cas.getPmin();
	Pmax = cas.getPmax();

	// unleash powe
	Pmin.set(0, 0, -POWERLIMIT);
	Pmax.set(_nAgentTrue, 0, POWERLIMIT);
	Pmin.set(_nAgentTrue, 0, -POWERLIMIT);
	Cp = cas.getb();

	MatrixCPU Lb(cas.getLb());
	MatrixCPU Ub(cas.getUb());

	int indice = 0;

	for (int idAgent = 0; idAgent < _nAgent; idAgent++) {
		int Nvoisinmax = (int) nVoisin.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			matLb.set(indice, 0, Lb.get(idAgent, 0));
			matUb.set(indice, 0, Ub.get(idAgent, 0));
			indice = indice + 1;
		}
	}

	float Ploss = Pn.get(0, 0);
	float Qloss = Pn.get(_nAgentTrue, 0);

	// pour essayer que cela marche
	Pn.add(&Pmin, &Pmax);
	Pn.divide(2);
	Pn.set(0, 0, Ploss);
	Pn.set(_nAgentTrue, 0, Qloss);

	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);
	Cp.multiplyT(&nVoisin);
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION


}

void MarEndoCons::init(const Simparam& sim, const StudyCase& cas)
{
	DELETEB(OPF);
	// initilisation des matrixs et variables 
	isAC = true;
	initSize(cas);


	initSimParam(sim);
	
	_rhoSO = _rhog;//_rhoSO = _rho1;
	
	//std::cout << "precision demandee " << epsG << " " << epsGC << " ratio " << _ratioEps << std::endl;
	paramOPF = sim;
	paramOPF.setItG(_iterIntern);
	paramOPF.setEpsG(_epsIntern);
	

	//std::cout << "Market" << std::endl;
	//std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	initCaseParam(sim, cas);
	Pmin.set(0, 0, -100000); // unleash power !!!
	Pmin.set(_nAgentTrue, 0, -100000); // unleash power !!!	
	Pmax.set(_nAgentTrue, 0, 100000); // unleash power !!!
	
	if (initWithMarketClear) {
		ADMMMarket market;
		Simparam res(sim);
		market.solve(&res, sim, cas);
		//res.display();
		LAMBDA = res.getLambda();
		trade = res.getTrade();
		Pn = res.getPn();
	}
	
	//std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	//std::cout << "time : " << (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / 1000000000 << std::endl;
	
	//Pn.display();
	paramOPF.setPn(&Pn);
	
	//std::cout << "mise sous forme lin�aire" << std::endl;
	initLinForm(cas);
	
	//std::cout << "donnees sur CPU pour le grid" << std::endl;
	
	PSO = MatrixCPU(_nAgent, 1);
	etaSO = MatrixCPU(_nAgent, 1);
	Bp3 = MatrixCPU(_nAgent, 1);
	_radial = cas.isRadial();

	if (_radial) {
		OPF = new OPFADMMCons;
	}
	else {
		throw std::invalid_argument("WIP : must be a radial grid");
	}
	
	OPF->initConsensus(paramOPF, cas, _rhoSO);

	//std::cout << "autres donnee sur CPU" << std::endl;
	initP2PMarket();
	Ap3 = nVoisin;
	Ap3.multiplyT(&nVoisin);
	Ap3.multiply(_rhoSO);
	// on en veut pas que l'agent des pertes consomme plus que n�cessaire !!!
	//a.set(0, 0, 1);
	//a.set(_nAgentTrue, 0, 1);

	
	Ap123 = Ap3;
	Ap123.add(&Ap12);
	

	/*std::cout << _at1 << " " << _at2 << std::endl;

	Ct.display();
	Ap1.display();
	Ap2.display();
	Ap3.display();
	Ap123.display();
	Cp.display();

	Pmin.display();
	Pmax.display();
	matLb.display();
	matUb.display();*/

	/*PSO.display();
	Pn.display();*/
	//std::cout << "******" << std::endl;
	//std::cout << "updateGlobal" << std::endl;
	updateGlobalProb();
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	//std::cout << "fin init " << std::endl;

}

void MarEndoCons::updateGlobalProb() {
	// FB 3a
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	std::chrono::high_resolution_clock::time_point t2;
#endif // INSTRUMENTATION
	
	float eps = MYMIN(_resG * _delta, _epsIntern);
	
	//std::cout << "SolveOPF" << std::endl;
	if (_iterGlobal % _stepIntern == 0) {
		OPF->solveConsensus(eps, &PSO);/**/
		/*float Ploss = OPF->getPLoss();
		float Qloss = OPF->getQLoss();
		PSO.set(0, 0, Ploss);
		PSO.set(_nAgentTrue, 0, Qloss);
		PSO.display();*/
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 4, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	// FB 3b
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	//std::cout << "update OPF" << std::endl;
		updatePn();
		//Pn.display();
		//PSO.display();
		OPF->updateConsensus(&Pn);
		/*for (int n = 0; n < _nAgent; n++) {
			PSO.set(n, 0, (PSO.get(n, 0) + Pn.get(n, 0)) / 2 + etaSO.get(n, 0));
		}*/
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 3, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		//Agent des pertes
		//std::cout << "update Market" << std::endl;
		/*Pmin.set(0, 0, Ploss / nVoisin.get(0,0));
		Pmax.set(0, 0, Ploss / nVoisin.get(0, 0));
		
		Pmin.set(_nAgentTrue, 0, Qloss / nVoisin.get(_nAgentTrue, 0));
		Pmax.set(_nAgentTrue, 0, Qloss / nVoisin.get(_nAgentTrue, 0));*/

		for (int j = 0; j < _nAgentTrue - 1; j++) {
			if (PSO.get(_nAgentTrue,0) > 0) {
				matUb.set(_nTradeP + j, 0, PSO.get(_nAgentTrue, 0));
				matLb.set(_nTradeP + j, 0, 0);
			}
			else {
				matLb.set(_nTradeP + j, 0, PSO.get(_nAgentTrue, 0));
				matUb.set(_nTradeP + j, 0, 0);
			}
		}
	// FB 3c
	updateEtaSO();
	updateBp3();
	}


	updateLambda();
	updateBt1();
	//Bp3.display();
	
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 5, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
	
}

void MarEndoCons::updateLocalProb() {
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
		int nVoisinLocal = (int) nVoisin.get(i, 0);
		int beginLocal = (int) CoresAgentLin.get(i, 0);
		int endLocal = beginLocal + nVoisinLocal;
		float m = 0;
		for (int j = beginLocal; j < endLocal; j++) {
			m += Tlocal.get(j, 0);
		}
		Tmoy.set(i, 0, m/nVoisinLocal);
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

void MarEndoCons::updateLambda()
{
	for (int t = 0; t < _nTrade; t++) {
		int k = (int) CoresLinTrans.get(t, 0);
		float lamb = 0.5f * _rhog * (tradeLin.get(t, 0) + tradeLin.get(k, 0));
		LAMBDALin.set(t, 0, LAMBDALin.get(t,0) + lamb);
	}
}

void MarEndoCons::updateEtaSO()
{
	for (int n = 0; n < _nAgent; n++) { 
		float eta = 0.5f * (Pn.get(n, 0) - PSO.get(n, 0));
		etaSO.set(n, 0, etaSO.get(n, 0) + eta);
	}
}



void MarEndoCons::updateBp3()
{
	Bp3.add(&PSO, &Pn);
	Bp3.multiply(0.5);
	Bp3.subtract(&etaSO);
	Bp3.divideT(&nVoisin);
}


void MarEndoCons::updateBt1()
{
	
	// subtractTrans
	for (int t = 0; t < _nTrade; t++) {
		int k = CoresLinTrans.get(t,0);
		Bt1.set(t, 0, tradeLin.get(t, 0) - tradeLin.get(k, 0));
	}
	Bt1.multiply(0.5f*_rhog); 
	Bt1.subtract(&LAMBDALin);
	Bt1.divide(_rhog);
}

void MarEndoCons::updateBt2()
{
	for (int i = 0; i < _nAgent; i++) {
		int nVoisinLocal = (int) nVoisin.get(i,0);
		int beginLocal = (int) CoresAgentLin.get(i,0);
		int endLocal = beginLocal + nVoisinLocal; 
		for (int j = beginLocal; j < endLocal; j++) {
			float m = Tlocal_pre.get(j,0) - Tmoy.get(i, 0) + P.get(i, 0) - MU.get(i, 0); 
			Bt2.set(j, 0, m);
		}
	}
}

void MarEndoCons::updateBp1()
{
	Bp1.add(&MU, &Tmoy);
}

void MarEndoCons::updateTl()
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

float MarEndoCons::calcRes()
{
	float d1 = Tlocal.max2(&Tlocal_pre);
	float d2 = Tmoy.max2(&P);

	return MYMAX(d1, d2);
}

float MarEndoCons::updateResBis(MatrixCPU* res, int iter, MatrixCPU* tempNN)
{
	
	//tradeLin.display();
	for (int t = 0; t < _nTrade; t++) {
		int k = (int) CoresLinTrans.get(t, 0);
		tempNN->set(t, 0, tradeLin.get(t, 0) + tradeLin.get(k, 0));
	}
	float resR = tempNN->max2();
	float resS = tradeLin.max2(&Tlocal);
	float resXf = PSO.max2(&Pn);
	
	/*for (int i = 1; i < _nAgentTrue; i++) {
		resXf = MYMAX(abs(PSO.get(i,0) - Pn.get(i,0)), resXf);
		resXf = MYMAX(abs(PSO.get(i + _nAgentTrue, 0) - Pn.get(i + _nAgentTrue, 0)), resXf);
	}*/

	res->set(0, iter, resR);
	res->set(1, iter, resS);
	res->set(2, iter, resXf);
	return MYMAX(MYMAX(resXf * _ratioEps, resS), resR);
}

void MarEndoCons::updateP()
{
	P.multiplyT(&Ap1, &Bp1);
	tempN1.multiplyT(&Bp3, &Ap3);
	P.add(&tempN1);
	P.subtract(&Cp);
	
	P.divideT(&Ap123);
	P.project(&Pmin, &Pmax);
}

void MarEndoCons::updateMU()
{
	MU.add(&Tmoy);
	MU.subtract(&P);
}



void MarEndoCons::display() {

	if (_iterGlobal == 0) {
		std::cout << "algorithm not launch" << std::endl;
	}
	else if (_iterGlobal < _iterG) {
		std::cout << "method " << _name << " converged in " << _iterGlobal << " iterations." << std::endl;
		std::cout << "Converged in " << (float) tMarket / CLOCKS_PER_SEC << " seconds" << std::endl;

	}
	else {
		std::cout << "method " << _name << " not converged in " << _iterGlobal << " iterations." << std::endl;
		std::cout << "time taken " << (float) tMarket / CLOCKS_PER_SEC << " seconds" << std::endl;
	}
	std::cout << "The power error of this state is (constraint) " << resF.get(0, _iterGlobal / _stepG) << " and convergence " << resF.get(1, _iterGlobal / _stepG) << std::endl;
	std::cout << "===============================================================|" << std::endl;
	std::cout << "      System Summary                                           |" << std::endl;
	std::cout << "===============================================================|" << std::endl;
	std::cout << "Agent            " << _nAgentTrue << std::endl;
	


	std::cout << std::endl << std::endl;
	std::cout << std::endl << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "      Agent Data                                                                                        |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Agent  |  Cost    |  Cost    |          Power Injection           |           Power Injection          |" << std::endl;
	std::cout << "  #     |   a (pu) |   b (pu) |  P (pu)  | Pmin (pu)  | Pmax (pu)  |  Q (pu)   | Qmin (pu)  | Qmax (pu) |" << std::endl;
	std::cout << "--------|----------|----------|----------|------------|------------|-----------|------------|-----------|" << std::endl;

	for (int n = 0; n < _nAgentTrue; n++) {
		
		std::cout << std::setw(8) << n << "|" << std::setw(9) << a.get(n, 0) << " |" << std::setw(10)
			<< b.get(n, 0) << "|" << std::setw(10) << Pn.get(n, 0) << "|" << std::setw(12)
			<< Pmin.get(n, 0) * nVoisin.get(n, 0) << "|" << std::setw(12) << Pmax.get(n, 0) * nVoisin.get(n, 0)
			<< "|" << std::setw(11) << Pn.get(n + _nAgentTrue, 0) << "|" << std::setw(12) << Pmin.get(n + _nAgentTrue, 0) * nVoisin.get(n + _nAgentTrue, 0)
			<< "|" << std::setw(11) << Pmax.get(n + _nAgentTrue, 0) * nVoisin.get(n + _nAgentTrue, 0) << "|" << std::endl;
	}


	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "                      END PRINT                                                                         |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
}
