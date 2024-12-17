#include "../head/MarEndoCons.h"
#define MAX(X, Y) X * (X >= Y) + Y * (Y > X)




MarEndoCons::MarEndoCons() : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " MarEndoCons Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilisé pendant la simu
}


MarEndoCons::MarEndoCons(float rho) : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default MarEndoCons Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
	timePerBlock = MatrixCPU(1, 8, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 8, 0); //nb de fois utilisé pendant la simu
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
	
	clock_t tall = clock();
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
		timePerBlock.set(0, 0, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		occurencePerBlock.set(0, 0, 1);
#endif // INSTRUMENTATION
	}
	_rhog = sim.getRho();
	_at1 = _rhog;
	_iterG = sim.getIterG();
	int iterL = sim.getIterL();
	int stepL = sim.getStepL();

	
	float epsL = sim.getEpsL();
	float epsG = sim.getEpsG();
	

	float fc = 0;

	int iterLocal = 0;
	_resG = 2 * epsG;
	float resL = 2 * epsL;
	_iterGlobal = 0;
	while ((_iterGlobal < _iterG) && (_resG>epsG) || (_iterGlobal <=_stepG)) { // || (_iterGlobal <= _stepG)
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
		
		
		resL = 2 * epsL;
		iterLocal = 0;

#ifdef INSTRUMENTATION
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		while (iterLocal< iterL && resL>epsL) {
			
			updateLocalProb();
			// FB 3
			if (!(iterLocal % stepL)) {

				resL = calcRes();
			
			}
			Tlocal_pre.swap(&Tlocal); 
			iterLocal++;
		}

#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 1, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
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
			timePerBlock.increment(0, 6, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
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
	//trade.display();
	/*std::cout << "PSO, Pn" << std::endl;
	PSO.display();
	Pn.display();*/

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
	timeMarketEndo = clock() - tall;
	result->setTime((float) timeMarketEndo / CLOCKS_PER_SEC);


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
		int Nvoisinmax = nVoisin.get(idAgent, 0);
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
	timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION


}

void MarEndoCons::init(const Simparam& sim, const StudyCase& cas)
{
	DELETEB(OPF);
	// intitilisation des matrixs et variables 
	
	
	_rhog = sim.getRho();
	_rhoSO = _rhog;
	//_rhoSO = sim.getRho1();
	_iterG = sim.getIterG();
	_stepG = sim.getStepG();
	float epsG = sim.getEpsG();
	float epsGC = sim.getEpsGC();
	_stepL = sim.getStepL();
	_ratioEps = epsG / epsGC;
	_nAgentTrue = sim.getNAgent();
	_nAgent = 2 * _nAgentTrue;

	paramOPF = sim;
	paramOPF.setItG(sim.getIterL());

	_rhol = _rho; //*nAgent
	//std::cout << "rho " << _rho << std::endl;
	if (_rho == 0) {
		_rhol = _rhog;
	}

	nVoisin = cas.getNvoi();
	_nTrade = nVoisin.sum();
	_nTradeP = 0;
	
	for (int n = 0; n < _nAgentTrue; n++) {
		_nTradeP += nVoisin.get(n, 0);
	}
	_nTradeQ = _nTrade - _nTradeP;
	if (_nTradeQ != (_nAgentTrue * (_nAgentTrue - 1))) {
		std::cout << "err ADMM : " << _nAgent << " " << _nAgentTrue << " " << _nTrade << " " << _nTradeP << " " << _nTradeQ << std::endl;
		throw std::invalid_argument("Agent must be fully conected for the Q echanges, WIP");
	}
	if (initWithMarketClear) {
		ADMMMarket market;
		Simparam res(sim);
		market.solve(&res, sim, cas);
		//res.display();
		LAMBDA = res.getLambda();
		trade = res.getTrade();
		Pn = res.getPn();
	}
	else {
		LAMBDA = sim.getLambda();
		trade = sim.getTrade();
		Pn = sim.getPn(); // somme des trades
	}
	//Pn.display();
	paramOPF.setPn(&Pn);
	
	_at1 = _rhog; 
	_at2 = _rhol;

	resF = MatrixCPU(3, (_iterG / _stepG) + 1);
	

	MatrixCPU BETA(cas.getBeta());
	
	MatrixCPU Ub(cas.getUb());
	MatrixCPU Lb(cas.getLb());
	

	//std::cout << "mise sous forme linéaire" << std::endl;
	
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
	Bt2 = MatrixCPU(_nTrade, 1);

	int indice = 0;

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

	

	//std::cout << "donnees sur CPU pour le grid" << std::endl;

	Ap3 = nVoisin;
	Ap3.multiplyT(&nVoisin);
	Ap3.multiply(_rhoSO);
	PSO = MatrixCPU(_nAgent, 1);
	etaSO = MatrixCPU(_nAgent, 1);
	Bp3 = MatrixCPU(_nAgent, 1);
	_radial = cas.isRadial();

	if (_radial) {
		OPF = new OPFADMMCons;
	}
	else {
		OPF = new OPFPDIPM;
	}
	

	OPF->initConsensus(paramOPF, cas, _rhoSO);



	//std::cout << "autres donnée sur CPU" << std::endl;
	tempNN = MatrixCPU(_nTrade, 1, 0);
	tempN1 = MatrixCPU(_nAgent, 1, 0); // plutôt que de re-allouer de la mémoire à chaque utilisation
	

	Tlocal = MatrixCPU(_nTrade, 1, 0);
	P = Pn; // moyenne des trades
	P.divideT(&nVoisin);

	a = cas.geta();
	b = cas.getb();

	// on enn veut pas que l'agent des pertes consomme plus que nécessaire !!!
	//a.set(0, 0, 1);
	//a.set(_nAgentTrue, 0, 1);


	Ap2 = a;
	Ap1 = nVoisin;
	Ap123 = Ap3;

	Bt1 = MatrixCPU(_nTrade, 1, 0);
	Cp = b;
	
	Bp1 = MatrixCPU(_nAgent, 1, 0);

	Pmin = cas.getPmin();
	Pmin.set(0, 0, -100000); // unleash power !!!
	Pmin.set(_nAgentTrue, 0, -100000); // unleash power !!!	
	Pmax = cas.getPmax();
	Pmax.set(_nAgentTrue, 0, 100000); // unleash power !!!
	
	MU = sim.getMU(); // facteur reduit i.e lambda_l/_rho
	Tmoy = sim.getPn();


	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);
	Ap1.multiply(_rhol);
	Cp.multiplyT(&nVoisin);
	Tmoy.divideT(&nVoisin);
	
	Ap2.multiplyT(&nVoisin);
	Ap2.multiplyT(&nVoisin);
	Ap123.add(&Ap1);
	Ap123.add(&Ap2);

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
	
	float eps = min(_resG / 10, 0.01);
	

	//std::cout << "SolveOPF" << std::endl;
	if (_iterGlobal % _stepL == 0) {
		OPF->solveConsensus(eps, &PSO);/**/
	
	/*float Ploss = OPF->getPLoss();
	float Qloss = OPF->getQLoss();
	PSO.set(0, 0, Ploss);
	PSO.set(_nAgentTrue, 0, Qloss);
	PSO.display();*/

#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();

	timePerBlock.increment(0, 4, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 3b
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	//std::cout << "update OPF" << std::endl;
	
	

		updatePn(&Pn,&P,&nVoisin);
		//Pn.display();
		//PSO.display();
		OPF->updateConsensus(&Pn);

		/*for (int n = 0; n < _nAgent; n++) {
			PSO.set(n, 0, (PSO.get(n, 0) + Pn.get(n, 0)) / 2 + etaSO.get(n, 0));
		}*/
	
	

#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 3, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	
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
	/**/
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
	timePerBlock.increment(0, 5, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
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
		Tmoy.set(i, 0, m/nVoisinLocal);
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

void MarEndoCons::updateLambda()
{
	for (int t = 0; t < _nTrade; t++) {
		int k = CoresLinTrans.get(t, 0);
		float lamb = 0.5 * _rhog * (tradeLin.get(t, 0) + tradeLin.get(k, 0));
		LAMBDALin.set(t, 0, LAMBDALin.get(t,0)+lamb);
	}
}

void MarEndoCons::updateEtaSO()
{
	for (int n = 0; n < _nAgent; n++) { 
		float eta = 0.5 * (Pn.get(n, 0) - PSO.get(n, 0));
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
	Bt1.multiply(0.5*_rhog); 
	Bt1.subtract(&LAMBDALin);
	Bt1.divide(_rhog);
}

void MarEndoCons::updateBt2()
{
	for (int i = 0; i < _nAgent; i++) {
		int nVoisinLocal = nVoisin.get(i,0);
		int beginLocal = CoresAgentLin.get(i,0);
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

	return MAX(d1, d2);
}

float MarEndoCons::updateResBis(MatrixCPU* res, int iter, MatrixCPU* tempNN)
{
	
	//tradeLin.display();
	for (int t = 0; t < _nTrade; t++) {
		int k = CoresLinTrans.get(t, 0);
		tempNN->set(t, 0, tradeLin.get(t, 0) + tradeLin.get(k, 0));
	}
	float resR = tempNN->max2();

	float resS = tradeLin.distance2(&Tlocal);
	float resXf = PSO.max2(&Pn) * _ratioEps;
	
	/*for (int i = 1; i < _nAgentTrue; i++) {
		resXf = MAX(abs(PSO.get(i,0) - Pn.get(i,0)), resXf);
		resXf = MAX(abs(PSO.get(i + _nAgentTrue, 0) - Pn.get(i + _nAgentTrue, 0)), resXf);
	}*/

	
	
	res->set(0, iter, resR);
	res->set(1, iter, resS);
	res->set(2, iter, resXf);
	return MAX(MAX(resXf, resS), resR);
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
		std::cout << "Converged in " << (float) timeMarketEndo / CLOCKS_PER_SEC << " seconds" << std::endl;

	}
	else {
		std::cout << "method " << _name << " not converged in " << _iterGlobal << " iterations." << std::endl;
		std::cout << "time taken " << (float) timeMarketEndo / CLOCKS_PER_SEC << " seconds" << std::endl;
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
