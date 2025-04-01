#include "../head/ADMMMarketOpenMP.h"
 


ADMMMarketOpenMP::ADMMMarketOpenMP() : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " ADMMMarketOpenMP Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
}


ADMMMarketOpenMP::ADMMMarketOpenMP(float rho) : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default ADMMMarketOpenMP Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
}

ADMMMarketOpenMP::~ADMMMarketOpenMP()
{
}
void ADMMMarketOpenMP::setParam(float rho)
{
	_rho = rho;
}

void ADMMMarketOpenMP::setTau(float tau)
{
	if (tau < 1) {
		throw std::invalid_argument("tau must be greater than 1");
	}
	_tau = tau;
}

void ADMMMarketOpenMP::solveWithMinPower(Simparam* result, const Simparam& sim, const StudyCase& cas)
{
	init(sim, cas);
	int nCons = cas.getNCons();

#pragma omp parallel for
	for (int n = 0; n < nCons; n++) {
		Pmax.set(n, 0, Pmin.get(n, 0));
	}

	solve(result, sim, cas);
}

void ADMMMarketOpenMP::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
		timePerBlock.increment(0, 0, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
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
	float resG = 2 * epsG;
	float resL = 2 * epsL;
	_iterGlobal = 0;

	while ((_iterGlobal < _iterG) && (resG>epsG)) {
		resL = 2 * epsL;
		iterLocal = 0;
		while (iterLocal< iterL && resL>epsL) {
			// FB 1a
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

			#pragma omp parallel for
			for (int i = 0; i < _nAgent; i++) {
				int nVoisinLocal = nVoisin.get(i, 0);
				int beginLocal = CoresAgentLin.get(i, 0);
				int endLocal = beginLocal + nVoisinLocal;
				float moy = 0;
				for (int j = beginLocal; j < endLocal; j++) {
					float m = Tlocal_pre.get(j, 0) - Tmoy.get(i, 0) + P.get(i, 0) - MU.get(i, 0);
					float t = (_at1 * Bt1.get(j, 0) + m * _at2 - Ct.get(j, 0))/ ( _at1 + _at2);
					float ub = matUb.get(j, 0);
					float lb = matLb.get(j, 0);
					t = t + (ub - t) * (t > ub) + (lb - t) * (lb > t);
					Tlocal.set(j, 0, t);
					moy += t;
					
				}
				moy /= nVoisinLocal;
				Tmoy.set(i, 0, moy);
				float bp1 = MU.get(i, 0) + moy;//Bp1.add(&MU, &Tmoy);
				float p = (Ap1.get(i, 0) * bp1 - Cp.get(i, 0)) / Ap12.get(i, 0);
				float ub = Pmax.get(i, 0);
				float lb = Pmin.get(i, 0);
				p = p + (ub - p) * (ub < p) + (lb - p) * (lb > p);
				P.set(i, 0, p);
				float mu = MU.get(i, 0) - p + moy;
				MU.set(i, 0, mu);

			}
			
#ifdef INSTRUMENTATION
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 1, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION	

		
			// FB 2
			if (!(iterLocal % stepL)) {
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
		if (iterLocal == iterL) {
			//std::cout << _iterGlobal << " " << iterLocal << " " << resL << " " << resG << std::endl;
		}
#ifdef INSTRUMENTATION
		occurencePerBlock.increment(0, 1, iterLocal);
		occurencePerBlock.increment(0, 4, iterLocal/stepL);
#endif // INSTRUMENTATION
		Tlocal_pre.swap(&Tlocal);
		tradeLin.swap(&Tlocal);
		
#ifdef INSTRUMENTATION
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		#pragma omp parallel for
		for (int t = 0; t < _nTrade; t++) {
			int k = CoresLinTrans.get(t, 0);
			float lamb = LAMBDALin.get(t, 0) + 0.5 * _rhog * (tradeLin.get(t, 0) + tradeLin.get(k, 0));
			float bt = 0.5 * (tradeLin.get(t, 0) - tradeLin.get(k, 0)) - lamb / _rhog;
			LAMBDALin.set(t, 0, lamb);
			Bt1.set(t, 0, bt);
		}
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 5, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		// FB 4
		if (!(_iterGlobal % _stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateRes(_iterGlobal / _stepG);
#ifdef INSTRUMENTATION
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		}
		//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;

		_iterGlobal++;
	}
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
	/*std::cout << "Trade" << std::endl;
	tradeLin.display();
	std::cout << "power" << std::endl;
	Tmoy.display();
	P.display();
	Pn.display();*/

	fc = calcFc();
	// FB 5
	//std::cout << _iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, (_iterGlobal - 1) / _stepG) << " " << resF.get(1, (_iterGlobal - 1) / _stepG) << " " << resF.get(2, (_iterGlobal - 1) / _stepG) << " " << fc << std::endl;

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

void ADMMMarketOpenMP::updateP0(const StudyCase& cas)
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

#pragma omp parallel for
	for (int idAgent = 0; idAgent < _nAgent; idAgent++) {
		Pmin.set(idAgent, 0, Pmin.get(idAgent, 0) / nVoisin.get(idAgent, 0));
		Pmax.set(idAgent, 0, Pmax.get(idAgent, 0) / nVoisin.get(idAgent, 0));
		Cp.set(idAgent, 0, Cp.get(idAgent, 0) * nVoisin.get(idAgent, 0));

	}


#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION

}

void ADMMMarketOpenMP::init(const Simparam& sim, const StudyCase& cas)
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
	_nTrade = nVoisin.sum();
	_nTradeP = 0;
	int testReduc = 0;
	if (isAC) {	
		#pragma omp parallel for reduction(+ : testReduc)
		for (int n = 0; n < _nAgentTrue; n++) {
			testReduc += nVoisin.get(n, 0);
		}
		_nTradeP = testReduc;
		_nTradeQ = _nTrade - _nTradeP;
		if (_nTradeQ != (_nAgentTrue * (_nAgentTrue - 1))) {
			std::cout << "err OpenMp : " << _nAgent << " " << _nAgentTrue << " " << _nTrade << " " << _nTradeP << " " << _nTradeQ << std::endl;

			throw std::invalid_argument("Agent must be fully conected for the Q echanges, WIP");
		}
	}
	else {
		_nTradeP = _nTrade;
	}
	//std::cout << isAC << " " <<  _nAgentTrue << " " << _nAgent << " " << _nTrade << " " << _nTradeP << " " << _nTradeQ << std::endl;
	


	_at1 = _rhog; // represente en fait 2*a
	_at2 = _rhol;

	resF = MatrixCPU(3, (_iterG / _stepG) + 1);

	MatrixCPU BETA(cas.getBeta());
	MatrixCPU Ub(cas.getUb());
	MatrixCPU Lb(cas.getLb());
	LAMBDA = sim.getLambda();
	trade = sim.getTrade();

	//std::cout << "mise sous forme lin�aire" << std::endl;
	
	
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
				matLb.set(indice, 0, Lb.get(idAgent, 0));
				matUb.set(indice, 0, Ub.get(idAgent, 0));
				//Ct.set(indice, 0, BETA.get(idAgent, idVoisin));
				tradeLin.set(indice, 0, trade.get(idAgent, idVoisin));
				Tlocal_pre.set(indice, 0, trade.get(idAgent, idVoisin));
				LAMBDALin.set(indice, 0, LAMBDA.get(idAgent, idVoisin));
				CoresLinAgent.set(indice, 0, idAgent );
				CoresLinVoisin.set(indice, 0, idVoisin + _nAgentTrue);
				CoresMatLin.set(idAgent, idVoisin, indice);
				indice = indice + 1;
			}
		}
		CoresAgentLin.set(idAgent + 1, 0, indice);
	}
	#pragma omp parallel for
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

	

	//std::cout << "autres donn�e sur CPU" << std::endl;
	tempNN = MatrixCPU(_nTrade, 1, 0);
	tempN1 = MatrixCPU(_nAgent, 1, 0); // plut�t que de re-allouer de la m�moire � chaque utilisation
	//MatrixCPU temp1N(1, _nAgent, 0, 1);

	Tlocal = MatrixCPU(_nTrade, 1, 0);
	
	
	
	
	Pn = sim.getPn(); // somme des trades
	P = Pn; // moyenne des trades, ici c'est juste pour qu'il ait la m�me taille sans avoir besoin de se poser de question

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


#pragma omp parallel for
	for (int t = 0; t < _nTrade; t++) {
		int k = CoresLinTrans.get(t, 0);
		float lamb = LAMBDALin.get(t, 0) + 0.5 * _rhog * (tradeLin.get(t, 0) + tradeLin.get(k, 0));
		float bt = 0.5 * (tradeLin.get(t, 0) - tradeLin.get(k, 0)) - lamb / _rhog;
		LAMBDALin.set(t, 0, lamb);
		Bt1.set(t, 0, bt);
	}
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	

}


void ADMMMarketOpenMP::display() {

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
	std::cout << "Nombre d'�change " << _nTrade << std::endl;

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
