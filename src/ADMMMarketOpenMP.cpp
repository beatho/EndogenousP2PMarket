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

	float fc = 0;

	int iterLocal = 0;
	_resG = 2 * _epsG;
	float resL = 2 * _epsL;
	_iterGlobal = 0;

	while ((_iterGlobal < _iterG) && (_resG > _epsG)) {
		resL = 2 * _epsL;
		iterLocal = 0;
		while (iterLocal< _iterL && resL> _epsL) {
			// FB 1a
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

			#pragma omp parallel for
			for (int i = 0; i < _nAgent; i++) {
				int nVoisinLocal = (int) nVoisin.get(i, 0);
				int beginLocal = (int) CoresAgentLin.get(i, 0);
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
			//std::cout << _iterGlobal << " " << iterLocal << " " << resL << " " << _resG << std::endl;
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
			int k = (int) CoresLinTrans.get(t, 0);
			float lamb = LAMBDALin.get(t, 0) + 0.5f * _rhog * (tradeLin.get(t, 0) + tradeLin.get(k, 0));
			float bt = 0.5f * (tradeLin.get(t, 0) - tradeLin.get(k, 0)) - lamb / _rhog;
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
			_resG = updateRes(_iterGlobal / _stepG);
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
	
	/*std::cout << "Trade" << std::endl;
	tradeLin.display();
	std::cout << "power" << std::endl;
	Tmoy.display();
	P.display();
	Pn.display();*/

	
	// FB 5
	//std::cout << _iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, (_iterGlobal - 1) / _stepG) << " " << resF.get(1, (_iterGlobal - 1) / _stepG) << " " << resF.get(2, (_iterGlobal - 1) / _stepG) << " " << fc << std::endl;

	setResult(result, cas.isAC());

#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 7, 1);

	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
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
		int Nvoisinmax = (int) nVoisin.get(idAgent, 0);
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
	initSize(cas);
	initSimParam(sim);

	//std::cout << isAC << " " <<  _nAgentTrue << " " << _nAgent << " " << _nTrade << " " << _nTradeP << " " << _nTradeQ << std::endl;
	
	_at1 = _rhog; // represente en fait 2*a
	_at2 = _rhol;

	
	initCaseParam(sim, cas);
	//std::cout << "mise sous forme lineaire" << std::endl;
	initLinForm(cas);
	//std::cout << "autres donnee sur CPU" << std::endl;
	
	initP2PMarket();

	#pragma omp parallel for
	for (int t = 0; t < _nTrade; t++) {
		int k = (int) CoresLinTrans.get(t, 0);
		float lamb = LAMBDALin.get(t, 0) + 0.5f * _rhog * (tradeLin.get(t, 0) + tradeLin.get(k, 0));
		float bt = 0.5f * (tradeLin.get(t, 0) - tradeLin.get(k, 0)) - lamb / _rhog;
		LAMBDALin.set(t, 0, lamb);
		Bt1.set(t, 0, bt);
	}
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	

}



