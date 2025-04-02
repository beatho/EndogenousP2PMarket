#include "../head/ADMMConst.h"
 


ADMMConst::ADMMConst() : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " ADMMConst Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu
}


ADMMConst::ADMMConst(float rho) : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default ADMMConst Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu

}

ADMMConst::~ADMMConst()
{
}
void ADMMConst::setParam(float rho)
{
	_rho = rho;
}

void ADMMConst::setTau(float tau)
{
	if (tau < 1) {
		throw std::invalid_argument("tau must be greater than 1");
	}
	_tau = tau;
}



void ADMMConst::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
{
#ifdef DEBUG_SOLVE
	cas.display();
	sim.display(1);
#endif // DEBUG_SOLVE
	
	clock_t tall =clock();
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
		occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION
	}
	_rhog = sim.getRho();
	_at1 = _rhog;
		

	float fc = 0;

	int iterLocal = 0;
	float resG = 2 * _epsG;
	float resL = 2 * _epsL;
	int iterGlobal = 0;
	while ((iterGlobal < _iterG) && (resG>_epsG)) {
		resL = 2 * _epsL;
		iterLocal = 0;
#ifdef INSTRUMENTATION
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		while (iterLocal< _iterL && resL> _epsL) {
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
		//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;
		Tlocal_pre.swap(&Tlocal);
		tradeLin.swap(&Tlocal);
		
		updateGlobalProb();

		// FB 4
		if (!(iterGlobal % _stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateResEndo(iterGlobal / _stepG);
			//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG)
			//	<< " " << resF.get(1, iterGlobal / stepG) << " " << resF.get(2, iterGlobal / stepG) << std::endl;
#ifdef INSTRUMENTATION
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION		
		}
		iterGlobal++;
	}
	//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, (iterGlobal - 1) / stepG) << " " << resF.get(1, (iterGlobal - 1) / stepG) << " " << resF.get(2, (iterGlobal - 1) / stepG) << std::endl;
#ifdef INSTRUMENTATION
	occurencePerBlock.increment(0, 1, iterGlobal);
	occurencePerBlock.increment(0, 3, iterGlobal);
	occurencePerBlock.increment(0, 4, iterGlobal);
	occurencePerBlock.increment(0, 5, iterGlobal);
	occurencePerBlock.increment(0, 6, iterGlobal / stepG);

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	Kappa1.projectNeg(); //delta1
	Kappa2.projectNeg(); // delta2


	int indice = 0;
	for (int idAgent = 0; idAgent < _nAgent; idAgent++) {
		MatrixCPU omega(cas.getVoisin(idAgent));
		int Nvoisinmax = nVoisin.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			int idVoisin = omega.get(voisin, 0);
			trade.set(idAgent, idVoisin, tradeLin.get(indice, 0));
			LAMBDA.set(idAgent, idVoisin, LAMBDALin.get(indice, 0));
			indice = indice + 1;
		}
	}
	
	fc = calcFc();
	// FB 5
	if (cas.isAC()) {
		MatrixCPU tradeTot(2 * _nAgent, _nAgent);
		MatrixCPU LAMBDATot(2 * _nAgent, _nAgent);
		MatrixCPU PnTot(2 * _nAgent, 1);
		MatrixCPU MUTot(2 * _nAgent, 1);

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
	else {
		
		result->setLAMBDA(&LAMBDA);
		result->setTrade(&trade);
		result->setMU(&MU);
		result->setPn(&Pn);
	}


	result->setResF(&resF);
	if (_nLine) {
		result->setDelta1(&Kappa1);
		result->setDelta2(&Kappa2);
	}
	result->setIter(iterGlobal);
	
	
	
	result->setFc(fc);
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 7, 1);
	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
	tall = clock() - tall;
	result->setTime((float)tall / CLOCKS_PER_SEC);
	
}

void ADMMConst::init(const Simparam& sim, const StudyCase& cas)
{
	// intitilisation des matrixs et variables 
	
	//std::cout << "init " << std::endl;
	isAC = false;
	initSimParam(sim);
	initSize(cas);

	//std::cout << "mise sous forme lin�aire" << std::endl;
	initLinForm(sim, cas);
	


	// transfert des mises lineaire
	

	//std::cout << "donnees sur CPU pour le grid" << std::endl;
	
	//std::cout << "autres donn�e sur CPU" << std::endl;
		
	initCaseParam(sim, cas);
	initDCEndoMarket();


	updateGlobalProb();
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	//std::cout << "fin init " << std::endl;

}

void ADMMConst::updateGlobalProb() {
	
	// FB 3a
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	

	updatePn();
	

#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	timePerBlock.increment(0, 3, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 3b
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	alpha.multiplyTVector(&G, &Pn, 0);
	updateQ();
	
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 4, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	
	Kappa1_pre.set(&Kappa1);
	Kappa2_pre.set(&Kappa2);
	updateKappa();
	updateCp2();
	Cp.add(&Cp1, &Cp2);
	// FB 3c
	updateLambda();
	updateBt1();
	
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 5, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
	
}

void ADMMConst::updateLocalProb() {
	// FB 1a

	updateBt2();
	updateTl();

	
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

	updateBp1();
	updateP();
	updateMU();

	
}




void ADMMConst::updateBt1()
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

void ADMMConst::updateBt2()
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

void ADMMConst::updateBp1()
{
	Bp1.add(&MU, &Tmoy);
}

void ADMMConst::updateTl()
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


void ADMMConst::updateP()
{
	P.multiplyT(&Ap1, &Bp1);
	P.subtract(&Cp);
	
	P.divideT(&Ap12);
	P.project(&Pmin, &Pmax);
}

void ADMMConst::updateMU()
{
	MU.add(&Tmoy);
	MU.subtract(&P);
}

void ADMMConst::updateQ()
{
	for (int l = 0; l < _nLine; l++) {
		float qt = 0;
		for (int n = _nAgent - 1; n >= 0; n--) {
			qt += alpha.get(l, n);
			if (n > 0) {
				Qpart.set(l, n - 1, qt);
			}
		}
		Qtot.set(l, 0, qt);
	}
}
