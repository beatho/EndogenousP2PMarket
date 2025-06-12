#include "../head/EndoPF.h"
 


EndoPF::EndoPF() : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " EndoPF Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu
}


EndoPF::EndoPF(float rho) : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default EndoPF Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu
}

EndoPF::~EndoPF()
{
	DELETEB(pf);
}
void EndoPF::setParam(float rho)
{
	_rho = rho;
}

void EndoPF::setTau(float tau)
{
	if (tau < 1) {
		throw std::invalid_argument("tau must be greater than 1");
	}
	_tau = tau;
}



void EndoPF::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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

	std::cout << "solve de EndoPF " <<std::endl;
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
	
	int iterLocal = 0;
	float resG = 2 * _epsG;
	float resL = 2 * _epsL;
	_iterGlobal = 0;
	//Pn.display();
	//std::cout << "*******" << std::endl;
	while ((_iterGlobal < _iterG) && (resG>_epsG) || (_iterGlobal <= _stepG)) {
		
		resL = 2 * _epsL;
		iterLocal = 0;
		/*std::cout << "probl�me local" << std::endl;
		std::cout << _at1 << " " << _at2 << std::endl;
		Bt1.display();
		Ct.display();
		matLb.display();
		matUb.display();
		Ap1.display();
		Ap12.display();
		Bp1.display();
		Cp.display();
		Pmin.display();
		Pmax.display();*/

#ifdef INSTRUMENTATION
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		while (iterLocal< _iterL && resL>_epsL) {
			updateLocalProb();
			//FB 3
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
		/*td::cout << "********power**************" << std::endl;
		P.display();
		Pmin.display();
		Pmax.display();
		std::cout << _iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, _iterGlobal / _stepG) << " " << resF.get(1, _iterGlobal / _stepG) << " " << resF.get(2, (_iterGlobal - 1) / _stepG) << std::endl;
		*/
		Tlocal_pre.swap(&Tlocal);
		tradeLin.swap(&Tlocal);
		
		updateGlobalProb();

		// FB 4
		if (!(_iterGlobal % _stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateResEndo(_iterGlobal / _stepG);
			//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG)
			//	<< " " << resF.get(1, iterGlobal / stepG) << " " << resF.get(2, iterGlobal / stepG) << std::endl;
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
	
	setResult(result, cas.isAC());

	/*pf->display2(true);
	std::cout << " Y " << std::endl;
	Y.display();

	trade.display();
	std::cout << "Trade" << std::endl;
	tradeLin.display();
	matLb.display();
	matUb.display();
	std::cout << "power" << std::endl;	
	Tmoy.display();
	P.display();
	Pn.display();
	Pmin.display();
	Pmax.display();
	*/


#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 7, 1);
	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
}

void EndoPF::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	Pmin = cas.getPmin();
	Pmax = cas.getPmax();


	MatrixCPU Lb(cas.getLb());

	b = cas.getb();
	Cp1 = b;
	int indice = 0;

	for (int idAgent = 0; idAgent < _nAgent; idAgent++) {
		int Nvoisinmax = (int) nVoisin.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			matLb.set(indice, 0, Lb.get(idAgent, 0));
			indice = indice + 1;
		}
	}

	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);
	Cp1.multiplyT(&nVoisin);
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION


}

void EndoPF::init(const Simparam& sim, const StudyCase& cas)
{
	DELETEB(pf);
	// intitilisation des matrixs et variables 
	
	//std::cout << "init " << std::endl;
	if (!cas.isAC()) {
		throw std::invalid_argument("Wrong studyCase must be AC");
	}
	isAC = true;
	initSize(cas);
	isRadial = cas.isRadial();
	//_nLine = cas.getNLine(true);
	_nVarPF = _nLine + 2 * _nBus;

	initSimParam(sim);

	initCaseParam(sim, cas);
	if (initWithMarketClear) {
		ADMMMarket market;
		Simparam res(sim);
		market.solve(&res, sim, cas);
		//res.display();
		LAMBDA = res.getLambda();
		trade = res.getTrade();
		Pn = res.getPn();

	}
	Pnpre = Pn;
	//Pnpre.display();
	Tmoy = Pnpre;
	
	//std::cout << "*******" << std::endl;
	//std::cout << "mise sous forme lineaire" << std::endl;
	
	initLinForm(cas);

	//std::cout << "donnees sur CPU pour le grid" << std::endl;
	delta1 = MatrixCPU(_nVarPF, 1, 0);
	delta2 = MatrixCPU(_nVarPF, 1, 0);
	Z1 = MatrixCPU(_nVarPF, 1, 0);
	Z2 = MatrixCPU(_nVarPF, 1, 0);
	Y = MatrixCPU(_nVarPF, 1, 0);
	Ypre = MatrixCPU(_nVarPF, 1, 0);
	dY = MatrixCPU(_nVarPF, 1, 0);

	Ylimit = MatrixCPU(_nVarPF, 1, 0); // angle, amplitude, flux
	YOffset = MatrixCPU(_nVarPF, 1, 0); // angle, amplitude, flux
	G = MatrixCPU(_nVarPF, _nAgent);
	SensiBis = MatrixCPU(_nVarPF, 1);
	

	MatrixCPU LimitsUb = cas.getUpperBound(); // angle, amplitude, flux
	MatrixCPU LimitsLb = cas.getLowerBound();
	
	
	for (int i = 0; i < _nVarPF; i++) {
		float ub = LimitsUb.get(i, 0);
		float lb = LimitsLb.get(i, 0);
		float mid = (ub + lb) / 2; // 0 pour angle, 0, pour line, V0 pour tension 
		float lim = ub - mid; // 
		Ylimit.set(i, 0, lim);
		YOffset.set(i, 0, mid);
	}
	if (isRadial) {
		pf = new CPUPFdistPQ;
	}
	else {
		pf = new CPUPF;
	}
	//Ylimit.display();
	//YOffset.display();
	
	pf->init(cas, &Pnpre);


	//std::cout << "autres donn�e sur CPU" << std::endl;
	
	dP = MatrixCPU(_nAgent, 1, 0);

	initP2PMarket();
	
	Cp2 = MatrixCPU(_nAgent, 1, 0);
	Cp1 = b;
	Cp1.multiplyT(&nVoisin);
	
	updateGlobalProb();
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	//std::cout << "fin init " << std::endl;

}

void EndoPF::updateGlobalProb() {
	

#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	//std::cout << "*************" << std::endl;
	Pn.swap(&Pnpre);
	updatePn();
	//std::cout << " Tmoy " << std::endl;
	//Tmoy.display();
	
	/*std::cout << " Pn " << std::endl;
	Pn.display();
	tradeLin.display();
	*/
	pf->updatePQ(&Pn);

#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	timePerBlock.increment(0, 3, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	// FB 3
	pf->solve();
	tempL1 = pf->getY();
	Ypre.swap(&Y);
	Y.subtract(&tempL1, &YOffset);
	
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 4, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	//std::cout << "Y " << std::endl;
	//Y.display();

	float Ploss = - pf->getPloss();
	float Qloss = - pf->getQloss();	
	//std::cout << " Ploss " << Ploss << " Qloss " << Qloss << std::endl;

	Pmin.set(0, 0, Ploss / nVoisin.get(0, 0));
	Pmax.set(0, 0, Ploss / nVoisin.get(0, 0));
	Pmin.set(_nAgentTrue, 0, Qloss / nVoisin.get(_nAgentTrue, 0));
	Pmax.set(_nAgentTrue, 0, Qloss / nVoisin.get(_nAgentTrue, 0));

	for (int j = 0; j < _nAgentTrue - 1; j++) {
		if (Qloss > 0) {
			matUb.set(_nTradeP + j, 0, Qloss);
		}
		else {
			matLb.set(_nTradeP + j, 0, Qloss);
		}
	}
	
	 

	updateZ();

	updateDelta();
	updateSensi();
	/*
	std::cout << " Z " << std::endl;
	Z1.display();
	Z2.display();
	std::cout << " Delta " << std::endl;
	delta1.display();
	delta2.display();
	
	std::cout << " Sensi " << std::endl;
	SensiBis.display();
	G.display();
	std::cout << " Cp " << std::endl;
	Cp2.display();
	Cp.display();

	*/
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

void EndoPF::updateLocalProb() {
	// FB 1a
	updateBt2();
	updateTl();

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
	updateBp1();
	updateP();
	updateMU();
	
}


void EndoPF::updateBt1()
{
	// subtractTrans
	for (int t = 0; t < _nTrade; t++) {
		int k = (int) CoresLinTrans.get(t,0);
		Bt1.set(t, 0, tradeLin.get(t, 0) - tradeLin.get(k, 0));
	}
	Bt1.multiply(0.5f*_rhog); 
	Bt1.subtract(&LAMBDALin);
	Bt1.divide(_rhog);
}

void EndoPF::updateBt2()
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

void EndoPF::updateBp1()
{
	Bp1.add(&MU, &Tmoy);
}

void EndoPF::updateTl()
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

void EndoPF::updateSensi()
{
	SensiBis.set(&Y);
	SensiBis.multiply(2);
	SensiBis.add(&Z1);
	SensiBis.subtract(&Z2);
	SensiBis.add(&delta2);
	SensiBis.subtract(&delta1);

	dY.subtract(&Y, &Ypre);
	dP.subtract(&Pn, &Pnpre);
	
	//dY.display();
	//dP.display();

	for (int i = 0; i < _nVarPF; i++) {
		for (int n = 1; n < _nAgent; n++) {
			if (abs(dP.get(n,0)) > 0.01) {
				G.set(i, n, dY.get(i, 0) / dP.get(n, 0));
			}
		}
	}
	//G.display();
}

void EndoPF::updateZ()
{
	Z1.add(&Ylimit, &delta1);
	Z1.subtract(&Y);
	Z1.projectPos();

	Z2.add(&Ylimit, &delta2);
	Z2.add(&Y);
	Z2.projectPos();
}

void EndoPF::updateDelta()
{

	delta1.add(&Ylimit, &delta1);
	delta1.subtract(&Y);
	delta1.projectNeg();

	delta2.add(&Ylimit, &delta2);
	delta2.add(&Y);
	delta2.projectNeg();
	/*delta1.add(&Ylimit);
	delta1.subtract(&Z1);
	delta1.subtract(&Y);

	delta2.add(&Ylimit);
	delta2.subtract(&Z2);
	delta2.add(&Y);*/
}


float EndoPF::updateResEndo(int iter)
{
	for (int t = 0; t < _nTrade; t++) {
		int k = CoresLinTrans.get(t, 0);
		tempNN.set(t, 0, tradeLin.get(t, 0) + tradeLin.get(k, 0));
	}
	float resR = tempNN.max2();
	float resS = Tlocal.max2(&tradeLin);

	// Residus reseau
	
	tempL1.subtract(&Ylimit, &Y);
	tempL1.projectNeg();

	float resXf = _ratioEps * tempL1.max2();
	resF.set(0, iter, resR);
	resF.set(1, iter, resS);
	resF.set(2, iter, resXf);
	return MYMAX(MYMAX(resXf, resS), resR);
}

void EndoPF::updateP()
{
	P.multiplyT(&Ap1, &Bp1);
	P.subtract(&Cp);
	
	P.divideT(&Ap12);
	P.project(&Pmin, &Pmax);
	/*std::cout << "********power**************" << std::endl;
	P.display();
	Pmin.display();
	Pmax.display();
	P.project(&Pmin, &Pmax);
	P.display();*/
}

void EndoPF::updateMU()
{
	MU.add(&Tmoy);
	MU.subtract(&P);
}



void EndoPF::updateCp2()
{
	for (int n = 0; n < _nAgent; n++) {
		float sum = 0;
		for (int i = 0; i < _nVarPF; i++) {
			sum += SensiBis.get(i, 0) * G.get(i, n);
		}
		Cp2.set(n, 0, sum * _rho1 * nVoisin.get(n, 0));
	}



}

void EndoPF::display() {

	std::cout << _name << std::endl;
}
