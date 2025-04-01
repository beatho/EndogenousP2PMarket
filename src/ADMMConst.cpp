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
	timePerBlock = MatrixCPU(1, 8, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 8, 0); //nb de fois utilis� pendant la simu

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
	int iterG = sim.getIterG();
	int iterL = sim.getIterL();
	int stepL = sim.getStepL();
	int stepG = sim.getStepG();
	
	float epsL = sim.getEpsL();
	float epsG = sim.getEpsG();
	

	float fc = 0;

	int iterLocal = 0;
	float resG = 2 * epsG;
	float resL = 2 * epsL;
	int iterGlobal = 0;
	while ((iterGlobal < iterG) && (resG>epsG) || (iterGlobal <= stepG)) {
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
		timePerBlock.increment(0, 1, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION			
		//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;
		Tlocal_pre.swap(&Tlocal);
		tradeLin.swap(&Tlocal);
		
		updateGlobalProb();

		// FB 4
		if (!(iterGlobal % stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateResEndo(iterGlobal / stepG);
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

void ADMMConst::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	if (cas.isAC()) {
		
		MatrixCPU bT = cas.getb();
		MatrixCPU PminT = cas.getPmin();
		MatrixCPU PmaxT = cas.getPmax();
	

		for (int n = 0; n < _nAgent; n++) {
						b.set(n, 0, bT.get(n, 0));
			Pmin.set(n, 0, PminT.get(n, 0));
			Pmax.set(n, 0, PmaxT.get(n, 0));
		}

	}
	else {
		b = cas.getb();
		Pmin = cas.getPmin();
		Pmax = cas.getPmax();
	}
	


	MatrixCPU Lb(cas.getLb());

	
	Cp1 = b;
	int indice = 0;

	for (int idAgent = 0; idAgent < _nAgent; idAgent++) {
		int Nvoisinmax = nVoisin.get(idAgent, 0);
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

void ADMMConst::init(const Simparam& sim, const StudyCase& cas)
{
	// intitilisation des matrixs et variables 
	
	//std::cout << "init " << std::endl;
	_rhog = sim.getRho();
	_rho1 = sim.getRho1();
	const int iterG = sim.getIterG();
	const int stepG = sim.getStepG();
	float epsG = sim.getEpsG();
	float epsGC = sim.getEpsGC();
	_ratioEps = epsG / epsGC;
	_nAgent = sim.getNAgent();
	

	_rhol = _rho; //*nAgent
	//std::cout << "rho " << _rho << std::endl;
	if (_rho == 0) {
		_rhol = _rhog;
	}
	if (cas.isAC()) {
		MatrixCPU nVoisinT = cas.getNvoi();
		nVoisin = MatrixCPU(_nAgent, 1);
		for (int n = 0; n < _nAgent; n++) {
			nVoisin.set(n, 0, nVoisinT.get(n, 0));
		}
	}
	else {
		nVoisin = cas.getNvoi();
	}
	

	_nLine = cas.getNLine();
	//std::cout << "_nLine " << _nLine << std::endl;
	_nBus = cas.getNBus();

	_nTrade = nVoisin.sum();
	_at1 = _rhog; // represente en fait 2*a
	_at2 = _rhol;

	resF = MatrixCPU(3, (iterG / stepG) + 1);
	resX = MatrixCPU(4, (iterG / stepG) + 1);

	MatrixCPU BETA(cas.getBeta());
	
	MatrixCPU Ub(cas.getUb());
	MatrixCPU Lb(cas.getLb());
	LAMBDA = sim.getLambda();
	trade = sim.getTrade();

	//std::cout << "mise sous forme lin�aire" << std::endl;
	
	CoresMatLin = MatrixCPU(_nAgent, _nAgent, -1);
	CoresLinAgent = MatrixCPU(_nTrade, 1);
	CoresAgentLin = MatrixCPU(_nAgent + 1, 1);
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

	for (int idAgent = 0; idAgent < _nAgent; idAgent++) {
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
	for (int lin = 0; lin < _nTrade; lin++) {
		int i = CoresLinAgent.get(lin, 0);
		int j = CoresLinVoisin.get(lin, 0);
		int k = CoresMatLin.get(j, i);
		CoresLinTrans.set(lin, 0, k);
	}


	// transfert des mises lineaire
	

	//std::cout << "donnees sur CPU pour le grid" << std::endl;
	Kappa1 = MatrixCPU(_nLine, 1, 0);
	Kappa2 = MatrixCPU(_nLine, 1, 0);
	Kappa1_pre = MatrixCPU(_nLine, 1, 0);
	Kappa2_pre = MatrixCPU(_nLine, 1, 0);
	Qpart = MatrixCPU(_nLine, _nAgent, 0);
	Qtot = MatrixCPU(_nLine, 1, 0);
	alpha = MatrixCPU(_nLine, _nAgent, 0);

	G = MatrixCPU(cas.getPowerSensi());
	
	

	lLimit = MatrixCPU(cas.getLineLimit());

	GTrans = MatrixCPU(_nAgent, _nLine);
	GTrans.setTrans(&G);

	G2 = GTrans;
	G2.multiplyT(&GTrans);


	//std::cout << "autres donn�e sur CPU" << std::endl;
	tempNN = MatrixCPU(_nTrade, 1, 0);
	tempN1 = MatrixCPU(_nAgent, 1, 0); // plut�t que de re-allouer de la m�moire � chaque utilisation
	tempL1 = MatrixCPU(_nLine, 1, 0);
	//MatrixCPU temp1N(1, _nAgent, 0, 1);

	Tlocal = MatrixCPU(_nTrade, 1, 0);
	P = MatrixCPU(_nAgent, 1, 0); // moyenne des trades
	Pn = MatrixCPU(_nAgent, 1, 0); // somme des trades

	// si cas AC, a, b ,Nvoisin, Pmin,Pma n'ont pas la bonne taille !!!
	if (cas.isAC()) {
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

		for (int n = 0; n < _nAgent; n++) {
			a.set(n, 0, aT.get(n, 0));
			b.set(n, 0, bT.get(n, 0));
			Pmin.set(n, 0, PminT.get(n, 0));
			Pmax.set(n, 0, PmaxT.get(n, 0));
			MU.set(n, 0, MUT.get(n, 0));
			Tmoy.set(n, 0, TmoyT.get(n, 0));
		}

	}
	else {
		a = MatrixCPU(cas.geta());
		b = MatrixCPU(cas.getb());
		
		
		
		Pmin = MatrixCPU(cas.getPmin());
		Pmax = MatrixCPU(cas.getPmax());
		MU = MatrixCPU(sim.getMU()); // facteur reduit i.e lambda_l/_rho
		Tmoy = MatrixCPU(sim.getPn());
	}
	Ap1 = nVoisin;
	Ap2 = a;
	Cp1 = b;
	
	Ap12 = MatrixCPU(_nAgent, 1, 0);

	Bt1 = MatrixCPU(_nTrade, 1, 0);
	Cp = MatrixCPU(_nAgent, 1, 0);
	Cp2 = MatrixCPU(_nAgent, 1, 0);
	Bp1 = MatrixCPU(_nAgent, 1, 0);

	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);
	Ap1.multiply(_rhol);
	Cp1.multiplyT(&nVoisin);
	Tmoy.divideT(&nVoisin);
	tempN1.sum(&G2);
	
	tempN1.multiply(2 * _rho1);
	Ap2.add(&tempN1);
	Ap2.multiplyT(&nVoisin);
	Ap2.multiplyT(&nVoisin);
	Ap12.add(&Ap1, &Ap2);


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

void ADMMConst::display() {

	std::cout << _name << std::endl;
}
