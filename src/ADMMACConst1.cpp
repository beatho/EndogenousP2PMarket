#include "../head/ADMMACConst1.h"
#define MAX(X, Y) X * (X >= Y) + Y * (Y > X)


ADMMACConst1::ADMMACConst1() : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " ADMMACConst1 Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;

}


ADMMACConst1::ADMMACConst1(float rho) : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default ADMMACConst1 Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
}

ADMMACConst1::~ADMMACConst1()
{
}
void ADMMACConst1::setParam(float rho)
{
	_rho = rho;
}

void ADMMACConst1::setTau(float tau)
{
	if (tau < 1) {
		throw std::invalid_argument("tau must be greater than 1");
	}
	_tau = tau;
}



void ADMMACConst1::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
		timePerBlock.set(0, 0, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		occurencePerBlock.set(0, 0, 1);
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
	iterGlobal = 0;
	while ((iterGlobal < iterG) && (resG>epsG)) {
		resL = 2 * epsL;
		iterLocal = 0;
		while (iterLocal< iterL && resL>epsL) {
			updateLocalProb();
			// FB 3
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
#ifdef INSTRUMENTATION
			occurencePerBlock.increment(0, 1, iterLocal);
			occurencePerBlock.increment(0, 2, iterLocal);
			occurencePerBlock.increment(0, 3, iterLocal);
			occurencePerBlock.increment(0, 4, iterLocal / stepL);
#endif // INSTRUMENTATION
		
		/*Bt2.display();
		std::cout << "----trade--upper--lower--result--" << std::endl;
		matUb.display();
		matLb.display();
		Tlocal_pre.display();*/
		//std::cout << "----obj--upper--lower--result--" << std::endl;
		//b.display();
		//Pmax.display();
		//Pmin.display();
		//P.display();
		//Tmoy.display();
		
		/*std::cout << "----upper--lower--result------" << std::endl;
		UpperBound.display();
		LowerBound.display();
		Constraint->display();*/
		Tlocal_pre.swap(&Tlocal);
		tradeLin.swap(&Tlocal);
		
		updateGlobalProb();

		// FB 4
		if (!(iterGlobal % stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateResBis(&resF, (iterGlobal / stepG), &tempNN);
			//P.display();
			//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, (iterGlobal - 1) / stepG) << " " << resF.get(1, (iterGlobal - 1) / stepG) << " " << resF.get(2, (iterGlobal - 1) / stepG) << std::endl;

#ifdef INSTRUMENTATION
			cudaDeviceSynchronize();
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION		
		}
		iterGlobal++;
	}
	#ifdef INSTRUMENTATION
	occurencePerBlock.increment(0, 5, iterGlobal);
	occurencePerBlock.increment(0, 6, iterGlobal);
	occurencePerBlock.increment(0, 7, iterGlobal);
	occurencePerBlock.increment(0, 8, iterGlobal / stepG);

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, (iterGlobal - 1) / stepG) << " " << resF.get(1, (iterGlobal - 1) / stepG) << " " << resF.get(2, (iterGlobal - 1) / stepG) << std::endl;

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

	//fc = calcFc(&a, &b, &tradeLin, &Pn, &Ct, &tempN1, &tempNN);
	std::cout << "----trade--upper--lower--result--" << std::endl;
	matUb.display();
	matLb.display();
	Tlocal_pre.display();
	LAMBDALin.display();
	Bt1.display();
	std::cout << "---P--upper--lower--result--" << std::endl;
	
	Pmax.display();
	Pmin.display();
	P.display();

	std::cout << "--Const--upper--lower--result------" << std::endl;
	UpperBound.display();
	LowerBound.display();
	Constraint->display();

	std::cout << "kappa" << std::endl;
	Kappa1.display();
	Kappa2.display();
	// FB 5
	/*result->setResF(&resF);
	result->setLAMBDA(&LAMBDA);
	result->setTrade(&trade);
	result->setDelta1(&Kappa1);
	result->setDelta2(&Kappa2);
	result->setIter(iterGlobal);
	result->setMU(&MU);
	
	result->setPn(&Pn);*/
	
	result->setFc(fc);
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 9, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 9, 1);
	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
	tall = clock() - tall;
	result->setTime((float)tall / CLOCKS_PER_SEC);
	
}

void ADMMACConst1::updateP0(const StudyCase& cas)
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
	timePerBlock.increment(0, 10, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 10, 1);
#endif // INSTRUMENTATION


}

void ADMMACConst1::init(const Simparam& sim, const StudyCase& cas)
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
	_nAgent2 = 2 * _nAgent;

	_rhol = _rho; //*nAgent
	//std::cout << "rho " << _rho << std::endl;
	if (_rho == 0) {
		_rhol = _rhog;
	}

	nVoisin = cas.getNvoi();

	_nLine = cas.getNLine();
	//std::cout << _nLine << std::endl;
	_nBus = cas.getNBus();

	_nTrade = nVoisin.sum();
	_nVoisinRef = nVoisin.get(0, 0);
	_nTradeQ = _nAgent * (_nAgent - 1);
	_nConstraint = _nLine + 2 * _nBus;
	_nTradeP = _nTrade - _nTradeQ;

	_at1 = _rhog; // represente en fait 2*a
	_at2 = _rhol;

	resF = MatrixCPU(3, (iterG / stepG) + 1);
	resX = MatrixCPU(4, (iterG / stepG) + 1);

	MatrixCPU BETA(cas.getBeta());
	MatrixCPU Ub(cas.getUb());
	MatrixCPU Lb(cas.getLb());
	LAMBDA = sim.getLambda(); 
	trade = sim.getTrade();

	//std::cout << "mise sous forme linéaire" << std::endl;
	
	CoresMatLin = MatrixCPU(_nAgent, _nAgent, -1);
	CoresLinAgent = MatrixCPU(_nTrade, 1);
	CoresAgentLin = MatrixCPU(_nAgent2 + 1, 1);
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
	int indice2 = _nTradeP;
	for (int idAgent = 0; idAgent < _nAgent; idAgent++) {
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
		for (int i = 0; i < _nAgent; i++) {
			if (i != idAgent) {
				matLb.set(indice2, 0, Lb.get(idAgent + _nAgent, 0));
				matUb.set(indice2, 0, Ub.get(idAgent + _nAgent, 0));
				CoresLinAgent.set(indice2, 0, idAgent);
				CoresLinVoisin.set(indice2, 0, i);
				indice2++;
				// reste à 0
			}
		}
		CoresAgentLin.set(_nAgent+ idAgent + 1, 0, indice2);

		CoresAgentLin.set(idAgent + 1, 0, indice);
	}
	for (int lin = 0; lin < _nTradeP; lin++) {
		int i = CoresLinAgent.get(lin, 0);
		int j = CoresLinVoisin.get(lin, 0);
		int k = CoresMatLin.get(j, i);
		CoresLinTrans.set(lin, 0, k);
	}
	for (int lin = _nTradeP; lin < _nTrade; lin++) {
		int i = CoresLinAgent.get(lin, 0);
		int j = CoresLinVoisin.get(lin, 0);
		
		if(i > j) 
		{
			i = i - 1;
		}
		int k = _nTradeP + j*(_nAgent-1) + i;
		CoresLinTrans.set(lin, 0, k);
	}
	for (int n = 0; n < _nAgent; n++) {
		for (int i = 0; i < _nAgent; i++) {
			if (i != n) {
				matLb.set(indice, 0, Lb.get(n + _nAgent, 0));
				matUb.set(indice, 0, Ub.get(n + _nAgent, 0));
				// reste à 0
				indice = indice + 1;
			}
		}
	}
	
	/*CoresLinAgent.display();
	CoresLinVoisin.display();
	CoresLinTrans.display();*/
	//std::cout << "donnees sur CPU pour le grid" << std::endl;
	
	Kappa1 = MatrixCPU(_nConstraint, 1, 0);
	Kappa2 = MatrixCPU(_nConstraint, 1, 0);
	Kappa1_pre = MatrixCPU(_nConstraint, 1, 0);
	Kappa2_pre = MatrixCPU(_nConstraint, 1, 0);
	Qpart = MatrixCPU(_nConstraint, _nAgent2, 0);
	Qtot = MatrixCPU(_nConstraint, 1, 0);
	rhoMn = MatrixCPU(nVoisin);
	rhoMn.multiply(_rho1);
	rhoMn2 = MatrixCPU(nVoisin);
	rhoMn2.multiplyT(&rhoMn);

	G2 = MatrixCPU(_nConstraint, _nAgent2);

	UpperBound = cas.getUpperBound();
	LowerBound = cas.getLowerBound();
	
	DiffBound = UpperBound;
	DiffBound.add(&LowerBound);


	//std::cout << "autres donnée sur CPU" << std::endl;
	tempNN = MatrixCPU(_nTrade, 1, 0);
	tempN1 = MatrixCPU(_nAgent, 1, 0); // plutôt que de re-allouer de la mémoire à chaque utilisation
	tempL1 = MatrixCPU(_nLine, 1, 0);
	tempC1 = MatrixCPU(_nConstraint, 1, 0);
	tempCN = MatrixCPU(_nConstraint, _nAgent2);
	//MatrixCPU temp1N(1, _nAgent, 0, 1);

	Tlocal = MatrixCPU(_nTrade, 1, 0);
	P = MatrixCPU(_nAgent2, 1, 0); // moyenne des trades
	Pn = MatrixCPU(sim.getPn()); // somme des trades
	PQ = MatrixCPU(_nAgent2, 1);
	PQ.setBloc(0, _nAgent, 0, 1, &Pn);
	PF.init(cas, &PQ);

	a = MatrixCPU(cas.geta());
	b = MatrixCPU(cas.getb());
	Ap2 = a;
    Ap2.multiplyT(&nVoisin);
	Ap2.multiplyT(&nVoisin);

	Ap1 = nVoisin;
	Ap1.multiply(_rhol);

	Ap3 = MatrixCPU(_nAgent2, 1);

	Ap12 = MatrixCPU(_nAgent2, 1, 0);
    Ap12.add(&Ap1, &Ap2);
	Ap123 = MatrixCPU(_nAgent2, 1, 0);


	Cp1 = b;
	Cp1.multiplyT(&nVoisin);

	Bt1 = MatrixCPU(_nTrade, 1, 0);
	Cp = MatrixCPU(_nAgent2, 1, 0);
	Cp2 = MatrixCPU(_nAgent2, 1, 0);
	
	Bp1 = MatrixCPU(_nAgent2, 1, 0);

	Pmin = MatrixCPU(cas.getPmin());
	Pmax = MatrixCPU(cas.getPmax());

	MU = MatrixCPU(_nAgent2, 1);// facteur reduit i.e lambda_l/_rho
	MU.setBloc(0, _nAgent, 0, 1, &sim.getMU());
	
	Tmoy = MatrixCPU(_nAgent2, 1, 0);
	Tmoy.setBloc(0, _nAgent, 0, 1, &Pn);

	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);
	
	
	Tmoy.divideT(&nVoisin);
	


	updateGlobalProb();
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	std::cout << "fin init" << std::endl;
}

void ADMMACConst1::updateGlobalProb() {
	
	// FB 3a
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	

	updatePn(&PQ,&Tmoy,&nVoisin);
	//if (!(iterGlobal % 50)) {
		PF.updatePQ(&PQ);
		PF.solve();
		Ploss = PF.getPloss();
		Qloss = PF.getQloss();
		if (Ploss < 0) {
			std::cout << "WARNING : negative 'loss', WTF ???" << std::endl;
			Ploss = 0;
		}
		else {
			//std::cout << "Power on the line : " << Ploss << " " << Qloss << std::endl;

		}
		Pmax.set(0, 0, -Ploss / _nVoisinRef);
		Pmax.set(_nAgent, 0, -Qloss / (_nAgent - 1));
		Pmin.set(0, 0, -Ploss / _nVoisinRef);
		Pmin.set(_nAgent, 0, -Qloss / (_nAgent - 1));

		//std::cout << "update phi" << std::endl;
		PF.calcPhi();
		//std::cout << "update G" << std::endl;
		G = PF.calcG();
		//std::cout << "update Constraint" << std::endl;
		Constraint = PF.calcY();

		updateQ();
	//}
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	timePerBlock.increment(0, 5, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 3b
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	

	Kappa1_pre.set(&Kappa1);
	Kappa2_pre.set(&Kappa2);
	
	//std::cout << "update Kappa" << std::endl;
	Kappa1.projectNeg();
	Kappa1.add(&UpperBound);
	Kappa1.subtract(Constraint);

	Kappa2.projectNeg();
	Kappa2.subtract(&LowerBound);
	Kappa2.add(Constraint);
	
	//std::cout << "updateCp" << std::endl;
	//std::cout << "updateC1" << std::endl;
	tempC1.subtractAbs(&Kappa1, &Kappa2);
	tempC1.multiply(0.5);
	tempC1.add(Constraint);
	tempC1.subtract(&Qtot);
	tempC1.multiply(2);
	tempC1.subtract(&DiffBound);

	//std::cout << "updateCN" << std::endl;
	tempCN.set(&Qpart);
	tempCN.multiply(2);
	tempCN.addVector(&tempC1);
	tempCN.multiplyT(G);

	
	//std::cout << "updateCp2" << std::endl;
	Cp2.sumT(&tempCN, 1);
	Cp2.multiplyT(&rhoMn);
	
	Cp.add(&Cp1, &Cp2);

	// updateAp3
	//std::cout << "updateAp3" << std::endl;
	G2.set(G);
	G2.multiplyT(G);

	Ap3.sumT(&G2, 1);
	Ap3.multiplyT(&rhoMn2);

	Ap123.add(&Ap12, &Ap3);

#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 6, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	
	//std::cout << "lambda" << std::endl;
	// FB 3c
	updateLambda();
	updateBt1();
	
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 7, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
	
}

void ADMMACConst1::updateLocalProb() {
	// FB 2a
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	
	//std::cout << "update Bt2" << std::endl;
	updateBt2();
	//std::cout << "update Tl" << std::endl;
	updateTl();
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	timePerBlock.increment(0, 1, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 2b
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	
	
	for (int i = 0; i < _nAgent2; i++) {
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
	
	// FB 2c
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

void ADMMACConst1::updateLambda()
{
	for (int t = 0; t < _nTrade; t++) {
		int k = CoresLinTrans.get(t, 0);
		float lamb = 0.5 * _rhog * (tradeLin.get(t, 0) + tradeLin.get(k, 0));
		LAMBDALin.set(t, 0, LAMBDALin.get(t,0)+lamb);
	}
}

void ADMMACConst1::updateBt1()
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

void ADMMACConst1::updateBt2()
{
	

	for (int i = 0; i < _nAgent2; i++) {
		int nVoisinLocal = nVoisin.get(i,0);
		int beginLocal = CoresAgentLin.get(i,0);
		int endLocal = beginLocal + nVoisinLocal; 
		for (int j = beginLocal; j < endLocal; j++) {
			float m = Tlocal_pre.get(j,0) - Tmoy.get(i, 0) + P.get(i, 0) - MU.get(i, 0); 
			Bt2.set(j, 0, m);
		}
	}
	/*Tlocal_pre.display();
	Tmoy.display();
	P.display();
	MU.display();
	std::cout << "----" << std::endl;
	Bt2.display();*/
}

void ADMMACConst1::updateBp1()
{
	Bp1.add(&MU, &Tmoy);
}

void ADMMACConst1::updateTl()
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

float ADMMACConst1::calcRes()
{
	MatrixCPU temp(Tlocal);
	temp.subtract(&Tlocal_pre);

	MatrixCPU temp2(Tmoy);
	temp2.subtract(&P);
	float d1 = temp.max2();
	float d2 = temp2.max2();
	

	return d1 * (d1 > d2) + d2 * (d2 >= d1);
}

float ADMMACConst1::updateResBis(MatrixCPU* res, int iter, MatrixCPU* tempNN)
{
	for (int t = 0; t < _nTrade; t++) {
		int k = CoresLinTrans.get(t, 0);
		tempNN->set(t, 0, tradeLin.get(t, 0) + tradeLin.get(k, 0));
	}
	float resR = tempNN->max2();

	MatrixCPU temp2(Tlocal);
	
	temp2.subtract(&tradeLin);
	float resS = temp2.max2();

	// version de l'article
	MatrixCPU tempL(Kappa1);
	MatrixCPU tempL2(Kappa2);
	Kappa1_pre.projectNeg();
	Kappa2_pre.projectNeg();
	tempL.projectNeg();
	tempL2.projectNeg();
	tempL.subtract(&Kappa1_pre);
	tempL2.subtract(&Kappa2_pre);
	tempL.multiplyT(&tempL);
	tempL2.multiplyT(&tempL2);
	tempL.add(&tempL2);

	float resXf = _ratioEps * sqrt(tempL.max2());
	res->set(0, iter, resR);
	res->set(1, iter, resS);
	res->set(2, iter, resXf);
	return MAX(MAX(resXf, resS), resR);
}

void ADMMACConst1::updateP()
{
	P.multiplyT(&Ap1, &Bp1);
	P.subtract(&Cp);
	
	P.divideT(&Ap123);
	P.project(&Pmin, &Pmax);
}

void ADMMACConst1::updateMU()
{
	MU.add(&Tmoy);
	MU.subtract(&P);
}

void ADMMACConst1::updateQ()
{
	for (int l = 0; l < _nLine; l++) {
		float qt = 0;
		for (int n = _nAgent2 - 1; n >= 0; n--) {
			qt += G->get(l, n) * PQ.get(n, 0);
			if (n > 0) {
				Qpart.set(l, n - 1, qt);
			}
		}
		Qtot.set(l, 0, qt);
	}
}

void ADMMACConst1::display() {

	std::cout << _name << std::endl;
}
