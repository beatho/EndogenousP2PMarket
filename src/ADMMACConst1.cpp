#include "../head/ADMMACConst1.h"



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
	int iterG = sim.getIterG();
	int iterL = sim.getIterL();
	int stepL = sim.getStepL();
	int stepG = sim.getStepG();
	
	float epsL = sim.getEpsL();
	float epsG = sim.getEpsG();
	

	float fc = 0;

	int iterLocal = 0;
	_resG = 2 * epsG;
	float resL = 2 * epsL;
	_iterGlobal = 0;
	while ((_iterGlobal < iterG) && (_resG>epsG)) {
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
				timePerBlock.increment(0, 4, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
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
		if (!(_iterGlobal % stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			_resG = updateResBis(&resF, (_iterGlobal / stepG), &tempNN);
			//P.display();
			//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, (iterGlobal - 1) / stepG) << " " << resF.get(1, (iterGlobal - 1) / stepG) << " " << resF.get(2, (iterGlobal - 1) / stepG) << std::endl;

#ifdef INSTRUMENTATION
			cudaDeviceSynchronize();
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION		
		}
		_iterGlobal++;
	}
	#ifdef INSTRUMENTATION
	occurencePerBlock.increment(0, 5, iterGlobal);
	occurencePerBlock.increment(0, 6, iterGlobal);
	occurencePerBlock.increment(0, 7, iterGlobal);
	occurencePerBlock.increment(0, 8, iterGlobal / stepG);

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	std::cout << _iterGlobal << " " << iterLocal << " " <<
	 resL << " " << resF.get(0, (_iterGlobal - 1) / stepG) << " " <<
	  resF.get(1, (_iterGlobal - 1) / stepG) << " " << resF.get(2, (_iterGlobal - 1) / stepG) << std::endl;


	//fc = calcFc(&a, &b, &tradeLin, &Pn, &Ct, &tempN1, &tempNN);
	/*std::cout << "----trade--upper--lower--result--" << std::endl;
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
	Kappa2.display();*/
	// FB 5
	
	setResult(result, cas.isAC());

#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 9, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 9, 1);
	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
	
	
}

void ADMMACConst1::init(const Simparam& sim, const StudyCase& cas)
{
	// intitilisation des matrixs et variables 
	//std::cout << "init " << std::endl;
	isAC = true;
	initSize(cas);
	initSimParam(sim);
	
	
	_nVoisinRef = (int) nVoisin.get(0, 0);
	_nConstraint = _nLine + 2 * _nBus;
	
	initCaseParam(sim, cas);
	//std::cout << "mise sous forme lineaire" << std::endl;
	initLinForm(cas);
	//std::cout << "donnees sur CPU pour le grid" << std::endl;
	
	Kappa1     = MatrixCPU(_nConstraint, 1, 0);
	Kappa2     = MatrixCPU(_nConstraint, 1, 0);
	Kappa1_pre = MatrixCPU(_nConstraint, 1, 0);
	Kappa2_pre = MatrixCPU(_nConstraint, 1, 0);
	Qpart      = MatrixCPU(_nConstraint, _nAgent, 0);
	Qtot       = MatrixCPU(_nConstraint, 1, 0);
	rhoMn      = MatrixCPU(nVoisin);
	rhoMn.multiply(_rho1);
	rhoMn2     = MatrixCPU(nVoisin);
	rhoMn2.multiplyT(&rhoMn);

	G2 = MatrixCPU(_nConstraint, _nAgent);

	UpperBound = cas.getUpperBound();
	LowerBound = cas.getLowerBound();
	
	DiffBound = UpperBound;
	DiffBound.add(&LowerBound);


	//std::cout << "autres donn�e sur CPU" << std::endl;

	tempC1 = MatrixCPU(_nConstraint, 1, 0);
	tempCN = MatrixCPU(_nConstraint, _nAgent);
	//MatrixCPU temp1N(1, _nAgent, 0, 1);


	PQ = MatrixCPU(_nAgent, 1);
	PQ.setBloc(0, _nAgentTrue, 0, 1, &Pn);
	PF.init(cas, &PQ);

	initP2PMarket();
	
	Ap3 = MatrixCPU(_nAgent, 1);
	Ap123 = MatrixCPU(_nAgent, 1, 0);


	Cp1 = b;
	Cp1.multiplyT(&nVoisin);

	Cp = MatrixCPU(_nAgent, 1, 0);
	Cp2 = MatrixCPU(_nAgent, 1, 0);


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
	

	updatePn();
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
		Pmax.set(_nAgentTrue, 0, -Qloss / (_nAgentTrue - 1));
		Pmin.set(0, 0, -Ploss / _nVoisinRef);
		Pmin.set(_nAgentTrue, 0, -Qloss / (_nAgentTrue - 1));

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

	timePerBlock.increment(0, 5, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

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
	tempCN.multiplyT(&G);

	
	//std::cout << "updateCp2" << std::endl;
	Cp2.sumT(&tempCN, 1);
	Cp2.multiplyT(&rhoMn);
	
	Cp.add(&Cp1, &Cp2);

	// updateAp3
	//std::cout << "updateAp3" << std::endl;
	G2.set(&G);
	G2.multiplyT(&G);

	Ap3.sumT(&G2, 1);
	Ap3.multiplyT(&rhoMn2);

	Ap123.add(&Ap12, &Ap3);

#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	
	//std::cout << "lambda" << std::endl;
	// FB 3c
	updateLambda();
	updateBt1();
	
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
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

	timePerBlock.increment(0, 1, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 2b
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
	
	// FB 2c
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

void ADMMACConst1::updateLambda()
{
	for (int t = 0; t < _nTrade; t++) {
		int k = (int) CoresLinTrans.get(t, 0);
		float lamb = 0.5f * _rhog * (tradeLin.get(t, 0) + tradeLin.get(k, 0));
		LAMBDALin.set(t, 0, LAMBDALin.get(t,0)+lamb);
	}
}

void ADMMACConst1::updateBt1()
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

void ADMMACConst1::updateBt2()
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
		int k = (int) CoresLinTrans.get(t, 0);
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
	return MYMAX(MYMAX(resXf, resS), resR);
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
		for (int n = _nAgent - 1; n >= 0; n--) {
			qt += G.get(l, n) * PQ.get(n, 0);
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
