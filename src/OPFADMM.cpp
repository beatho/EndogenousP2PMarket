#include "../head/OPFADMM.h"
 


OPFADMM::OPFADMM() : MethodOPF()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " OPFADMM Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	timePerBlock = MatrixCPU(1, 12, 0); // Fb0, Fb11abcd, FB12, Fb2, Fb3, Fb4, Fb5,FB6, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 12, 0); //nb de fois utilis� pendant la simu
}


OPFADMM::OPFADMM(float rho) : MethodOPF()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default OPFADMM Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
	timePerBlock = MatrixCPU(1, 12, 0); // Fb0, Fb11, FB12, Fb2, Fb3, Fb4, Fb5, FB6, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 12, 0); //nb de fois utilis� pendant la simu
}

OPFADMM::~OPFADMM()
{
	DELETEA(tempM1);
	DELETEA(tempM);

	DELETEA(X);
	DELETEA(Ypre);
	DELETEA(Y);
	DELETEA(YTrans);
	DELETEA(Mu);

	DELETEA(Hinv);
	DELETEA(A);
	DELETEA(Q);

	DELETEA(Childs);
	
}
void OPFADMM::setParam(float rho)
{
	_rho = rho;
}

bool OPFADMM::chekcase()
{
	if (_nBus != (_nLine + 1)) {
		std::cout << "wrong number of line " << _nLine << "against " << _nBus << std::endl;
		return false;
	}
	
	for (int i = 0; i < _nLine; i++) {
		if (CoresLineBus.get(i, 1) != (i + 1)) {
			std::cout << "wrong numerotation of line " << CoresLineBus.get(i, 1) << "against " << (i + 1) << std::endl;
			return false;
		}
		if (CoresLineBus.get(i, 0) > CoresLineBus.get(i, 1)) {
			std::cout << "wrong numeoration of bus " << CoresLineBus.get(i, 0) << "against " << CoresLineBus.get(i, 1) << std::endl;
			return false;
		}
	}
	if (ZsRe.getNLin() == 0  || ZsIm.getNLin() == 0) {
		std::cout << "matrice non defined, ZsRe, Zs Im, Yd" << std::endl;
		ZsRe.display();
		ZsIm.display();
		return false;
	}

	return true;
}

void OPFADMM::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
		timePerBlock.increment(0, 0, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION
	}
	
	_iterG = sim.getIterG();
	int iterL = sim.getIterL();
	_stepG = sim.getStepG();
	int stepL = sim.getStepL();
	
	float epsG = sim.getEpsG();
	float epsL = sim.getEpsL();
	
	
	float fc = 0;
	float resG = 2 * epsG;
	float resL = 2 * epsL;
	_iterGlobal = 0;
	
	/*CoresSoloBusAgent.display();

	Apt1.display();
	Apt2.display();
	Bpt2.display();
	Cost1.display();
	Cost2.display();
	Pmin.display();
	Pmax.display();

	PnPre.display();
	PnMoy.display();
	PnTilde.display();
	MuL.display();
	_nAgentByBus.display();
	std::cout << _rhol << " " << epsL << " " << iterL << " " << _nAgent << " " << _nBus;
	std::cout << "------" << std::endl;*/

	while ((_iterGlobal < _iterG) && (resG>epsG)) {
		/*std::cout << "--------" << std::endl;
		Pn.display();
		for (int i = 0; i < 2; i++) {
			std::cout << " X " << i << std::endl;
			X[i].display();
			//std::cout << " Y " << i << std::endl;
			//Y[i].display();
			//std::cout << " Q " << i << std::endl;
			//Q[i].display();
			//std::cout << " Mu " << i << std::endl;
			//Mu[i].display();
		}
		Chat.display();
		Bpt2.display();
		*/
		_iterLocal = 0;
		resL = 2 * epsL;
		while (_iterLocal< iterL && resL>epsL) {
			updateLocalProb();
			// FB 3
			if (!(_iterLocal % stepL)) {
#ifdef INSTRUMENTATION
				t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
				resL = calcRes();
#ifdef INSTRUMENTATION
				t2 = std::chrono::high_resolution_clock::now();
				timePerBlock.increment(0, 4, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION			
			}
			_iterLocal++;
		}
		/*std::cout << " Pn " << std::endl;
		Pn.display();
		std::cout << " PnMoy " << std::endl;
		PnMoy.display();
		std::cout << " PnTilde " << std::endl;
		PnTilde.display();
		std::cout << " MuL " << std::endl;
		MuL.display();*/
		/*if (_iterLocal == iterL) {
			std::cout << _iterGlobal << " " << _iterLocal << " " << resL << " " << resG << std::endl;
		}*/
#ifdef INSTRUMENTATION
		occurencePerBlock.increment(0, 1, _iterLocal);
		occurencePerBlock.increment(0, 2, _iterLocal);
		occurencePerBlock.increment(0, 3, _iterLocal);
		occurencePerBlock.increment(0, 4, _iterLocal / stepL);
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

		updateX();

#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 5, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		//std::cout << "---------" << std::endl;
		CommunicationX();
		
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		updateGlobalProb();
		updateMu();
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		updateChat();
		
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		// FB 4
		if (!(_iterGlobal % _stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateRes(_iterGlobal / _stepG);
			//std::cout << _iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, _iterGlobal / _stepG) << " " << resF.get(1, _iterGlobal / _stepG) << std::endl;
#ifdef INSTRUMENTATION
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 9, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		}
		//std::cout << iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;

		_iterGlobal++;
	}
	//std::cout << iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, (iterGlobal - 1) / stepG) << " " << resF.get(1, (iterGlobal - 1) / stepG) << " " << resF.get(2, (iterGlobal - 1) / stepG) << std::endl;


#ifdef INSTRUMENTATION	
	occurencePerBlock.increment(0, 5, _iterGlobal);
	occurencePerBlock.increment(0, 6, _iterGlobal);
	occurencePerBlock.increment(0, 7, _iterGlobal);
	occurencePerBlock.increment(0, 8, _iterGlobal);
	occurencePerBlock.increment(0, 9, _iterGlobal / _stepG);

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	
	
	fc = calcFc(&Cost1, &Cost2, &Pn, &tempN2);
	// FB 5
	
	Pn.set(0, 0, getPLoss());
	Pn.set(_nAgent, 0, getQLoss());
	result->setResF(&resF);
	

	/*std::cout << " Pn " << std::endl;
	Pn.display();
	std::cout << " PnTilde " << std::endl;
	PnTilde.display();*/
	
	MatrixCPU Pb(getPb());
	MatrixCPU Phi(getPhi());
	//Phi.display();
	MatrixCPU E(getE());
	
	result->setE(&E);
	result->setPhi(&Phi);
	result->setPb(&Pb);


	result->setIter(_iterGlobal);
	result->setPn(&Pn);
	result->setFc(fc);

#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 10, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 10, 1);

	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
	tall = clock() - tall;
	timeOPF = tall;

	result->setTime((float)tall / CLOCKS_PER_SEC);
	
}


void OPFADMM::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	

	Pmin = cas.getPmin();
	Pmax = cas.getPmax();
	Cost2 = cas.getb();

	// pour essayer que cela marche
	Pn.add(&Pmin, &Pmax);
	Pn.divide(2);
	PnTilde.set(0.0);

	for (int n = 1; n < _nAgent; n++) {
		int b = _CoresBusAgent.get(n, 0);
		PnTilde.increment(b, 0, Pn.get(n, 0));
		PnTilde.increment(b + _nBus, 0, Pn.get(n + _nAgent, 0));
	}

	for (int b = 0; b < _nBus; b++) {
		int nb = _nAgentByBus.get(b, 0);
		if (nb > 0) {
			PnTilde.set(b, 0, PnTilde.get(b, 0) / nb);
			PnTilde.set(b + _nBus, 0, PnTilde.get(b + _nBus, 0) / nb);
		}
		X[b].set(4, 0, nb * PnTilde.get(b, 0)); // pi = sum(pn)
		X[b].set(5, 0, nb * PnTilde.get(b + _nBus, 0)); // qi = sum(pn)
	}
	PnMoy.set(&PnTilde);
	DFSP(0); // Pi
	X[0].set(0, 0, 0);
	DFSQ(0); // Qi
	X[0].set(1, 0, 0);
	for (int i = 0; i < _nBus; i++) {
		//X[i].set(2, 0, 1);
		float Si = X[i].get(0, 0) * X[i].get(0, 0) + X[i].get(1, 0) * X[i].get(1, 0);
		X[i].set(3, 0, Si / X[i].get(2, 0)); // li = Si^2/vi
	}
	for (int i = 0; i < _nBus; i++) {
		int m = nChild.get(i, 0);
		for (int j = 0; j < m; j++) {// (Pci, Qci, lci) for all child Ci
			int c = Childs[i].get(j, 0);
			X[i].set(7 + 3 * j, 0, X[c].get(0, 0));
			X[i].set(8 + 3 * j, 0, X[c].get(1, 0));
			X[i].set(9 + 3 * j, 0, X[c].get(3, 0));
		}
		Y[i].set(&X[i]);
		//Mu[i].set(0.0);
	}/**/
	updateChat();
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 11, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 11, 1);
#endif // INSTRUMENTATION

}

void OPFADMM::init(const Simparam& sim, const StudyCase& cas)
{
	// intitilisation des matrixs et variables 
	//cas.display();
	//sim.display(1);
	clock_t t = clock();
	
	_rho = sim.getRho();
	
	if (_rhol == 0) {
		_rhol = _rho;
	}
	
	
	const int iterG = sim.getIterG();
	const int stepG = sim.getStepG();
	
	_nAgent = cas.getNagent();
	
	_nBus = cas.getNBus();
	_nLine = cas.getNLine(true); // ne doit pas etre reduit ici !!!

	//std::cout << _nAgent << " " << _nBus << " " << _nLine << std::endl;
	_nAgentByBus = cas.getNagentByBus();
	nChild = MatrixCPU(_nBus, 1);
	CoresLineBus = cas.getCoresLineBus(true);
	_CoresBusAgent = cas.getCoresBusAgentLin(); // Cores[n] = b
	Ancestor = MatrixCPU(_nBus, 1, 0); // A_i = bus ant�c�dent de i
	PosChild = MatrixCPU(_nBus, 1, 0); // indice du bus i dans Child[Ai]
	Ancestor.set(0, 0, -1); // the slack bus has no ancestor
	ZsRe = cas.getZsRe();
	ZsIm = cas.getZsImag();
	ZsNorm = MatrixCPU(_nLine, 1);
	
	if (!chekcase()) {
		throw std::invalid_argument("not a radial case");
	}
	
	for (int lold = 0; lold < _nLine; lold++) {
		int l = lold + 1;
		int busTo = l ;
		if(busTo != CoresLineBus.get(lold, 1)){
			throw std::invalid_argument(" Problem on the link between buses, line l must go to bus l+1");
		}
		int busFrom = CoresLineBus.get(lold, 0);
		Ancestor.set(busTo, 0, busFrom);
		nChild.increment(busFrom, 0, 1);
		ZsNorm.set(lold, 0, ZsRe.get(lold, 0) * ZsRe.get(lold, 0) + ZsIm.get(lold, 0) * ZsIm.get(lold, 0));
	}
	//ZsNorm.display();
	//Ancestor.display();
	//nChild.display();

	_rhoInv = 1 / _rho;
	resF = MatrixCPU(3, (iterG / stepG) + 1);

	
	MatrixCPU lowerBound(cas.getLowerBound()); //voltage angle, voltage, line...
	MatrixCPU upperBound(cas.getUpperBound()); //voltage angle, voltage, line...
	
	//lowerBound.display();
	//upperBound.display();

	//std::cout << " local resolution " << std::endl;
	// local resolution
	tempN2 = MatrixCPU(2 * _nAgent, 1);
	tempB2 = MatrixCPU(2 * _nBus, 1);
	CoresSoloBusAgent = MatrixCPU(_nBus, 1, -1);
	Pn = sim.getPn();
	//Pn.display();
	Pmin = cas.getPmin();
	Pmax = cas.getPmax();
	//Pmin.display();
	//Pmax.display();

	PnTmin = MatrixCPU(2 * _nBus, 1);
	PnTmax = MatrixCPU(2 * _nBus, 1);

	Cost1 = MatrixCPU(cas.geta());
	Cost2 = MatrixCPU(cas.getb());

	//Cost1.display();
	//Cost2.display();

	PnMoy = MatrixCPU(2 * _nBus, 1);
	PnPre = sim.getPn();
	MuL = MatrixCPU(2 * _nBus, 1);
	PnTilde = MatrixCPU(2 * _nBus, 1);
	Ap12 = Cost1;

	Ap12.add(_rhol);
	
	Bp1 = MatrixCPU(2 * _nAgent, 1);
	Bpt1 = MatrixCPU(2 * _nBus, 1);
	Bpt2 = MatrixCPU(2 * _nBus, 1);
	Apt1 = MatrixCPU(2 * _nBus, 1);
	Apt2 = MatrixCPU(2 * _nBus, 1);
	Apt12 = MatrixCPU(2 * _nBus, 1);
	
	if (Pn.max2() == 0) {
		Pn.add(&Pmin, &Pmax);;
		Pn.divide(2);
	}
	
	_nAgentByBus.increment(0, 0, -1); // don't want the grid agent in this resolution
	//_nAgentByBus.display();
	for (int n = 1; n < _nAgent; n++) {
		int b = _CoresBusAgent.get(n, 0);
		PnTilde.increment(b, 0, Pn.get(n, 0));
		PnTilde.increment(b + _nBus, 0, Pn.get(n + _nAgent, 0));
		
		PnTmax.increment(b, 0, Pmax.get(n, 0));
		PnTmax.increment(b + _nBus, 0, Pmax.get(n + _nAgent, 0));
		PnTmin.increment(b, 0, Pmin.get(n, 0));
		PnTmin.increment(b + _nBus, 0, Pmin.get(n + _nAgent, 0));

		if (_nAgentByBus.get(b, 0) == 1) {
			CoresSoloBusAgent.set(b, 0, n);
		}
	}
	
	for (int b = 0; b < _nBus; b++) {
		int nb = _nAgentByBus.get(b, 0);
		if (nb > 0) {
			Apt1.set(b, 0, nb * _rhol);
			Apt1.set(b + _nBus, 0, nb * _rhol);

			Apt2.set(b, 0, nb * nb * _rho);
			Apt2.set(b + _nBus, 0, nb * nb * _rho);

			PnTilde.set(b, 0, PnTilde.get(b, 0) / nb);
			PnTilde.set(b + _nBus, 0, PnTilde.get(b + _nBus, 0) / nb);

			PnTmax.set(b, 0, PnTmax.get(b, 0) / nb);
			PnTmax.set(b + _nBus, 0, PnTmax.get(b + _nBus, 0) / nb);

			PnTmin.set(b, 0, PnTmin.get(b, 0) / nb);
			PnTmin.set(b + _nBus, 0, PnTmin.get(b + _nBus, 0) / nb);
		}
		
	}
	PnMoy.set(&PnTilde);
	
	Apt12.add(&Apt1, &Apt2);
	
	/*std::cout << "Apt1 " << std::endl;
	Apt1.display();
	std::cout << "Apt2 " << std::endl;
	Apt2.display();
	std::cout << " Pmin " << std::endl;
	Pmin.display();
	std::cout << " Pmax " << std::endl;
	Pmax.display();

	std::cout << " Pn  limits " << std::endl;
	PnTmin.display();
	PnTmax.display();*/



	//std::cout << " creation " << std::endl;
	X = new MatrixCPU[_nBus];
	Ypre = new MatrixCPU[_nBus];
	Y = new MatrixCPU[_nBus];
	YTrans = new MatrixCPU[_nBus];
	Mu = new MatrixCPU[_nBus];
	
	tempN1 = MatrixCPU(_nAgent, 1);
	tempNN = MatrixCPU(_nAgent, _nAgent);
	//tempM1 = new MatrixCPU[_nAgent];
	tempM = new MatrixCPU[_nBus];
	

	Hinv = new MatrixCPU[_nBus];
	A = new MatrixCPU[_nBus];
	Q = new MatrixCPU[_nBus];
	
	Childs = new MatrixCPU[_nBus];

	Chat = MatrixCPU(_nBus, 4);
	VoltageLimit = MatrixCPU(_nBus, 2); // min, max
	VoltageLimitReal = MatrixCPU(_nBus, 2); // min, max
	sizeOPFADMM = MatrixCPU(_nBus, 1);
	
	
	int indice = 0;
	MatrixCPU nChildTemp(_nBus, 1, 0);
	for (int i = 0; i < _nBus; i++) {
		VoltageLimitReal.set(i, 0, lowerBound.get(_nBus + i, 0));
		VoltageLimitReal.set(i, 1, upperBound.get(_nBus + i, 0));
		VoltageLimit.set(i, 0, lowerBound.get(_nBus + i, 0) * lowerBound.get(_nBus + i, 0) * sqrt((nChild.get(i, 0) + 1) / 2));
		VoltageLimit.set(i, 1, upperBound.get(_nBus + i, 0) * upperBound.get(_nBus + i, 0) * sqrt((nChild.get(i, 0) + 1) / 2));
		
		Childs[i] = MatrixCPU(nChild.get(i, 0), 1);
		
		int sizeOPF =  3 * nChild.get(i, 0) + 7;
		int sizeOPF2 = 1 * nChild.get(i, 0) + 9;
		sizeOPFADMM.set(i, 0, sizeOPF);

		X[i] = MatrixCPU(sizeOPF, 1); // {Pi, Qi, vi, li, pi, qi, vAi, (Pci, Qci, lci) for all child Ci}
		Ypre[i] = MatrixCPU(sizeOPF, 1);// Y[i][j] not� dans l'article Yji est ce que i connait sur j
		Y[i] = MatrixCPU(sizeOPF, 1); //Y[i] = {Pi, Qi, vi, li, pi, qi, vAi, (Pci, Qci, lci) for all child Ci}
		YTrans[i] = MatrixCPU(sizeOPF2, 1);   // (Pi,Qi,vi,li,pi,qi,pai,qai,lai,(vai) for all child Ci) !!!
		Mu[i] = MatrixCPU(sizeOPF, 1);
		A[i] = MatrixCPU(2 + 1*(i>0), sizeOPF);
		Hinv[i] = MatrixCPU(sizeOPF, sizeOPF);
		Q[i] = MatrixCPU(sizeOPF, 1, 0);
		tempM[i] = MatrixCPU(sizeOPF, 1);
	}
	for (int i = 0; i < _nBus; i++) {
		if (i > 0) {
			int Ai = Ancestor.get(i, 0);
			Childs[Ai].set(nChildTemp.get(Ai, 0), 0, i);
			PosChild.set(i, 0, nChildTemp.get(Ai, 0));
			nChildTemp.increment(Ai, 0, 1);
		}
	}

	

	//std::cout << " init valeur " << std::endl;
	for (int i = 0; i < _nBus; i++) {
		// vi = 1; vai = 1
		X[i].set(2, 0, 1); 
		X[i].set(6, 0, 1);

		// pi & qi
		int Nb = _nAgentByBus.get(i, 0);
		X[i].set(4, 0, Nb * PnTilde.get(i,0)); // pi = sum(pn)
		X[i].set(5, 0, Nb * PnTilde.get(i + _nBus, 0)); // qi = sum(pn)
	}
	DFSP(0); // Pi
	X[0].set(0, 0, 0);
	DFSQ(0); // Qi
	X[0].set(1, 0, 0); 
	
	for (int i = 0; i < _nBus; i++) {
		float Si = X[i].get(0, 0) * X[i].get(0, 0) + X[i].get(1, 0) * X[i].get(1, 0);
		X[i].set(3, 0, Si / X[i].get(2, 0)); // li = Si^2/vi
	}
	for (int i = 0; i < _nBus; i++) {
		int m = nChild.get(i, 0);
		for (int j = 0; j < m; j++) {
			// (Pci, Qci, lci) for all child Ci
			
			int c = Childs[i].get(j, 0);
			X[i].set(7 + 3 * j, 0, X[c].get(0, 0));
			X[i].set(8 + 3 * j, 0, X[c].get(1, 0));
			X[i].set(9 + 3 * j, 0, X[c].get(3, 0));
		}
		Y[i].set(&X[i]);
	}


	//std::cout << " Hinv " << std::endl;
	for (int i = 0; i < _nBus; i++) {
		if (i > 0) {
			A[i].set(2, 0, 2 * ZsRe.get(i - 1, 0));
			A[i].set(2, 1, 2 * ZsIm.get(i - 1, 0));
			A[i].set(2, 3, -ZsNorm.get(i - 1, 0));
			A[i].set(2, 2, -1);
			A[i].set(2, 6, 1);
			A[i].set(0, 0, -1);
			A[i].set(1, 1, -1);
		}
		A[i].set(0, 4, 1);
		A[i].set(1, 5, 1);
		for (int j = 0; j < nChild.get(i, 0); j++) {
			int c = Childs[i].get(j, 0);
			A[i].set(0, 7 + 3 * j, 1); // Pci
			A[i].set(1, 8 + 3 * j, 1); // Qci
			A[i].set(0, 9 + 3 * j, -ZsRe.get(c - 1, 0)); // -R l
			A[i].set(1, 9 + 3 * j, -ZsIm.get(c - 1, 0)); // -X l
		}
		//A[i].display();
		MatrixCPU temp33(2 + 1 * (i > 0), 2 + 1 * (i > 0));
		MatrixCPU temp3M(2 + 1 * (i > 0), sizeOPFADMM.get(i, 0));
		MatrixCPU tempMM(sizeOPFADMM.get(i, 0), sizeOPFADMM.get(i, 0));
		temp33.multiplyTrans(&A[i], &A[i]); // (3*o_b) * (o_b*3) -> 9 * o_b
		temp33.invertGaussJordan(&temp33); // 3^3 = 27 (fixe !!!)
		temp3M.MultiplyMatMat(&temp33, &A[i]); // (3*3) * (3*o_b) -> 27 *o_b
		Hinv[i].multiplyTrans(&A[i], &temp3M, 0); // (o_b*3) * (3*o_b) -> 9 * o_b

		tempMM.setEyes(1);
		Hinv[i].subtract(&tempMM);
		Hinv[i].divide(_rho);
		//Hinv[i].display();
	}

	/*std::cout << "--------" << std::endl;
	std::cout << " Pn " << std::endl;
	Pn.display();
	PnMoy.display();
	std::cout << " PnTilde " << std::endl;
	PnTilde.display();
	for (int i = 0; i < 10; i++) {
		std::cout << " X " << i << std::endl;
		X[i].display();
		std::cout << " Y " << i << std::endl;
		Y[i].display();
		std::cout << " Q " << i << std::endl;
		Q[i].display();
		std::cout << " Mu " << i << std::endl;
		Mu[i].display();
	}*/
	updateChat();
	/*std::cout << " Chat " << std::endl;
	Chat.display();
	std::cout << " Bpt2 " << std::endl;
	Bpt2.display();
	std::cout << " Cp " << std::endl;
	Cost2.display();
	std::cout << " Ap2 " << std::endl;
	Cost1.display();
	std::cout << " Pmin " << std::endl;
	Pmin.display();
	std::cout << " Pma " << std::endl;
	Pmax.display();*/
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	//std::cout << "---------------------------------------------------------------------------------------" << std::endl;
}

void OPFADMM::solveConsensus(float eps, MatrixCPU* PSO)
{
	
	float epsG = eps;
	float epsL = eps/10;


	float fc = 0;
	float resG = 2 * epsG;
	float resL = 2 * epsL;
	_iterGlobal = 0;


	while ((_iterGlobal < _iterG) && (resG > epsG)) {
		_iterLocal = 0;
		resL = 2 * epsL;
		while (_iterLocal< _iterL && resL>epsL) {
			updateLocalProb();
			// FB 3
			if (!(_iterLocal % _stepL)) {
				resL = calcRes();		
			}
			_iterLocal++;
		}
		if (_iterLocal == _iterL) {
			std::cout << _iterGlobal << " " << _iterLocal << " " << resL << " " << resG << std::endl;
		}

		updateX();
		CommunicationX();

		updateGlobalProb();
		updateMu();
		updateChat();

		// FB 4
		if (!(_iterGlobal % _stepG)) {
			resG = updateRes(_iterGlobal / _stepG);
		}
		//std::cout << iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;

		_iterGlobal++;
	}
	std::cout << "OPF : " << _iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, (_iterGlobal - 1) / _stepG) << " " << resF.get(1, (_iterGlobal - 1) / _stepG) << " " << resF.get(2, (_iterGlobal - 1) / _stepG) << std::endl;
	PSO->set(&Pn);
}

void OPFADMM::initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO)
{
	// intitilisation des matrixs et variables 

	clock_t t = clock();
	std::cout << "init " << std::endl;
	_rho = sim.getRho();
	_rhoSO = rhoSO;

	if (_rhol == 0) {
		_rhol = _rho;
	}


	_iterG = sim.getIterG();
	_iterL = sim.getIterL();

	_stepG = sim.getStepG();
	_stepL = sim.getStepL();

	_nAgent = cas.getNagent();

	_nBus = cas.getNBus();
	_nLine = cas.getNLine(true); // ne doit pas �tre r�duit ici !!!

	std::cout << _nAgent << " " << _nBus << " " << _nLine << std::endl;
	_nAgentByBus = cas.getNagentByBus();
	nChild = MatrixCPU(_nBus, 1);
	CoresLineBus = cas.getCoresLineBus(true);
	_CoresBusAgent = cas.getCoresBusAgentLin(); // Cores[n] = b
	Ancestor = MatrixCPU(_nBus, 1, 0); // A_i = bus ant�c�dent de i
	PosChild = MatrixCPU(_nBus, 1, 0); // indice du bus i dans Child[Ai]
	Ancestor.set(0, 0, -1); // the slack bus has no ancestor
	ZsRe = cas.getZsRe();
	ZsIm = cas.getZsImag();
	ZsNorm = MatrixCPU(_nLine, 1);

	if (!chekcase()) {
		throw std::invalid_argument("not a radial case");
	}
	std::cout << "a" << std::endl;
	for (int lold = 0; lold < _nLine; lold++) {
		int l = lold + 1;
		int busTo = l;
		int busFrom = CoresLineBus.get(lold, 0);
		Ancestor.set(busTo, 0, busFrom);
		nChild.increment(busFrom, 0, 1);
		ZsNorm.set(lold, 0, ZsRe.get(lold, 0) * ZsRe.get(lold, 0) + ZsIm.get(lold, 0) * ZsIm.get(lold, 0));
	}



	_rhoInv = 1 / _rho;
	resF = MatrixCPU(3, (_iterG / _stepG) + 1);


	MatrixCPU lowerBound(cas.getLowerBound()); //voltage angle, voltage, line...
	MatrixCPU upperBound(cas.getUpperBound()); //voltage angle, voltage, line...


	std::cout << " local resolution " << std::endl;
	// local resolution
	tempN2 = MatrixCPU(2 * _nAgent, 1);
	tempB2 = MatrixCPU(2 * _nBus, 1);
	CoresSoloBusAgent = MatrixCPU(_nBus, 1, -1);
	Pn = MatrixCPU(2 * _nAgent, 1); // not the real agent
	Pmin = MatrixCPU(2 * _nAgent, 1, -1000000); // must not be the real one
	Pmax = MatrixCPU(2 * _nAgent, 1, 1000000); // idem
	

	PnTmin = MatrixCPU(2 * _nBus, 1, -1000000);
	PnTmax = MatrixCPU(2 * _nBus, 1, 1000000);

	Cost1 = MatrixCPU(2 * _nAgent, 1, rhoSO); //must be the consensus rhoPSO, pas important pour l'instant
	Cost2 = MatrixCPU(2 * _nAgent, 1);


	PnMoy = MatrixCPU(2 * _nBus, 1);
	PnPre = sim.getPn();
	MuL = MatrixCPU(2 * _nBus, 1);
	PnTilde = MatrixCPU(2 * _nBus, 1);
	Ap12 = Cost1;
	Ap12.add(_rhol);
	Bp1 = MatrixCPU(2 * _nAgent, 1);
	Bpt1 = MatrixCPU(2 * _nBus, 1);
	Bpt2 = MatrixCPU(2 * _nBus, 1);
	Apt1 = MatrixCPU(2 * _nBus, 1);
	Apt2 = MatrixCPU(2 * _nBus, 1);
	Apt12 = MatrixCPU(2 * _nBus, 1);

	if (Pn.max2() == 0) {
		Pn.add(&Pmin, &Pmax);
		Pn.divide(2);
	}

	_nAgentByBus.increment(0, 0, -1); // don't want the grid agent in this resolution
	//_nAgentByBus.display();
	for (int n = 1; n < _nAgent; n++) {
		int b = _CoresBusAgent.get(n, 0);
		PnTilde.increment(b, 0, Pn.get(n, 0));
		PnTilde.increment(b + _nBus, 0, Pn.get(n + _nAgent, 0));

		if (_nAgentByBus.get(b, 0) == 1) {
			CoresSoloBusAgent.set(b, 0, n);
		}
	}

	for (int b = 0; b < _nBus; b++) {
		int nb = _nAgentByBus.get(b, 0);
		if (nb > 0) {
			Apt1.set(b, 0, nb * _rhol);
			Apt1.set(b + _nBus, 0, nb * _rhol);

			Apt2.set(b, 0, nb * nb * _rho);
			Apt2.set(b + _nBus, 0, nb * nb * _rho);

			PnTilde.set(b, 0, PnTilde.get(b, 0) / nb);
			PnTilde.set(b + _nBus, 0, PnTilde.get(b + _nBus, 0) / nb);

			PnTmax.set(b, 0, PnTmax.get(b, 0) / nb);
			PnTmax.set(b + _nBus, 0, PnTmax.get(b + _nBus, 0) / nb);

			PnTmin.set(b, 0, PnTmin.get(b, 0) / nb);
			PnTmin.set(b + _nBus, 0, PnTmin.get(b + _nBus, 0) / nb);
		}

	}
	PnMoy.set(&PnTilde);
	Apt12.add(&Apt1, &Apt2);

	//std::cout << " creation " << std::endl;
	X = new MatrixCPU[_nBus];
	Ypre = new MatrixCPU[_nBus];
	Y = new MatrixCPU[_nBus];
	YTrans = new MatrixCPU[_nBus];
	Mu = new MatrixCPU[_nBus];

	tempN1 = MatrixCPU(_nAgent, 1);
	tempNN = MatrixCPU(_nAgent, _nAgent);
	//tempM1 = new MatrixCPU[_nAgent];
	tempM = new MatrixCPU[_nBus];


	Hinv = new MatrixCPU[_nBus];
	A = new MatrixCPU[_nBus];
	Q = new MatrixCPU[_nBus];

	Childs = new MatrixCPU[_nBus];

	Chat = MatrixCPU(_nBus, 4);
	VoltageLimit = MatrixCPU(_nBus, 2); // min, max
	VoltageLimitReal = MatrixCPU(_nBus, 2); // min, max
	sizeOPFADMM = MatrixCPU(_nBus, 1);


	int indice = 0;
	MatrixCPU nChildTemp(_nBus, 1, 0);
	for (int i = 0; i < _nBus; i++) {
		VoltageLimitReal.set(i, 0, lowerBound.get(_nBus + i, 0));
		VoltageLimitReal.set(i, 1, upperBound.get(_nBus + i, 0));
		VoltageLimit.set(i, 0, lowerBound.get(_nBus + i, 0) * lowerBound.get(_nBus + i, 0) * sqrt((nChild.get(i, 0) + 1) / 2));
		VoltageLimit.set(i, 1, upperBound.get(_nBus + i, 0) * upperBound.get(_nBus + i, 0) * sqrt((nChild.get(i, 0) + 1) / 2));

		Childs[i] = MatrixCPU(nChild.get(i, 0), 1);

		int sizeOPF = 3 * nChild.get(i, 0) + 7;
		int sizeOPF2 = 1 * nChild.get(i, 0) + 9;
		sizeOPFADMM.set(i, 0, sizeOPF);

		X[i] = MatrixCPU(sizeOPF, 1); // {Pi, Qi, vi, li, pi, qi, vAi, (Pci, Qci, lci) for all child Ci}
		Ypre[i] = MatrixCPU(sizeOPF, 1);// Y[i][j] not� dans l'article Yji est ce que i connait sur j
		Y[i] = MatrixCPU(sizeOPF, 1); //Y[i] = {Pi, Qi, vi, li, pi, qi, vAi, (Pci, Qci, lci) for all child Ci}
		YTrans[i] = MatrixCPU(sizeOPF2, 1);   // (Pi,Qi,vi,li,pi,qi,pai,qai,lai,(vai) for all child Ci) !!!
		Mu[i] = MatrixCPU(sizeOPF, 1);
		A[i] = MatrixCPU(2 + 1 * (i > 0), sizeOPF);
		Hinv[i] = MatrixCPU(sizeOPF, sizeOPF);
		Q[i] = MatrixCPU(sizeOPF, 1, 0);
		tempM[i] = MatrixCPU(sizeOPF, 1);
	}
	for (int i = 0; i < _nBus; i++) {
		if (i > 0) {
			int Ai = Ancestor.get(i, 0);
			Childs[Ai].set(nChildTemp.get(Ai, 0), 0, i);
			PosChild.set(i, 0, nChildTemp.get(Ai, 0));
			nChildTemp.increment(Ai, 0, 1);
		}
	}



	std::cout << " init valeur " << std::endl;
	for (int i = 0; i < _nBus; i++) {
		// vi = 1; vai = 1
		X[i].set(2, 0, 1);
		X[i].set(6, 0, 1);

		// pi & qi
		int Nb = _nAgentByBus.get(i, 0);
		X[i].set(4, 0, Nb * PnTilde.get(i, 0)); // pi = sum(pn)
		X[i].set(5, 0, Nb * PnTilde.get(i + _nBus, 0)); // qi = sum(pn)
	}
	DFSP(0); // Pi
	X[0].set(0, 0, 0);
	DFSQ(0); // Qi
	X[0].set(1, 0, 0);

	for (int i = 0; i < _nBus; i++) {
		float Si = X[i].get(0, 0) * X[i].get(0, 0) + X[i].get(1, 0) * X[i].get(1, 0);
		X[i].set(3, 0, Si / X[i].get(2, 0)); // li = Si^2/vi
	}
	for (int i = 0; i < _nBus; i++) {
		int m = nChild.get(i, 0);
		for (int j = 0; j < m; j++) {
			// (Pci, Qci, lci) for all child Ci

			int c = Childs[i].get(j, 0);
			X[i].set(7 + 3 * j, 0, X[c].get(0, 0));
			X[i].set(8 + 3 * j, 0, X[c].get(1, 0));
			X[i].set(9 + 3 * j, 0, X[c].get(3, 0));
		}
		Y[i].set(&X[i]);
	}


	std::cout << " Hinv " << std::endl;
	for (int i = 0; i < _nBus; i++) {
		if (i > 0) {
			A[i].set(2, 0, 2 * ZsRe.get(i - 1, 0));
			A[i].set(2, 1, 2 * ZsIm.get(i - 1, 0));
			A[i].set(2, 3, -ZsNorm.get(i - 1, 0));
			A[i].set(2, 2, -1);
			A[i].set(2, 6, 1);
			A[i].set(0, 0, -1);
			A[i].set(1, 1, -1);
		}
		A[i].set(0, 4, 1);
		A[i].set(1, 5, 1);
		for (int j = 0; j < nChild.get(i, 0); j++) {
			int c = Childs[i].get(j, 0);
			A[i].set(0, 7 + 3 * j, 1); // Pci
			A[i].set(1, 8 + 3 * j, 1); // Qci
			A[i].set(0, 9 + 3 * j, -ZsRe.get(c - 1, 0)); // -R l
			A[i].set(1, 9 + 3 * j, -ZsIm.get(c - 1, 0)); // -X l
		}
		//A[i].display();
		MatrixCPU temp33(2 + 1 * (i > 0), 2 + 1 * (i > 0));
		MatrixCPU temp3M(2 + 1 * (i > 0), sizeOPFADMM.get(i, 0));
		MatrixCPU tempMM(sizeOPFADMM.get(i, 0), sizeOPFADMM.get(i, 0));
		temp33.multiplyTrans(&A[i], &A[i]); // (3*o_b) * (o_b*3) -> 9 * o_b^2
		temp33.invertGaussJordan(&temp33); // 3^3 = 27 (fixe !!!)
		temp3M.MultiplyMatMat(&temp33, &A[i]); // (3*3) * (3*o_b) -> 27 *o_b
		Hinv[i].multiplyTrans(&A[i], &temp3M, 0); // (o_b*3) * (3*o_b) -> 9 * o_b

		tempMM.setEyes(1);
		Hinv[i].subtract(&tempMM);
		Hinv[i].divide(_rho);
		//Hinv[i].display();
	}

	
	updateChat();
	
	std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	std::cout << "---------------------------------------------------------------------------------------" << std::endl;


}

void OPFADMM::updateConsensus(MatrixCPU* Pmarket)
{
	for (int n = 1; n < _nAgent; n++) { // pas l'agent des pertes
		float eta = 0.5 * (Pn.get(n, 0) - Pmarket->get(n, 0));
		etaSO.set(n, 0, etaSO.get(n, 0) + eta);
		eta = 0.5 * (Pn.get(n + _nAgent, 0) - Pmarket->get(n + _nAgent, 0));
		etaSO.set(n + _nAgent, 0, etaSO.get(n + _nAgent, 0) + eta);
	}



	Cost2.add(&Pn, Pmarket);
	Cost2.set(0, 0, 0);
	Cost2.set(_nAgent, 0, 0);
	Cost2.multiply(-0.5);
	Cost2.add(&etaSO);
	Cost2.multiply(_rhoSO);

}

void OPFADMM::updateGlobalProb() {
	
	
	for (int i = 0; i < _nBus; i++) {
		//std::cout << " Q " << std::endl;
		//Q[i].display();
		Ypre[i].swap(&Y[i]);
		Y[i].MultiplyMatVec(&Hinv[i], &Q[i]); // solve system by using the inverse
	}
	Y[0].set(2, 0, 1);
	Y[0].set(6, 0, 1);


	// communication of y, mu

}

void OPFADMM::updateLocalProb() {
	// FB 11a
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	updateBp1();
	updatepl();

#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 1, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 11b
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	PnMoy.set(0.0);
	for (int n = 0; n < _nAgent; n++) {
		int bus = _CoresBusAgent.get(n, 0);
		PnMoy.increment(bus, 0, Pn.get(n, 0));
		PnMoy.increment(bus + _nBus, 0, Pn.get(n + _nAgent, 0));

	}
	for (int b = 0; b < _nBus; b++) {
		int nb = _nAgentByBus.get(b, 0);
		if (nb > 0) {
			PnMoy.set(b, 0, PnMoy.get(b, 0) / nb );
			PnMoy.set(b + _nBus, 0, PnMoy.get(b + _nBus, 0) / nb);
		}
		
	}

#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 2, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 11c
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	for (int b = 0; b < _nBus; b++) {
		int nb = _nAgentByBus.get(b, 0);
		if (nb > 0) {
			Bpt1.set(b, 0, MuL.get(b, 0) + PnMoy.get(b, 0));
			Bpt1.set(b +_nBus, 0, MuL.get(b + _nBus, 0) + PnMoy.get(b + _nBus, 0));
		}
	}
	
	updatePTilde();
	float mu = 0;
	for (int b = 0; b < _nBus; b++) {
		int nb = _nAgentByBus.get(b, 0);
		if (nb > 0) {
			mu = MuL.get(b, 0) + PnMoy.get(b, 0) - PnTilde.get(b, 0);
			MuL.set(b, 0, mu);
			mu = MuL.get(b + _nBus, 0) + PnMoy.get(b + _nBus, 0) - PnTilde.get(b + _nBus, 0);
			MuL.set(b + _nBus, 0, mu);
		}
	}
	
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 3, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
	
}

void OPFADMM::updateX()
{
	double x1, x2, x3, x4, c1, c2, c3, c4, lambdaLo, lambdaUp, x3min, x3max, gamma, k2;
	double c1122;
	int nSol = 0;
	int typeSol = 0;
	int BestRoot = 0;
	double bestGamma = -1;
	double p = 0;
	int nRoot = 0;

	for (int i = 0; i < _nBus; i++) {

		bool goodSol = false;
		k2 = sqrt(2.0 / (nChild.get(i, 0) + 1));
		if (i == 0) { // slack bus
			goodSol = true;
			c3 = -2 * Chat.get(i, 2) / k2;
			x1 = 0;
			x2 = 0;
			x4 = 0;
			x3 = 1 / k2;
			gamma = 0;
		}
		else {
			c1 = -2 * Chat.get(i, 0);
			c2 = -2 * Chat.get(i, 1);
			c3 = -2 * Chat.get(i, 2) / k2;
			c4 = -2 * Chat.get(i, 3);
			c1122 = c1 * c1 + c2 * c2;
			
			
			x3min = VoltageLimit.get(i, 0);
			x3max = VoltageLimit.get(i, 1);
			
			// case without constraint
			
			x1 = -c1 / 2;
			x2 = -c2 / 2;
			x3 = -c3 / 2;
			x4 = -c4 / 2;
			lambdaUp = 0;
			lambdaLo = 0;
			
			if (x3 < x3min) {
				x3 = x3min;
				lambdaLo = (2 * x3 + c3);
			}
			else if (x3 > x3max) {
				x3 = x3max;
				lambdaUp = -(2 * x3 + c3);
			}
			gamma = k2 * x4 - (x1 * x1 + x2 * x2) / x3; // ce n'est pas vraiment gamma, doit �tre positif
			//std::cout << "x 1 : " << x1 << " " << x2 << " " << x3 * k2 << " " << x4 << " " << (x1 * x1 + x2 * x2) / x3  - k2 * x4 << std::endl;

			if (gamma >= 0) {
				// the solution is good !
				goodSol = true;
			}
			else {
				if (gamma > bestGamma) {
					typeSol = 1;
					bestGamma = gamma;
				}
			}
		}
		if (!goodSol) { // cas d�g�n�r�
			if (c1122 == 0) {
				std::cout << " bus " << i << " : c1= " << c1 << " c2=" << c2 << " c4=" << c4 << " gamma= " << gamma << std::endl;

				x4 = 0;
				goodSol = true;
			}
		}
		// case x3 = x3max lambdaLo = 0
		if (!goodSol) {
			x3 = x3max;

			coefPoly2[0] = 2 * (c4 / (k2 * x3) + 1);
			coefPoly2[1] = 1 / x3;
			coefPoly2[0] = coefPoly2[0] * k2 * k2 / (4 * c1122);
			coefPoly2[1] = coefPoly2[1] * k2 * k2 / (4 * c1122);

			nRoot = resolveRealPolynome3without2term(root2, coefPoly2);
			
			for (int n = 0; n < nRoot; n++) {
				p = root2[n];
								
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;
				lambdaUp = -(2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));
				//std::cout << "x2 : " << x1 << " " << x2 << " " << x3 * k2 << " " << x4 << " " << gamma << " " << lambdaUp << std::endl;
				if (gamma >=  0 && lambdaUp >=  0) {
					// the solution is good 
					goodSol = true;
					//nSol = n;
					break;
				}
				if (gamma > bestGamma && lambdaUp > bestGamma) {
					typeSol = 2;
					bestGamma = MYMIN(gamma, lambdaLo);
					BestRoot = n;
				}

			}
			// case x3 = x3min lambdaUp = 0
			if (!goodSol) {
				x3 = x3min;
				 
				coefPoly2[0] = 2 * (c4 / (k2 * x3) + 1);
				coefPoly2[1] = 1 / x3;
				coefPoly2[0] = coefPoly2[0] * k2 * k2 / (4 * c1122);
				coefPoly2[1] = coefPoly2[1] * k2 * k2 / (4 * c1122);

				nRoot = resolveRealPolynome3without2term(root3, coefPoly2);

				for (int n = 0; n < nRoot; n++) {
					p = root3[n];
					//std::cout << "poly " << coefPoly2[0] * p + coefPoly2[1] + p * p * p << std::endl;
					x1 = p * c1 * x3;
					x2 = p * c2 * x3;
					x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
					gamma = (2 * x4 + c4) / k2;
					lambdaLo = (2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));
					
					//std::cout << "x 3: " << x1 << " " << x2 << " " << x3 * k2 << " " << x4 << " " << gamma << " " << lambdaLo << std::endl;

					if (gamma >= 0 && lambdaLo >= 0) {
						// the solution is good !
						goodSol = true;
						break;
					}
					if (gamma > bestGamma && lambdaLo > bestGamma) {
						typeSol = 3;
						bestGamma = MYMIN(gamma, lambdaLo);
						BestRoot = n;
					}
				}
			}
			// case xmin<x3<xmax lambdaLo = 0 lambdaUp = 0
			if (!goodSol) {
				 
				coefPoly3[0] = c1122 / k2 * (2 * c3 / k2 - c4);
				coefPoly3[1] = (c3 - 2 * c4 / k2);
				coefPoly3[2] = -1;
				coefPoly3[0] = coefPoly3[0] * k2 * k2 / (c1122 * c1122);
				coefPoly3[1] = coefPoly3[1] * k2 * k2 / (c1122 * c1122);
				coefPoly3[2] = coefPoly3[2] * k2 * k2 / (c1122 * c1122);

				nRoot = resvolveRealPolynome4without2term(root4, coefPoly3, Lagrange);

				for (int n = 0; n < nRoot; n++) {
					p = root4[n];
					//std::cout << "poly " <<p * p * p * p + coefPoly3[0] * p*p*p + coefPoly3[1]*p + coefPoly3[2] << std::endl;
					x3 = -(c1122 * p + 2 * c3) / (2 * (c1122 * p * p + 2));
					x1 = p * c1 * x3;
					x2 = p * c2 * x3;
					x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
					gamma = (2 * x4 + c4) / k2;
					//std::cout << "x 4: " << x1 << " " << x2 << " " << x3 * k2 << " " << x4 << " " << gamma << " " << std::endl;

					if (gamma >= 0 && x3 <= x3max && x3 >= x3min) {
						// the solution is good !
						goodSol = true;
						break;
					}if (gamma > bestGamma && (x3max - x3) > bestGamma && (x3 - x3min) > bestGamma) {
						typeSol = 4;
						bestGamma = MYMIN(MYMIN(gamma, (x3max - x3)), (x3 - x3min));
						BestRoot = n;
					}
				}
			}
		}
		if (!goodSol) {
			if (typeSol == 1) {
				// case without constraint
				x1 = -c1 / 2;
				x2 = -c2 / 2;
				x3 = -c3 / 2;
				x3 = (x3max - x3) * (x3 > x3max) + (x3min - x3) * (x3min > x3) + x3;
				x4 = -c4 / 2; // ou  (x1 * x1 + x2 * x2) / (k2 * x3)
			}
			else {
				if (typeSol == 2) {
					x3 = x3max;
					p = root2[BestRoot];
				}
				else if (typeSol == 3) {
					x3 = x3min;
					p = root3[BestRoot];
				}
				else if (typeSol == 4) {
					p = root4[BestRoot];
					x3 = -(c1122 * p + 2 * c3) / (2 * (c1122 * p * p + 2));
					x3 = (x3max - x3) * (x3 > x3max) + (x3min - x3) * (x3min > x3) + x3;
				}
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
			}
		}
		
		// X =  {Pi, Qi, vi, li, pi, qi, vAi, (Pci, Qci, lci) for all child Ci}
		X[i].set(0, 0, x1);
		X[i].set(1, 0, x2);
		X[i].set(2, 0, x3* k2);
		X[i].set(3, 0, x4);
		X[i].set(4, 0, PnTilde.get(i, 0) * _nAgentByBus.get(i, 0));
		X[i].set(5, 0, PnTilde.get(i + _nBus, 0) * _nAgentByBus.get(i, 0));

		//std::cout << "x F : " << x1 << " " << x2 << " " << x3*k2 << " " << x4 << " " << gamma << std::endl;
		
	}

}



void OPFADMM::updateMu()
{
	for (int i = 0; i < _nBus; i++) {
		tempM[i].subtract(&X[i], &Y[i]);
		tempM[i].multiply(_rho);
		if (consensus) {
			tempM[i].divide(2);
		}
		Mu[i].add(&tempM[i]);
	}
}

void OPFADMM::updateBp1()
{
	for (int i = 0; i < _nAgent; i++) {
		int bus = _CoresBusAgent.get(i, 0);
		int Nb = _nAgentByBus.get(bus, 0);
		if (Nb > 1) {
			float m = Pn.get(i, 0) - PnMoy.get(bus, 0) + PnTilde.get(bus, 0) - MuL.get(bus, 0);
			Bp1.set(i, 0, m);
			m = Pn.get(i + _nAgent, 0) - PnMoy.get(bus + _nBus, 0) + PnTilde.get(bus + _nBus, 0) - MuL.get(bus + _nBus, 0);
			Bp1.set(i + _nAgent, 0, m);
		}
		
	}
}



void OPFADMM::updatepl()
{
	float pn = 0;
	PnPre.swap(&Pn);
	for (int n = 0; n < _nAgent; n++) {
		int b = _CoresBusAgent.get(n, 0);
		int Nb = _nAgentByBus.get(b, 0);
		if (Nb == 1) { // � ne faire qu'� la premier it�ration...
			if (_iterLocal == 0) {
				float ub = Pmax.get(n, 0);
				float lb = Pmin.get(n, 0);
				pn = (_rho * Bpt2.get(b, 0) - Cost2.get(n, 0)) / ((_rho + Cost1.get(n, 0)));
				pn = (ub - pn) * (pn > ub) + (lb - pn) * (pn < lb) + pn;
				Pn.set(n, 0, pn);
				PnPre.set(n, 0, pn);

				ub = Pmax.get(n + _nAgent, 0);
				lb = Pmin.get(n + _nAgent, 0);
				pn = (_rho * Bpt2.get(b + _nBus, 0) - Cost2.get(n + _nAgent, 0)) / ((_rho + Cost1.get(n + _nAgent, 0)));
				pn = (ub - pn) * (pn > ub) + (lb - pn) * (pn < lb) + pn;
				//std::cout << _rho << " "<< Bpt2.get(b + _nBus, 0)<< " "<< Cost2.get(n + _nAgent, 0) << " "<< Cost1.get(n + _nAgent, 0)<<" "<< lb << " "<< ub << " " <<   pn << std::endl;
				Pn.set(n + _nAgent, 0, pn);
				PnPre.set(n + _nAgent, 0, pn);
			}
		}
		else {
			float ub = Pmax.get(n, 0);
			float lb = Pmin.get(n, 0);
			pn = (Bp1.get(n, 0) * _rhol - Cost2.get(n, 0)) / Ap12.get(n, 0);
			pn = (ub - pn) * (pn > ub) + (lb - pn) * (pn < lb) + pn;
			Pn.set(n, 0, pn);

			ub = Pmax.get(n + _nAgent, 0);
			lb = Pmin.get(n + _nAgent, 0);
			pn = (Bp1.get(n + _nAgent, 0) * _rhol - Cost2.get(n + _nAgent, 0)) / Ap12.get(n + _nAgent, 0);
			pn = (ub - pn) * (pn > ub) + (lb - pn) * (pn < lb) + pn;
			Pn.set(n + _nAgent, 0, pn);
		}
	}
	/*
	Pn.set(&Bp1);
	Pn.multiply(_rhol);

	Pn.subtract(&Cost2);
	Pn.divideT(&Ap12);
	Pn.project(&Pmin, &Pmax);*/
}

void OPFADMM::updatePTilde()
{
	float pn = 0;
	for (int b = 0; b < _nBus; b++) {
		int n = CoresSoloBusAgent.get(b, 0);
		if (n != -1) {
			if (_iterLocal == 0) {
				PnTilde.set(b, 0, PnMoy.get(b, 0));
				PnTilde.set(b + _nBus, 0, PnMoy.get(b + _nBus, 0));
			}
		}
		else {
			int nb = _nAgentByBus.get(b, 0);
			if (nb > 0) {
				pn = (Bpt1.get(b, 0) * Apt1.get(b, 0) + Bpt2.get(b, 0) * Apt2.get(b, 0)) / Apt12.get(b, 0);
				PnTilde.set(b, 0, pn);
				pn = (Bpt1.get(b + _nBus, 0) * Apt1.get(b + _nBus, 0) + Bpt2.get(b + _nBus, 0) * Apt2.get(b + _nBus, 0)) / Apt12.get(b + _nBus, 0);
				PnTilde.set(b + _nBus, 0, pn);
			}
		}
	}

	/*PnTilde.set(&Bpt1);
	PnTilde.multiplyT(&Apt1);
	tempB2.set(&Bpt2);
	tempB2.multiplyT(&Apt2);
	PnTilde.add(&tempB2);

	PnTilde.divideT(&Apt12);

	PnTilde.project(&PnTmin, &PnTmax);*/
}

float OPFADMM::calcRes()
{
	float d1 = Pn.max2(&PnPre);
	float d2 = PnMoy.max2(&PnTilde);

	return d1 * (d1 > d2) + d2 * (d2 >= d1);
}



void OPFADMM::updateChat()
{
	//Ytrans (Pi, Qi, vi, li, pi, qi, Pai, Qai, lai, vij) !!!
	//Y      (Pi, Qi, vi, li, pi, qi, vai, Pji, Qji, lji)
	for (int i = 0; i < _nBus; i++) {
		
		float Phat, Qhat, lhat, phat, qhat;
		float vhat = 0;
		float muhat = 0;
		float m = nChild.get(i, 0);
		// est ce que c'est y - MU/rho ou y + MU/rho (comme dans le papier mais le signe change je ne sais pas pourquoi)
		if (consensus) {
			Phat = (X[i].get(0, 0) + Y[i].get(0, 0)) / 4 - Mu[i].get(0, 0) / (2 * _rho);
			Qhat = (X[i].get(1, 0) + Y[i].get(1, 0)) / 4 - Mu[i].get(1, 0) / (2 * _rho);
			lhat = (X[i].get(3, 0) + Y[i].get(3, 0)) / 4 - Mu[i].get(3, 0) / (2 * _rho);
			if (i > 0) {
				int Ai = Ancestor.get(i, 0);
				int c = PosChild.get(i, 0);
				Phat += (X[Ai].get(7 + 3 * c, 0) + Y[Ai].get(7 + 3 * c, 0)) / 4 - Mu[Ai].get(7 + 3 * c, 0) / (2 * _rho); // or Y[Ai].get(7 + 3 * c, 0)
				Qhat += (X[Ai].get(8 + 3 * c, 0) + Y[Ai].get(8 + 3 * c, 0)) / 4 - Mu[Ai].get(8 + 3 * c, 0) / (2 * _rho);
				lhat += (X[Ai].get(9 + 3 * c, 0) + Y[Ai].get(9 + 3 * c, 0)) / 4 - Mu[Ai].get(9 + 3 * c, 0) / (2 * _rho);
			}
			phat = (X[i].get(4, 0) + Y[i].get(4, 0)) / 2 - Mu[i].get(4, 0) / _rho;
			qhat = (X[i].get(5, 0) + Y[i].get(5, 0)) / 2 - Mu[i].get(5, 0) / _rho;
			for (int j = 0; j < m; j++) { // Vai of all childs
				int c = Childs[i].get(j, 0);
				vhat += Y[c].get(6, 0) + X[c].get(6, 0);
				muhat += Mu[c].get(6, 0);
			}
			vhat = (Y[i].get(2, 0) + X[i].get(2, 0) + vhat) / (2 * m + 2) - (Mu[i].get(2, 0) + muhat) / (_rho * (m + 1));
		}
		else {
			Phat = Y[i].get(0, 0) / 2 - Mu[i].get(0, 0) / (2 * _rho);
			Qhat = Y[i].get(1, 0) / 2 - Mu[i].get(1, 0) / (2 * _rho);
			lhat = Y[i].get(3, 0) / 2 - Mu[i].get(3, 0) / (2 * _rho);
			if (i > 0) {
				int Ai = Ancestor.get(i, 0);
				int c = PosChild.get(i, 0);
				Phat += Y[Ai].get(7 + 3 * c, 0) / 2 - Mu[Ai].get(7 + 3 * c, 0) / (2 * _rho); // or Y[Ai].get(7 + 3 * c, 0)
				Qhat += Y[Ai].get(8 + 3 * c, 0) / 2 - Mu[Ai].get(8 + 3 * c, 0) / (2 * _rho);
				lhat += Y[Ai].get(9 + 3 * c, 0) / 2 - Mu[Ai].get(9 + 3 * c, 0) / (2 * _rho);
			}
			phat = Y[i].get(4, 0) - Mu[i].get(4, 0) / _rho;
			qhat = Y[i].get(5, 0) - Mu[i].get(5, 0) / _rho;
			for (int j = 0; j < m; j++) { // Vai of all childs
				int c = Childs[i].get(j, 0);
				muhat += Mu[c].get(6, 0);
				vhat += Y[c].get(6, 0);
			}
			vhat = (Y[i].get(2, 0) + vhat) / (m + 1) - (Mu[i].get(2, 0) + muhat) / (_rho * (m + 1));
		}
				

		Chat.set(i, 0, Phat);
		Chat.set(i, 1, Qhat);
		Chat.set(i, 2, vhat);
		Chat.set(i, 3, lhat);
		int nA = _nAgentByBus.get(i, 0);
		if (nA > 0) {
			Bpt2.set(i, 0, phat / nA);
			Bpt2.set(i + _nBus, 0, qhat / nA);
		}
		
		//Chat.set(i, 4, phat);
		//Chat.set(i, 5, qhat);
	}


}

void OPFADMM::CommunicationX()
{
/**/ // X = { Pi, Qi, vi, li, pi, qi, vAi, (Pci, Qci, lci) for all child Ci }
	
	for (int i = 0; i < _nBus; i++) {

		if (i > 0) {
			int Ai = Ancestor.get(i, 0);
			X[i].set(6, 0, X[Ai].get(2, 0));
		}

		int m = nChild.get(i, 0);
		for (int j = 0; j < m; j++) {
			int c = Childs[i].get(j, 0);
			X[i].set(7 + 3 * j, 0, X[c].get(0, 0));
			X[i].set(8 + 3 * j, 0, X[c].get(1, 0));
			X[i].set(9 + 3 * j, 0, X[c].get(3, 0));

		}
	}
	//Ytrans (Pi, Qi, vi, li, pi, qi, pai, qai, lai, vij) !!!
	//Y      (Pi, Qi, vi, li, pi, qi, vai, Pji, Qji, lji)
	// Q udate in argmin 0.5yHy + Qy

	for (int i = 0; i < _nBus; i++) {
		for (int j = 0; j < sizeOPFADMM.get(i, 0); j++) {
			if (consensus) {
				Q[i].set(j,0, -(Mu[i].get(j, 0) + _rho * (X[i].get(j, 0) + Y[i].get(j, 0)) / 2));
			}
			else {
				Q[i].set(j, 0, -(Mu[i].get(j, 0) + _rho * X[i].get(j, 0)));
			}
			
			
		}
	}

}

void OPFADMM::CommunicationY()
{
	//Ytrans (Pi, Qi, vi, li, pi, qi, pai, qai, lai, vij) !!!
	//Y      (Pi, Qi, vi, li, pi, qi, vai, Pji, Qji, lji)
	for (int i = 0; i < _nBus; i++) {
		for (int j = 0; j < 6; j++) {
			YTrans[i].set(j, 0, Y[i].get(j, 0));
		}

		if (i > 0) { // recieve yiAi from ancestor
			int Ai = Ancestor.get(i, 0);
			int m = PosChild.get(i, 0);
			YTrans[i].set(6, 0, Y[Ai].get(7 + m * 3, 0)); //pai = pij
			YTrans[i].set(7, 0, Y[Ai].get(8 + m * 3, 0)); //qai = qij
			YTrans[i].set(8, 0, Y[Ai].get(9 + m * 3, 0)); //lai = lij

		}
		for (int n = 0; n < nChild.get(i, 0); n++) {// recieve yij from child
			int c = Childs[i].get(n, 0);
			YTrans[i].set(9 + n, 0, Y[c].get(6, 0)); // vij = vAi
		}
	}
}

float OPFADMM::updateRes(int indice) 
{
	float resS = 0;
	float resR = 0;
	float resV = 0;
	for (int i = 0; i < _nBus; i++) {
		
		float resTempS = Y[i].max2(&Ypre[i]);
		float resTempR = Y[i].max2(&X[i]);

		if (resTempS > resS) {
			resS = resTempS;
		}
		if (resTempR > resR) {
			resR = resTempR;
		}

	}
	float oldrho = _rho;
	resF.set(0, indice, resR);
	resF.set(1, indice, oldrho * resS);
	resF.set(2, indice, resV);

	if (_tau > 1) {
		if (resR > _mu * resS) {
			_rho = _tau * _rho;
			Apt2.multiply(_tau);
			for (int i = 0; i < _nBus; i++) {
				Hinv[i].divide(_tau);
			}
			Apt12.add(&Apt1, &Apt2);
			//std::cout << _iterGlobal << "rho augmente " << _rho << std::endl;
		}
		else if (resS > _mu * resR) {// rho = rho / tau_inc;
			_rho = _rho / _tau;
			Apt2.divide(_tau);
			for (int i = 0; i < _nBus; i++) {
				Hinv[i].multiply(_tau);
			}
			Apt12.add(&Apt1, &Apt2);
			//std::cout << _iterGlobal << "rho diminue " << _rho << std::endl;
		}/**/
	}
	


	return MYMAX(MYMAX(resV, oldrho *resS),  resR);
}

int OPFADMM::feasiblePoint()
{

	MatrixCPU test(_nBus, 1, -1);
	int counter = 0;
	for (int bus = 0; bus < _nBus; bus++) {
		
		float Si = X[bus].get(0, 0) * X[bus].get(0, 0) + X[bus].get(1, 0) * X[bus].get(1, 0);
		float li = X[bus].get(3, 0);
		float vi = X[bus].get(2, 0);
		float err = Si - li * vi;
		test.set(bus, 0, err);
		if (abs(err) > 0.0001) {
			counter++;
		}
	}
	//std::cout << " erreur sur la relaXation " << test.max2() << " " << counter << std::endl;
	//test.display();
	resF.set(2, (_iterGlobal - 1) / _stepG, test.max2());
	return counter;
}

// (Pi,Qi,vi,li,pi,qi,vai,pji,qji,lji)
 MatrixCPU OPFADMM::getPb(){
	MatrixCPU Pb(2*_nBus, 1);
	for (int i = 0; i < _nBus; i++)
	{
		Pb.set(i, 0, Y[i].get(4,0));
		Pb.set(i + _nBus, 0, Y[i].get(5,0));
	}
	return Pb;
	
}
 MatrixCPU OPFADMM::getPhi(){
	MatrixCPU Phi(2*_nLine, 1);
	for (int i = 0; i < _nLine; i++)
	{
		Phi.set(i, 0, Y[i + 1].get(0,0));
		Phi.set(i + _nLine, 0, Y[i + 1].get(1,0));
	}
	return Phi;
}
 MatrixCPU OPFADMM::getE(){
	MatrixCPU E(2*_nBus, 1);
	for (int i = 0; i < _nBus; i++)
	{
		E.set(i, 0, Y[i].get(3,0)); //l_i
		E.set(i + _nBus, 0, Y[i].get(2,0)); // v_i
	}
	return E;
}

void OPFADMM::display() {

	std::cout.precision(3);

	if (_iterGlobal == 0) {
		std::cout << "algorithm not launch" << std::endl;
	}
	else if (_iterGlobal < _iterG) {
		std::cout << "method " << _name << " converged in " << _iterGlobal << " iterations." << std::endl;
		std::cout << "Converged in " << (float) timeOPF / CLOCKS_PER_SEC << " seconds" << std::endl;

	}
	else {
		std::cout << "method " << _name << " not converged in " << _iterGlobal << " iterations." << std::endl;
		std::cout << "time taken " << (float) timeOPF / CLOCKS_PER_SEC << " seconds" << std::endl;
	}
	std::cout << "The power error of this state is (constraint) " << resF.get(0, _iterGlobal / _stepG) << " and convergence " << resF.get(1, _iterGlobal / _stepG) << std::endl;
	std::cout << "===============================================================|" << std::endl;
	std::cout << "      System Summary                                           |" << std::endl;
	std::cout << "===============================================================|" << std::endl;
	std::cout << "Buses            " << _nBus << std::endl;
	std::cout << "Branches         " << _nLine << std::endl;
	std::cout << "Agent            " << _nAgent << std::endl;
	std::cout << "Ploss            " << getPLoss() << std::endl;
	std::cout << "Qloss            " << getQLoss() << std::endl;


	std::cout << std::endl << std::endl;
	
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "      Bus Data                                                                                          |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Bus |    Voltage  |   Power = Generation  + Load    |                Mu voltage and power              |" << std::endl;
	std::cout << "  #  |     Mag(pu) |    P (pu)      |      Q (pu)    |     V (pu)     |      P (pu)    |      Q (pu)    |" << std::endl;
	std::cout << "-----|-------------|----------------|----------------|----------------|----------------|----------------|" << std::endl;

		
	float seuil = 0.0001;
		
	for (int b = 0; b < _nBus; b++) {
	std::cout << std::setw(5) << b << "|" << std::setw(12) << sqrt(X[b].get(2,0)) << " |" << std::setw(16)
			<< (abs(X[b].get(4, 0)) > seuil) * X[b].get(4, 0) << "|" << std::setw(16) << (abs(X[b].get(5, 0)) > seuil) * X[b].get(5, 0)
			<< "|" << std::setw(16) << Mu[b].get(2, 0) << "|" << std::setw(16)
			<< Mu[b].get(4, 0) << "|" << std::setw(16) << Mu[b].get(5, 0) << "|" << std::endl;

	}
	std::cout << std::endl << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "      Line Data                                                                                         |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Line |    From     |    To      |                           Upstream flow                              |" << std::endl;
	std::cout << "  #   |    Bus      |    Bus     |    P (pu)      |    Q (pu)      |     l (pu)     |     Loss (pu)     |" << std::endl;
	std::cout << "------|-------------|------------|----------------|----------------|----------------|-------------------|" << std::endl;

	for (int l = 0; l < _nLine; l++) {
		int b = l + 1;
		std::cout << std::setw(6) << l << "|" << std::setw(12) << CoresLineBus.get(l, 0) << " |" << std::setw(12)
			<< CoresLineBus.get(l, 1) << "|" << std::setw(16) << X[b].get(0, 0)
			<< "|" << std::setw(16) << X[b].get(1, 0) << "|" << std::setw(16)
			<< X[b].get(3, 0) << "|" << std::setw(19) << X[b].get(3, 0) * ZsRe.get(l, 0) << "|" << std::endl;
	}
	std::cout << std::endl << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "     Constraints                                                                                        |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Bus | Voltage | Voltage | Voltage |        Power Injection          |          Power Injection         |" << std::endl;
	std::cout << "  #  | Mag(pu) | MIN(pu) |  MYMAX(pu)|  P (pu) | Pmin (pu) | Pmax (pu) |  Q (pu)  | Qmin (pu) | Qmax (pu) |" << std::endl;
	std::cout << "-----|---------|---------|---------|---------|-----------|-----------|----------|-----------|-----------|" << std::endl;
	

	for (int b = 0; b < _nBus; b++) {
		int nb = _nAgentByBus.get(b, 0);
		std::cout << std::setw(5) << b << "|" << std::setw(8) << sqrt(Y[b].get(2, 0)) << " |" << std::setw(9)
			<< VoltageLimitReal.get(b, 0) << "|" << std::setw(9) << VoltageLimitReal.get(b, 1)
			<< "|" << std::setw(9) << Y[b].get(4, 0) << "|" << std::setw(11)
			<< PnTmin.get(b, 0) * nb << "|" << std::setw(11) << PnTmax.get(b, 0) * nb << "|" << std::setw(10) << Y[b].get(5, 0) 
			<< "|" << std::setw(11) << PnTmin.get(b + _nBus, 0) * nb << "|" << std::setw(11) << PnTmax.get(b + _nBus, 0) * nb << "|" << std::endl;

	}
	std::cout << std::endl << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "      Agent Data                                                                                        |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Agent |  Bus  |  Cost   |  Cost   |        Power Injection          |          Power Injection         |" << std::endl;
	std::cout << "  #    |   #   |  a (pu) |  b (pu) |  P (pu) | Pmin (pu) | Pmax (pu) |  Q (pu)  | Qmin (pu) | Qmax (pu) |" << std::endl;
	std::cout << "-------|-------|---------|---------|---------|-----------|-----------|----------|-----------|-----------|" << std::endl;

	for (int n = 0; n < _nAgent; n++) {
		int b = _CoresBusAgent.get(n, 0);
		std::cout << std::setw(7) << n << "|" << std::setw(7) << b << "|" << std::setw(8) << Cost1.get(n,0) << " |" << std::setw(9)
			<< Cost2.get(n, 0) << "|" << std::setw(9) << Pn.get(n,0) << "|" << std::setw(11)
			<< Pmin.get(n, 0) << "|" << std::setw(11) << Pmax.get(n, 0) << "|" << std::setw(10) << Pn.get(n + _nAgent, 0)
			<< "|" << std::setw(11) << Pmin.get(n + _nAgent, 0) << "|" << std::setw(11) << Pmax.get(n + _nAgent, 0) << "|" << std::endl;
	}


	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "                      END PRINT                                                                         |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;

}

float OPFADMM::DFSP(int j)
{
	//std::cout << "DFSP " << j << std::endl;
	float p = X[j].get(4, 0);
	for (int i = 0; i < nChild.get(j, 0); i++) {
		int c = Childs[j].get(i, 0);
		p += DFSP(c);
	}
	X[j].set(0, 0, p);
	return p;
}
float OPFADMM::DFSQ(int j)
{
	float q = X[j].get(5, 0);
	for (int i = 0; i < nChild.get(j, 0); i++) {
		int c = Childs[j].get(i, 0);
		q += DFSQ(c);
	}
	X[j].set(1, 0, q);
	return q;
}