#include "../head/OPFADMMCons.h"
#define MAX(X, Y) X * (X >= Y) + Y * (Y > X)


OPFADMMCons::OPFADMMCons() : MethodOPF()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " OPFADMMCons Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	timePerBlock = MatrixCPU(1, 12, 0); // Fb0, Fb11abcd, FB12, Fb2, Fb3, Fb4, Fb5,FB6, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 12, 0); //nb de fois utilisé pendant la simu
}


OPFADMMCons::OPFADMMCons(float rho) : MethodOPF()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default OPFADMMCons Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
	timePerBlock = MatrixCPU(1, 12, 0); // Fb0, Fb11, FB12, Fb2, Fb3, Fb4, Fb5, FB6, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 12, 0); //nb de fois utilisé pendant la simu
}

OPFADMMCons::~OPFADMMCons()
{
	DELETEA(tempM1);
	DELETEA(tempM);

	DELETEA(X);
	DELETEA(Ypre);
	DELETEA(Y);
	DELETEA(Mu);

	DELETEA(Hinv);
	DELETEA(A);
	DELETEA(Q);

	DELETEA(Childs);
	DELETEA(Chat);
}
void OPFADMMCons::setParam(float rho)
{
	_rho = rho;
}

bool OPFADMMCons::chekcase()
{
	if (_nBus != (_nLine + 1)) {
		std::cout << "wrong number of line " << _nLine << "against " << _nBus << std::endl;
		return false;
	}
	//CoresLineBus.display();
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

void OPFADMMCons::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
		timePerBlock.increment(0, 0, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
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
	_iterGlobal = 0;
	

	//Chat.display();
	//Bpt2.display();

	while ((_iterGlobal < _iterG) && (resG>epsG)) {
		/*std::cout << "---------------------------------" << std::endl;
		for (int i = 0; i < 3; i++) {
			std::cout << " X " << i << std::endl;
			X[i].display();
			std::cout << " Chat " << i << std::endl;
			Chat[i].display();
			std::cout << " Q " << i << std::endl;
			Q[i].display();
			std::cout << " Y " << i << std::endl;
			Y[i].display();
			std::cout << " Mu " << i << std::endl;
			Mu[i].display();
			
		}*/
		
		
		updateXWOCurrent();
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 5, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		
		CommunicationX();
		
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 6, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		updateGlobalProb();
		updateMu();
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 7, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		updateChat();
		
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		// FB 4
		if (!(_iterGlobal % _stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateRes(_iterGlobal / _stepG);
			//std::cout << _iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, _iterGlobal / _stepG) << " " << resF.get(1, _iterGlobal / _stepG) << std::endl;
			//resG = 1;
#ifdef INSTRUMENTATION
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 9, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
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
	for (int b = 0; b < _nBus; b++) {
		// pi & qi
		int Nb = _nAgentByBus.get(b, 0);
		int begin = _CoresAgentBusBegin.get(b, 0);
		for (int In = 0; In < Nb; In++) {
			int n = _CoresAgentBus.get(In + begin, 0);
			// pn & qn
			Pn.set(n, 0, X[b].get(5 + 2 * In, 0));
			Pn.set(n + _nAgent, 0, X[b].get(6 + 2 * In, 0));
		}
	}


	
	fc = calcFc(&Cost1, &Cost2, &Pn, &tempN2);
	// FB 5
	
	result->setResF(&resF);
	std::cout << "--------" << std::endl;
	std::cout << " Pn " << std::endl;
	Pn.display();
	for (int i = 0; i < 3; i++) {
		std::cout << " X " << i << std::endl;
		X[i].display();
	}
	/**/
	
	

	result->setIter(_iterGlobal);
	

	result->setPn(&Pn);
	
	result->setFc(fc);

#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 10, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 10, 1);

	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
	tall = clock() - tall;
	timeOPF = tall;

	result->setTime((float)tall / CLOCKS_PER_SEC);
	
}

void OPFADMMCons::updateP0(const StudyCase& cas)
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
	
	for (int b = 0; b < _nBus; b++) {
		X[b].set(4, 0, 0);
		X[b].set(5, 0, 0);
		int Nb = _nAgentByBus.get(b, 0);
		int begin = _CoresAgentBusBegin.get(b, 0);
		for (int In = 0; In < Nb; In++) {
			int n = _CoresAgentBus.get(In + begin, 0);
			X[b].increment(4, 0, Pn.get(n, 0));
			X[b].increment(5, 0, Pn.get(n + _nAgent, 0));
	
			// pn & qn
			X[b].set(7 + 2 * In, 0, Pn.get(n, 0));
			X[b].set(8 + 2 * In, 0, Pn.get(n + _nAgent, 0));
	}
	}
	
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
		int Nb = _nAgentByBus.get(i, 0);
		for (int j = 0; j < m; j++) {// (Pci, Qci, lci) for all child Ci
			int c = Childs[i].get(j, 0);
			X[i].set(7 + 2 * Nb + 3 * j, 0, X[c].get(0, 0));
			X[i].set(8 + 2 * Nb + 3 * j, 0, X[c].get(1, 0));
			X[i].set(9 + 2 * Nb + 3 * j, 0, X[c].get(3, 0));
		}
		Y[i].set(&X[i]);
		//Mu[i].set(0.0);
	}/**/
	updateChat();
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 11, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 11, 1);
#endif // INSTRUMENTATION

}

void OPFADMMCons::init(const Simparam& sim, const StudyCase& cas)
{
	// intitilisation des matrixs et variables 
	
	clock_t t = clock();
	std::cout << "init " << std::endl;
	_rho = sim.getRho();
	
	
	const int iterG = sim.getIterG();
	const int stepG = sim.getStepG();
	
	_nAgent = cas.getNagent();
	
	_nBus = cas.getNBus();
	_nBusWLoss = _nBus + 1;
	_nLine = cas.getNLine(true); // ne doit pas être réduit ici !!!

	std::cout << _nAgent << " " << _nBus << " " << _nLine << std::endl;
	_CoresAgentBus = cas.getCoresAgentBusLin();
	_CoresAgentBusBegin = cas.getCoresAgentBusLinBegin();
	_nAgentByBus = cas.getNagentByBus();
	nChild = MatrixCPU(_nBus, 1);
	CoresLineBus = cas.getCoresLineBus(true);
	_CoresBusAgent = cas.getCoresBusAgentLin(); // Cores[n] = b
	Ancestor = MatrixCPU(_nBus, 1, 0); // A_i = bus antécédent de i
	PosChild = MatrixCPU(_nBus, 1, 0); // indice du bus i dans Child[Ai]
	PosAgent = MatrixCPU(_nAgent, 1, -1); // indice de l'agent i dans _CoresAgentBus[CoresAgentBegin]
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
		int busFrom = CoresLineBus.get(lold, 0);
		Ancestor.set(busTo, 0, busFrom);
		nChild.increment(busFrom, 0, 1);
		ZsNorm.set(lold, 0, ZsRe.get(lold, 0) * ZsRe.get(lold, 0) + ZsIm.get(lold, 0) * ZsIm.get(lold, 0));
	}
	


	_rhoInv = 1 / _rho;
	resF = MatrixCPU(3, (iterG / stepG) + 1);

	
	MatrixCPU lowerBound(cas.getLowerBound()); //voltage angle, voltage, line...
	MatrixCPU upperBound(cas.getUpperBound()); //voltage angle, voltage, line...
	

	//std::cout << " local resolution " << std::endl;
	// local resolution
	tempN2 = MatrixCPU(2 * _nAgent, 1);
	tempB2 = MatrixCPU(2 * _nBus, 1);
	CoresSoloBusAgent = MatrixCPU(_nBus, 1, -1);
	Pn = sim.getPn();
	Pmin = cas.getPmin();
	Pmax = cas.getPmax();
	Pbmax = MatrixCPU(2 * _nBus, 1);
	Pbmin = MatrixCPU(2 * _nBus, 1);
	Pb = MatrixCPU(2 * _nBus, 1);
	//Pmin.display();
	//Pn.display();


	Cost1 = MatrixCPU(cas.geta());
	Cost2 = MatrixCPU(cas.getb());

	
	if (Pn.max2() == 0) {
		Pn.add(&Pmin, &Pmax);
		Pn.divide(2);
	}
	
	_nAgentByBus.increment(0, 0, -1); // don't want the grid agent in this resolution
	_CoresAgentBusBegin.increment(0, 0, 1); // idem
	//_nAgentByBus.display();
	
	DELETEA(tempM);
	DELETEA(X);
	DELETEA(Ypre);
	DELETEA(Y);
	DELETEA(Mu);
	DELETEA(Hinv);
	DELETEA(A);
	DELETEA(Q);
	DELETEA(Childs);
	DELETEA(Chat);

	//std::cout << " creation " << std::endl;
	X = new MatrixCPU[_nBusWLoss];
	Ypre = new MatrixCPU[_nBusWLoss];
	Y = new MatrixCPU[_nBusWLoss];
	
	Mu = new MatrixCPU[_nBusWLoss];
	
	tempN1 = MatrixCPU(_nAgent, 1);
	tempNN = MatrixCPU(_nAgent, _nAgent);
	//tempM1 = new MatrixCPU[_nAgent];
	tempM = new MatrixCPU[_nBusWLoss];
	

	Hinv = new MatrixCPU[_nBusWLoss];
	A = new MatrixCPU[_nBusWLoss];
	Q = new MatrixCPU[_nBusWLoss];
	
	Childs = new MatrixCPU[_nBus];

	Chat = new MatrixCPU[_nBusWLoss];
	
	
	VoltageLimit = MatrixCPU(_nBus, 2); // min, max
	VoltageLimitReal = MatrixCPU(_nBus, 2); // min, max
	if (cas.isCurrentLimit()) {
		FluxLimit = cas.getCurrentLimit();
	}
	else {
		FluxLimit = MatrixCPU(_nBus, 1, 1000); // max
	}/**/
	sizeOPFADMMCons = MatrixCPU(_nBusWLoss, 1);
	
	
	int indice = 0;
	MatrixCPU nChildTemp(_nBus, 1, 0);
	for (int i = 0; i < _nBus; i++) {
		int Nb = _nAgentByBus.get(i, 0);
		VoltageLimitReal.set(i, 0, lowerBound.get(_nBus + i, 0));
		VoltageLimitReal.set(i, 1, upperBound.get(_nBus + i, 0));
		VoltageLimit.set(i, 0, lowerBound.get(_nBus + i, 0) * lowerBound.get(_nBus + i, 0) * sqrt((nChild.get(i, 0) + 1) / 2));
		VoltageLimit.set(i, 1, upperBound.get(_nBus + i, 0) * upperBound.get(_nBus + i, 0) * sqrt((nChild.get(i, 0) + 1) / 2));
		
		Childs[i] = MatrixCPU(nChild.get(i, 0), 1);
		
		int sizeOPF = 3 * nChild.get(i, 0) + 5 + 2 * Nb;
		sizeOPFADMMCons.set(i, 0, sizeOPF);

		X[i] = MatrixCPU(sizeOPF, 1); // {Pi, Qi, vi, li, vAi, (pn, qn), (Pci, Qci, lci) for all child Ci}
		Ypre[i] = MatrixCPU(sizeOPF, 1);// Y[i][j] noté dans l'article Yji est ce que i connait sur j
		Y[i] = MatrixCPU(sizeOPF, 1); //Y[i] = {Pi, Qi, vi, li, vAi, (pn, qn),  (Pci, Qci, lci) for all child Ci}
		Mu[i] = MatrixCPU(sizeOPF, 1);
		A[i] = MatrixCPU(2 + 1*(i>0), sizeOPF);
		Hinv[i] = MatrixCPU(sizeOPF, sizeOPF);
		Q[i] = MatrixCPU(sizeOPF, 1, 0);
		tempM[i] = MatrixCPU(sizeOPF, 1);

		Chat[i] = MatrixCPU(4 + 2 * Nb, 1);

	}
	for (int i = 0; i < _nBus; i++) {
		if (i > 0) {
			int Ai = Ancestor.get(i, 0);
			Childs[Ai].set(nChildTemp.get(Ai, 0), 0, i);
			PosChild.set(i, 0, nChildTemp.get(Ai, 0));
			nChildTemp.increment(Ai, 0, 1);
		}
		
	}
	// bus fictif
	int sizeOPF = 0;
	switch (losstype)
	{
	case LossType::POWER:
		sizeOPF = 2 * _nAgent;
		break;
	case LossType::CURRENT:
		sizeOPF = 2 + _nLine;
		break;
	}


	X[_nBus] = MatrixCPU(sizeOPF, 1, 0); // {Ploss, Pn for all agent}
	Ypre[_nBus] = MatrixCPU(sizeOPF, 1, 0);// Y[i][j] noté dans l'article Yji est ce que i connait sur j
	Y[_nBus] = MatrixCPU(sizeOPF, 1, 0); // {Ploss, Pn for all agent}
	Mu[_nBus] = MatrixCPU(sizeOPF, 1, 0);
	A[_nBus] = MatrixCPU(2, sizeOPF, 0);
	Hinv[_nBus] = MatrixCPU(sizeOPF, sizeOPF, 1.0/_nAgent);
	Q[_nBus] = MatrixCPU(sizeOPF, 1);
	tempM[_nBus] = MatrixCPU(sizeOPF, 1, 0);
	Chat[_nBus] = MatrixCPU(2, 1); // que ploss et qloss à gerer 
	sizeOPFADMMCons.set(_nBus, 0, sizeOPF);


	//std::cout << " init valeur " << std::endl;
	for (int i = 0; i < _nBus; i++) {
		// vi = 1; vai = 1
		X[i].set(2, 0, 1); 
		X[i].set(4, 0, 1);

		// pi & qi
		int Nb = _nAgentByBus.get(i, 0);
		int begin = _CoresAgentBusBegin.get(i, 0);
		for (int In = 0; In < Nb; In++) {
			int n = _CoresAgentBus.get(In + begin, 0);
			PosAgent.set(n, 0, In);
			Pb.increment(i, 0, Pn.get(n, 0));
			Pb.increment(i + _nBus, 0, Pn.get(n + _nAgent, 0));
			//X[i].increment(4, 0, Pn.get(n, 0));
			//X[i].increment(5, 0, Pn.get(n + _nAgent, 0));
			
			Pbmax.increment(i, 0, Pmax.get(n, 0));
			Pbmin.increment(i, 0, Pmin.get(n, 0));

			Pbmax.increment(i + _nBus, 0, Pmax.get(n + _nAgent, 0));
			Pbmin.increment(i + _nBus, 0, Pmin.get(n + _nAgent, 0));

			if (Nb == 1) {
				CoresSoloBusAgent.set(i, 0, n);
			}
			// pn & qn
			X[i].set(5 + 2 * In, 0, Pn.get(n, 0));
			X[i].set(6 + 2 * In, 0, Pn.get(n + _nAgent, 0));
		}
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
		int Ni = _nAgentByBus.get(i, 0);
		for (int j = 0; j < m; j++) {
			// (Pci, Qci, lci) for all child Ci
			
			int c = Childs[i].get(j, 0);
			X[i].set(5 + 2 * Ni + 3 * j, 0, X[c].get(0, 0));
			X[i].set(6 + 2 * Ni + 3 * j, 0, X[c].get(1, 0));
			X[i].set(7 + 2 * Ni + 3 * j, 0, X[c].get(3, 0));
		}
		Y[i].set(&X[i]);
	}

	// bus factice
	float sumP = 0;
	float sumQ = 0;
	switch (losstype)
	{
	case LossType::POWER:
		for (int n = 1; n < _nAgent; n++) {
			int bus = _CoresBusAgent.get(n, 0);
			int In = PosAgent.get(n, 0);

			X[_nBus].set(n, 0, X[bus].get(5 + 2 * In, 0));
			X[_nBus].set(n + _nAgent, 0, X[bus].get(6 + 2 * In, 0));

			sumP += X[bus].get(5 + 2 * In, 0);
			sumQ += X[bus].get(6 + 2 * In, 0);
		}
		X[_nBus].set(0, 0, -sumP);
		X[_nBus].set(_nAgent, 0, -sumQ);

		for (int i = 0; i < _nAgent; i++) {
			A[_nBus].set(0, i, 1); // sum(p) + Ploss = 0
			A[_nBus].set(1, i + _nAgent, 1); // Qloss + sum(q) = 0
		}
		break;
	case LossType::CURRENT:
		for (int i = 1; i < _nBus; i++) {
			X[_nBus].set(i + 1, 0, X[i].get(3, 0));
			sumP += X[i].get(3, 0) * ZsRe.get(i - 1, 0);
			sumQ += X[i].get(3, 0) * ZsIm.get(i - 1, 0);
		}
		for (int i = 0; i < _nLine; i++) {
			A[_nBus].set(0, i + 2, ZsRe.get(i, 0)); // sum(l R)= Ploss
			A[_nBus].set(1, i + 2, ZsIm.get(i, 0)); // Qloss=sum(l L)
		}
		X[_nBus].set(0, 0, -sumP);
		X[_nBus].set(1, 0, -sumQ);

		A[_nBus].set(0, 0, 1); // ploss
		A[_nBus].set(1, 1, 1); // qloss
		break;
	}


	Y[_nBus].set(&X[_nBus]);

	//std::cout << " Hinv " << std::endl;
	for (int i = 0; i < _nBus; i++) {
		int Ni = _nAgentByBus.get(i, 0);
		if (i > 0) {
			A[i].set(2, 0, 2 * ZsRe.get(i - 1, 0));
			A[i].set(2, 1, 2 * ZsIm.get(i - 1, 0));
			A[i].set(2, 3, -ZsNorm.get(i - 1, 0));
			A[i].set(2, 2, -1);
			A[i].set(2, 4, 1);
			A[i].set(0, 0, -1);
			A[i].set(1, 1, -1);
		}
		
		// pi = sum(pn) & qi = sum(qn)
	
		for (int In = 0; In < Ni; In++) {
			A[i].set(0, 5 + 2 * In, 1);
			A[i].set(1, 6 + 2 * In, 1);
		}
		
		
		for (int j = 0; j < nChild.get(i, 0); j++) {
			int c = Childs[i].get(j, 0);
			A[i].set(0, 5 + 2 * Ni + 3 * j, 1); // Pci
			A[i].set(1, 6 + 2 * Ni + 3 * j, 1); // Qci
			A[i].set(0, 7 + 2 * Ni + 3 * j, -ZsRe.get(c - 1, 0)); // -R l
			A[i].set(1, 7 + 2 * Ni + 3 * j, -ZsIm.get(c - 1, 0)); // -X l
		}
		//A[i].display();
		MatrixCPU temp33(2 + 1 * (i > 0), 2 + 1 * (i > 0));
		MatrixCPU temp3M(2 + 1 * (i > 0), sizeOPFADMMCons.get(i, 0));
		MatrixCPU tempMM(sizeOPFADMMCons.get(i, 0), sizeOPFADMMCons.get(i, 0));
		temp33.multiplyTrans(&A[i], &A[i]); // (3*o_b) * (o_b*3) -> 9 * o_b^2
		temp33.invertEigen(&temp33); // 3^3 = 27 (fixe !!!)
		temp3M.MultiplyMatMat(&temp33, &A[i]); // (3*3) * (3*o_b) -> 27 *o_b
		Hinv[i].multiplyTrans(&A[i], &temp3M, 0); // (o_b*3) * (3*o_b) -> 9 * o_b

		tempMM.setEyes(1);
		Hinv[i].subtract(&tempMM);
		Hinv[i].divide(_rho);
		//Hinv[i].display();
	}
	
	//A[_nBus].display();
	MatrixCPU temp22(2, 2);
	MatrixCPU temp2M(2, sizeOPFADMMCons.get(_nBus, 0));
	MatrixCPU tempMM(sizeOPFADMMCons.get(_nBus, 0), sizeOPFADMMCons.get(_nBus, 0));
	temp22.multiplyTrans(&A[_nBus], &A[_nBus]);
	temp22.invertEigen(&temp22); // 3^3 = 27 (fixe !!!)
	temp2M.MultiplyMatMat(&temp22, &A[_nBus]); // (3*3) * (3*o_b) -> 27 *o_b
	Hinv[_nBus].multiplyTrans(&A[_nBus], &temp2M, 0); // (o_b*3) * (3*o_b) -> 9 * o_b

	tempMM.setEyes(1);
	Hinv[_nBus].subtract(&tempMM);
	Hinv[_nBus].divide(_rho);
	//Hinv[_nBus].display();
	
	//std::cout << "updateChat" << std::endl;
	updateChat();
	
	
	/*std::cout << "--------" << std::endl;
	std::cout << " Pn " << std::endl;
	Pn.display();*/
	/*for (int i = 0; i < _nBus; i++) {
		std::cout << " X " << i << std::endl;
		X[i].display();
		std::cout << " Y " << i << std::endl;
		Y[i].display();
		std::cout << " Q " << i << std::endl;
		Q[i].display();
		std::cout << " Mu " << i << std::endl;
		Mu[i].display();
		std::cout << " Chat " << i << std::endl;
		Chat[i].display();
	}*/
	

	
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	std::cout << "---------------------------------------------------------------------------------------" << std::endl;
}

void OPFADMMCons::solveConsensus(float eps, MatrixCPU* PSO)
{
	float epsG = eps;
	
	float fc = 0;
	float resG = 2 * epsG;
	_iterGlobal = 0;
	timeOPF = clock();

	/*std::cout << "****** Solve Consensus OPF part ************" << std::endl;
	Chat[1].display();
	Cost2.display();*/

	while ((_iterGlobal < _iterG) && (resG > epsG)) {

		
		updateXWOCurrent();
		CommunicationX();

		
		updateGlobalProb();
		updateMu();
		
		updateChat();

		// FB 4
		if (!(_iterGlobal % _stepG)) {
			resG = updateResRhoFixe(_iterGlobal / _stepG);
		}
		//std::cout << iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;

		_iterGlobal++;
	}
	if (_iterG == _iterGlobal) {
		//std::cout << "OPF 2 : " << _iterGlobal << " " << resF.get(0, (_iterGlobal - 1) / _stepG) << " " << resF.get(1, (_iterGlobal - 1) / _stepG) << " " << resF.get(2, (_iterGlobal - 1) / _stepG) << std::endl;
	}
	//X[1].display();
	
	
	for (int b = 0; b < _nBus; b++) {
		// pi & qi
		int Nb = _nAgentByBus.get(b, 0);
		int begin = _CoresAgentBusBegin.get(b, 0);
		for (int In = 0; In < Nb; In++) {
			int n = _CoresAgentBus.get(In + begin, 0);
			// pn & qn
			Pn.set(n, 0, X[b].get(5 + 2 * In, 0));
			Pn.set(n + _nAgent, 0, X[b].get(6 + 2 * In, 0));
		}
	}
	
	//X[_nBus].display();
	//Y[_nBus].display();
	//Pn.display();
	switch (losstype)
	{
	case LossType::POWER:
		_Ploss = X[_nBus].get(0, 0);
		_Qloss = X[_nBus].get(_nAgent, 0);
		break;
	case LossType::CURRENT:
		_Ploss = X[_nBus].get(0, 0);
		_Qloss = X[_nBus].get(1, 0);
		break;
	}

	//Pn.set(0, 0, _Ploss);
	//Pn.set(_nAgent, 0, _Qloss);
	
	PSO->set(&Pn);
	PSO->set(0, 0, _Ploss);
	PSO->set(_nAgent, 0, _Qloss);

	//PSO->display();

	timeOPF = clock() - timeOPF;

}


void OPFADMMCons::initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO)
{
	// intitilisation des matrixs et variables 

	clock_t t = clock();
	//std::cout << "init " << std::endl;
	_rho = sim.getRho();
	_rhoSO = rhoSO;

	_iterG = sim.getIterG();
	_stepG = sim.getStepG();

	_nAgent = cas.getNagent();

	_nBus = cas.getNBus();
	_nBusWLoss = _nBus + 1;
	_nLine = cas.getNLine(true); // ne doit pas être réduit ici !!!

	//std::cout << _nAgent << " " << _nBus << " " << _nLine << std::endl;
	_CoresAgentBus = cas.getCoresAgentBusLin();
	_CoresAgentBusBegin = cas.getCoresAgentBusLinBegin();
	_nAgentByBus = cas.getNagentByBus();
	nChild = MatrixCPU(_nBus, 1);
	CoresLineBus = cas.getCoresLineBus(true);
	_CoresBusAgent = cas.getCoresBusAgentLin(); // Cores[n] = b
	Ancestor = MatrixCPU(_nBus, 1, 0); // A_i = bus antécédent de i
	PosChild = MatrixCPU(_nBus, 1, 0); // indice du bus i dans Child[Ai]
	PosAgent = MatrixCPU(_nAgent, 1, -1); // indice de l'agent i dans _CoresAgentBus[CoresAgentBegin]
	Ancestor.set(0, 0, -1); // the slack bus has no ancestor
	ZsRe = cas.getZsRe();
	ZsIm = cas.getZsImag();	
	ZsNorm = MatrixCPU(_nLine, 1);

	if (!chekcase()) {
		throw std::invalid_argument("not a radial case");
	}
	if (losstype != LossType::CURRENT && losstype != LossType::POWER ) {
		throw std::invalid_argument("unkown lossType");
	}
	

	
	for (int lold = 0; lold < _nLine; lold++) {
		int l = lold + 1;
		int busTo = l;
		int busFrom = CoresLineBus.get(lold, 0);
		Ancestor.set(busTo, 0, busFrom);
		nChild.increment(busFrom, 0, 1);
		ZsNorm.set(lold, 0, ZsRe.get(lold, 0) * ZsRe.get(lold, 0) + ZsIm.get(lold, 0) * ZsIm.get(lold, 0));
	}
	
	/*
	ZsNorm.display();*/


	_rhoInv = 1 / _rho;
	resF = MatrixCPU(3, (_iterG / _stepG) + 1);


	MatrixCPU lowerBound(cas.getLowerBound()); //voltage angle, voltage, line...
	MatrixCPU upperBound(cas.getUpperBound()); //voltage angle, voltage, line...


	//std::cout << " local resolution " << std::endl;
	// local resolution
	tempN2 = MatrixCPU(2 * _nAgent, 1);
	tempB2 = MatrixCPU(2 * _nBus, 1);
	CoresSoloBusAgent = MatrixCPU(_nBus, 1, -1);
	Pn = sim.getPn(); // not the real agent
	Pmin = MatrixCPU(2 * _nAgent, 1, -1000000); // must not be the real one
	Pmax = MatrixCPU(2 * _nAgent, 1, 1000000); // idem

	// the loss provider
	/*Pmin.set(0, 0, 0);
	Pmax.set(0, 0, 0);
	Pmin.set(_nAgent, 0, 0);
	Pmax.set(_nAgent, 0, 0);*/

	Pbmax = MatrixCPU(2 * _nBus, 1);
	Pbmin = MatrixCPU(2 * _nBus, 1);
	Pb = MatrixCPU(2 * _nBus, 1);
	//Pmin.display();
	//Pn.display();


	Cost1 = MatrixCPU(2 * _nAgent, 1, _rhoSO); //
	Cost1.set(0, 0, 0);
	Cost1.set(_nAgent, 0, 0);
	Cost2 = MatrixCPU(2 * _nAgent, 1);
	etaSO = MatrixCPU(2 * _nAgent, 1);

	
	_nAgentByBus.increment(0, 0, -1); // don't want the grid agent in this resolution
	_CoresAgentBusBegin.increment(0, 0, 1); // idem
	//_nAgentByBus.display();

	DELETEA(tempM);
	DELETEA(X);
	DELETEA(Ypre);
	DELETEA(Y);
	DELETEA(Mu);
	DELETEA(Hinv);
	DELETEA(A);
	DELETEA(Q);
	DELETEA(Childs);
	DELETEA(Chat);


	//std::cout << " creation " << std::endl;
	X = new MatrixCPU[_nBusWLoss];
	Ypre = new MatrixCPU[_nBusWLoss];
	Y = new MatrixCPU[_nBusWLoss];

	Mu = new MatrixCPU[_nBusWLoss];

	tempN1 = MatrixCPU(_nAgent, 1);
	tempNN = MatrixCPU(_nAgent, _nAgent);
	//tempM1 = new MatrixCPU[_nAgent];
	tempM = new MatrixCPU[_nBusWLoss];


	Hinv = new MatrixCPU[_nBusWLoss];
	A = new MatrixCPU[_nBusWLoss];
	Q = new MatrixCPU[_nBusWLoss];

	Childs = new MatrixCPU[_nBus];

	Chat = new MatrixCPU[_nBusWLoss];
	VoltageLimit = MatrixCPU(_nBus, 2); // min, max
	VoltageLimitReal = MatrixCPU(_nBus, 2); // min, max
	if (cas.isCurrentLimit()) {
		FluxLimit = cas.getCurrentLimit();
		isCurrentLimited = true;
	}
	else {
		FluxLimit = MatrixCPU(_nLine, 1, 1000); // max
	}
	//FluxLimit.display();
	sizeOPFADMMCons = MatrixCPU(_nBusWLoss, 1);


	int indice = 0;
	MatrixCPU nChildTemp(_nBus, 1, 0);
	for (int i = 0; i < _nBus; i++) {
		int Nb = _nAgentByBus.get(i, 0);
		VoltageLimitReal.set(i, 0, lowerBound.get(_nBus + i, 0));
		VoltageLimitReal.set(i, 1, upperBound.get(_nBus + i, 0));
		VoltageLimit.set(i, 0, lowerBound.get(_nBus + i, 0) * lowerBound.get(_nBus + i, 0) * sqrt((nChild.get(i, 0) + 1) / 2));
		VoltageLimit.set(i, 1, upperBound.get(_nBus + i, 0) * upperBound.get(_nBus + i, 0) * sqrt((nChild.get(i, 0) + 1) / 2));
		
		Childs[i] = MatrixCPU(nChild.get(i, 0), 1);

		int sizeOPF = 3 * nChild.get(i, 0) + 5 + 2 * Nb;
		sizeOPFADMMCons.set(i, 0, sizeOPF);

		X[i] = MatrixCPU(sizeOPF, 1); // {Pi, Qi, vi, li, vAi, (pn, qn), (Pci, Qci, lci) for all child Ci}
		Ypre[i] = MatrixCPU(sizeOPF, 1);// Y[i][j] noté dans l'article Yji est ce que i connait sur j
		Y[i] = MatrixCPU(sizeOPF, 1); //Y[i] = {Pi, Qi, vi, li, vAi, (pn, qn),  (Pci, Qci, lci) for all child Ci}
		Mu[i] = MatrixCPU(sizeOPF, 1);
		A[i] = MatrixCPU(2 + 1 * (i > 0), sizeOPF);
		Hinv[i] = MatrixCPU(sizeOPF, sizeOPF);
		Q[i] = MatrixCPU(sizeOPF, 1, 0);
		tempM[i] = MatrixCPU(sizeOPF, 1);

		Chat[i] = MatrixCPU(4 + 2 * Nb, 1);

	}
	for (int i = 0; i < _nBus; i++) {
		if (i > 0) {
			int Ai = Ancestor.get(i, 0);
			Childs[Ai].set(nChildTemp.get(Ai, 0), 0, i);
			PosChild.set(i, 0, nChildTemp.get(Ai, 0));
			nChildTemp.increment(Ai, 0, 1);
		}

	}
	

	// bus fictif
	int sizeOPF = 2;
	switch (losstype)
	{
	case LossType::POWER:
		sizeOPF = 2 * _nAgent;
		break;
	case LossType::CURRENT:
		sizeOPF = 2 + _nLine;
		break;
	}/**/
	
	
	X[_nBus] = MatrixCPU(sizeOPF, 1, 0); // {Ploss, Qloss, lj for all line}
	Ypre[_nBus] = MatrixCPU(sizeOPF, 1, 0);// Y[i][j] noté dans l'article Yji est ce que i connait sur j
	Y[_nBus] = MatrixCPU(sizeOPF, 1, 0); // {Ploss, Qloss, lj for all line}
	Mu[_nBus] = MatrixCPU(sizeOPF, 1, 0);
	//A[_nBus] = MatrixCPU(2, sizeOPF, 0);
	Hinv[_nBus] = MatrixCPU(sizeOPF, sizeOPF);
	Q[_nBus] = MatrixCPU(sizeOPF, 1);
	tempM[_nBus] = MatrixCPU(sizeOPF, 1, 0);
	Chat[_nBus] = MatrixCPU(2, 1); // que ploss et qloss à gerer 
	sizeOPFADMMCons.set(_nBus, 0, sizeOPF);


	//std::cout << " init valeur " << std::endl;
	for (int i = 0; i < _nBus; i++) {
		// vi = 1; vai = 1
		X[i].set(2, 0, 1);
		X[i].set(4, 0, 1);

		// pi & qi
		int Nb = _nAgentByBus.get(i, 0);
		int begin = _CoresAgentBusBegin.get(i, 0);
		for (int In = 0; In < Nb; In++) {
			int n = _CoresAgentBus.get(In + begin, 0);
			PosAgent.set(n, 0, In);
			Pb.increment(i, 0, Pn.get(n, 0));
			Pb.increment(i + _nBus, 0, Pn.get(n + _nAgent, 0));
			//X[i].increment(4, 0, Pn.get(n, 0));
			//X[i].increment(5, 0, Pn.get(n + _nAgent, 0));

			Pbmax.increment(i, 0, Pmax.get(n, 0));
			Pbmin.increment(i, 0, Pmin.get(n, 0));

			Pbmax.increment(i + _nBus, 0, Pmax.get(n + _nAgent, 0));
			Pbmin.increment(i + _nBus, 0, Pmin.get(n + _nAgent, 0));

			if (Nb == 1) {
				CoresSoloBusAgent.set(i, 0, n);
			}
			// pn & qn
			X[i].set(5 + 2 * In, 0, Pn.get(n, 0));
			X[i].set(6 + 2 * In, 0, Pn.get(n + _nAgent, 0));
		}
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
		int Ni = _nAgentByBus.get(i, 0);
		for (int j = 0; j < m; j++) {
			// (Pci, Qci, lci) for all child Ci

			int c = Childs[i].get(j, 0);
			X[i].set(5 + 2 * Ni + 3 * j, 0, X[c].get(0, 0));
			X[i].set(6 + 2 * Ni + 3 * j, 0, X[c].get(1, 0));
			X[i].set(7 + 2 * Ni + 3 * j, 0, X[c].get(3, 0));
		}
		Y[i].set(&X[i]);
	}
	// bus factice
	float sumP = 0;
	float sumQ = 0;
	switch (losstype)
	{
	case LossType::POWER:
		for (int n = 1; n < _nAgent; n++) {
			int bus = _CoresBusAgent.get(n, 0);
			int In = PosAgent.get(n, 0);

			X[_nBus].set(n, 0, X[bus].get(5 + 2 * In, 0));
			X[_nBus].set(n + _nAgent, 0, X[bus].get(6 + 2 * In, 0));

			sumP += X[bus].get(5 + 2 * In, 0);
			sumQ += X[bus].get(6 + 2 * In, 0);
		}
		X[_nBus].set(0		, 0, -sumP);
		X[_nBus].set(_nAgent, 0, -sumQ);
		
		/*for (int i = 0; i < _nAgent; i++) {
			A[_nBus].set(0, i, 1); // sum(p) + Ploss = 0
			A[_nBus].set(1, i + _nAgent, 1); // Qloss + sum(q) = 0
		}*/
		Hinv[_nBus].setEyes(-1);
		Hinv[_nBus].set(0, 0, 0);
		Hinv[_nBus].set(_nAgent, _nAgent, 0);
		for (int i = 1; i < _nAgent; i++) {
			Hinv[_nBus].set(0, i, 1); // sum(p) + Ploss = 0
			Hinv[_nBus].set(_nAgent, i + _nAgent, 1); // Qloss + sum(q) = 0
		}
		break;
	case LossType::CURRENT:
		for (int i = 1; i < _nBus; i++) {
			//X[_nBus].set(i + 1, 0, X[i].get(3, 0));
			sumP += X[i].get(3, 0) * ZsRe.get(i - 1, 0);
			sumQ += X[i].get(3, 0) * ZsIm.get(i - 1, 0);
		}
		X[_nBus].set(0, 0, -sumP);
		X[_nBus].set(1, 0, -sumQ);
		
		/*for (int i = 0; i < _nLine; i++) {
			A[_nBus].set(0, i + 2, ZsRe.get(i,0)); // sum(l R)= Ploss
			A[_nBus].set(1, i + 2, ZsIm.get(i,0)); // Qloss=sum(l L)
		}
		A[_nBus].set(0, 0, 1); // ploss
		A[_nBus].set(1, 1, 1); // qloss
		*/
		Hinv[_nBus].setEyes(-1);
		Hinv[_nBus].set(0, 0, 0);
		Hinv[_nBus].set(1, 1, 0);
		for (int i = 0; i < _nLine; i++) {
			Hinv[_nBus].set(0, i + 2, ZsRe.get(i, 0)); // sum(p) + Ploss = 0
			Hinv[_nBus].set(1, i + 2, ZsIm.get(i, 0)); // Qloss + sum(q) = 0
		}
		
		break;
	}
	
	//X[_nBus].set(0, 0, getPLoss());
	//X[_nBus].set(1, 0, getQLoss());
	Hinv[_nBus].divide(_rho);
	Y[_nBus].set(&X[_nBus]);

	//std::cout << " Hinv " << std::endl;
	for (int i = 0; i < _nBus; i++) {
		int Ni = _nAgentByBus.get(i, 0);
		if (i > 0) {
			A[i].set(2, 0, 2 * ZsRe.get(i - 1, 0));
			A[i].set(2, 1, 2 * ZsIm.get(i - 1, 0));
			A[i].set(2, 3, -ZsNorm.get(i - 1, 0));
			A[i].set(2, 2, -1);
			A[i].set(2, 4, 1);
			A[i].set(0, 0, -1);
			A[i].set(1, 1, -1);
		}

		// pi = sum(pn) & qi = sum(qn)

		for (int In = 0; In < Ni; In++) {
			A[i].set(0, 5 + 2 * In, 1);
			A[i].set(1, 6 + 2 * In, 1);
		}


		for (int j = 0; j < nChild.get(i, 0); j++) {
			int c = Childs[i].get(j, 0);
			A[i].set(0, 5 + 2 * Ni + 3 * j, 1); // Pci
			A[i].set(1, 6 + 2 * Ni + 3 * j, 1); // Qci
			A[i].set(0, 7 + 2 * Ni + 3 * j, -ZsRe.get(c - 1, 0)); // -R l
			A[i].set(1, 7 + 2 * Ni + 3 * j, -ZsIm.get(c - 1, 0)); // -X l
		}
		//A[i].display();
		MatrixCPU temp33(2 + 1 * (i > 0), 2 + 1 * (i > 0));
		MatrixCPU temp3M(2 + 1 * (i > 0), sizeOPFADMMCons.get(i, 0));
		MatrixCPU tempMM(sizeOPFADMMCons.get(i, 0), sizeOPFADMMCons.get(i, 0));
		temp33.multiplyTrans(&A[i], &A[i]); // (3*o_b) * (o_b*3) -> 9 * o_b^2
		temp33.invertEigen(&temp33); // 3^3 = 27 (fixe !!!)
		temp3M.MultiplyMatMat(&temp33, &A[i]); // (3*3) * (3*o_b) -> 27 *o_b
		Hinv[i].multiplyTrans(&A[i], &temp3M, 0); // (o_b*3) * (3*o_b) -> 9 * o_b

		tempMM.setEyes(1);
		Hinv[i].subtract(&tempMM);
		Hinv[i].divide(_rho);
		//Hinv[i].display();
	}

	// bus fictif

	
	
	//A[_nBus].display();
	/*MatrixCPU temp22(2, 2);
	MatrixCPU temp2M(2, sizeOPFADMMCons.get(_nBus, 0));
	MatrixCPU tempMM(sizeOPFADMMCons.get(_nBus, 0), sizeOPFADMMCons.get(_nBus, 0));
	temp22.multiplyTrans(&A[_nBus], &A[_nBus]);
	temp22.invertEigen(&temp22); // 3^3 = 27 (fixe !!!)
	temp2M.MultiplyMatMat(&temp22, &A[_nBus]); // (3*3) * (3*o_b) -> 27 *o_b
	Hinv[_nBus].multiplyTrans(&A[_nBus], &temp2M, 0); // (o_b*3) * (3*o_b) -> 9 * o_b

	tempMM.setEyes(1);
	Hinv[_nBus].subtract(&tempMM);
	Hinv[_nBus].divide(_rho);*/
	//Hinv[_nBus].display();

	//std::cout << "updateChat" << std::endl;
	updateChat();
	

	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	//std::cout << "---------------------------------------------------------------------------------------" << std::endl;
}


void OPFADMMCons::updateConsensus(MatrixCPU* Pmarket)
{
	// z_loss = (P_0 + Ploss)/2
	// lambda_loss += 0.5 (Ploss - P_0)  

	/*std::cout << "Pmarket, Pn" << std::endl;
	Pmarket->display();
	Pn.display();
	std::cout << " Ploss, Qloss " << _Ploss << " " << _Qloss << std::endl;**/

	/*float eta = 0.5 * (_Ploss - Pmarket->get(0, 0));
	etaSO.set(0, 0, etaSO.get(0, 0) + eta);
	eta = 0.5 * (_Qloss - Pmarket->get( _nAgent, 0));
	etaSO.set(_nAgent, 0, etaSO.get(_nAgent, 0) + eta);
	int omega = _nAgent - 1;*/
	//Pmarket->display();
	
	for (int n = 1; n < _nAgent; n++) { 
		float eta = etaSO.get(n, 0) +  0.5 * (Pn.get(n, 0) - Pmarket->get(n, 0));
		etaSO.set(n, 0, eta);
		eta = etaSO.get(n + _nAgent, 0) + 0.5 * (Pn.get(n + _nAgent, 0) - Pmarket->get(n + _nAgent, 0));
		etaSO.set(n + _nAgent, 0,  eta);
	}
	
	Cost2.add(&Pn, Pmarket);
	Cost2.set(0, 0, 0);
	Cost2.set(_nAgent, 0, 0);
	//Cost2.set(0, 0, _Ploss + Pmarket->get(0, 0));
	//Cost2.set(_nAgent, 0, _Qloss + Pmarket->get(_nAgent, 0));
	Cost2.multiply(-0.5);
	Cost2.add(&etaSO);
	Cost2.multiply(_rhoSO);

	/*std::cout << "Cost 2" << std::endl;
	Cost2.display();
	std::cout << "*********" << std::endl;*/

}



void OPFADMMCons::updateGlobalProb() {
	
	
	for (int i = 0; i < _nBusWLoss; i++) {
		//std::cout << " Q " << std::endl;
		//Q[i].display();
		Ypre[i].swap(&Y[i]);
		Y[i].MultiplyMatVec(&Hinv[i], &Q[i]); // solve system by using the inverse
	}
	/*Ypre[_nBus].swap(&Y[_nBus]);
	Y[_nBus].set(0, 0, getPLoss());
	Y[_nBus].set(1, 0, getQLoss());*/

	Y[0].set(2, 0, 1);
	Y[0].set(4, 0, 1);


	// communication of y, mu

}

void OPFADMMCons::updateX()
{
	double x1, x2, x3, x4, c1, c2, c3, c4, lambdaLo, lambdaUp, delta, x3min, x3max, x4max, gamma, k2;
	double c1122;
	int nSol = 0;
	int typeSol = 0;
	int BestRoot = 0;
	double bestGamma = -1;
	double p = 0;
	int nRoot = 0;
	bool neg = false;


	for (int i = 0; i < _nBus; i++) {

		bool goodSol = false;
		k2 = sqrt(2.0 / (nChild.get(i, 0) + 1));
		typeSol = 0;
		if (i == 0) { // slack bus
			goodSol = true;
			c3 = -2 * Chat[i].get(2, 0) / k2;
			x1 = 0;
			x2 = 0;
			x4 = 0;
			x4max = 0;
			x3 = 1 / k2;
			gamma = 0;
			typeSol = 1;
		}
		else {
			c1 = -2 * Chat[i].get(0, 0);
			c2 = -2 * Chat[i].get(1, 0);
			c3 = -2 * Chat[i].get(2, 0) / k2;
			c4 = -2 * Chat[i].get(3, 0);
			c1122 = c1 * c1 + c2 * c2;


			x3min = VoltageLimit.get(i, 0);
			x3max = VoltageLimit.get(i, 1);
			x4max = FluxLimit.get(i - 1, 0);

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

			if (x4 > x4max) {
				x4 = x4max;
			}

			gamma = k2 * x4 - (x1 * x1 + x2 * x2) / x3; // ce n'est pas vraiment gamma, doit être positif
			//std::cout << "x 1 : " << x1 << " " << x2 << " " << x3 * k2 << " " << x4 << " " << (x1 * x1 + x2 * x2) / x3  - k2 * x4 << std::endl;

			if (gamma >= 0) {
				// the solution is good !
				typeSol = 1;
				goodSol = true;
			}
			else {
				if (c1122 == 0) { // cas dégénéré
					std::cout << " bus " << i << " : c1= " << c1 << " c2=" << c2 << " c4=" << c4 << " gamma= " << gamma << std::endl;
					x4 = 0;
					goodSol = true;
				}
				else if (gamma > bestGamma) {
					typeSol = 1;
					bestGamma = gamma;
				}
			}
		}

		// cas x4 = x4 max 

		/*	//x3 = x3max
		if (!goodSol) {
			x4 = x4max;
			x3 = x3max;
			p = sqrt((k2 * x4) / (c1122 * x3)); // plus ou mois ce truc !!!
			x1 = p * c1 * x3;
			x2 = p * c2 * x3;
			if (abs(c1) > 0) {
				gamma = -(2 * x1 + c1) * x3 / (2 * x1);
			}
			else {
				gamma = -(2 * x2 + c2) * x3 / (2 * x2);
			}
			lambdaUp = -(2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));
			delta = k2 * gamma - 2 * x4 - c4;
			if (gamma >= 0 && lambdaUp >= 0 && delta > 0) {
				typeSol = 5;
				goodSol = true;
			}
			else if (gamma > bestGamma && lambdaUp > bestGamma && delta > bestGamma) {
				typeSol = 5;
				bestGamma = min(min(gamma, lambdaUp), delta);
			}
		}
		if (!goodSol) {
			p = -p;
			x1 = p * c1 * x3;
			x2 = p * c2 * x3;
			if (abs(c1) > 0) {
				gamma = -(2 * x1 + c1) * x3 / (2 * x1);
			}
			else {
				gamma = -(2 * x2 + c2) * x3 / (2 * x2);
			}
			lambdaUp = -(2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));
			delta = k2 * gamma - 2 * x4 - c4;
			if (gamma >= 0 && lambdaUp >= 0 && delta >= 0) {
				// the solution is good 
				goodSol = true;
				typeSol = 5;
				//nSol = n;
			}
			if (gamma > bestGamma && lambdaUp > bestGamma && delta > bestGamma) {
				typeSol = 5;
				neg = true;
				bestGamma = min(min(gamma, lambdaUp), delta);
			}
		}
		//x3 = x3min
		if (!goodSol) {
			// cas x3 = xmin
			x3 = x3min;
			p = sqrt((k2 * x4) / (c1122 * x3));
			x1 = p * c1 * x3;
			x2 = p * c2 * x3;
			if (abs(c1) > 0) {
				gamma = -(2 * x1 + c1) * x3 / (2 * x1);
			}
			else {
				gamma = -(2 * x2 + c2) * x3 / (2 * x2);
			}

			lambdaLo = -(2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));
			delta = k2 * gamma - 2 * x4 - c4;
			if (gamma >= 0 && lambdaLo >= 0 && delta >= 0) {
				// the solution is good 
				typeSol = 6;
				goodSol = true;
				//nSol = n;

			}
			else if (gamma > bestGamma && lambdaLo > bestGamma && delta > bestGamma) {
				typeSol = 6;
				bestGamma = min(min(gamma, lambdaLo), delta);
			}
		}
		if (!goodSol) {
			p = -p;
			x1 = p * c1 * x3;
			x2 = p * c2 * x3;
			if (abs(c1) > 0) {
				gamma = -(2 * x1 + c1) * x3 / (2 * x1);
			}
			else {
				gamma = -(2 * x2 + c2) * x3 / (2 * x2);
			}

			lambdaLo = -(2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));
			delta = k2 * gamma - 2 * x4 - c4;
			if (gamma >= 0 && lambdaLo >= 0 && delta >= 0) {
				// the solution is good 
				goodSol = true;
				typeSol = 6;
				//nSol = n;

			}
			else if (gamma > bestGamma && lambdaLo > bestGamma && delta > bestGamma) {
				typeSol = 6;
				bestGamma = min(min(gamma, lambdaLo), delta);
				neg = true;
			}
		}
		// x3min <x3 < x3max
		if (!goodSol) {
			// cas tension libre
			coefPoly2[0] = (c3 + k2 * x4) / 2;
			coefPoly2[1] = sqrt(k2 * x4 * c1122) / 4;
			//std::cout << " polynome " << coefPoly2[0] << " " << coefPoly2[1] << std::endl;

			nRoot = resolveRealPolynome3without2term(root5, coefPoly2);
			for (int n = 0; n < nRoot; n++) {
				double sqrtX3 = root5[n];
				//std::cout << "root5 " << root5[n] << std::endl;
				if (sqrtX3 >= 0) {
					x3 = sqrtX3 * sqrtX3;
					p = sqrt((k2 * x4) / (c1122 * x3));

					x1 = p * c1 * x3;
					x2 = p * c2 * x3;


					if (abs(c1) > 0) {
						gamma = -(2 * x1 + c1) * x3 / (2 * x1);
					}
					else {
						gamma = -(2 * x2 + c2) * x3 / (2 * x2);
					}
					delta = k2 * gamma - 2 * x4 - c4;
					//std::cout << "x : " << x1 << " " << x2 << " " << x3 * k2 << " " << x4 << " " << gamma << " " << delta << std::endl;
					if (gamma >= 0 && delta >= 0 && x3 <= x3max && x3 >= x3min) {
						// the solution is good 
						typeSol = 7;
						goodSol = true;
						//nSol = n;
						break;
					}
					if (gamma > bestGamma && delta > bestGamma && (x3max - x3) > bestGamma && (x3 - x3min) > bestGamma) {
						typeSol = 7;
						bestGamma = min(min(min(gamma, (x3max - x3)), (x3 - x3min)), delta);
						BestRoot = n;
					}
				}

			}
		}
		if (!goodSol) {
			coefPoly2[0] = (c3 + k2 * x4) / 2;
			coefPoly2[1] = -sqrt(k2 * x4 * c1122) / 4;
			//std::cout << " polynome " << coefPoly2[0] << " " << coefPoly2[1] << std::endl;

			nRoot = resolveRealPolynome3without2term(root6, coefPoly2);
			for (int n = 0; n < nRoot; n++) {
				double sqrtX3 = root6[n];
				//std::cout << "root6 " << root6[n] << std::endl;
				if (sqrtX3 > 0) {
					x3 = sqrtX3 * sqrtX3;
					p = -sqrt((k2 * x4) / (c1122 * x3));

					x1 = p * c1 * x3;
					x2 = p * c2 * x3;

					if (abs(c1) > 0) {
						gamma = -(2 * x1 + c1) * x3 / (2 * x1);
					}
					else {
						gamma = -(2 * x2 + c2) * x3 / (2 * x2);
					}
					delta = k2 * gamma - 2 * x4 - c4;
					//std::cout << "x : " << x1 << " " << x2 << " " << x3 * k2 << " " << x4 << " " << gamma << " " << delta << std::endl;
					if (gamma >= 0 && delta >= 0 && x3 <= x3max && x3 >= x3min) {
						// the solution is good 
						typeSol = 8;
						goodSol = true;
						//nSol = n;
						break;
					}
					if (gamma > bestGamma && delta > bestGamma && (x3max - x3) > bestGamma && (x3 - x3min) > bestGamma) {
						typeSol = 8;
						bestGamma = min(min(min(gamma, (x3max - x3)), (x3 - x3min)), delta);
						BestRoot = n;
						neg = true;
					}
				}

			}

		}
		*/

		// cas x4 < x4 max 
			// case x3 = x3max lambdaLo = 0 delta = 0
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
				if (gamma >= 0 && lambdaUp >= 0 && x4 <= x4max) {
					// the solution is good 
					goodSol = true;
					typeSol = 2;
					//nSol = n;
					break;
				}
				if (gamma > bestGamma && lambdaUp > bestGamma && (x4max - x4) > bestGamma) {
					typeSol = 2;
					bestGamma = min((x4max - x4), min(gamma, lambdaUp));
					BestRoot = n;
				}

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

				if (gamma >= 0 && lambdaLo >= 0 && x4 <= x4max) {
					// the solution is good !
					typeSol = 3;
					goodSol = true;
					break;
				}
				if (gamma > bestGamma && lambdaLo > bestGamma && (x4max - x4) > bestGamma) {
					typeSol = 3;
					bestGamma = min((x4max - x4), min(gamma, lambdaLo));
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

				if (gamma >= 0 && x3 <= x3max && x3 >= x3min && x4 <= x4max) {
					// the solution is good !
					typeSol = 4;
					goodSol = true;
					break;
				}if (gamma > bestGamma && (x3max - x3) > bestGamma && (x3 - x3min) > bestGamma && (x4max - x4) > bestGamma) {
					typeSol = 4;
					bestGamma = min((x4max - x4), min(min(gamma, (x3max - x3)), (x3 - x3min)));
					BestRoot = n;
				}
			}
		}



		if (!goodSol) {
			std::cout << "*|*" << bestGamma << " " << typeSol << std::endl;	
			Chat[1].display();
			std::cout << "x : " << x1 << " " << x2 << " " << x3 * k2 << " " << x4 << " " << std::endl;
			if (typeSol == 1) {
				// case without constraint
				x1 = -c1 / 2;
				x2 = -c2 / 2;
				x3 = -c3 / 2;
				x3 = (x3max - x3) * (x3 > x3max) + (x3min - x3) * (x3min > x3) + x3;
				x4 = -c4 / 2; // ou  (x1 * x1 + x2 * x2) / (k2 * x3)
				x4 = (x4max - x4) * (x4 > x4max) + x4;
			}
			else if (typeSol > 4) {
				x4 = x4max;
				if (typeSol == 5) {
					x3 = x3max;
				}
				else if (typeSol == 6) {
					x3 = x3min;
				}
				else if (typeSol == 7) {
					x3 = root5[BestRoot];
				}
				else if (typeSol == 8) {
					x3 = root6[BestRoot];
				}
				p = sqrt((k2 * x4) / (c1122 * x3));
				if (neg) {
					p = -p;
				}

				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
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
				x4 = (x4max - x4) * (x4 > x4max) + x4;
			}
		}

		// X =  {Pi, Qi, vi, li, vAi, (pn, qn), (Pci, Qci, lci) for all child Ci}


		if (typeSol) {
			if (x4 > x4max) {
				std::cout << "probleme bus " << i << " " << x4max << " " << goodSol << " " << typeSol << std::endl;
				std::cout << "x : " << x1 << " " << x2 << " " << x3 * k2 << " " << x4 << " " << std::endl;
			}

			X[i].set(0, 0, x1);
			X[i].set(1, 0, x2);
			X[i].set(2, 0, x3 * k2);
			X[i].set(3, 0, x4);
		}
		else {
			std::cout << "pas de solution, pas de changement" << std::endl;
			std::cout << "probleme bus " << i << " " << x4max << " " << goodSol << " " << typeSol << std::endl;
			std::cout << "x : " << x1 << " " << x2 << " " << x3 * k2 << " " << x4 << " " << std::endl;

		}
		
		
		//std::cout << "x F : " << x1 << " " << x2 << " " << x3*k2 << " " << x4 << " " << gamma << std::endl;
		
		
		// pn & qn
		int Nb = _nAgentByBus.get(i, 0);
		int begin = _CoresAgentBusBegin.get(i, 0);
		for (int In = 0; In < Nb; In++) {
			int n = _CoresAgentBus.get(In + begin, 0);
			float ub = Pmax.get(n, 0);
			float lb = Pmin.get(n, 0);
			float pn = (_rho * Chat[i].get(4 + 2 * In, 0) - Cost2.get(n, 0)) / (Cost1.get(n, 0) + _rho);
			pn = ub * (ub < pn) + lb * (lb > pn) + pn * (pn >= lb) * (pn <= ub);

			ub = Pmax.get(n + _nAgent, 0);
			lb = Pmin.get(n + _nAgent, 0);
			float qn = (_rho * Chat[i].get(5 + 2 * In, 0) - Cost2.get(n + _nAgent, 0)) / (Cost1.get(n + _nAgent, 0) + _rho);
			qn = ub * (ub < qn) + lb * (lb > qn) + qn * (qn >= lb) * (qn <= ub);
			
			// pn & qn
			X[i].set(5 + 2 * In, 0, pn);
			X[i].set(6 + 2 * In, 0, qn);
		}

	}

	//bus fictif
	float pn = (_rho * Chat[_nBus].get(0, 0) - Cost2.get(0, 0)) / (Cost1.get(0, 0) + _rho);
	float qn = (_rho * Chat[_nBus].get(1, 0) - Cost2.get(_nAgent, 0)) / (Cost1.get(_nAgent, 0) + _rho);


	switch (losstype)
	{
	case LossType::CURRENT:
		X[_nBus].set(0, 0, pn);
		X[_nBus].set(1, 0, qn);
		break;
	case LossType::POWER:
		X[_nBus].set(0, 0, pn);
		X[_nBus].set(_nAgent, 0, qn);
		break;
	}

	
	

}

void OPFADMMCons::updateXWOCurrent()
{
	double x1, x2, x3, x4, c1, c2, c3, c4, lambdaLo, lambdaUp, x3min, x3max, gamma, k2;
	double c1122;
	int nSol = 0;
	int typeSol = 0;
	int BestRoot = 0;
	double bestGamma = -1;
	double p = 0;
	int nRoot = 0;

	for (int i = 1; i < _nBus; i++) {

		bool goodSol = false;
		k2 = sqrt(2.0 / (nChild.get(i, 0) + 1));


		c1 = -2 * Chat[i].get(0, 0);
		c2 = -2 * Chat[i].get(1, 0);
		c3 = -2 * Chat[i].get(2, 0) / k2;
		c4 = -2 * Chat[i].get(3, 0);
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
		gamma = k2 * x4 - (x1 * x1 + x2 * x2) / x3; // ce n'est pas vraiment gamma, doit être positif
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

		if (!goodSol) { // cas dégénéré
			if (c1122 == 0) {
				//std::cout << " bus " << i << " : c1= " << c1 << " c2=" << c2 << " c4=" << c4 << " gamma= " << gamma << std::endl;

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
				if (gamma >= 0 && lambdaUp >= 0) {
					// the solution is good 
					goodSol = true;
					//nSol = n;
				}
				if (gamma > bestGamma && lambdaUp > bestGamma) {
					typeSol = 2;
					bestGamma = min(gamma, lambdaUp);
					BestRoot = n;
				}

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
					bestGamma = min(gamma, lambdaLo);
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

			nRoot = resvolveRealPolynome4without2term(root4, coefPoly3);

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
					bestGamma = min(min(gamma, (x3max - x3)), (x3 - x3min));
					BestRoot = n;
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

		// X =  {Pi, Qi, vi, li, vAi, (pn, qn), (Pci, Qci, lci) for all child Ci}
		X[i].set(0, 0, x1);
		X[i].set(1, 0, x2);
		X[i].set(2, 0, x3 * k2);
		X[i].set(3, 0, x4);

		//std::cout << "x F : " << x1 << " " << x2 << " " << x3*k2 << " " << x4 << " " << gamma << std::endl;

	}
	for (int i = 0; i < _nBus; i++) {
		// pn & qn
		int Nb = _nAgentByBus.get(i, 0);
		int begin = _CoresAgentBusBegin.get(i, 0);
		for (int In = 0; In < Nb; In++) {
			int n = _CoresAgentBus.get(In + begin, 0);
			float ub = Pmax.get(n, 0);
			float lb = Pmin.get(n, 0);
			float pn = (_rho * Chat[i].get(4 + 2 * In, 0) - Cost2.get(n, 0)) / (Cost1.get(n, 0) + _rho);
			pn = (ub - pn) * (pn > ub) + (lb - pn) * (pn < lb) + pn;


			ub = Pmax.get(n + _nAgent, 0);
			lb = Pmin.get(n + _nAgent, 0);
			float qn = (_rho * Chat[i].get(5 + 2 * In, 0) - Cost2.get(n + _nAgent, 0)) / (Cost1.get(n + _nAgent, 0) + _rho);
			qn = (ub - qn) * (qn > ub) + (lb - qn) * (qn < lb) + qn;

			// pn & qn
			X[i].set(5 + 2 * In, 0, pn);
			X[i].set(6 + 2 * In, 0, qn);
		}
	}

	//bus fictif
	float pn = (_rho * Chat[_nBus].get(0, 0) - Cost2.get(0, 0)) / (Cost1.get(0, 0) + _rho);
	float qn = (_rho * Chat[_nBus].get(1, 0) - Cost2.get(_nAgent, 0)) / (Cost1.get(_nAgent, 0) + _rho);

	//std::cout << " bus fictif : " << _rho << " " << Chat[_nBus].get(0, 0) << " " << Cost1.get(0, 0) << " " << Cost2.get(0, 0) << std::endl;
	//std::cout << " pn : " << pn << std::endl;


	//X[_nBus].set(0, 0, pn);
	//X[_nBus].set(1, 0, qn);

	switch (losstype)
	{
	case LossType::CURRENT:
		X[_nBus].set(0, 0, pn);
		X[_nBus].set(1, 0, qn);
		break;
	case LossType::POWER:
		X[_nBus].set(0, 0, pn);
		X[_nBus].set(_nAgent, 0, qn);
		break;
	}/**/
}

void OPFADMMCons::updateMu()
{
	for (int i = 0; i < _nBusWLoss; i++) {
		tempM[i].subtract(&X[i], &Y[i]);
		tempM[i].multiply(_rho);
		Mu[i].add(&tempM[i]);
	}
}



float OPFADMMCons::getPLoss()
{
	_Ploss = 0;
	switch (losstype)
	{
	case LossType::POWER:
		for (int n = 1; n < _nAgent; n++) {
			int bus = _CoresBusAgent.get(n, 0);
			int In = PosAgent.get(n, 0);

			_Ploss -= X[bus].get(5 + 2 * In, 0);
		}
		break;
	case LossType::CURRENT:
		for (int i = 1; i < _nBus; i++) {
			_Ploss -= X[i].get(3, 0) * ZsRe.get(i - 1, 0);
		}
		break;
	}
	return _Ploss;
}

float OPFADMMCons::getQLoss()
{
	_Qloss = 0;

	switch (losstype)
	{
	case LossType::POWER:
		for (int n = 1; n < _nAgent; n++) {
			int bus = _CoresBusAgent.get(n, 0);
			int In = PosAgent.get(n, 0);

			_Qloss -= X[bus].get(6 + 2 * In, 0);
		}
		break;
	case LossType::CURRENT:
		for (int i = 1; i < _nBus; i++) {
			_Qloss -= X[i].get(3, 0) * ZsIm.get(i - 1, 0);
		}
		break;
	}
	return _Qloss;
}

void OPFADMMCons::updateChat()
{
	//Y      (Pi, Qi, vi, li, pi, qi, vai,(pn, qn), (Pji, Qji, lji))
	for (int i = 0; i < _nBus; i++) {
		int Nb = _nAgentByBus.get(i, 0);
		float Phat, Qhat, lhat, phat, qhat;
		float vhat = 0;
		float muhat = 0;
		float m = nChild.get(i, 0);
		
		int divideV = m + 1;
		int divideL = 2;
		int divideP = 1;
		
		if (losstype == LossType::CURRENT) {
			divideL += 1;
		}
		else if (losstype == LossType::POWER) { // POWER
			divideP += 1;
		}/**/

		// est ce que c'est y - MU/rho ou y + MU/rho (comme dans le papier mais le signe change je ne sais pas pourquoi)
		
		
		Phat = Y[i].get(0, 0) / 2 - Mu[i].get(0, 0) / (2 * _rho);
		Qhat = Y[i].get(1, 0) / 2 - Mu[i].get(1, 0) / (2 * _rho);
		lhat = Y[i].get(3, 0)     - Mu[i].get(3, 0) / _rho;
		if (i > 0) {
			int Ai = Ancestor.get(i, 0);
			int c = PosChild.get(i, 0);
			int NAi = _nAgentByBus.get(Ai, 0);
			Phat += Y[Ai].get(5 + 2 * NAi + 3 * c, 0) / 2 - Mu[Ai].get(5 + 2 * NAi + 3 * c, 0) / (2 * _rho);
			Qhat += Y[Ai].get(6 + 2 * NAi + 3 * c, 0) / 2 - Mu[Ai].get(6 + 2 * NAi + 3 * c, 0) / (2 * _rho);
			lhat += Y[Ai].get(7 + 2 * NAi + 3 * c, 0)     - Mu[Ai].get(7 + 2 * NAi + 3 * c, 0) / _rho;
			if (losstype == LossType::CURRENT) {
				lhat += Y[_nBus].get(i + 1, 0) + Mu[_nBus].get(i + 1, 0) / _rho;
			}/**/
		
		}
		for (int j = 0; j < m; j++) { // Vai of all childs
			int c = Childs[i].get(j, 0);
			muhat += Mu[c].get(4, 0);
			vhat += Y[c].get(4, 0);
		}
		
		vhat = (Y[i].get(2, 0) + vhat)  - (Mu[i].get(2, 0) + muhat) / _rho;
		
				

		Chat[i].set(0, 0, Phat);
		Chat[i].set(1, 0, Qhat);
		Chat[i].set(2, 0, vhat / divideV);
		Chat[i].set(3, 0, lhat / divideL);
		int begin = _CoresAgentBusBegin.get(i, 0);

		for (int In = 0; In < Nb; In++) {
			int n = _CoresAgentBus.get(In + begin, 0);
			phat = Y[i].get(5 + 2 * In, 0) - Mu[i].get(5 + 2 * In, 0) / _rho;
			qhat = Y[i].get(6 + 2 * In, 0) - Mu[i].get(6 + 2 * In, 0) / _rho;
			if (losstype == LossType::POWER) {
				phat +=  Y[_nBus].get(n, 0)			  - Mu[_nBus].get(n, 0) / (_rho);
				qhat +=  Y[_nBus].get(n + _nAgent, 0) - Mu[_nBus].get(n + _nAgent, 0) / (_rho);
			}
		
			
			Chat[i].set(4 + 2 * In, 0, phat / divideP);
			Chat[i].set(5 + 2 * In, 0, qhat / divideP);
		}
		
	}
	float phatLoss = 0;// Y[_nBus].get(0, 0) - Mu[_nBus].get(0, 0) / _rho;
	float qhatLoss = 0;// Y[_nBus].get(1, 0) - Mu[_nBus].get(1, 0) / _rho;
	switch (losstype)
	{
	case LossType::POWER:
		phatLoss = Y[_nBus].get(0, 0) - Mu[_nBus].get(0, 0) / _rho;
		qhatLoss = Y[_nBus].get(_nAgent, 0) - Mu[_nBus].get(_nAgent, 0) / _rho;
		break;
	case LossType::CURRENT:
		phatLoss = Y[_nBus].get(0, 0) - Mu[_nBus].get(0, 0) / _rho;
		qhatLoss = Y[_nBus].get(1, 0) - Mu[_nBus].get(1, 0) / _rho;
		break;
	}/**/

	// Y     (Ploss, Pn, Qn)
	Chat[_nBus].set(0, 0, phatLoss);
	Chat[_nBus].set(1, 0, qhatLoss);

	// le reste c'est que des copies ?

}

void OPFADMMCons::CommunicationX()
{
/**/ // X = { Pi, Qi, vi, li, vAi, (pn, qn) (Pci, Qci, lci) for all child Ci }
	
	for (int i = 0; i < _nBus; i++) {
		int Ni = _nAgentByBus.get(i, 0);
		if (i > 0) {
			int Ai = Ancestor.get(i, 0);
			X[i].set(4, 0, X[Ai].get(2, 0));
		}

		int m = nChild.get(i, 0);
		for (int j = 0; j < m; j++) {
			int c = Childs[i].get(j, 0);
			X[i].set(5 + 2 * Ni + 3 * j, 0, X[c].get(0, 0));
			X[i].set(6 + 2 * Ni + 3 * j, 0, X[c].get(1, 0));
			X[i].set(7 + 2 * Ni + 3 * j, 0, X[c].get(3, 0));

		}
	}
	// bus fictif
	switch (losstype)
	{
	case LossType::POWER:
		for (int n = 1; n < _nAgent; n++) {
			int bus = _CoresBusAgent.get(n, 0);
			int In = PosAgent.get(n, 0);

			X[_nBus].set(n, 0, X[bus].get(5 + 2 * In, 0));
			X[_nBus].set(n + _nAgent, 0, X[bus].get(6 + 2 * In, 0));
		}
		break;
	case LossType::CURRENT:
		for (int i = 1; i < _nBus; i++) {
			X[_nBus].set(1 + i, 0, X[i].get(3, 0));
		}
		break;
	}/**/
	
	


	
	//Y      (Pi, Qi, vi, li, pi, qi, vai, (pn, qn) Pji, Qji, lji)
	// Q udate in argmin 0.5yHy + Qy

	for (int i = 0; i < _nBusWLoss; i++) {
		for (int j = 0; j < sizeOPFADMMCons.get(i, 0); j++) {	
			Q[i].set(j, 0, -(Mu[i].get(j, 0) + _rho * X[i].get(j, 0)));
		}
	}   

}


float OPFADMMCons::updateRes(int indice) 
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

	/*std::cout << resS << " " << resR << std::endl;

	for (int i = 0; i < _nBus; i++) {
		std::cout << " Y " << std::endl;
		Y[i].display();
	}

	for (int i = 0; i < _nBus; i++) {
		std::cout << " X " << std::endl;
		X[i].display();
	}*/
	if (_tau > 1) {
		if (resR > _mu * resS) {
			_rho = _tau * _rho;
			for (int i = 0; i < _nBusWLoss; i++) {
				Hinv[i].divide(_tau);
			}
		 
			//std::cout << _iterGlobal << "rho augmente " << _rho << std::endl;
		}
		else if (resS > _mu * resR) {// rho = rho / tau_inc;
			_rho = _rho / _tau;
			for (int i = 0; i < _nBusWLoss; i++) {
				Hinv[i].multiply(_tau);
			}
		
		 
			//std::cout << _iterGlobal << "rho diminue " << _rho << std::endl;
		}
	}
	


	return MAX(MAX(resV, oldrho * resS), resR);
}

float OPFADMMCons::updateResRhoFixe(int indice)
{
	float resS = 0;
	float resR = 0;
	float resV = 0;
	for (int i = 0; i < _nBusWLoss; i++) {

		float resTempS = _rho * Y[i].max2(&Ypre[i]);
		float resTempR = Y[i].max2(&X[i]);

		if (resTempS > resS) {
			resS = resTempS;
		}
		if (resTempR > resR) {
			resR = resTempR;
		}

	}
	resF.set(0, indice, resR);
	resF.set(1, indice, resS);
	resF.set(2, indice, resV);

	return MAX(MAX(resV, resS), resR);
}

int OPFADMMCons::feasiblePoint()
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




void OPFADMMCons::display() {

	std::cout.precision(3);
	Pb.set(0.0);
	for (int i = 0; i < _nBus; i++) {
		int Nb = _nAgentByBus.get(i, 0);
		int begin = _CoresAgentBusBegin.get(i, 0);
		for (int In = 0; In < Nb; In++) {
			int n = _CoresAgentBus.get(In + begin, 0);
			Pb.increment(i, 0, Pn.get(n, 0));
			Pb.increment(i + _nBus, 0, Pn.get(n + _nAgent, 0));
		}
	}



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
	std::cout << "Ploss            " << X[_nBus].get(0, 0) << std::endl;
	std::cout << "Qloss            " << X[_nBus].get(1, 0) << std::endl;


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
			<< (abs(Pb.get(b, 0)) > seuil) * Pb.get(b,0) << "|" << std::setw(16) << (abs(Pb.get(b + _nBus, 0)) > seuil) * Pb.get(b + _nBus, 0)
			<< "|" << std::setw(16) << Mu[b].get(2, 0) << "|" << std::setw(16)
			<< Mu[b].get(0, 0) << "|" << std::setw(16) << Mu[b].get(1, 0) << "|" << std::endl;

	}
	std::cout << std::endl << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "      Line Data                                                                                         |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Line |    From     |    To      |                Upstream flow                      |    Constraint    |" << std::endl;
	std::cout << "  #   |    Bus      |    Bus     |   P (pu)   |   Q (pu)   |   l (pu)   |  Loss (pu) |      lmax        |" << std::endl;
	std::cout << "------|-------------|------------|------------|------------|------------|------------|------------------|" << std::endl;

	for (int l = 0; l < _nLine; l++) {
		int b = l + 1;
		std::cout << std::setw(6) << l << "|" << std::setw(12) << CoresLineBus.get(l, 0) << " |" << std::setw(12)
			<< CoresLineBus.get(l, 1) << "|" << std::setw(12) << X[b].get(0, 0)
			<< "|" << std::setw(12) << X[b].get(1, 0) << "|" << std::setw(12)
			<< X[b].get(3, 0) << "|" << std::setw(12) << X[b].get(3, 0) * ZsRe.get(l, 0) << "|" <<
			std::setw(18) << FluxLimit.get(l,0) << "|" << std::endl;
	}
	std::cout << std::endl << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "     Constraints                                                                                        |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Bus | Voltage | Voltage | Voltage |        Power Injection          |          Power Injection         |" << std::endl;
	std::cout << "  #  | Mag(pu) | MIN(pu) |  MAX(pu)|  P (pu) | Pmin (pu) | Pmax (pu) |  Q (pu)  | Qmin (pu) | Qmax (pu) |" << std::endl;
	std::cout << "-----|---------|---------|---------|---------|-----------|-----------|----------|-----------|-----------|" << std::endl;
	

	for (int b = 0; b < _nBus; b++) {
		int nb = _nAgentByBus.get(b, 0);
		std::cout << std::setw(5) << b << "|" << std::setw(8) << sqrt(Y[b].get(2, 0)) << " |" << std::setw(9)
			<< VoltageLimitReal.get(b, 0) << "|" << std::setw(9) << VoltageLimitReal.get(b, 1)
			<< "|" << std::setw(9) << Pb.get(b, 0) << "|" << std::setw(11)
			<< Pbmin.get(b, 0) << "|" << std::setw(11) << Pbmax.get(b,0) << "|" << std::setw(10) << Pb.get(b + _nBus, 0)
			<< "|" << std::setw(11) << Pbmin.get(b + _nBus, 0) << "|" << std::setw(11) << Pbmax.get(b + _nBus, 0) << "|" << std::endl;

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

float OPFADMMCons::DFSP(int j)
{
	//std::cout << "DFSP " << j << std::endl;
	float p = Pb.get(j,0);
	for (int i = 0; i < nChild.get(j, 0); i++) {
		int c = Childs[j].get(i, 0);
		p += DFSP(c);
	}
	X[j].set(0, 0, p);
	return p;
}
float OPFADMMCons::DFSQ(int j)
{
	float q = Pb.get(j + _nBus, 0);
	for (int i = 0; i < nChild.get(j, 0); i++) {
		int c = Childs[j].get(i, 0);
		q += DFSQ(c);
	}
	X[j].set(1, 0, q);
	return q;
}




















