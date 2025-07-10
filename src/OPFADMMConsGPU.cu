#include "../head/OPFADMMConsGPU.cuh"
 


OPFADMMConsGPU::OPFADMMConsGPU() : MethodOPFGPU()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " OPFADMMConsGPU Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	timePerBlock = MatrixCPU(1, 12, 0); // Fb0, Fb11abcd, FB12, Fb2, Fb3, Fb4, Fb5,FB6, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 12, 0); //nb de fois utilisé pendant la simu
}


OPFADMMConsGPU::OPFADMMConsGPU(float rho) : MethodOPFGPU()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default OPFADMMConsGPU Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
	timePerBlock = MatrixCPU(1, 12, 0); // Fb0, Fb11, FB12, Fb2, Fb3, Fb4, Fb5, FB6, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 12, 0); //nb de fois utilisé pendant la simu
}

OPFADMMConsGPU::~OPFADMMConsGPU()
{
	 
}
void OPFADMMConsGPU::setParam(float rho)
{
	_rho = rho;
}

bool OPFADMMConsGPU::chekcase()
{
	if (_nBus != (_nLine + 1)) {
		std::cout << "wrong number of line " << _nLine << "against " << _nBus << std::endl;
		return false;
	}
	//CoresLineBus.display();
	for (int i = 0; i < _nLine; i++) {
		if (CoresLineBusCPU.get(i, 1) != (i + 1)) {
			std::cout << "wrong numerotation of line " << CoresLineBusCPU.get(i, 1) << "against " << (i + 1) << std::endl;
			return false;
		}
		if (CoresLineBusCPU.get(i, 0) > CoresLineBusCPU.get(i, 1)) {
			std::cout << "wrong numeoration of bus " << CoresLineBusCPU.get(i, 0) << "against " << CoresLineBusCPU.get(i, 1) << std::endl;
			return false;
		}
	}
	if (ZsRe.getNLin() == 0  || ZsIm.getNLin() == 0) {
		std::cout << "matrice non defined, ZsRe, Zs Im, Yd" << std::endl;
		ZsRe.display(true);
		ZsIm.display(true);
		return false;
	}

	return true;
}

void OPFADMMConsGPU::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
	_iterGlobal = 0;
	

	//Chat.display();
	//Bpt2.display();
	
	while ((_iterGlobal < _iterG) && (resG>epsG)) {
		
		
		updateXWOCurrent();


		// 
		//updateXWOCurrentOnCPU();
		//updateXWOCurrentOnCPUBis();
		CHECK_LAST_CUDA_ERROR();
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 5, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		
		CommunicationX();
		
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		updateGlobalProb();
		CHECK_LAST_CUDA_ERROR();
		updateMu();
		CHECK_LAST_CUDA_ERROR();
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		updateChat();

		CHECK_LAST_CUDA_ERROR();
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
			//resG = 1;
#ifdef INSTRUMENTATION
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 9, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		}
		//std::cout << iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;

		_iterGlobal++;
	}
	std::cout << "---------------------------------" << std::endl;
		/*std::cout << " X " << std::endl;
		X.display(true);
		std::cout << " Chat " << std::endl;
		Chat.display(true);
		std::cout << " Q " << std::endl;
		Q.display(true);
		std::cout << "Y " << std::endl;
		Y.display(true);
		std::cout << "Mu " << std::endl;
		Mu.display(true);
		*/
	std::cout << _iterGlobal << " " << resF.get(0, (_iterGlobal - 1) / _stepG) << " " << resF.get(1, (_iterGlobal - 1) / _stepG) << " " << resF.get(2, (_iterGlobal - 1) / _stepG) << std::endl;


#ifdef INSTRUMENTATION	
	occurencePerBlock.increment(0, 5, _iterGlobal);
	occurencePerBlock.increment(0, 6, _iterGlobal);
	occurencePerBlock.increment(0, 7, _iterGlobal);
	occurencePerBlock.increment(0, 8, _iterGlobal);
	occurencePerBlock.increment(0, 9, _iterGlobal / _stepG);

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	//setPnFromX << < _nBus, _blockSizeSmall >> > (Pn._matrixGPU, X._matrixGPU, _indiceBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent);
	////CHECK_LAST_CUDA_ERROR();
	MatrixCPU PnCPU;
	Pn.toMatCPU(PnCPU);
	PnCPU.set(0, 0, getPLoss());
	PnCPU.set(_nAgent, 0, getQLoss());
	
	fc = calcFc(&Cost1, &Cost2, &Pn, &tempN2);
	// FB 5
	
	result->setResF(&resF);
	
	/*std::cout << "--------" << std::endl;
	std::cout << " Pn " << std::endl;
	Pn.display();
	for (int i = 0; i < 3; i++) {
		std::cout << " X " << i << std::endl;
		X[i].display();
	}
	*/
	MatrixCPU Pb(getPb());
	MatrixCPU Phi(getPhi());
	MatrixCPU E(getE());

	result->setE(&E);
	result->setPhi(&Phi);
	result->setPb(&Pb);
	

	result->setIter(_iterGlobal);
	

	result->setPn(&PnCPU);
	
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

void OPFADMMConsGPU::updateP0(const StudyCase& cas)
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
	
	// remove loss agent
	Pn.set(0, 0, 0, 1);
	Pmin.set(0, 0, 0, 1);
	Pmax.set(0, 0, 0, 1);
	Pn.set(_nAgent, 0, 0, 1);
	Pmin.set(_nAgent, 0, 0, 1);
	Pmax.set(_nAgent, 0, 0, 1);
	
	
	ComputePFromAgentToBus();

	initPQAgent << < _nBus, _blockSize >> > (X._matrixGPU, _indiceBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, Pn._matrixGPU, _nAgent);

	//_global__ void initDFSPQ(float* X, float* Pb, float* nChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, int nBus)
	initDFSPQ << <1, _nBus, _nBus* (sizeof(bool) + sizeof(int)) >> > (X._matrixGPU, Pb._matrixGPU, nChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nBus);
	communicateX << <_nBus, _blockSize >> > (X._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nAgentByBus._matrixGPU, _nBus);


	Y.set(&X);

	updateChat();
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 11, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 11, 1);
#endif // INSTRUMENTATION

}

void OPFADMMConsGPU::init(const Simparam& sim, const StudyCase& cas)
{

	if (_CoresChatBegin.getPos()) {

		_CoresChatBegin.transferCPU();
		_indiceBusBegin.transferCPU();

		Ancestor.transferCPU();
		PosChild.transferCPU();

		_indiceChildBegin.transferCPU();
		Childs.transferCPU();

		ZsIm.transferCPU();
		ZsRe.transferCPU();
	}


	// intitilisation des matrixs et variables 

	clock_t t = clock();
	
	std::cout << "init " << std::endl;
	_rho = sim.getRho();
	
	_iterG = sim.getIterG();
	_stepG = sim.getStepG();

	_nAgent = cas.getNagent();

	_nBus = cas.getNBus();
	_nBusWLoss = _nBus + 1;
	_nLine = cas.getNLine(true); // ne doit pas être réduit ici !!!

	_debutloss = 3 * _nLine + 5 * _nBus + 2 * (_nAgent - 1); // L = nChild.sum()
	_sizeOPFADMMConsTotal = _debutloss;
	_sizeChat = 4 * _nBus + 2 * _nAgent;

	if (losstype == LossType::CURRENT) {
		_sizeOPFADMMConsTotal += (_nBus + 2); // pertes et courants sauf premier bus ou + 2
	}
	else if (losstype == LossType::POWER) {
		_sizeOPFADMMConsTotal += _nAgent;
	}

	_numBlocksB = ceil((_nBus + _blockSize - 1) / _blockSize);
	_numBlocksH = ceil((_sizeOPFADMMConsTotal + _blockSize - 1) / _blockSize);
	_numBlocksN = ceil((_nAgent + _blockSize - 1) / _blockSize);

	


	//std::cout << _nAgent << " " << _nBus << " " << _nLine << std::endl;
	_CoresAgentBus = MatrixGPU(cas.getCoresAgentBusLin(), 1);
	_CoresAgentBusBegin = MatrixGPU(cas.getCoresAgentBusLinBegin(), 1);
	_nAgentByBus = MatrixGPU(cas.getNagentByBus(), 1);
	_nAgentByBusCPU = cas.getNagentByBus();
	PosAgent = MatrixGPU(_nAgent, 1, 0, 1);

	initPosAgent << <_nBus, _blockSizeSmall >> > (PosAgent._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU);
	
	nChildCPU = MatrixCPU(_nBus, 1);

	CoresLineBusCPU = cas.getCoresLineBus(true);
	CoresLineBus = MatrixGPU(CoresLineBusCPU, 1);

	_CoresBusAgent = MatrixGPU(cas.getCoresBusAgentLin(), 1); // Cores[n] = b

	Ancestor = MatrixGPU(_nBus, 1, 0); // A_i = bus antécédent de i
	PosChild = MatrixGPU(_nBus, 1, 0); // indice du bus i dans Child[Ai]
	Ancestor.set(0, 0, -1); // the slack bus has no ancestor

	ZsRe = cas.getZsRe();
	ZsIm = cas.getZsImag();
	ZsNorm = MatrixCPU(_nLine, 1);

	if (!chekcase()) {
		throw std::invalid_argument("not a radial case");
	}
	if (losstype != LossType::CURRENT && losstype != LossType::POWER) {
		throw std::invalid_argument("unkown lossType");
	}



	for (int lold = 0; lold < _nLine; lold++) {
		int l = lold + 1;
		int busTo = l;
		int busFrom = CoresLineBusCPU.get(lold, 0);
		Ancestor.set(busTo, 0, busFrom);
		nChildCPU.increment(busFrom, 0, 1);
		ZsNorm.set(lold, 0, ZsRe.get(lold, 0) * ZsRe.get(lold, 0) + ZsIm.get(lold, 0) * ZsIm.get(lold, 0));
	}
	nChild = MatrixGPU(nChildCPU, 1);

	/*
	ZsNorm.display();*/


	_rhoInv = 1 / _rho;
	resF = MatrixCPU(3, (_iterG / _stepG) + 1);


	std::cout << " local resolution " << std::endl;
	// local resolution
	tempN2 = MatrixGPU(2 * _nAgent, 1, 0, 1);
	tempB2 = MatrixGPU(2 * _nBus, 1, 0, 1);

	CoresSoloBusAgent = MatrixGPU(_nBus, 1, -1, 1);
	Pn = MatrixGPU(sim.getPn(), 1); // not the real agent
	Pmin = MatrixGPU(cas.getPmin(), 1); 
	Pmax = MatrixGPU(cas.getPmax(), 1); // idem

	// the loss provider
	/*Pmin.set(0, 0, 0);
	Pmax.set(0, 0, 0);
	Pmin.set(_nAgent, 0, 0);
	Pmax.set(_nAgent, 0, 0);*/

	Pbmax = MatrixGPU(2 * _nBus, 1, 0, 1);
	Pbmin = MatrixGPU(2 * _nBus, 1, 0, 1);
	Pb = MatrixGPU(2 * _nBus, 1, 0, 1);
	//Pmin.display();
	//Pn.display();


	Cost1 = MatrixGPU(cas.geta(), 1);
	Cost2 = MatrixGPU(cas.getb(), 1);

	Pn.preallocateReduction();
	if (Pn.max2() < 0.00001) {
		Pn.add(&Pmin, &Pmax);
		Pn.divide(2);
	}
	Pn.set(0, 0, 0, 1);
	Pmin.set(0, 0, 0, 1);
	Pmax.set(0, 0, 0, 1);
	Pn.set(_nAgent, 0, 0, 1);
	Pmin.set(_nAgent, 0, 0, 1);
	Pmax.set(_nAgent, 0, 0, 1);

	_nAgentByBusCPU.increment(0, 0, -1);
	_nAgentOn0 = _nAgentByBusCPU.get(0, 0);

	removeLossAgent << <1, 1 >> > (_nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU);
	 
	ComputePFromAgentToBus();
	 
	//_nAgentByBus.display();


	std::cout << " creation " << std::endl;
	X = MatrixGPU(_sizeOPFADMMConsTotal, 1, 0, 1);
	Ypre = MatrixGPU(_sizeOPFADMMConsTotal, 1, 0, 1);
	Y = MatrixGPU(_sizeOPFADMMConsTotal, 1, 0, 1);
	Y.preallocateReduction();

	Mu = MatrixGPU(_sizeOPFADMMConsTotal, 1, 0, 1);

	tempNN = MatrixGPU(_nAgent, _nAgent, 0, 1);
	//tempM1 = new MatrixGPU[_nAgent];
	tempM = MatrixGPU(_sizeOPFADMMConsTotal, 1, 0, 1);

	sizeOPFADMMConsGPU = MatrixGPU(_nBusWLoss, 1, 0, 1);
	sizeOPFADMMConsGPU.preallocateReduction();
	sizeOPFADMMConsGPUBig = MatrixGPU(_sizeOPFADMMConsTotal, 1, 0, 1);

	indiceBusBeginCPU = MatrixCPU(_nBusWLoss, 1);
	_indiceBusBeginBig = MatrixGPU(_sizeOPFADMMConsTotal, 1, 0, 1);
	CoresChatBeginCPU = MatrixCPU(_nBusWLoss, 1);


	int debut = 0;
	int debutChat = 0;
	for (int i = 0; i < _nBus; i++) {
		int m = nChildCPU.get(i, 0);
		int nB = _nAgentByBusCPU.get(i, 0);
		indiceBusBeginCPU.set(i, 0, debut);
		CoresChatBeginCPU.set(i, 0, debutChat);
		int sizeA = m * 3 + 5 + 2 * nB;
		debut += sizeA;
		debutChat += 4 + 2 * nB;
	}
	indiceBusBeginCPU.set(_nBus, 0, debut);
	CoresChatBeginCPU.set(_nBus, 0, debutChat);


	_CoresChatBegin = MatrixGPU(CoresChatBeginCPU, 1);
	_indiceBusBegin = MatrixGPU(indiceBusBeginCPU, 1);
	defineSizeBig << <_nBusWLoss, _blockSize >> > (sizeOPFADMMConsGPUBig._matrixGPU, nChild._matrixGPU, _indiceBusBegin._matrixGPU, sizeOPFADMMConsGPU._matrixGPU, _indiceBusBeginBig._matrixGPU, _nAgentByBus._matrixGPU, losstype, _nBus, _nAgent);
	 

	_sizeOPFADMMConsMax = sizeOPFADMMConsGPU.max2();
	Hinv = MatrixGPU(_sizeOPFADMMConsTotal, _sizeOPFADMMConsMax, 0, 1);
	Q = MatrixGPU(_sizeOPFADMMConsTotal, 1, 0, 1);

	Childs = MatrixGPU(_nLine, 1);
	PosChild = MatrixGPU(_nBus, 1, -1);
	_indiceChildBegin = MatrixGPU(_nBus, 1);
	Chat = MatrixGPU(_sizeChat, 1, 0, 1);

	MatrixGPU lowerBound(cas.getLowerBound(), 1); //voltage angle, voltage, line...
	MatrixGPU upperBound(cas.getUpperBound(), 1); //voltage angle, voltage, line...

	VoltageLimit = MatrixGPU(2, _nBus, 0, 1); // min, max
	VoltageLimitCPU = MatrixCPU(2, _nBus, 0);
	VoltageLimitReal = MatrixGPU(2, _nBus, 0, 1); // min, max

	if (cas.isCurrentLimit()) {
		FluxLimit = MatrixGPU(cas.getCurrentLimit(), 1);
		isCurrentLimited = true;
	}
	else {
		FluxLimit = MatrixGPU(_nLine, 1, 1000, 1); // max
	}
	//FluxLimit.display();
	initVoltageBound << < _numBlocksB, _blockSize >> > (VoltageLimitReal._matrixGPU, VoltageLimit._matrixGPU, lowerBound._matrixGPU, upperBound._matrixGPU, nChild._matrixGPU, _nBus);
	VoltageLimit.toMatCPU(VoltageLimitCPU);


	int debutChild = 0;
	MatrixCPU nChildTemp(_nBus, 1, 0);

	for (int i = 0; i < _nBus; i++) {
		if (i > 0) {
			_indiceChildBegin.set(i - 1, 0, debutChild);

			int Ai = Ancestor.get(i, 0);
			Childs.set(_indiceChildBegin.get(Ai, 0) + nChildTemp.get(Ai, 0), 0, i);
			PosChild.set(i, 0, nChildTemp.get(Ai, 0));
			nChildTemp.increment(Ai, 0, 1);
			debutChild += nChildCPU.get(i - 1, 0);
		}

	}
	Ancestor.transferGPU();
	PosChild.transferGPU();

	debut = 0;
	//std::cout << " Hinv " << std::endl;
	for (int i = 0; i < _nBus; i++) {
		// (Pi, Qi, li, vi, pn..., qn..., vai, Pci ..., Qci... , lci...) !!!!!
		int m = nChildCPU.get(i, 0);
		int nB = _nAgentByBusCPU.get(i, 0);
		int sizeA = m * 3 + 5 + 2 * nB;
		MatrixCPU A(2 + 1 * (i > 0), sizeA);

		if (i > 0) {
			A.set(2, 0, 2 * ZsRe.get(i - 1, 0));
			A.set(2, 1, 2 * ZsIm.get(i - 1, 0));
			A.set(2, 2, -ZsNorm.get(i - 1, 0));
			A.set(2, 3, -1);
			A.set(2, 4 + 2 * nB, 1); // vai
			A.set(0, 0, -1);
			A.set(1, 1, -1);
		}
		for (int In = 0; In < nB; In++) {
			A.set(0, 4 + In, 1);
			A.set(1, 4 + nB + In, 1);
		}

		for (int j = 0; j < m; j++) {
			int c = Childs.get(_indiceChildBegin.get(i, 0) + j, 0);
			A.set(0, 5 + 2 * nB + j, 1); // Pci
			A.set(1, 5 + 2 * nB + m + j, 1); // Qci
			A.set(0, 5 + 2 * nB + 2 * m + j, -ZsRe.get(c - 1, 0)); // -R l
			A.set(1, 5 + 2 * nB + 2 * m + j, -ZsIm.get(c - 1, 0)); // -X l
		}

		//A.display();

		MatrixCPU temp33(2 + 1 * (i > 0), 2 + 1 * (i > 0));
		MatrixCPU temp3M(2 + 1 * (i > 0), sizeA);
		MatrixCPU tempMM(sizeA, sizeA);
		MatrixCPU tempMMbis(sizeA, sizeA);

		temp33.multiplyTrans(&A, &A);
		temp33.invertGaussJordan(&temp33);
		temp3M.MultiplyMatMat(&temp33, &A);
		tempMM.multiplyTrans(&A, &temp3M, 0);

		tempMMbis.setEyes(-1);
		tempMMbis.add(&tempMM);
		MatrixGPU tempMMGPU = MatrixGPU(tempMMbis, 1);

		Hinv.setBloc(debut, debut + sizeA, 0, sizeA, &tempMMGPU);
		debut += sizeA;
	}

	// bus factice
	int sizeA = 0;
	MatrixGPU A;
	switch (losstype)
	{
	case LossType::POWER:
		sizeA = 2 * _nAgent;
		A = MatrixGPU(sizeA, sizeA);
		A.setEyes(-1);
		A.set(0, 0, 0);
		A.set(_nAgent, _nAgent, 0);
		for (int i = 1; i < _nAgent; i++) {
			A.set(0, i, 1); // sum(p) + Ploss = 0
			A.set(_nAgent, i + _nAgent, 1); // Qloss + sum(q) = 0
		}
		A.transferGPU();
		Hinv.setBloc(debut, debut + sizeA, 0, sizeA, &A);
		/*for (int i = 0; i < _nAgentTrue; i++) {
			A[_nBus].set(0, i, 1); // sum(p) + Ploss = 0
			A[_nBus].set(1, i + _nAgentTrue, 1); // Qloss + sum(q) = 0
		}*/
		break;
	case LossType::CURRENT:
		sizeA = 2 + _nBus;
		A = MatrixGPU(sizeA, sizeA);
		A.setEyes(-1);
		A.set(0, 0, 0);
		A.set(1, 1, 0);
		for (int i = 0; i < _nLine; i++) {
			A.set(0, i + 3, ZsRe.get(i, 0)); // sum(p) + Ploss = 0
			A.set(1, i + 3, ZsIm.get(i, 0)); // Qloss + sum(q) = 0
		}
		A.transferGPU();
		Hinv.setBloc(debut, debut + sizeA, 0, sizeA, &A);
		//A[_nBus].set(0, 0, 1); // ploss
		//A[_nBus].set(1, 1, 1); // qloss
		break;
	}
	Hinv.divide(_rho);
	ZsRe.transferGPU();
	ZsIm.transferGPU();
	_indiceChildBegin.transferGPU();
	Childs.transferGPU();

	//	std::cout << " init valeur " << std::endl;

	initPQAgentV << < _nBus, _blockSizeSmall >> > (X._matrixGPU, _indiceBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, Pn._matrixGPU, _nAgent);
	

	initDFSPQ << <1, _nBus, _nBus* (8*sizeof(bool) + sizeof(int)) >> > (X._matrixGPU, Pb._matrixGPU, nChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nBus);
	//std::cout << " X " << std::endl;
	//X.display(true);

	communicateX << <_nBusWLoss, _blockSize >> > (X._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nAgentByBus._matrixGPU, _CoresBusAgent._matrixGPU, PosAgent._matrixGPU, losstype, _nBus, _nAgent);
	 

	computeLoss();
	 

	Y.set(&X);

	// bus factice




	//std::cout << "updateChat" << std::endl;
	updateChat();
	 

	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	//std::cout << "---------------------------------------------------------------------------------------" << std::endl;
}

void OPFADMMConsGPU::solveConsensus(float eps, MatrixGPU* PSO)
{
	float epsG = eps;
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
	
	//if (_iterG == _iterGlobal) {
		//std::cout << "OPF GPU : " << _iterGlobal << " " << resF.get(0, (_iterGlobal - 1) / _stepG) << " " << resF.get(1, (_iterGlobal - 1) / _stepG) << " " << resF.get(2, (_iterGlobal - 1) / _stepG) << std::endl;
	//}
	//X[1].display();
	
	//setPnFromX << < _nBus, _blockSizeSmall >> > (Pn._matrixGPU, X._matrixGPU, _indiceBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent);

	//MatrixCPU PnCPU;
	//Pn.toMatCPU(PnCPU);
	
	
	//X.display(true);
	//Y[_nBus].display();
	//Pn.display(true);

	//Pn.set(0, 0, _Ploss);
	//Pn.set(_nAgent, 0, _Qloss);
	
	PSO->set(&Pn);
	PSO->set(0, 0, getPLoss(), true);
	PSO->set(_nAgent, 0, getQLoss(), true);

	//PSO->display();

	timeOPF = clock() - timeOPF;

}


void OPFADMMConsGPU::solveConsensus(float eps, MatrixCPU* PSO) {
	
	solveConsensus(eps, &PSOGPU);

	PSOGPU.toMatCPU(*PSO);

}


void OPFADMMConsGPU::initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO)
{
	if (_CoresChatBegin.getPos()) {

		_CoresChatBegin.transferCPU();
		_indiceBusBegin.transferCPU();

		Ancestor.transferCPU();
		PosChild.transferCPU();

		_indiceChildBegin.transferCPU();
		Childs.transferCPU();

		ZsIm.transferCPU();
		ZsRe.transferCPU();
	}
	// intitilisation des matrixs et variables 

	clock_t t = clock();
	//std::cout << "init OPF " << std::endl;
	_rho = sim.getRho();
	_rhoSO = rhoSO;

	_iterG = sim.getIterG();
	_stepG = sim.getStepG();

	_nAgent = cas.getNagent();

	_nBus = cas.getNBus();
	_nBusWLoss = _nBus + 1;
	_nLine = cas.getNLine(true); // ne doit pas être réduit ici !!!

	_debutloss = 3 * _nLine + 5 * _nBus + 2 * (_nAgent - 1); // L = nChild.sum()
	_sizeOPFADMMConsTotal = _debutloss;
	_sizeChat = 4 * _nBus + 2 * _nAgent;

	if (losstype == LossType::CURRENT) {
		_sizeOPFADMMConsTotal += (_nBus + 2); // pertes et courants sauf premier bus ou + 2
	}
	else if (losstype == LossType::POWER) {
		_sizeOPFADMMConsTotal += _nAgent;
	}

	_numBlocksB = ceil((_nBus + _blockSize - 1) / _blockSize);
	_numBlocksH = ceil((_sizeOPFADMMConsTotal + _blockSize - 1) / _blockSize);
	_numBlocksN = ceil((_nAgent + _blockSize - 1) / _blockSize);
	
	
	//std::cout << _nAgent << " " << _nBus << " " << _nLine << std::endl;

	_CoresAgentBus = MatrixGPU(cas.getCoresAgentBusLin(), 1);
	_CoresAgentBusBegin = MatrixGPU(cas.getCoresAgentBusLinBegin(), 1);
	_nAgentByBus = MatrixGPU(cas.getNagentByBus(), 1);
	_nAgentByBusCPU = cas.getNagentByBus();
	PosAgent = MatrixGPU(_nAgent, 1, 0, 1);

	////CHECK_LAST_CUDA_ERROR();
	initPosAgent << <_nBus, _blockSizeSmall >> > (PosAgent._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU);
	////CHECK_LAST_CUDA_ERROR();


	nChildCPU = MatrixCPU(_nBus, 1);

	CoresLineBusCPU = cas.getCoresLineBus(true);
	CoresLineBus = MatrixGPU(CoresLineBusCPU, 1);

	_CoresBusAgent = MatrixGPU(cas.getCoresBusAgentLin(), 1); // Cores[n] = b

	Ancestor = MatrixGPU(_nBus, 1, 0); // A_i = bus antécédent de i
	PosChild = MatrixGPU(_nBus, 1, 0); // indice du bus i dans Child[Ai]
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
		int busFrom = CoresLineBusCPU.get(lold, 0);
		Ancestor.set(busTo, 0, busFrom);
		nChildCPU.increment(busFrom, 0, 1);
		ZsNorm.set(lold, 0, ZsRe.get(lold, 0) * ZsRe.get(lold, 0) + ZsIm.get(lold, 0) * ZsIm.get(lold, 0));
	}
	nChild = MatrixGPU(nChildCPU, 1);

	/*
	ZsNorm.display();*/


	_rhoInv = 1 / _rho;
	resF = MatrixCPU(3, (_iterG / _stepG) + 1);


	//std::cout << " local resolution " << std::endl;
	// local resolution
	tempN2 = MatrixGPU(2 * _nAgent, 1, 0, 1);
	tempB2 = MatrixGPU(2 * _nBus, 1, 0, 1);

	CoresSoloBusAgent = MatrixGPU(_nBus, 1, -1, 1);

	Pn = MatrixGPU(sim.getPn(), 1); // not the real agent
	PSOGPU = MatrixGPU(2 * _nAgent, 1, 0, 1);
	Pmin = MatrixGPU(2 * _nAgent, 1, -1000000, 1); // must not be the real one
	Pmax = MatrixGPU(2 * _nAgent, 1, 1000000, 1); // idem

	// the loss provider
	/*Pmin.set(0, 0, 0);
	Pmax.set(0, 0, 0);
	Pmin.set(_nAgent, 0, 0);
	Pmax.set(_nAgent, 0, 0);*/

	Pbmax = MatrixGPU(2 * _nBus, 1, 0, 1);
	Pbmin = MatrixGPU(2 * _nBus, 1, 0, 1);
	Pb = MatrixGPU(2 * _nBus, 1, 0, 1);
	//Pmin.display();
	//Pn.display();


	Cost1 = MatrixGPU(2 * _nAgent, 1, _rhoSO, 1); //
	Cost1.set(0, 0, 0, true);
	Cost1.set(_nAgent, 0, 0, true);
	Cost2 = MatrixGPU(2 * _nAgent, 1, 0, 1);
	etaSO = MatrixGPU(2 * _nAgent, 1, 0, 1);

	_nAgentByBusCPU.increment(0, 0, -1);
	_nAgentOn0 = _nAgentByBusCPU.get(0, 0);

	////CHECK_LAST_CUDA_ERROR();
	removeLossAgent << <1, 1 >> > (_nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU);

	////CHECK_LAST_CUDA_ERROR();
	ComputePFromAgentToBus();
	////CHECK_LAST_CUDA_ERROR();
	//_nAgentByBus.display();


	//std::cout << " creation " << std::endl;
	X = MatrixGPU(_sizeOPFADMMConsTotal, 1, 0, 1);
	Ypre = MatrixGPU(_sizeOPFADMMConsTotal, 1, 0, 1);
	Y = MatrixGPU(_sizeOPFADMMConsTotal, 1, 0, 1);
	Y.preallocateReduction();

	Mu = MatrixGPU(_sizeOPFADMMConsTotal, 1, 0, 1);

	tempNN = MatrixGPU(_nAgent, _nAgent, 0, 1);
	//tempM1 = new MatrixGPU[_nAgent];
	tempM = MatrixGPU(_sizeOPFADMMConsTotal, 1, 0, 1);

	sizeOPFADMMConsGPU = MatrixGPU(_nBusWLoss, 1, 0, 1);
	sizeOPFADMMConsGPU.preallocateReduction();
	sizeOPFADMMConsGPUBig = MatrixGPU(_sizeOPFADMMConsTotal, 1, 0, 1);

	_indiceBusBegin = MatrixGPU(_nBusWLoss, 1);
	_indiceBusBeginBig = MatrixGPU(_sizeOPFADMMConsTotal, 1, 0, 1);
	_CoresChatBegin = MatrixGPU(_nBusWLoss, 1);


	int debut = 0;
	int debutChat = 0;
	for (int i = 0; i < _nBus; i++) {
		int m = nChildCPU.get(i, 0);
		int nB = _nAgentByBusCPU.get(i, 0);
		_indiceBusBegin.set(i, 0, debut);
		_CoresChatBegin.set(i, 0, debutChat);
		int sizeA = m * 3 + 5 + 2 * nB;
		debut += sizeA;
		debutChat += 4 + 2 * nB;
	}
	_indiceBusBegin.set(_nBus, 0, debut);
	_CoresChatBegin.set(_nBus, 0, debutChat);


	_CoresChatBegin.transferGPU();
	_indiceBusBegin.transferGPU();
	////CHECK_LAST_CUDA_ERROR();
	defineSizeBig << <_nBusWLoss, _blockSize >> > (sizeOPFADMMConsGPUBig._matrixGPU, nChild._matrixGPU, _indiceBusBegin._matrixGPU, sizeOPFADMMConsGPU._matrixGPU, _indiceBusBeginBig._matrixGPU, _nAgentByBus._matrixGPU, losstype, _nBus, _nAgent);
	////CHECK_LAST_CUDA_ERROR();

	_sizeOPFADMMConsMax = sizeOPFADMMConsGPU.max2();
	Hinv = MatrixGPU(_sizeOPFADMMConsTotal, _sizeOPFADMMConsMax, 0, 1);
	Q = MatrixGPU(_sizeOPFADMMConsTotal, 1, 0, 1);

	Childs = MatrixGPU(_nLine, 1);
	PosChild = MatrixGPU(_nBus, 1, -1);
	_indiceChildBegin = MatrixGPU(_nBus, 1);

	Chat = MatrixGPU(_sizeChat, 1, 0, 1);

	MatrixGPU lowerBound(cas.getLowerBound(), 1); //voltage angle, voltage, line...
	MatrixGPU upperBound(cas.getUpperBound(), 1); //voltage angle, voltage, line...

	VoltageLimit = MatrixGPU(2, _nBus, 0, 1); // min, max
	VoltageLimitReal = MatrixGPU(2, _nBus, 0, 1); // min, max


	if (cas.isCurrentLimit()) {
		FluxLimit = MatrixGPU(cas.getCurrentLimit(), 1);
		isCurrentLimited = true;
	}
	else {
		FluxLimit = MatrixGPU(_nLine, 1, 1000, 1); // max
	}
	//FluxLimit.display();
	////CHECK_LAST_CUDA_ERROR();
	initVoltageBound << < _numBlocksB, _blockSize >> > (VoltageLimitReal._matrixGPU, VoltageLimit._matrixGPU, lowerBound._matrixGPU, upperBound._matrixGPU, nChild._matrixGPU, _nBus);
	////CHECK_LAST_CUDA_ERROR();

	int debutChild = 0;
	MatrixCPU nChildTemp(_nBus, 1, 0);
	

	for (int i = 0; i < _nBus; i++) {
		if (i > 0) {
			_indiceChildBegin.set(i - 1, 0, debutChild);

			int Ai = Ancestor.get(i, 0);
			Childs.set(_indiceChildBegin.get(Ai, 0) + nChildTemp.get(Ai, 0), 0, i);
			PosChild.set(i, 0, nChildTemp.get(Ai, 0));
			nChildTemp.increment(Ai, 0, 1);
			debutChild += nChildCPU.get(i - 1, 0);
		}

	}
	Ancestor.transferGPU();
	PosChild.transferGPU();

	debut = 0;
	for (int i = 0; i < _nBus; i++) {
		// (Pi, Qi, li, vi, pn..., qn..., vai, Pci ..., Qci... , lci...) !!!!!
		int m = nChildCPU.get(i, 0);
		int nB = _nAgentByBusCPU.get(i, 0);
		int sizeA = m * 3 + 5 + 2 * nB;
		MatrixCPU A(2 + 1 * (i > 0), sizeA);

		if (i > 0) {
			A.set(2, 0, 2 * ZsRe.get(i - 1, 0));
			A.set(2, 1, 2 * ZsIm.get(i - 1, 0));
			A.set(2, 2, -ZsNorm.get(i - 1, 0));
			A.set(2, 3, -1);
			A.set(2, 4 + 2 * nB, 1); // vai
			A.set(0, 0, -1);
			A.set(1, 1, -1);
		}
		for (int In = 0; In < nB; In++) {
			A.set(0, 4 + In, 1);
			A.set(1, 4 + nB + In, 1);
		}

		for (int j = 0; j < m; j++) {
			int c = Childs.get(_indiceChildBegin.get(i, 0) + j, 0);
			A.set(0, 5 + 2 * nB + j, 1); // Pci
			A.set(1, 5 + 2 * nB + m + j, 1); // Qci
			A.set(0, 5 + 2 * nB + 2 * m + j, -ZsRe.get(c - 1, 0)); // -R l
			A.set(1, 5 + 2 * nB + 2 * m + j, -ZsIm.get(c - 1, 0)); // -X l
		}

		//A.display();

		MatrixCPU temp33(2 + 1 * (i > 0), 2 + 1 * (i > 0));
		MatrixCPU temp3M(2 + 1 * (i > 0), sizeA);
		MatrixCPU tempMM(sizeA, sizeA);
		MatrixCPU tempMMbis(sizeA, sizeA);

		temp33.multiplyTrans(&A, &A);
		temp33.invertGaussJordan(&temp33);
		temp3M.MultiplyMatMat(&temp33, &A);
		tempMM.multiplyTrans(&A, &temp3M, 0);

		tempMMbis.setEyes(-1);
		tempMMbis.add(&tempMM);
		MatrixGPU tempMMGPU = MatrixGPU(tempMMbis, 1);

		Hinv.setBloc(debut, debut + sizeA, 0, sizeA, &tempMMGPU);
		debut += sizeA;
	}

	// bus factice
	int sizeA = 0;
	MatrixGPU A;
	switch (losstype)
	{
	case LossType::POWER:
		sizeA = 2 *_nAgent;
		A = MatrixGPU(sizeA, sizeA);
		A.setEyes(-1);
		A.set(0, 0, 0);
		A.set(_nAgent, _nAgent, 0);
		for (int i = 1; i < _nAgent; i++) {
			A.set(0, i, 1); // sum(p) + Ploss = 0
			A.set(_nAgent, i + _nAgent, 1); // Qloss + sum(q) = 0
		}
		A.transferGPU();
		Hinv.setBloc(debut, debut + sizeA, 0, sizeA, &A);
		/*for (int i = 0; i < _nAgentTrue; i++) {
			A[_nBus].set(0, i, 1); // sum(p) + Ploss = 0
			A[_nBus].set(1, i + _nAgentTrue, 1); // Qloss + sum(q) = 0
		}*/
		break;
	case LossType::CURRENT:
		sizeA = 2 + _nBus;
		A = MatrixGPU(sizeA, sizeA);
		A.setEyes(-1);
		A.set(0, 0, 0);
		A.set(1, 1, 0);
		for (int i = 0; i < _nLine; i++) {
			A.set(0, i + 3, ZsRe.get(i, 0)); // sum(p) + Ploss = 0
			A.set(1, i + 3, ZsIm.get(i, 0)); // Qloss + sum(q) = 0
		}
		A.transferGPU();
		Hinv.setBloc(debut, debut + sizeA, 0, sizeA, &A);
		//A[_nBus].set(0, 0, 1); // ploss
		//A[_nBus].set(1, 1, 1); // qloss
		break;
	}

	Hinv.divide(_rho);

	ZsRe.transferGPU();
	ZsIm.transferGPU();
	_indiceChildBegin.transferGPU();
	Childs.transferGPU();


	//std::cout << " init valeur " << std::endl;
	////CHECK_LAST_CUDA_ERROR();
	initPQAgentV << < _nBus, _blockSizeSmall >> > (X._matrixGPU, _indiceBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, Pn._matrixGPU, _nAgent);
	////CHECK_LAST_CUDA_ERROR();

	////CHECK_LAST_CUDA_ERROR();
	initDFSPQ << <1, _nBus, _nBus* (8 * sizeof(bool) + sizeof(int)) >> > (X._matrixGPU, Pb._matrixGPU, nChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nBus);
	////CHECK_LAST_CUDA_ERROR();

	communicateX << <_nBusWLoss, _blockSize >> > (X._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nAgentByBus._matrixGPU, _CoresBusAgent._matrixGPU, PosAgent._matrixGPU, losstype, _nBus, _nAgent);
	////CHECK_LAST_CUDA_ERROR();

	computeLoss();
	////CHECK_LAST_CUDA_ERROR();

	Y.set(&X);

	// bus factice

	
	

	//std::cout << "updateChat" << std::endl;
	updateChat();
	
	cudaDeviceSynchronize();
	//CHECK_LAST_CUDA_ERROR();
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	//std::cout << "---------------------------------------------------------------------------------------" << std::endl;
}


void OPFADMMConsGPU::updateConsensus(MatrixGPU* Pmarket)
{
	
	//Pmarket->display(true);
	CHECK_LAST_CUDA_ERROR();
	updateConsensusGPU << <_numBlocksN, _blockSize >> > (Cost2._matrixGPU, etaSO._matrixGPU, Pn._matrixGPU, Pmarket->_matrixGPU, _rhoSO, _nAgent);
	
	CHECK_LAST_CUDA_ERROR();
	/*std::cout << "Cost 2" << std::endl;
	Cost2.display(true);
	std::cout << "*********" << std::endl;*/

}

void OPFADMMConsGPU::updateConsensus(MatrixCPU* Pmarket) {

	PSOGPU.set(Pmarket);
	updateConsensus(&PSOGPU);

}



void OPFADMMConsGPU::updateGlobalProb() {
	
	
	Ypre.swap(&Y);
	int numBlock = _sizeOPFADMMConsTotal;
	switch (_blockSize) {
	case 512:
		updateY<512> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMConsGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFADMMConsMax);
		break;
	case 256:
		updateY<256> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMConsGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFADMMConsMax);
		break;
	case 128:
		updateY<128> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMConsGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFADMMConsMax);
		break;
	case 64:
		updateY< 64> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMConsGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFADMMConsMax);
		break;
	case 32:
		updateY< 32> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMConsGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFADMMConsMax);
		break;
	case 16:
		updateY< 16> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMConsGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFADMMConsMax);
		break;
	case  8:
		updateY<  8> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMConsGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFADMMConsMax);
		break;
	case  4:
		updateY<  4> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMConsGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFADMMConsMax);
		break;
	case  2:
		updateY<  2> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMConsGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFADMMConsMax);
		break;
	case  1:
		updateY<  1> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMConsGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFADMMConsMax);
		break;
	}

	Y.set(3, 0, 1, 1);
	Y.set(4 + 2 * _nAgentOn0, 0, 1, 1);

	// communication of y, mu

}

void OPFADMMConsGPU::updateX()
{
	/*double x1, x2, x3, x4, c1, c2, c3, c4, lambdaLo, lambdaUp, delta, x3min, x3max, x4max, gamma, k2;
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

			//x3 = x3max
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
*/
	
}

void OPFADMMConsGPU::updateXWOCurrent()
{
	updateXOPFADMMCons << <_nBus, _blockSizeSmall >> > (X._matrixGPU, Pn._matrixGPU, Chat._matrixGPU, VoltageLimit._matrixGPU, _nAgentByBus._matrixGPU, nChild._matrixGPU, _indiceBusBegin._matrixGPU, _CoresChatBegin._matrixGPU,
		_CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, Cost1._matrixGPU, Cost2._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _rho, losstype, _nBus, _nAgent, Lagrange);

}

void OPFADMMConsGPU::updateXWOCurrentOnCPU()
{
	X.transferCPU();
	Chat.transferCPU();
	double x1, x2, x3, x4, c1, c2, c3, c4, lambdaLo, lambdaUp, x3min, x3max, gamma, k2; // double peut �tre necessaire
	double c1122; // c3 : voltage -> indice + 3, c4 : current -> indice + 2;
	
	


	for (int bus = 1; bus < _nBus; bus++) {

		int typeSol = 0;
		int BestRoot = 0;
		double bestGamma = -1;
		double p = 0;

		int nRoot = 0;

		int begining = indiceBusBeginCPU.get(bus, 0);
		int nC = nChildCPU.get(bus,0);
		int beginChat = CoresChatBeginCPU.get(bus, 0);
		bool goodSol = false;
		k2 = sqrt(2.0 / (nC + 1));

		c1 = -2 * Chat.get(beginChat, 0);
		c2 = -2 * Chat.get(beginChat + 1, 0);
		c4 = -2 * Chat.get(beginChat + 2, 0);
		c3 = -2 * Chat.get(beginChat + 3, 0) / k2;

		c1122 = c1 * c1 + c2 * c2;
		x3min = VoltageLimitCPU.get(0, bus);
		x3max = VoltageLimitCPU.get(1, bus);

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

		if (gamma >= 0) {
			// the solution is good !
			goodSol = true;
		}
		else {
			if (c1122 == 0) {
				x4 = 0;
				goodSol = true;
			}
			if (gamma > bestGamma) {
				typeSol = 1;
				bestGamma = gamma;
			}
		}
		//}
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

				if (gamma >= 0 && lambdaUp >= 0) {
					goodSol = true;
					break;
				}
				if (gamma > bestGamma && lambdaUp > bestGamma) {
					typeSol = 2;
					bestGamma = MYMIN(gamma, lambdaUp);
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
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;
				lambdaLo = (2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));

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

			nRoot = resvolveRealPolynome4without2term(root4, coefPoly3);

			for (int n = 0; n < nRoot; n++) {
				p = root4[n];
				x3 = -(c1122 * p + 2 * c3) / (2 * (c1122 * p * p + 2));
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;

				if (gamma >= 0 && x3 <= x3max && x3 >= x3min) {
					// the solution is good !
					goodSol = true;
					break;
				}
				if (gamma > bestGamma && (x3max - x3) > bestGamma && (x3 - x3min) > bestGamma) {
					typeSol = 4;
					bestGamma = MYMIN(MYMIN(gamma, (x3max - x3)), (x3 - x3min));
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

		X.set(begining, 0, x1);
		X.set(begining + 1, 0, x2);
		X.set(begining + 2, 0, x4);
		X.set(begining + 3, 0, x3 *k2);
	}



	X.transferGPU();
	Chat.transferGPU();


	updateXPnOPFADMMCons << <_nBusWLoss, _blockSizeSmall >> > (X._matrixGPU, Pn._matrixGPU, Chat._matrixGPU, _nAgentByBus._matrixGPU, _indiceBusBegin._matrixGPU, _CoresChatBegin._matrixGPU,
		_CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, Cost1._matrixGPU, Cost2._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _rho, losstype, _nBus, _nAgent);


}


void OPFADMMConsGPU::updateXWOCurrentOnCPUBis()
{
	X.transferCPU();
	Chat.transferCPU();
	double x1, x2, x3, x4, c1, c2, c3, c4, lambdaLo, lambdaUp, x3min, x3max, gamma, k2; // double peut �tre necessaire
	double c1122; // c3 : voltage -> indice + 3, c4 : current -> indice + 2;

	// résolution du polynome 4 sur GPU
	int nPoly = _nBus - 1;
	MatrixGPUD coefsGPU(4, nPoly);
	MatrixGPUD rootsGPU(4, nPoly, 0, 1);
	MatrixGPUD nRootGPU(nPoly, 1, 0, 1);

	for (int bus = 1; bus < _nBus; bus++) {
		int nC = nChildCPU.get(bus, 0);
		int beginChat = CoresChatBeginCPU.get(bus, 0);
		k2 = sqrt(2.0 / (nC + 1));
		c1 = -2 * Chat.get(beginChat, 0);
		c2 = -2 * Chat.get(beginChat + 1, 0);
		c4 = -2 * Chat.get(beginChat + 2, 0);
		c3 = -2 * Chat.get(beginChat + 3, 0) / k2;
		c1122 = c1 * c1 + c2 * c2;

		double aInv = k2 * k2 / (c1122 * c1122);
		double b = c1122 / k2 * (2 * c3 / k2 - c4);
		double d = (c3 - 2 * c4 / k2);
		double e = -1;

		b *= aInv;
		d *= aInv;
		e *= aInv;

		coefsGPU.set(0, bus - 1, b);
		coefsGPU.set(1, bus - 1, 0);
		coefsGPU.set(2, bus - 1, d);
		coefsGPU.set(3, bus - 1, e);
	}
	coefsGPU.transferGPU();
	resolveSeveralRealPolynome4WO2termGPULagrange << <_numBlocksB, _blockSizeSmall >> > (nRootGPU._matrixGPU, rootsGPU._matrixGPU, coefsGPU._matrixGPU, nPoly);
	nRootGPU.transferCPU();
	rootsGPU.transferCPU();


	for (int bus = 1; bus < _nBus; bus++) {

		int typeSol = 0;
		int BestRoot = 0;
		double bestGamma = -1;
		double p = 0;

		int nRoot = 0;

		int begining = indiceBusBeginCPU.get(bus, 0);
		int nC = nChildCPU.get(bus, 0);
		int beginChat = CoresChatBeginCPU.get(bus, 0);
		bool goodSol = false;
		k2 = sqrt(2.0 / (nC + 1));

		c1 = -2 * Chat.get(beginChat, 0);
		c2 = -2 * Chat.get(beginChat + 1, 0);
		c4 = -2 * Chat.get(beginChat + 2, 0);
		c3 = -2 * Chat.get(beginChat + 3, 0) / k2;

		c1122 = c1 * c1 + c2 * c2;
		x3min = VoltageLimitCPU.get(0, bus);
		x3max = VoltageLimitCPU.get(1, bus);

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

		if (gamma >= 0) {
			// the solution is good !
			goodSol = true;
		}
		else {
			if (c1122 == 0) {
				x4 = 0;
				goodSol = true;
			}
			if (gamma > bestGamma) {
				typeSol = 1;
				bestGamma = gamma;
			}
		}
		//}
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

				if (gamma >= 0 && lambdaUp >= 0) {
					goodSol = true;
					break;
				}
				if (gamma > bestGamma && lambdaUp > bestGamma) {
					typeSol = 2;
					bestGamma = MYMIN(gamma, lambdaUp);
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
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;
				lambdaLo = (2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));

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

			/*coefPoly3[0] = c1122 / k2 * (2 * c3 / k2 - c4);
			coefPoly3[1] = (c3 - 2 * c4 / k2);
			coefPoly3[2] = -1;
			coefPoly3[0] = coefPoly3[0] * k2 * k2 / (c1122 * c1122);
			coefPoly3[1] = coefPoly3[1] * k2 * k2 / (c1122 * c1122);
			coefPoly3[2] = coefPoly3[2] * k2 * k2 / (c1122 * c1122);

			nRoot = resvolveRealPolynome4without2term(root4, coefPoly3);*/

			nRoot = nRootGPU.get(bus - 1, 0);

			for (int n = 0; n < nRoot; n++) {
				root4[n] = rootsGPU.get(n, bus - 1);
				p = root4[n];
				x3 = -(c1122 * p + 2 * c3) / (2 * (c1122 * p * p + 2));
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;

				if (gamma >= 0 && x3 <= x3max && x3 >= x3min) {
					// the solution is good !
					goodSol = true;
					break;
				}
				if (gamma > bestGamma && (x3max - x3) > bestGamma && (x3 - x3min) > bestGamma) {
					typeSol = 4;
					bestGamma = MYMIN(MYMIN(gamma, (x3max - x3)), (x3 - x3min));
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

		X.set(begining, 0, x1);
		X.set(begining + 1, 0, x2);
		X.set(begining + 2, 0, x4);
		X.set(begining + 3, 0, x3 * k2);
	}



	X.transferGPU();
	Chat.transferGPU();


	updateXPnOPFADMMCons << <_nBusWLoss, _blockSizeSmall >> > (X._matrixGPU, Pn._matrixGPU, Chat._matrixGPU, _nAgentByBus._matrixGPU, _indiceBusBegin._matrixGPU, _CoresChatBegin._matrixGPU,
		_CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, Cost1._matrixGPU, Cost2._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _rho, losstype, _nBus, _nAgent);


}

void OPFADMMConsGPU::updateXWOCurrentOnCPUBis(bool first)
{
	X.transferCPU();
	Chat.transferCPU();
	double x1, x2, x3, x4, c1, c2, c3, c4, lambdaLo, lambdaUp, x3min, x3max, gamma, k2; // double peut �tre necessaire
	double c1122; // c3 : voltage -> indice + 3, c4 : current -> indice + 2;

	// résolution du polynome 4 sur GPU
	int nPoly = _nBus - 1;
	MatrixGPUD coefsGPU(4, nPoly);
	MatrixGPUD rootsGPU(3, nPoly, 0, 1);
	MatrixGPUD nRootGPU(nPoly, 1, 0, 1);

	for (int bus = 1; bus < _nBus; bus++) {
		int nC = nChildCPU.get(bus, 0);
		int beginChat = CoresChatBeginCPU.get(bus, 0);
		k2 = sqrt(2.0 / (nC + 1));
		c1 = -2 * Chat.get(beginChat, 0);
		c2 = -2 * Chat.get(beginChat + 1, 0);
		c4 = -2 * Chat.get(beginChat + 2, 0);
		c3 = -2 * Chat.get(beginChat + 3, 0) / k2;
		c1122 = c1 * c1 + c2 * c2;
		x3min = VoltageLimitCPU.get(0, bus);
		x3max = VoltageLimitCPU.get(1, bus);

		if (first){
			x3 = x3max;
		}
		else {
			x3 = x3min;
		}
		
		double aInv = k2 * k2 / (4 * c1122);
		double p = 2 * (c4 / (k2 * x3) + 1);
		double q = 1 / x3;

		

		p *= aInv;
		q *= aInv;
	
		coefsGPU.set(0, bus - 1, 1);
		coefsGPU.set(1, bus - 1, 0);
		coefsGPU.set(2, bus - 1, p);
		coefsGPU.set(3, bus - 1, q);
	}
	coefsGPU.transferGPU();
	resolveSeveralRealPolynome3termGPU << <_numBlocksB, _blockSizeSmall >> > (nRootGPU._matrixGPU, rootsGPU._matrixGPU, coefsGPU._matrixGPU, nPoly);
	nRootGPU.transferCPU();
	rootsGPU.transferCPU();


	for (int bus = 1; bus < _nBus; bus++) {

		int typeSol = 0;
		int BestRoot = 0;
		double bestGamma = -1;
		double p = 0;

		int nRoot = 0;

		int begining = indiceBusBeginCPU.get(bus, 0);
		int nC = nChildCPU.get(bus, 0);
		int beginChat = CoresChatBeginCPU.get(bus, 0);
		bool goodSol = false;
		k2 = sqrt(2.0 / (nC + 1));

		c1 = -2 * Chat.get(beginChat, 0);
		c2 = -2 * Chat.get(beginChat + 1, 0);
		c4 = -2 * Chat.get(beginChat + 2, 0);
		c3 = -2 * Chat.get(beginChat + 3, 0) / k2;

		c1122 = c1 * c1 + c2 * c2;
		x3min = VoltageLimitCPU.get(0, bus);
		x3max = VoltageLimitCPU.get(1, bus);

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

		if (gamma >= 0) {
			// the solution is good !
			goodSol = true;
		}
		else {
			if (c1122 == 0) {
				x4 = 0;
				goodSol = true;
			}
			if (gamma > bestGamma) {
				typeSol = 1;
				bestGamma = gamma;
			}
		}
		//}
		if (!goodSol) {
			x3 = x3max;

			if (first) {
				nRoot = nRootGPU.get(bus - 1, 0);
				for (int n = 0; n < nRoot; n++) {
					root2[n] = rootsGPU.get(n, bus - 1);
				}
			}
			else {
				coefPoly2[0] = 2 * (c4 / (k2 * x3) + 1);
				coefPoly2[1] = 1 / x3;
				coefPoly2[0] = coefPoly2[0] * k2 * k2 / (4 * c1122);
				coefPoly2[1] = coefPoly2[1] * k2 * k2 / (4 * c1122);
				nRoot = resolveRealPolynome3without2term(root2, coefPoly2);
			}
			

			for (int n = 0; n < nRoot; n++) {
				p = root2[n];

				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;
				lambdaUp = -(2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));

				if (gamma >= 0 && lambdaUp >= 0) {
					goodSol = true;
					break;
				}
				if (gamma > bestGamma && lambdaUp > bestGamma) {
					typeSol = 2;
					bestGamma = MYMIN(gamma, lambdaUp);
					BestRoot = n;
				}

			}
		}
		// case x3 = x3min lambdaUp = 0
		if (!goodSol) {
			x3 = x3min;

			if (!first) {
				nRoot = nRootGPU.get(bus - 1, 0);
				for (int n = 0; n < nRoot; n++) {
					root3[n] = rootsGPU.get(n, bus - 1);
				}
			}
			else {
				coefPoly2[0] = 2 * (c4 / (k2 * x3) + 1);
				coefPoly2[1] = 1 / x3;
				coefPoly2[0] = coefPoly2[0] * k2 * k2 / (4 * c1122);
				coefPoly2[1] = coefPoly2[1] * k2 * k2 / (4 * c1122);
				nRoot = resolveRealPolynome3without2term(root2, coefPoly2);
			}

			for (int n = 0; n < nRoot; n++) {
				p = root3[n];
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;
				lambdaLo = (2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));

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

			nRoot = resvolveRealPolynome4without2term(root4, coefPoly3, Lagrange);/**/

			

			for (int n = 0; n < nRoot; n++) {
				p = root4[n];
				x3 = -(c1122 * p + 2 * c3) / (2 * (c1122 * p * p + 2));
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;

				if (gamma >= 0 && x3 <= x3max && x3 >= x3min) {
					// the solution is good !
					goodSol = true;
					break;
				}
				if (gamma > bestGamma && (x3max - x3) > bestGamma && (x3 - x3min) > bestGamma) {
					typeSol = 4;
					bestGamma = MYMIN(MYMIN(gamma, (x3max - x3)), (x3 - x3min));
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

		X.set(begining, 0, x1);
		X.set(begining + 1, 0, x2);
		X.set(begining + 2, 0, x4);
		X.set(begining + 3, 0, x3 * k2);
	}



	X.transferGPU();
	Chat.transferGPU();


	updateXPnOPFADMMCons << <_nBusWLoss, _blockSizeSmall >> > (X._matrixGPU, Pn._matrixGPU, Chat._matrixGPU, _nAgentByBus._matrixGPU, _indiceBusBegin._matrixGPU, _CoresChatBegin._matrixGPU,
		_CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, Cost1._matrixGPU, Cost2._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _rho, losstype, _nBus, _nAgent);


}


void OPFADMMConsGPU::updateMu()
{
	updateMUGPU << <_numBlocksH, _blockSize >> > (Mu._matrixGPU, Y._matrixGPU, X._matrixGPU, _rho, _sizeOPFADMMConsTotal);

}



float OPFADMMConsGPU::getPLoss()
{
	_Ploss = Y.get(_debutloss, 0, false);

	return _Ploss;
}

float OPFADMMConsGPU::getQLoss()
{
	int indice = 1;
	if (losstype == LossType::POWER) {
		indice = _nAgent;
	}

	_Qloss = Y.get(_debutloss + indice, 0, false);

	return _Qloss;
}

void OPFADMMConsGPU::computeLoss()
{

	int numBlock = 1;
	switch (_blockSize) {
	case 512:
		ComputeLoss<512> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgent, _nBus);
		break;
	case 256:
		ComputeLoss<256> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgent, _nBus);
		break;
	case 128:
		ComputeLoss<128> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgent, _nBus);
		break;
	case 64:
		ComputeLoss< 64> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgent, _nBus);
		break;
	case 32:
		ComputeLoss< 32> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgent, _nBus);
		break;
	case 16:
		ComputeLoss< 16> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgent, _nBus);
		break;
	case  8:
		ComputeLoss<  8> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgent, _nBus);
		break;
	case  4:
		ComputeLoss<  4> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgent, _nBus);
		break;
	case  2:
		ComputeLoss<  2> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgent, _nBus);
		break;
	case  1:
		ComputeLoss<  1> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgent, _nBus);
		break;
	}
}

void OPFADMMConsGPU::updateChat()
{
	int numBlock = _nBusWLoss;
	switch (_blockSizeSmall) {
	case 512:
		updateChatGPU4<512> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case 256:
		updateChatGPU4<256> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case 128:
		updateChatGPU4<128> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case 64:
		updateChatGPU4< 64> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case 32:
		updateChatGPU4< 32> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case 16:
		updateChatGPU4< 16> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case  8:
		updateChatGPU4<  8> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case  4:
		updateChatGPU4<  4> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case  2:
		updateChatGPU4<  2> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case  1:
		updateChatGPU4<  1> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	}

}

void OPFADMMConsGPU::CommunicationX()
{
/**/ // X = { Pi, Qi, vi, li, vAi, (pn, qn) (Pci, Qci, lci) for all child Ci }
	
	communicateX << <_nBusWLoss, _blockSize >> > (X._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nAgentByBus._matrixGPU, _CoresBusAgent._matrixGPU, PosAgent._matrixGPU, losstype, _nBus, _nAgent);

	
	


	
	//Y      (Pi, Qi, vi, li, pi, qi, vai, (pn, qn) Pji, Qji, lji)
	// Q udate in argmin 0.5yHy + Qy

	updateQ << <_numBlocksH, _blockSize >> > (Q._matrixGPU, X._matrixGPU, Mu._matrixGPU, _rho, _sizeOPFADMMConsTotal);


}


float OPFADMMConsGPU::updateRes(int indice) 
{
	float resS = Y.max2(&Ypre);
	float resR = Y.max2(&X);
	float resV = 0;
	
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
		 
			Hinv.divide(_tau);
			//std::cout << _iterGlobal << "rho augmente " << _rho << std::endl;
		}
		else if (resS > _mu * resR) {// rho = rho / tau_inc;
		
			_rho = _rho / _tau;
		
			Hinv.multiply(_tau);
			//std::cout << _tau << " " << _mu << std::endl;
			//std::cout << _iterGlobal << "rho diminue " << _rho << std::endl;
		}
	}
	


	return MYMAX(MYMAX(resV, oldrho * resS), resR);
}

float OPFADMMConsGPU::updateResRhoFixe(int indice)
{
	float resS = _rho * Y.max2(&Ypre);
	float resR = Y.max2(&X);
	float resV = 0;

	resF.set(0, indice, resR);
	resF.set(1, indice, resS);
	resF.set(2, indice, resV);

	return MYMAX(MYMAX(resV, resS), resR);
}

int OPFADMMConsGPU::feasiblePoint()
{
	bool mustTrans = false;
	if (X.getPos()) {
		X.transferCPU();
		_indiceBusBegin.transferCPU();
		mustTrans = true;
	}
	// X  (Pi, Qi, li, vi, pn..., qn..., vai, Pci ..., Qci... , lci...) !!!!!

	MatrixCPU test(_nBus, 1, -1);
	int counter = 0;
	for (int bus = 0; bus < _nBus; bus++) {
		int begin = _indiceBusBegin.get(bus, 0);
		float Si = X.get(begin, 0) * X.get(begin, 0) + X.get(begin + 1, 0) * X.get(begin + 1, 0);
		float li = X.get(begin + 2, 0);
		float vi = X.get(begin + 3, 0);
		float err = Si - li * vi;
		test.set(bus, 0, err);
		if (abs(err) > 0.0001) {
			counter++;
		}
	}
	//std::cout << " erreur sur la relaXation " << test.max2() << " " << counter << std::endl;
	//test.display();

	if (mustTrans) {
		X.transferGPU();
		_indiceBusBegin.transferGPU();
	}
	resF.set(2, (_iterGlobal - 1) / _stepG, test.max2());
	return counter;
}

void OPFADMMConsGPU::ComputePFromAgentToBus()
{
	int numBlock = _nBus;
	switch (_blockSizeSmall) {
	case 512:
		ComputePFromAgentToBusGPU<512> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case 256:
		ComputePFromAgentToBusGPU<256> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case 128:
		ComputePFromAgentToBusGPU<128> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case 64:
		ComputePFromAgentToBusGPU< 64> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case 32:
		ComputePFromAgentToBusGPU< 32> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case 16:
		ComputePFromAgentToBusGPU< 16> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case  8:
		ComputePFromAgentToBusGPU<  8> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case  4:
		ComputePFromAgentToBusGPU<  4> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case  2:
		ComputePFromAgentToBusGPU<  2> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case  1:
		ComputePFromAgentToBusGPU<  1> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	}
}

MatrixCPU OPFADMMConsGPU::getPb(){
	MatrixCPU PbCPU;
	Pb.toMatCPU(PbCPU);
	return PbCPU;
}
MatrixCPU OPFADMMConsGPU::getPhi(){
	bool transferToDo = false;
	if(Y.getPos()){
		Y.transferCPU();
		_indiceBusBegin.transferCPU();
		transferToDo = true;
	}
	MatrixCPU Phi(2*_nLine, 1);
	
	for (int i = 0; i <_nLine; i++)
	{
		Phi.set(i,0, Y.get(_indiceBusBegin.get(i + 1,0) + 0, 0));
		Phi.set(i + _nLine,0, Y.get(_indiceBusBegin.get(i + 1,0) + 1, 0));
	}
	if(transferToDo){
		Y.transferGPU();
		_indiceBusBegin.transferGPU();
	}
	return Phi;
}
MatrixCPU OPFADMMConsGPU::getE(){
	bool transferToDo = false;
	if(Y.getPos()){
		Y.transferCPU();
		_indiceBusBegin.transferCPU();
		transferToDo = true;
	}
	MatrixCPU E(2*_nBus, 1);
	
	for (int i = 0; i <_nBus; i++)
	{
		E.set(i,0, Y.get(_indiceBusBegin.get(i, 0) + 2, 0));
		E.set(i + _nLine,0, Y.get(_indiceBusBegin.get(i, 0) + 3, 0));
	}
	if(transferToDo){
		Y.transferGPU();
		_indiceBusBegin.transferGPU();
	}
	return E;
}


void OPFADMMConsGPU::display() {

	
	X.transferCPU();
	Y.transferCPU();
	Mu.transferCPU();
	Pn.transferCPU();
	ZsRe.transferCPU();
	_indiceBusBegin.transferCPU();
	VoltageLimitReal.transferCPU();
	Pmin.transferCPU();
	Pmax.transferCPU();
	Pbmax.transferCPU();
	Pbmin.transferCPU();
	Pb.transferCPU();
	_CoresAgentBusBegin.transferCPU();
	_CoresAgentBus.transferCPU();
	_CoresBusAgent.transferCPU();
	Cost1.transferCPU();
	Cost2.transferCPU();
	


	if (_iterGlobal == 0) {
		std::cout << "algorithm not launch" << std::endl;
	}
	else if (_iterGlobal < _iterG) {
		std::cout << "method " << _name << " converged in " << _iterGlobal << " iterations." << std::endl;
		std::cout << "Converged in " << (float)timeOPF / CLOCKS_PER_SEC << " seconds" << std::endl;

	}
	else {
		std::cout << "method " << _name << " not converged in " << _iterGlobal << " iterations." << std::endl;
		std::cout << "time taken " << (float)timeOPF / CLOCKS_PER_SEC << " seconds" << std::endl;
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


	for (int b = 0; b < _nBus; b++) {
		int begining = _indiceBusBegin.get(b, 0);
		std::cout << std::setw(5) << b << "|" << std::setw(12) << sqrt(X.get(begining + 3, 0)) << " |" << std::setw(16)
			<< Pb.get(b, 0) << "|" << std::setw(16) << Pb.get(b, 0)
			<< "|" << std::setw(16) << Mu.get(begining + 3, 0) << "|" << std::setw(16)
			<< Mu.get(begining, 0) << "|" << std::setw(16) << Mu.get(begining + 1, 0) << "|" << std::endl;

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
		int begining = _indiceBusBegin.get(b, 0);
		std::cout << std::setw(6) << l << "|" << std::setw(12) << CoresLineBusCPU.get(l, 0) << " |" << std::setw(12)
			<< CoresLineBusCPU.get(l, 1) << "|" << std::setw(16) << X.get(begining + 0, 0)
			<< "|" << std::setw(16) << X.get(begining + 1, 0) << "|" << std::setw(16)
			<< X.get(begining + 2, 0) << "|" << std::setw(19) << X.get(begining + 2, 0) * ZsRe.get(l, 0) << "|" << std::endl;
	}
	std::cout << std::endl << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "     Constraints                                                                                        |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Bus | Voltage | Voltage | Voltage |        Power Injection          |          Power Injection         |" << std::endl;
	std::cout << "  #  | Mag(pu) | MIN(pu) |  MYMAX(pu)|  P (pu) | Pmin (pu) | Pmax (pu) |  Q (pu)  | Qmin (pu) | Qmax (pu) |" << std::endl;
	std::cout << "-----|---------|---------|---------|---------|-----------|-----------|----------|-----------|-----------|" << std::endl;


	for (int b = 0; b < _nBus; b++) {
		int begining = _indiceBusBegin.get(b, 0);
		std::cout << std::setw(5) << b << "|" << std::setw(8) << sqrt(Y.get(begining + 3, 0)) << " |" << std::setw(9)
			<< VoltageLimitReal.get(0, b) << "|" << std::setw(9) << VoltageLimitReal.get(1, b)
			<< "|" << std::setw(9) << Pb.get(b, 0) << "|" << std::setw(11)
			<< Pbmin.get(b, 0) << "|" << std::setw(11) << Pbmax.get(b, 0) << "|" << std::setw(10) << Pb.get(b + _nBus, 0)
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
		std::cout << std::setw(7) << n << "|" << std::setw(7) << b << "|" << std::setw(8) << Cost1.get(n, 0) << " |" << std::setw(9)
			<< Cost2.get(n, 0) << "|" << std::setw(9) << Pn.get(n, 0) << "|" << std::setw(11)
			<< Pmin.get(n, 0) << "|" << std::setw(11) << Pmax.get(n, 0) << "|" << std::setw(10) << Pn.get(n + _nAgent, 0)
			<< "|" << std::setw(11) << Pmin.get(n + _nAgent, 0) << "|" << std::setw(11) << Pmax.get(n + _nAgent, 0) << "|" << std::endl;
	}


	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "                      END PRINT                                                                         |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;

}

template <unsigned int _blockSizeSmall>
__global__ void updateChatGPU4(float* Chat, float* Y, float* MU, float* nChild, float* Ancestor, float* posChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, float* CoresChatBegin, float* indiceAgentBegin, float* CoresAgentBus,  float* nAgentByBus, float _rho, int losstype, int nBus) {

	int bus = blockIdx.x;
	int index = threadIdx.x;
	int step = blockDim.x;


	int beginChat = CoresChatBegin[bus];

	__shared__ float shArr[_blockSizeSmall]; // c'est grand pour pas grand chose...

	int beginBus = indiceBusBegin[bus];
	int beginChild = (bus < (nBus - 1)) ? indiceChildBegin[bus] : 0;
	int childCount = nChild[bus];
	int AncestorIndex = Ancestor[bus];
	int nAgent = nAgentByBus[bus];
	int c = posChild[bus];
	int beginAgent = indiceAgentBegin[bus];
	int beginLoss = indiceBusBegin[nBus];
	
	float var = 0;
	int borne = 4 + 2 * nAgent;
	int divideVar = 2 + ((losstype == 1) && (index == 2)) * 1 - ((losstype == 1) && (index > 3)) * 1;

	if (bus < nBus) {
		if (index < borne) {
			//float Phat, Qhat, lhat, vihat, pnhat..., qnhat...;
			var = Y[beginBus + index] / (divideVar) - MU[beginBus + index] / (divideVar * _rho);
			
			if (index < 3 && bus > 0) {
				int childCountAi = nChild[AncestorIndex];
				int nAgentAi = nAgentByBus[AncestorIndex];
				int indiceAncBus = indiceBusBegin[AncestorIndex] + 4 + 2 * nAgentAi + 1 + childCountAi * index + c;
				//var = indiceAncBus;
				var += Y[indiceAncBus] / divideVar - MU[indiceAncBus] / (divideVar * _rho);
				var += ((index == 2) && (losstype == 1)) ? (Y[beginLoss + 2 + bus] / divideVar - MU[beginLoss + 2 + bus] / (divideVar * _rho)) : 0;

			}
			if ((index > 3) && (losstype==0)) {
				int offset = index >= 4 + nAgent ? 4 + nAgent : 4;
				int n = CoresAgentBus[ beginAgent + index - offset];
				var += Y[beginLoss + n] / (divideVar) - MU[beginLoss + n] / (divideVar * _rho);
			}
		}
		
		float vhat = 0;
		float muhat = 0;
		for (int i = index; i < childCount; i += step) {
			int Bus2 = Childs[beginChild + i];
			int indiceBusChild = indiceBusBegin[Bus2];
			int nAgent2 = nAgentByBus[Bus2];
			muhat += MU[indiceBusChild + 4 + 2 * nAgent2]; // pas du tout coalescent
			vhat += Y[indiceBusChild + 4 + 2 * nAgent2]; // pas du tout coalescent
		}
		shArr[index] = vhat / (childCount + 1) - muhat / (_rho * (childCount + 1));
		__syncthreads();
		for (int size = _blockSizeSmall / 2; size > 0; size /= 2) { //uniform
			if (index < size) {
				shArr[index] += shArr[index + size];
			}
			__syncthreads();
		}

		if (index < borne) {
			if (index == 3) {
				var = shArr[0] + Y[beginBus + 3] / (childCount + 1) - MU[beginBus + 3] / (_rho * (childCount + 1)); //shArr[0];
			}
			Chat[beginChat + index] = var; // coalescent  !!!!
		}
	}
	else { //bus des pertes
		if (index == 0) {
			float phat;
			float qhat;
			if (losstype == 0) {
				phat = Y[beginBus] - MU[beginBus] / _rho;
				qhat = Y[beginBus + nAgent] - MU[beginBus + nAgent] / _rho;
			}
			else {
				phat = Y[beginBus] - MU[beginBus] / _rho;
				qhat = Y[beginBus + 1] - MU[beginBus + 1] / _rho;
			}
			Chat[beginChat] = phat;
			Chat[beginChat + 1] = qhat;
		}
	}

	
}




/*
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

*/

__global__ void updateConsensusGPU(float* Cost2, float* etaSO, float* Pn, float* Pmarket, float _rhoSO, int nAgent) {

	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int step = blockDim.x * gridDim.x;

	for (int agent = index; agent < nAgent; agent += step) {
		float eta = etaSO[agent] + 0.5 * (Pn[agent] - Pmarket[agent]);
		etaSO[agent] = eta;
		Cost2[agent] = _rhoSO * (eta - 0.5 * (Pn[agent] + Pmarket[agent]));
		eta = etaSO[agent + nAgent] + 0.5 * (Pn[agent + nAgent] - Pmarket[agent + nAgent]);
		etaSO[agent + nAgent] = eta;
		Cost2[agent + nAgent] = _rhoSO * (eta - 0.5 * (Pn[agent + nAgent] + Pmarket[agent + nAgent]));
	}


}