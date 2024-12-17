#include "../head/OPFADMMGPU2.cuh"
#define MAX(X, Y) X * (X >= Y) + Y * (Y > X)
#define NMAXAGENTPERTHREAD 5

OPFADMMGPU2::OPFADMMGPU2() : MethodOPFGPU()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " OPFADMMGPU2 Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	timePerBlock = MatrixCPU(1, 12, 0); // Fb0, Fb11abcd, FB12, Fb2, Fb3, Fb4, Fb5,FB6, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 12, 0); //nb de fois utilisé pendant la simu
}

OPFADMMGPU2::OPFADMMGPU2(float rho) : MethodOPFGPU()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default OPFADMMGPU2 Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
	timePerBlock = MatrixCPU(1, 12, 0); // Fb0, Fb11, FB12, Fb2, Fb3, Fb4, Fb5,FB6, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 12, 0); //nb de fois utilisé pendant la simu
}

OPFADMMGPU2::~OPFADMMGPU2()
{
}
void OPFADMMGPU2::setParam(float rho)
{
	_rho = rho;
}

bool OPFADMMGPU2::chekcase()
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

void OPFADMMGPU2::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
		cudaDeviceSynchronize();
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
	float rhoInit = sim.getRho();
	
	
	float fc = 0;
	float resG = 2 * epsG;
	float resL = 2 * epsL;
	_iterGlobal = 0;
	
	/*Chat.display(true);
	Bpt2.display(true);
	CoresSoloBusAgent.display();
	Cost1.display();
	Cost2.display();
	Pmin.display();
	Pmax.display();
	std::cout << "------" << std::endl;*/
	
	while ((_iterGlobal < _iterG) && (resG>epsG)) {
		
		
		/*std::cout << "--------" << std::endl;
		
		std::cout << " X " << std::endl;
		X.display(true);
		std::cout << " Q " << std::endl;
		Q.display(true);
		std::cout << " Y " << std::endl;
		Y.display(true);
		std::cout << " Mu " << std::endl;
		Mu.display(true);
		std::cout << " Chat " << std::endl;
		Chat.display(true);
		*/
		

#ifdef INSTRUMENTATION
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		// on pourrait faire en 2 appels pour mieu parraléliser
		// on pourrait serialiser la gestion des agents
		updateXOPFADMM <<<_nBus, _blockSizeSmall >> > (X._matrixGPU, Chat._matrixGPU, VoltageLimit._matrixGPU, _nAgentByBus._matrixGPU, nChild._matrixGPU, _indiceBusBegin._matrixGPU, _CoresChatBegin._matrixGPU, 
			_CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, Cost1._matrixGPU, Cost2._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _rho,  _nBus, _nAgent, Lagrange);

#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 1, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

		CommunicationX();

#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 6, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		
		updateGlobalProb();
		updateMu();
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 7, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		
		updateChat();
		
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
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
			cudaDeviceSynchronize();
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 9, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		}
		//std::cout << iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;

		_iterGlobal++;
	}
	//std::cout << iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, (iterGlobal - 1) / stepG) << " " << resF.get(1, (iterGlobal - 1) / stepG) << " " << resF.get(2, (iterGlobal - 1) / stepG) << std::endl;


#ifdef INSTRUMENTATION	
	occurencePerBlock.increment(0, 1, _iterGlobal);
	occurencePerBlock.increment(0, 6, _iterGlobal);
	occurencePerBlock.increment(0, 7, _iterGlobal);
	occurencePerBlock.increment(0, 8, _iterGlobal);
	occurencePerBlock.increment(0, 9, _iterGlobal / _stepG);

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	
	setPnFromX << < _nBus, _blockSizeSmall >> > (Pn._matrixGPU, X._matrixGPU, _indiceBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent);

	/*std::cout << "--------" << std::endl;
	std::cout << " X " << std::endl;
	X.display(true);
	std::cout << " Pn " << std::endl;
	Pn.display(true);
	*/
	
	fc = calcFc(&Cost1, &Cost2, &Pn, &tempN2);
	// FB 5
	
	result->setResF(&resF);
	result->setIter(_iterGlobal);
	MatrixCPU PnCPU;
	Pn.toMatCPU(PnCPU);

	result->setPn(&PnCPU);
	
	result->setFc(fc);

#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 10, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 10, 1);

	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
	tall = clock() - tall;
	timeOPF = tall;

	result->setTime((float)tall / CLOCKS_PER_SEC);
	
}

void OPFADMMGPU2::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif
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
	initDFSPQ << <1, _nBus, _nBus* (8*sizeof(bool) + sizeof(int)) >> > (X._matrixGPU, Pb._matrixGPU, nChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nBus);
	communicateX << <_nBus, _blockSize >> > (X._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nAgentByBus._matrixGPU, _nBus);



	
	Y.set(&X);
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 11, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 11, 1);
	t1 = std::chrono::high_resolution_clock::now();
#endif

	updateChat();
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);
	
#endif


}

void OPFADMMGPU2::init(const Simparam& sim, const StudyCase& cas)
{
	// intialisation des matrixs et variables 
	
	clock_t t = clock();
	//std::cout << "init " << std::endl;
	_rho = sim.getRho();
	
	if (_rhol == 0) {
		_rhol = _rho;
	}
	if (consensus) {
		std::cout << "pas coder pour update Q !!!" << std::endl;
		exit(-1);
	}
	
	const int iterG = sim.getIterG();
	const int stepG = sim.getStepG();
	
	_nAgent = cas.getNagent();
	
	_nBus = cas.getNBus();
	_nLine = cas.getNLine(true); // ne doit pas �tre r�duit ici !!!
	_sizeOPFTotal = 3 * _nLine + 5 * _nBus + 2 * (_nAgent - 1); // L = nChild.sum()
	_sizeChat = 4 * _nBus + 2 * (_nAgent - 1);
	_numBlocksN = ceil((_nAgent + _blockSize - 1) / _blockSize);
	_numBlocksM = ceil((_sizeOPFTotal + _blockSize - 1) / _blockSize);
	_numBlocksB = ceil((_nBus + _blockSize - 1) / _blockSize);

	
	
	// il faut remettre sur CPU ce qu'il faut !!!
	if (tempL.getPos()) {
		tempL.transferCPU();
		_CoresChatBegin.transferCPU();
		_indiceBusBegin.transferCPU();
		Ancestor.transferCPU();
		PosChild.transferCPU();

		_indiceChildBegin.transferCPU();
		Childs.transferCPU();
	}


	tempL = MatrixGPU(_nLine, 1);


	//std::cout << _nAgent << " " << _nBus << " " << _nLine << std::endl;
	
	nChildCPU = MatrixCPU(_nBus, 1);
	CoresLineBus = cas.getCoresLineBus(true);
	_CoresBusAgent = cas.getCoresBusAgentLin(); // Cores[n] = b
	Ancestor = MatrixGPU(_nBus, 1, 0); // A_i = bus ant�c�dent de i
	Ancestor.set(0, 0, -1); // the slack bus has no ancestor
	ZsRe = MatrixGPU(cas.getZsRe());
	ZsIm = MatrixGPU(cas.getZsImag());
	ZsNorm = MatrixGPU(_nLine, 1);
	ZsNorm.multiplyT(&ZsRe, &ZsRe);
	tempL.multiplyT(&ZsIm, &ZsIm);
	ZsNorm.add(&tempL);
	tempL.transferGPU();

	if (!chekcase()) {
		throw std::invalid_argument("not a radial case");
	}

	for (int lold = 0; lold < _nLine; lold++) {
		int l = lold + 1;
		int busTo = l ;
		int busFrom = CoresLineBus.get(lold, 0);
		Ancestor.set(busTo, 0, busFrom);
		nChildCPU.set(busFrom, 0, nChildCPU.get(busFrom, 0) + 1); // pas parallelisable -> reduction chelou
	}
	
	nChild = MatrixGPU(nChildCPU, 1);

	
	_rhoInv = 1 / _rho;
	resF = MatrixCPU(3, (iterG / stepG) + 1, 0);

	
	MatrixGPU lowerBound(cas.getLowerBound(), 1); //voltage angle, voltage, line...
	MatrixGPU upperBound(cas.getUpperBound(), 1); //voltage angle, voltage, line...
	

	//std::cout << " local resolution " << std::endl;
	// local resolution
	tempN2 = MatrixGPU(2 * _nAgent, 1, 0, 1);
	tempB2 = MatrixGPU(2 * _nBus, 1, 0, 1);
	CoresSoloBusAgent = MatrixGPU(_nBus, 1, -1, 1);
	Pn = MatrixGPU(sim.getPn(), 1);
	Pmin = MatrixGPU(cas.getPmin(), 1);
	Pmax = MatrixGPU(cas.getPmax(), 1);
	
	Pb = MatrixGPU(2 * _nBus, 1, 0, 1);
	Pbmin = MatrixGPU(2 * _nBus, 1, 0, 1);
	Pbmax = MatrixGPU(2 * _nBus, 1, 0, 1);

	Cost1 = MatrixGPU(cas.geta(), 1);
	Cost2 = MatrixGPU(cas.getb(), 1);

	Pn.preallocateReduction();
	if (Pn.max2() < 0.00001) {
		Pn.add(&Pmin, &Pmax);
		Pn.divide(2);
	}
	_CoresAgentBus = MatrixGPU(cas.getCoresAgentBusLin(), 1);
	_CoresAgentBusBegin = MatrixGPU(cas.getCoresAgentBusLinBegin(), 1);
	_nAgentByBus = MatrixGPU(cas.getNagentByBus(), 1);
	_nAgentByBusCPU = cas.getNagentByBus();
	// remove the grid agent


	Pn.set(0, 0, 0, 1);
	Pmin.set(0, 0, 0, 1);
	Pmax.set(0, 0, 0, 1);
	Pn.set(_nAgent, 0, 0, 1);
	Pmin.set(_nAgent, 0, 0, 1);
	Pmax.set(_nAgent, 0, 0, 1);

	_nAgentByBusCPU.increment(0, 0, -1);
	_nAgentOn0 = _nAgentByBusCPU.get(0, 0);
	//std::cout << "remove loss agent" << std::endl;
	removeLossAgent << <1, 1 >> > (_nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU);


	ComputePFromAgentToBus();


	
	//std::cout << " creation " << std::endl;
	X = MatrixGPU(_sizeOPFTotal, 1, 0, 1); // Changement d'ordre !!!!!!!!!!!!
	Ypre = MatrixGPU(_sizeOPFTotal, 1, 0, 1); // (Pi, Qi, li, vi, pn..., qn..., vai, Pci ..., Qci... , lci...) !!!!!
	Y = MatrixGPU(_sizeOPFTotal, 1, 0, 1);
	Y.preallocateReduction();
	//YTrans = MatrixGPU(_sizeOPFTotal, 1, 0, 1);
	Mu = MatrixGPU(_sizeOPFTotal, 1, 0, 1);
	
	tempN1 = MatrixGPU(_nAgent, 1, 0, 1);
	tempNN = MatrixGPU(_nAgent, _nAgent, 0, 1);
	//tempM1 = new MatrixGPU[_nAgent];
	tempM = MatrixGPU(_sizeOPFTotal, 1, 0, 1);
	
	sizeOPFADMMGPU2 = MatrixGPU(_nBus, 1, 0, 1);
	sizeOPFADMMGPU2.preallocateReduction();
	sizeOPFADMMGPU2Big = MatrixGPU(_sizeOPFTotal, 1, 0, 1);

	_indiceBusBegin = MatrixGPU(_nBus, 1);
	_indiceBusBeginBig = MatrixGPU(_sizeOPFTotal, 1, 0, 1);
	_CoresChatBegin = MatrixGPU(_nBus, 1);
	
	int debut = 0;
	int debutChat = 0;
	for (int i = 0; i < _nBus; i++) {
		int m = nChildCPU.get(i, 0);
		int nB = _nAgentByBusCPU.get(i, 0);
		_indiceBusBegin.set(i, 0, debut);
		_CoresChatBegin.set(i, 0, debutChat);
		int sizeA = m * 3 + 5 + 2 * nB;
		debut += sizeA;
		debutChat += (4 + 2 * nB);
	}
	//_CoresChatBegin.display();


	_CoresChatBegin.transferGPU();
	_indiceBusBegin.transferGPU();
	defineSizeBig <<<_nBus, _blockSize >> > (sizeOPFADMMGPU2Big._matrixGPU, nChild._matrixGPU, _indiceBusBegin._matrixGPU, sizeOPFADMMGPU2._matrixGPU, _indiceBusBeginBig._matrixGPU, _nAgentByBus._matrixGPU);
	
	//sizeOPFADMMGPU2.display(true);


	_sizeOPFMax = sizeOPFADMMGPU2.max2();
	Hinv = MatrixGPU(_sizeOPFTotal, _sizeOPFMax, 0, 1);
	Q = MatrixGPU(_sizeOPFTotal, 1, 0, 1);
	
	Childs = MatrixGPU(_nLine, 1);
	PosChild = MatrixGPU(_nBus, 1, -1);

	Chat = MatrixGPU( _sizeChat, 1, 0, 1);
	VoltageLimit = MatrixGPU(2, _nBus, 0, 1); // min, max
	VoltageLimitReal = MatrixGPU(2, _nBus, 0, 1); // min, max
	
	
	_indiceChildBegin = MatrixGPU(_nLine, 1);
	//int sizeOPF2 = 1 * nChild.get(i, 0) + 9;
	
	
	MatrixCPU nChildTemp(_nBus, 1, 0);
	//lowerBound.display(true);
	//upperBound.display(true);
	initVoltageBound <<< _numBlocksB, _blockSize >> > (VoltageLimitReal._matrixGPU, VoltageLimit._matrixGPU, lowerBound._matrixGPU, upperBound._matrixGPU, nChild._matrixGPU, _nBus);

	//nChild.display(true);
	//VoltageLimit.display(true);
	//VoltageLimitReal.display(true);
	//

	//nChild.display();
	
	int debutChild = 0;
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
	/*Childs.display();
	Ancestor.display();
	nChildCPU.display();
	PosChild.display();
	std::cout << " _indiceChildBegin " << std::endl;
	_indiceChildBegin.display(true);*/
	Ancestor.transferGPU();
	PosChild.transferGPU();
	debut = 0;
	//std::cout << " Hinv " << std::endl;
	for (int i = 0; i < _nBus; i++) {
		// (Pi, Qi, li, vi, pn..., qn..., vai, Pci ..., Qci... , lci...) !!!!!
		int m = nChildCPU.get(i, 0);
		int nB = _nAgentByBusCPU.get(i, 0);
		int sizeA = nChildCPU.get(i,0) * 3 + 5 + 2 * nB;
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
		temp33.invertEigen(&temp33);
		temp3M.MultiplyMatMat(&temp33, &A);
		tempMM.multiplyTrans(&A, &temp3M, 0);

		tempMMbis.setEyes(-1);
		tempMMbis.add(&tempMM);
		MatrixGPU tempMMGPU = MatrixGPU(tempMMbis, 1);
		
		Hinv.setBloc(debut, debut + sizeA, 0, sizeA, &tempMMGPU);
		debut += sizeA;
	}
	Hinv.divide(_rho);
	//Hinv.display(true);
	_indiceChildBegin.transferGPU();
	Childs.transferGPU();
	//std::cout << " Childs " << std::endl;
	//Childs.display(true);
    //std::cout << " init valeur " << std::endl;
	
	initPQAgentV <<< _nBus, _blockSizeSmall >> > (X._matrixGPU, _indiceBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, Pn._matrixGPU, _nAgent);
	
	
	/*std::cout << " X " << std::endl;
	X.display(true);
	std::cout << " _indiceBusBegin " << std::endl;
	_indiceBusBegin.display(true);
	std::cout << " _indiceChildBegin " << std::endl;
	_indiceChildBegin.display(true);
	std::cout << " Childs " << std::endl;
	Childs.display(true);
	std::cout << " nChild " << std::endl;
	nChild.display(true);
	std::cout << " posChild " << std::endl;
	PosChild.display(true);*/

	initDFSPQ << <1, _nBus, _nBus*(8*sizeof(bool) + sizeof(int)) >> > (X._matrixGPU, Pb._matrixGPU, nChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nBus);
	
	

	
	////CHECK_LAST_CUDA_ERROR();
	
	communicateX << <_nBus, _blockSize >> > (X._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nAgentByBus._matrixGPU, _nBus);

	
	Y.set(&X);
	/*std::cout << " X " << std::endl;
	X.display(true);
	std::cout << " Y " << std::endl;
	Y.display(true);
	std::cout << " Q " << std::endl;
	Q.display(true);
	std::cout << " Mu "<< std::endl;
	Mu.display(true);*/

	updateChat();
	/*std::cout << " Chat " << std::endl;
	Chat.display(true);*/
	//std::cout << " Bpt2 " << std::endl;
	//Bpt2.display(true);
	/*std::cout << " Nagent " << std::endl;
	_nAgentByBus.display(true);
	std::cout << " Bus Agent : agent->bus " << std::endl;
	_CoresBusAgent.display(true);
	std::cout << " Agent bus : bus->agent " << std::endl;
	_CoresAgentBus.display(true);
	std::cout << " Agent bus begin : bus->agent " << std::endl;
	_CoresAgentBusBegin.display(true);
	
	std::cout << "Apt2 " << std::endl;
	Apt2.display(true);
	std::cout << " Pmin " << std::endl;
	Pmin.display(true);
	std::cout << " Pma " << std::endl;
	Pmax.display(true);
	std::cout << " CoresSoloBusAgent " << std::endl;
	CoresSoloBusAgent.display(true);*/


	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	//std::cout << "---------------------------------------------------------------------------------------" << std::endl;
}

void OPFADMMGPU2::solveConsensus(float eps, MatrixCPU* PSO)
{
	throw std::invalid_argument("WIP !!");
}

void OPFADMMGPU2::initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO)
{
	throw std::invalid_argument("WIP !!");
}

void OPFADMMGPU2::updateConsensus(MatrixCPU* Pmarket)
{
	throw std::invalid_argument("WIP !!");
}

void OPFADMMGPU2::solveConsensus(float eps, MatrixGPU* PSO)
{
	throw std::invalid_argument("WIP !!");
}

void OPFADMMGPU2::updateConsensus(MatrixGPU* Pmarket)
{
	throw std::invalid_argument("WIP !!");
}



void OPFADMMGPU2::updateGlobalProb() {
	
	Ypre.swap(&Y);
	int numBlock = _sizeOPFTotal;
	switch (_blockSize) {
	case 512:
		updateY<512> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPU2Big._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case 256:
		updateY<256> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPU2Big._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case 128:
		updateY<128> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPU2Big._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case 64:
		updateY< 64> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPU2Big._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case 32:
		updateY< 32> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPU2Big._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case 16:
		updateY< 16> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPU2Big._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case  8:
		updateY<  8> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPU2Big._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case  4:
		updateY<  4> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPU2Big._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case  2:
		updateY<  2> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPU2Big._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case  1:
		updateY<  1> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPU2Big._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	}
	
	Y.set(3, 0, 1, 1);
	Y.set(4 + 2 * _nAgentOn0, 0, 1, 1);

}



void OPFADMMGPU2::updateMu()
{
	updateMUGPU << <_numBlocksM, _blockSize >> > (Mu._matrixGPU, Y._matrixGPU, X._matrixGPU, _rho, _sizeOPFTotal);	
}


float OPFADMMGPU2::getPLoss()
{
	
	return Pn.sum(0, _nAgent);
}

float OPFADMMGPU2::getQLoss()
{
	return Pn.sum(_nAgent, 2 * _nAgent);
}

void OPFADMMGPU2::ComputePFromAgentToBus()
{
	int numBlock = _nBus;
	switch (_blockSize) {
	case 512:
		ComputePFromAgentToBusGPU<512> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case 256:
		ComputePFromAgentToBusGPU<256> << <numBlock, _blockSizeSmall >>>  (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case 128:
		ComputePFromAgentToBusGPU<128> << <numBlock, _blockSizeSmall >>> (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case 64:
		ComputePFromAgentToBusGPU< 64> << <numBlock, _blockSizeSmall >>> (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case 32:
		ComputePFromAgentToBusGPU< 32> << <numBlock, _blockSizeSmall >>> (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case 16:
		ComputePFromAgentToBusGPU< 16> << <numBlock, _blockSizeSmall >>> (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case  8:
		ComputePFromAgentToBusGPU<  8> << <numBlock, _blockSizeSmall >>> (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case  4:
		ComputePFromAgentToBusGPU<  4> << <numBlock, _blockSizeSmall >>> (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case  2:
		ComputePFromAgentToBusGPU<  2> << <numBlock, _blockSizeSmall >>> (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case  1:
		ComputePFromAgentToBusGPU<  1> << <numBlock, _blockSizeSmall >>> (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	}
}

void OPFADMMGPU2::updateChat()
{
	int numBlock = _nBus;
	switch (_blockSizeSmall) {
	case 512:
		updateChatGPU2<512> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, _nBus);
		break;
	case 256:
		updateChatGPU2<256> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, _nBus);
		break;
	case 128:
		updateChatGPU2<128> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, _nBus);
		break;
	case 64:
		updateChatGPU2< 64> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, _nBus);
		break;
	case 32:
		updateChatGPU2< 32> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, _nBus);
		break;
	case 16:
		updateChatGPU2< 16> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, _nBus);
		break;
	case  8:
		updateChatGPU2<  8> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, _nBus);
		break;
	case  4:
		updateChatGPU2<  4> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, _nBus);
		break;
	case  2:
		updateChatGPU2<  2> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, _nBus);
		break;
	case  1:
		updateChatGPU2<  1> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, _nBus);
		break;
	}
}

void OPFADMMGPU2::CommunicationX()
{
	 // X = { Pi, Qi, li, vi, (pn ...), qn..., vAi,  Pci ... , Qci ... , lci ... for all child Ci }
	
	communicateX << <_nBus, _blockSize >> > (X._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nAgentByBus._matrixGPU, _nBus);

	// Y = { Pi, Qi, vi, li, vAi, (pn ...), qn...,  Pci ... , Qci ... , lci ... for all child Ci }

	updateQ << <_numBlocksM, _blockSize >> > (Q._matrixGPU, X._matrixGPU, Mu._matrixGPU, _rho, _sizeOPFTotal);
}



float OPFADMMGPU2::updateRes(int indice) 
{
	
	float resS = _rho * Y.max2(&Ypre);
	float resR = Y.max2(&X);
	float resV = 0;
	
	resF.set(0, indice, resR);
	resF.set(1, indice, resS);
	resF.set(2, indice, resV);

	/*std::cout << resS << " " << resR << std::endl;
	std::cout << " Y " << std::endl;
	Y.display(true);
	std::cout << " X " << std::endl;
	X.display(true);*/


	if (resR > _mu * resS) {
		_rho = _tau * _rho;
		
		Hinv.divide(_tau);
		//std::cout << _iterGlobal << "rho augmente " << _rho << std::endl;
	}
	else if (resS > _mu * resR) {// rho = rho / tau_inc;
		_rho = _rho / _tau;
		Hinv.multiply(_tau);
		
		//std::cout << _iterGlobal << "rho diminue " << _rho << std::endl;
	}/**/


	return MAX(MAX(resV, resS), resR);
}

int OPFADMMGPU2::feasiblePoint()
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




void OPFADMMGPU2::display() {

	std::cout.precision(3);

	X.transferCPU();
	Y.transferCPU();
	Mu.transferCPU();
	Pn.transferCPU();
	_indiceBusBegin.transferCPU();
	_nAgentByBus.transferCPU();
	VoltageLimitReal.transferCPU();
	Pmin.transferCPU();
	Pmax.transferCPU();
	Pbmax.transferCPU();
	Pbmin.transferCPU();
	Pb.transferCPU();
	Cost1.transferCPU();
	Cost2.transferCPU();

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

		
	for (int b = 0; b < _nBus; b++) {
		int begining = _indiceBusBegin.get(b, 0);
	std::cout << std::setw(5) << b << "|" << std::setw(12) << sqrt(X.get(begining + 3,0)) << " |" << std::setw(16)
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
		std::cout << std::setw(6) << l << "|" << std::setw(12) << CoresLineBus.get(l, 0) << " |" << std::setw(12)
			<< CoresLineBus.get(l, 1) << "|" << std::setw(16) << X.get(begining + 0, 0)
			<< "|" << std::setw(16) << X.get(begining + 1, 0) << "|" << std::setw(16)
			<< X.get(begining + 2, 0) << "|" << std::setw(19) << X.get(begining + 2, 0) * ZsRe.get(l, 0) << "|" << std::endl;
	}
	std::cout << std::endl << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "     Constraints                                                                                        |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Bus | Voltage | Voltage | Voltage |        Power Injection          |          Power Injection         |" << std::endl;
	std::cout << "  #  | Mag(pu) | MIN(pu) |  MAX(pu)|  P (pu) | Pmin (pu) | Pmax (pu) |  Q (pu)  | Qmin (pu) | Qmax (pu) |" << std::endl;
	std::cout << "-----|---------|---------|---------|---------|-----------|-----------|----------|-----------|-----------|" << std::endl;
	

	for (int b = 0; b < _nBus; b++) {
		int begining = _indiceBusBegin.get(b, 0);
		std::cout << std::setw(5) << b << "|" << std::setw(8) << sqrt(Y.get(begining + 3, 0)) << " |" << std::setw(9)
			<< VoltageLimitReal.get(0, b) << "|" << std::setw(9) << VoltageLimitReal.get(1, b)
			<< "|" << std::setw(9) << Pb.get(b, 0) << "|" << std::setw(11)
			<< Pbmin.get(b, 0) << "|" << std::setw(11) << Pbmax.get(b, 0)  << "|" << std::setw(10) << Pb.get(b + _nBus, 0)
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





template <unsigned int _blockSizeSmall>
__global__ void updateChatGPU2(float* Chat, float* Y, float* MU, float* nChild, float* Ancestor, float* posChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, float* CoresChatBegin, float* nAgentByBus, float _rho, int nBus) {

	int bus = blockIdx.x;
	int index = threadIdx.x;
	int step = blockDim.x;
	int begin = CoresChatBegin[bus];

	__shared__ float shArr[_blockSizeSmall]; // c'est grand pour pas grand chose...
	

	int indice = indiceBusBegin[bus];
	int indiceChild = (bus < (nBus - 1)) ? indiceChildBegin[bus] : 0;
	int nb = nChild[bus];
	int Ai = Ancestor[bus];
	int nAgent = nAgentByBus[bus];
	int c = posChild[bus];
	float var = 0;
	int borne = 4 + 2 * nAgent;

	if (index < borne) {
		//float Phat, Qhat, lhat, vihat, pnhat..., qnhat...;
		var =  Y[indice + index] / (1 + (index < 4)) - MU[indice + index] / ((1 + (index < 4)) * _rho);
		if (bus > 0) {
			if (index < 3) {
				int nAi = nChild[Ai];
				int nAgentAi = nAgentByBus[Ai];
				int indiceAncBus = indiceBusBegin[Ai] + 5 + 2 * nAgentAi +  nAi * index + c;
				//var = indiceAncBus;
				var += Y[indiceAncBus] / 2  - MU[indiceAncBus] / (2 * _rho); 
			}			
		}
	}
	float vhat = 0;
	float muhat = 0;
	for (int i = index; i < nb; i += step) {
		int Bus2 = Childs[indiceChild + i];
		int indiceBusChild = indiceBusBegin[Bus2];
		int nAgent2 = nAgentByBus[Bus2];
		muhat += MU[indiceBusChild + 4 + 2 * nAgent2]; // pas du tout coalescent
		vhat += Y[indiceBusChild + 4 + 2 * nAgent2]; // pas du tout coalescent
	}
	shArr[index] = vhat / (nb + 1) - muhat / (_rho * (nb + 1));
	__syncthreads();
	for (int size = _blockSizeSmall / 2; size > 0; size /= 2) { //uniform
		if (index < size) {
			shArr[index] += shArr[index + size];
		}
		__syncthreads();
	}

	if (index < borne) {
		if (index == 3) {
			var = shArr[0] + Y[indice + 3] / (nb + 1) - MU[indice + 3] / (_rho * (nb + 1)); //shArr[0];
		}
		Chat[begin + index] = var; // coalescent  !!!!
	}
}


