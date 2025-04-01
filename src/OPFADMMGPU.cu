#include "../head/OPFADMMGPU.cuh"
 
#define NMAXAGENTPERTHREAD 5

OPFADMMGPU::OPFADMMGPU() : MethodOPFGPU()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " OPFADMMGPU Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	timePerBlock = MatrixCPU(1, 12, 0); // Fb0, Fb11abcd, FB12, Fb2, Fb3, Fb4, Fb5,FB6, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 12, 0); //nb de fois utilisé pendant la simu
}




OPFADMMGPU::OPFADMMGPU(float rho) : MethodOPFGPU()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default OPFADMMGPU Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
	timePerBlock = MatrixCPU(1, 12, 0); // Fb0, Fb11, FB12, Fb2, Fb3, Fb4, Fb5,FB6, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 12, 0); //nb de fois utilisé pendant la simu
}

OPFADMMGPU::~OPFADMMGPU()
{
	/*DELETEA(tempM1);
	DELETEA(tempM);

	DELETEA(X);
	DELETEA(Ypre);
	DELETEA(Y);
	DELETEA(YTrans);
	DELETEA(Mu);

	DELETEA(Hinv);
	DELETEA(A);
	DELETEA(Q);

	DELETEA(Childs);*/

}
void OPFADMMGPU::setParam(float rho)
{
	_rho = rho;
}

bool OPFADMMGPU::chekcase()
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
	if (ZsRe.getNLin() == 0 || ZsIm.getNLin() == 0) {
		std::cout << "matrice non defined, ZsRe, Zs Im, Yd" << std::endl;
		ZsRe.display();
		ZsIm.display();
		return false;
	}

	return true;
}

void OPFADMMGPU::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
{
#ifdef DEBUG_SOLVE
	cas.display();
	sim.display(1);
#endif // DEBUG_SOLVE

	clock_t tall = clock();
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
		timePerBlock.increment(0, 0, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION
	}

	_iterG = sim.getIterG();
	int iterL = sim.getIterL();
	_stepG = sim.getStepG();
	int stepL = sim.getStepL();

	float epsG = sim.getEpsG();
	float epsL = MYMIN(sim.getEpsL(), epsG / 200);
	float rhoInit = sim.getRho();


	float fc = 0;
	float resG = 2 * epsG;
	float resL = 2 * epsL;
	_iterGlobal = 0;

	//Chat.display(true);
	//

	/*CoresSoloBusAgent.display(true);
	Apt1.display(true);
	Apt2.display(true);
	Bpt2.display(true);
	Cost1.display(true);
	Cost2.display(true);
	Pmin.display(true);
	Pmax.display(true);

	PnPre.display(true);
	PnMoy.display(true);
	PnTilde.display(true);
	MuL.display(true);
	_nAgentByBus.display(true);
	std::cout << _rhol << " " << epsL << " " << iterL << " " << _nAgent << " " << _nBus;
	_CoresAgentBus.display(true); 
	_CoresAgentBusBegin.display(true);

	std::cout << "------" << std::endl;*/

	while ((_iterGlobal < _iterG) && (resG > epsG)) {

		/*std::cout << "--------" << std::endl;
		std::cout << " Pn " << std::endl;
		Pn.display(true);
		std::cout << " PnTilde " << std::endl;
		PnTilde.display(true);
		std::cout << " X " << std::endl;
		X.display(true);



		std::cout << " Y " << std::endl;
		Y.display(true);
		std::cout << " Mu " << std::endl;
		Mu.display(true);
		Chat.display(true);

*/
#ifdef INSTRUMENTATION
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

		updateLocalProb(epsL, iterL);
		//CHECK_LAST_CUDA_ERROR();
		/*std::cout << " Pn " << std::endl;
		Pn.display(true);
		std::cout << " PnMoy " << std::endl;
		PnMoy.display(true);
		std::cout << " PnTilde " << std::endl;
		PnTilde.display(true);
		std::cout << " MuL " << std::endl;
		MuL.display(true);*/
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 1, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

		//std::cout << _numBlocksB << " " << _blockSize << std::endl;
		updateXOPFADMM << <_numBlocksB, _blockSize >> > (X._matrixGPU, Chat._matrixGPU, VoltageLimit._matrixGPU, PnTilde._matrixGPU, _nAgentByBus._matrixGPU, nChild._matrixGPU, _indiceBusBegin._matrixGPU, _nBus, Lagrange);
		//CHECK_LAST_CUDA_ERROR();
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 5, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

		CommunicationX();
		//CHECK_LAST_CUDA_ERROR();
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

		updateGlobalProb();
		updateMu();
		//CHECK_LAST_CUDA_ERROR();
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

		updateChat();
		CHECK_LAST_CUDA_ERROR();
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		// FB 4
		if (!(_iterGlobal % _stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateRes(_iterGlobal / _stepG);
			//std::cout << _iterGlobal << " " << _iterLocal << " " << _rho << " " << resL << " " << resF.get(0, _iterGlobal / _stepG) << " " << resF.get(1, _iterGlobal / _stepG) << std::endl;
			//resG = 1;
#ifdef INSTRUMENTATION
			cudaDeviceSynchronize();
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 9, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		}
		//std::cout << iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;

		_iterGlobal++;
	}
	//std::cout << iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, (iterGlobal - 1) / stepG) << " " << resF.get(1, (iterGlobal - 1) / stepG) << " " << resF.get(2, (iterGlobal - 1) / stepG) << std::endl;


#ifdef INSTRUMENTATION	
	occurencePerBlock.increment(0, 1, _iterGlobal);
	occurencePerBlock.increment(0, 5, _iterGlobal);
	occurencePerBlock.increment(0, 6, _iterGlobal);
	occurencePerBlock.increment(0, 7, _iterGlobal);
	occurencePerBlock.increment(0, 8, _iterGlobal);
	occurencePerBlock.increment(0, 9, _iterGlobal / _stepG);

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	//std::cout << "--------" << std::endl;
	/*std::cout << " Pn " << std::endl;
	Pn.display(true);
	std::cout << " PnTilde " << std::endl;
	PnTilde.display(true);*/

	fc = calcFc(&Cost1, &Cost2, &Pn, &tempN2);
	// FB 5

	result->setResF(&resF);

	/*MatrixGPU Param(1, 12);
	Param.set(0, 0, _nAgent);
	Param.set(0, 1, _nBus);
	Param.set(0, 2, _nLine);
	Param.set(0, 3, rhoInit);
	Param.set(0, 4, _rhol);
	Param.set(0, 5, _stepG);
	Param.set(0, 6, fc);
	Param.set(0, 7, Pn.get(1, 0));
	Param.set(0, 8, Pn.get(2, 0));
	Param.set(0, 9, Pn.get(4, 0));
	Param.set(0, 10, Pn.get(5, 0));
	Param.set(0, 11, _iterGlobal);

	std::string nameFile = "Residuals_RhoConst_W0agent.csv";
	Param.saveCSV(nameFile);
	resF.saveCSV(nameFile);*/

	result->setIter(_iterGlobal);
	MatrixCPU PnCPU;
	Pn.toMatCPU(PnCPU);
	result->setPn(&PnCPU);
	
	MatrixCPU Pb(getPb());
	MatrixCPU Phi(getPhi());
	MatrixCPU E(getE());
	
	result->setE(&E);
	result->setPhi(&Phi);
	result->setPb(&Pb);
	

	result->setFc(fc);

#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 10, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 10, 1);

	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
	tall = clock() - tall;
	timeOPF = tall;

	result->setTime((float)tall / CLOCKS_PER_SEC);

}

void OPFADMMGPU::updateP0(const StudyCase& cas)
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
	divideMultiplyByNagentByBus << <_numBlocksB, _blockSize >> > (Apt1._matrixGPU, Apt2._matrixGPU, PnTilde._matrixGPU, PnTmin._matrixGPU, PnTmax._matrixGPU, _nAgentByBus._matrixGPU, _rhol, _nBus);
	PnMoy.set(&PnTilde);
	initPQ << < _numBlocksB, _blockSize >> > (X._matrixGPU, _indiceBusBegin._matrixGPU, _nAgentByBus._matrixGPU, PnTilde._matrixGPU, _nBus);
	initDFSPQ << <1, _nBus, _nBus* (sizeof(bool) + sizeof(int)) >> > (X._matrixGPU, nChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nBus);
	communicateX << <_nBus, _blockSize >> > (X._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nBus);


	Y.set(&X);
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 11, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 11, 1);
	t1 = std::chrono::high_resolution_clock::now();
#endif

	updateChat();
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);

#endif


}

void OPFADMMGPU::init(const Simparam& sim, const StudyCase& cas)
{
	// intitilisation des matrixs et variables 

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
	_sizeOPFTotal = 3 * _nLine + 7 * _nBus; // L = nChild.sum()
	_numBlocksN = ceil((_nAgent + _blockSize - 1) / _blockSize);
	_numBlocksM = ceil((_sizeOPFTotal + _blockSize - 1) / _blockSize);
	_numBlocksB = ceil((_nBus + _blockSize - 1) / _blockSize);
	_nAgentByBus = MatrixGPU(cas.getNagentByBus(), 1);
	_nAgentByBus.preallocateReduction();
	// 2 sources d'erreurs il faut que 
	int nAgentMaByBus = _nAgentByBus.max2();
	if (nAgentMaByBus > _blockSizeSmall / 2) {
		throw std::invalid_argument("the number of agent by bus is too high , must change blocksize if possible");
	}

	// il faut remettre sur CPU ce qu'il faut !!!
	if (tempL.getPos()) {
		Ancestor.transferCPU();
		_indiceBusBegin.transferCPU();
		PosChild.transferCPU();
		_indiceChildBegin.transferCPU();
		Childs.transferCPU();
		tempL.transferCPU();
	}

	tempL = MatrixGPU(_nLine, 1);


	//std::cout << _nAgent << " " << _nBus << " " << _nLine << std::endl;

	nChildCPU = MatrixCPU(_nBus, 1);
	CoresLineBus = cas.getCoresLineBus();
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
		int busTo = l;
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
	Pn.preallocateReduction();
	Pmin = MatrixGPU(cas.getPmin(), 1);
	Pmax = MatrixGPU(cas.getPmax(), 1);


	PnTmin = MatrixGPU(2 * _nBus, 1, 0, 1);
	PnTmax = MatrixGPU(2 * _nBus, 1, 0, 1);

	Cost1 = MatrixGPU(cas.geta(), 1);
	Cost2 = MatrixGPU(cas.getb(), 1);


	PnMoy = MatrixGPU(2 * _nBus, 1, 0, 1);
	PnPre = MatrixGPU(sim.getPn(), 1);
	MuL = MatrixGPU(2 * _nBus, 1, 0, 1);
	PnTilde = MatrixGPU(2 * _nBus, 1, 0, 1);
	Bp1 = MatrixGPU(2 * _nAgent, 1, 0, 1);
	Bpt1 = MatrixGPU(2 * _nBus, 1, 0, 1);
	Bpt2 = MatrixGPU(2 * _nBus, 1, 0, 1);
	Apt1 = MatrixGPU(2 * _nBus, 1, 0, 1);
	Apt2 = MatrixGPU(2 * _nBus, 1, 0, 1);


	if (Pn.max2() < 0.00001) {
		Pn.add(&Pmin, &Pmax);
		Pn.divide(2);
	}
	_CoresAgentBus = MatrixGPU(cas.getCoresAgentBusLin(), 1);
	_CoresAgentBusBegin = MatrixGPU(cas.getCoresAgentBusLinBegin(), 1);

	// remove the grid agent


	Pn.set(0, 0, 0, 1);
	Pmin.set(0, 0, 0, 1);
	Pmax.set(0, 0, 0, 1);
	Pn.set(_nAgent, 0, 0, 1);
	Pmin.set(_nAgent, 0, 0, 1);
	Pmax.set(_nAgent, 0, 0, 1);



	//std::cout << "remove loss agent" << std::endl;
	removeLossAgent << <1, 1 >> > (_nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU);


	ComputePFromAgentToBus();


	divideMultiplyByNagentByBus << <_numBlocksB, _blockSize >> > (Apt1._matrixGPU, Apt2._matrixGPU, PnTilde._matrixGPU, PnTmin._matrixGPU, PnTmax._matrixGPU, _nAgentByBus._matrixGPU, _rhol, _nBus);


	PnMoy.set(&PnTilde);
	
	
	/*std::cout << "Apt1 " << std::endl;
	Apt1.display(true);
	std::cout << "Apt2 " << std::endl;
	Apt2.display(true);
	std::cout << " Pmin " << std::endl;
	Pmin.display(true);
	std::cout << " Pmax " << std::endl;
	Pmax.display(true);

	std::cout << " Pn  limits " << std::endl;
	PnTmin.display(true);
	PnTmax.display(true);*/
	


	//std::cout << " creation " << std::endl;
	X = MatrixGPU(_sizeOPFTotal, 1, 0, 1); // Changement d'ordre !!!!!!!!!!!!
	Ypre = MatrixGPU(_sizeOPFTotal, 1, 0, 1); // (Pi, Qi, li, vi, pi, qi, vai, Pci ..., Qci... , lci...) !!!!!
	Y = MatrixGPU(_sizeOPFTotal, 1, 0, 1);
	Y.preallocateReduction();
	//YTrans = MatrixGPU(_sizeOPFTotal, 1, 0, 1);
	Mu = MatrixGPU(_sizeOPFTotal, 1, 0, 1);

	tempN1 = MatrixGPU(_nAgent, 1);
	tempNN = MatrixGPU(_nAgent, _nAgent);
	//tempM1 = new MatrixGPU[_nAgent];
	tempM = MatrixGPU(_sizeOPFTotal, 1, 0, 1);
	sizeOPFADMMGPU = MatrixGPU(_nBus, 1, 0, 1);
	sizeOPFADMMGPU.preallocateReduction();
	sizeOPFADMMGPUBig = MatrixGPU(_sizeOPFTotal, 1, 0, 1);
	_indiceBusBegin = MatrixGPU(_nBus, 1);
	_indiceBusBeginBig = MatrixGPU(_sizeOPFTotal, 1, 0, 1);
	int debut = 0;
	for (int i = 0; i < _nBus; i++) {
		int m = nChildCPU.get(i, 0);
		_indiceBusBegin.set(i, 0, debut);
		int sizeA = m * 3 + 7;
		debut += sizeA;
	}
	_indiceBusBegin.transferGPU();
	defineSizeBig << <_nBus, _blockSize >> > (sizeOPFADMMGPUBig._matrixGPU, nChild._matrixGPU, _indiceBusBegin._matrixGPU, sizeOPFADMMGPU._matrixGPU, _indiceBusBeginBig._matrixGPU);
	/*sizeOPFADMMGPU = nChild;
	sizeOPFADMMGPU.multiply(3);
	sizeOPFADMMGPU.add(7);
	sizeOPFADMMGPU.display(true);
	sizeOPFADMMGPUBig.display(true);*/

	_sizeOPFMax = sizeOPFADMMGPU.max2();
	Hinv = MatrixGPU(_sizeOPFTotal, _sizeOPFMax, 0, 1);
	Q = MatrixGPU(_sizeOPFTotal, 1, 0, 1);

	Childs = MatrixGPU(_nLine, 1);
	PosChild = MatrixGPU(_nBus, 1, -1);
	Chat = MatrixGPU(6, _nBus, 0, 1);
	VoltageLimit = MatrixGPU(2, _nBus, 0, 1); // min, max
	VoltageLimitReal = MatrixGPU(2, _nBus, 0, 1); // min, max


	_indiceChildBegin = MatrixGPU(_nLine, 1);
	//int sizeOPF2 = 1 * nChild.get(i, 0) + 9;


	MatrixCPU nChildTemp(_nBus, 1, 0);
	//lowerBound.display(true);
	//upperBound.display(true);
	initVoltageBound << < _numBlocksB, _blockSize >> > (VoltageLimitReal._matrixGPU, VoltageLimit._matrixGPU, lowerBound._matrixGPU, upperBound._matrixGPU, nChild._matrixGPU, _nBus);

	//nChild.display(true);
	//VoltageLimit.display(true);
	//VoltageLimitReal.display(true);
	//

	//nChild.display();
	//std::cout << " Child " << std::endl;
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
		int m = nChildCPU.get(i, 0);
		int sizeA = nChildCPU.get(i, 0) * 3 + 7;
		MatrixCPU A(2 + (i > 0), sizeA);

		if (i > 0) {
			A.set(2, 0, 2 * ZsRe.get(i - 1, 0));
			A.set(2, 1, 2 * ZsIm.get(i - 1, 0));
			A.set(2, 2, -ZsNorm.get(i - 1, 0));
			A.set(2, 3, -1);
			A.set(2, 6, 1);
			A.set(0, 0, -1);
			A.set(1, 1, -1);
		}
		A.set(0, 4, 1);
		A.set(1, 5, 1);

		for (int j = 0; j < m; j++) {
			int c = Childs.get(_indiceChildBegin.get(i, 0) + j, 0);
			A.set(0, 7 + j, 1); // Pci
			A.set(1, 7 + m + j, 1); // Qci
			A.set(0, 7 + 2 * m + j, -ZsRe.get(c - 1, 0)); // -R l
			A.set(1, 7 + 2 * m + j, -ZsIm.get(c - 1, 0)); // -X l
		}

		//A[i].display();

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
	Hinv.divide(_rho);
	//Hinv.display(true);
	_indiceChildBegin.transferGPU();
	Childs.transferGPU();
	//std::cout << " Childs " << std::endl;
	//Childs.display(true);
	//std::cout << " init valeur " << std::endl;

	initPQV << < _numBlocksB, _blockSize >> > (X._matrixGPU, _indiceBusBegin._matrixGPU, _nAgentByBus._matrixGPU, PnTilde._matrixGPU, _nBus);
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
	//
	initDFSPQ << <1, _nBus, _nBus* (sizeof(bool) + sizeof(int)) >> > (X._matrixGPU, nChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nBus);




	//CHECK_LAST_CUDA_ERROR();

	communicateX << <_nBus, _blockSize >> > (X._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nBus);


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
	Chat.display(true);
	std::cout << " Bpt2 " << std::endl;
	Bpt2.display(true);
	std::cout << " Cp " << std::endl;
	Cost2.display(true);
	std::cout << " Ap2 " << std::endl;
	Cost1.display(true);
	std::cout << " Nagent " << std::endl;
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
	CHECK_LAST_CUDA_ERROR();
}

void OPFADMMGPU::updateGlobalProb() {

	Ypre.swap(&Y);
	int numBlock = _sizeOPFTotal;
	switch (_blockSize) {
	case 512:
		updateY<512> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case 256:
		updateY<256> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case 128:
		updateY<128> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case 64:
		updateY< 64> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case 32:
		updateY< 32> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case 16:
		updateY< 16> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case  8:
		updateY<  8> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case  4:
		updateY<  4> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case  2:
		updateY<  2> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	case  1:
		updateY<  1> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeOPFADMMGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeOPFMax);
		break;
	}

	Y.set(3, 0, 1, 1);
	Y.set(6, 0, 1, 1);

}


void OPFADMMGPU::solveConsensus(float eps, MatrixCPU* PSO)
{
	throw std::invalid_argument("WIP !!");
}

void OPFADMMGPU::initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO)
{
	throw std::invalid_argument("WIP !!");
}

void OPFADMMGPU::updateConsensus(MatrixCPU* Pmarket)
{
	throw std::invalid_argument("WIP !!");
}

void OPFADMMGPU::solveConsensus(float eps, MatrixGPU* PSO)
{
	throw std::invalid_argument("WIP !!");
}

void OPFADMMGPU::updateConsensus(MatrixGPU* Pmarket)
{
	throw std::invalid_argument("WIP !!");
}


void OPFADMMGPU::updateLocalProb(float epsL, int nIterL) {

	int numBlocks = _nBus;
	switch (_blockSizeSmall) {
	case 512:
		updatePnPGPUSharedResidual<512> << <numBlocks, _blockSizeSmall >> > (Pn._matrixGPU, PnPre._matrixGPU, PnMoy._matrixGPU, PnTilde._matrixGPU, MuL._matrixGPU, _nAgentByBus._matrixGPU, _rhol, Cost1._matrixGPU, Cost2._matrixGPU,
			Pmin._matrixGPU, Pmax._matrixGPU, Apt1._matrixGPU, Apt2._matrixGPU, Bpt2._matrixGPU, CoresSoloBusAgent._matrixGPU, _CoresAgentBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, epsL, nIterL, _nAgent, _nBus);
		break;
	case 256:
		updatePnPGPUSharedResidual<256> << <numBlocks, _blockSizeSmall >> > (Pn._matrixGPU, PnPre._matrixGPU, PnMoy._matrixGPU, PnTilde._matrixGPU, MuL._matrixGPU, _nAgentByBus._matrixGPU, _rhol, Cost1._matrixGPU, Cost2._matrixGPU,
			Pmin._matrixGPU, Pmax._matrixGPU, Apt1._matrixGPU, Apt2._matrixGPU, Bpt2._matrixGPU, CoresSoloBusAgent._matrixGPU, _CoresAgentBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, epsL, nIterL, _nAgent, _nBus);
		break;
	case 128:
		updatePnPGPUSharedResidual<128> << <numBlocks, _blockSizeSmall >> > (Pn._matrixGPU, PnPre._matrixGPU, PnMoy._matrixGPU, PnTilde._matrixGPU, MuL._matrixGPU, _nAgentByBus._matrixGPU, _rhol, Cost1._matrixGPU, Cost2._matrixGPU,
			Pmin._matrixGPU, Pmax._matrixGPU, Apt1._matrixGPU, Apt2._matrixGPU, Bpt2._matrixGPU, CoresSoloBusAgent._matrixGPU, _CoresAgentBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, epsL, nIterL, _nAgent, _nBus);
		break;
	case 64:
		updatePnPGPUSharedResidual< 64> << <numBlocks, _blockSizeSmall >> > (Pn._matrixGPU, PnPre._matrixGPU, PnMoy._matrixGPU, PnTilde._matrixGPU, MuL._matrixGPU, _nAgentByBus._matrixGPU, _rhol, Cost1._matrixGPU, Cost2._matrixGPU,
			Pmin._matrixGPU, Pmax._matrixGPU, Apt1._matrixGPU, Apt2._matrixGPU, Bpt2._matrixGPU, CoresSoloBusAgent._matrixGPU, _CoresAgentBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, epsL, nIterL, _nAgent, _nBus);
		break;
	case 32:
		updatePnPGPUSharedResidual< 32> << <numBlocks, _blockSizeSmall >> > (Pn._matrixGPU, PnPre._matrixGPU, PnMoy._matrixGPU, PnTilde._matrixGPU, MuL._matrixGPU, _nAgentByBus._matrixGPU, _rhol, Cost1._matrixGPU, Cost2._matrixGPU,
			Pmin._matrixGPU, Pmax._matrixGPU, Apt1._matrixGPU, Apt2._matrixGPU, Bpt2._matrixGPU, CoresSoloBusAgent._matrixGPU, _CoresAgentBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, epsL, nIterL, _nAgent, _nBus);
		break;
	case 16:
		updatePnPGPUSharedResidual< 16> << <numBlocks, _blockSizeSmall >> > (Pn._matrixGPU, PnPre._matrixGPU, PnMoy._matrixGPU, PnTilde._matrixGPU, MuL._matrixGPU, _nAgentByBus._matrixGPU, _rhol, Cost1._matrixGPU, Cost2._matrixGPU,
			Pmin._matrixGPU, Pmax._matrixGPU, Apt1._matrixGPU, Apt2._matrixGPU, Bpt2._matrixGPU, CoresSoloBusAgent._matrixGPU, _CoresAgentBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, epsL, nIterL, _nAgent, _nBus);
		break;
	case  8:
		updatePnPGPUSharedResidual<  8> << <numBlocks, _blockSizeSmall >> > (Pn._matrixGPU, PnPre._matrixGPU, PnMoy._matrixGPU, PnTilde._matrixGPU, MuL._matrixGPU, _nAgentByBus._matrixGPU, _rhol, Cost1._matrixGPU, Cost2._matrixGPU,
			Pmin._matrixGPU, Pmax._matrixGPU, Apt1._matrixGPU, Apt2._matrixGPU, Bpt2._matrixGPU, CoresSoloBusAgent._matrixGPU, _CoresAgentBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, epsL, nIterL, _nAgent, _nBus);
		break;
	case  4:
		updatePnPGPUSharedResidual<  4> << <numBlocks, _blockSizeSmall >> > (Pn._matrixGPU, PnPre._matrixGPU, PnMoy._matrixGPU, PnTilde._matrixGPU, MuL._matrixGPU, _nAgentByBus._matrixGPU, _rhol, Cost1._matrixGPU, Cost2._matrixGPU,
			Pmin._matrixGPU, Pmax._matrixGPU, Apt1._matrixGPU, Apt2._matrixGPU, Bpt2._matrixGPU, CoresSoloBusAgent._matrixGPU, _CoresAgentBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, epsL, nIterL, _nAgent, _nBus);
		break;
	case  2:
		updatePnPGPUSharedResidual<  2> << <numBlocks, _blockSizeSmall >> > (Pn._matrixGPU, PnPre._matrixGPU, PnMoy._matrixGPU, PnTilde._matrixGPU, MuL._matrixGPU, _nAgentByBus._matrixGPU, _rhol, Cost1._matrixGPU, Cost2._matrixGPU,
			Pmin._matrixGPU, Pmax._matrixGPU, Apt1._matrixGPU, Apt2._matrixGPU, Bpt2._matrixGPU, CoresSoloBusAgent._matrixGPU, _CoresAgentBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, epsL, nIterL, _nAgent, _nBus);
		break;
	case  1:
		updatePnPGPUSharedResidual<  1> << <numBlocks, _blockSizeSmall >> > (Pn._matrixGPU, PnPre._matrixGPU, PnMoy._matrixGPU, PnTilde._matrixGPU, MuL._matrixGPU, _nAgentByBus._matrixGPU, _rhol, Cost1._matrixGPU, Cost2._matrixGPU,
			Pmin._matrixGPU, Pmax._matrixGPU, Apt1._matrixGPU, Apt2._matrixGPU, Bpt2._matrixGPU, CoresSoloBusAgent._matrixGPU, _CoresAgentBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, epsL, nIterL, _nAgent, _nBus);
		break;
	}
}





void OPFADMMGPU::updateMu()
{
	updateMUGPU << <_numBlocksM, _blockSize >> > (Mu._matrixGPU, Y._matrixGPU, X._matrixGPU, _rho, _sizeOPFTotal);
	/*tempM.subtract(&X, &Y);
	tempM.multiply(_rho);
	if (consensus) {
		tempM.divide(2);
	}
	Mu.add(&tempM);*/

}


float OPFADMMGPU::getPLoss()
{
	float Ploss = 0;
	for (int i = 1; i < _nAgent; i++) {
		Ploss += Pn.get(i, 0);
	}
	return Ploss;
}

float OPFADMMGPU::getQLoss()
{
	float Qloss = 0;
	for (int i = 1; i < _nAgent; i++) {
		Qloss += Pn.get(i + _nAgent, 0);
	}
	return Qloss;
}

void OPFADMMGPU::ComputePFromAgentToBus()
{
	int numBlock = _nBus;
	switch (_blockSize) {
	case 512:
		ComputePFromAgentToBusGPU<512> << <numBlock, _blockSize >> > (PnTilde._matrixGPU, PnTmin._matrixGPU, PnTmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case 256:
		ComputePFromAgentToBusGPU<256> << <numBlock, _blockSize >> > (PnTilde._matrixGPU, PnTmin._matrixGPU, PnTmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case 128:
		ComputePFromAgentToBusGPU<128> << <numBlock, _blockSize >> > (PnTilde._matrixGPU, PnTmin._matrixGPU, PnTmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case 64:
		ComputePFromAgentToBusGPU< 64> << <numBlock, _blockSize >> > (PnTilde._matrixGPU, PnTmin._matrixGPU, PnTmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case 32:
		ComputePFromAgentToBusGPU< 32> << <numBlock, _blockSize >> > (PnTilde._matrixGPU, PnTmin._matrixGPU, PnTmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case 16:
		ComputePFromAgentToBusGPU< 16> << <numBlock, _blockSize >> > (PnTilde._matrixGPU, PnTmin._matrixGPU, PnTmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case  8:
		ComputePFromAgentToBusGPU<  8> << <numBlock, _blockSize >> > (PnTilde._matrixGPU, PnTmin._matrixGPU, PnTmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case  4:
		ComputePFromAgentToBusGPU<  4> << <numBlock, _blockSize >> > (PnTilde._matrixGPU, PnTmin._matrixGPU, PnTmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case  2:
		ComputePFromAgentToBusGPU<  2> << <numBlock, _blockSize >> > (PnTilde._matrixGPU, PnTmin._matrixGPU, PnTmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	case  1:
		ComputePFromAgentToBusGPU<  1> << <numBlock, _blockSize >> > (PnTilde._matrixGPU, PnTmin._matrixGPU, PnTmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent, _nBus);
		break;
	}
}

void OPFADMMGPU::updateChat()
{
	int numBlock = _nBus;
	switch (_blockSizeSmall) {
	case 512:
		updateChatGPU<512> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _rho, _nBus);
		break;
	case 256:
		updateChatGPU<256> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _rho, _nBus);
		break;
	case 128:
		updateChatGPU<128> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _rho, _nBus);
		break;
	case 64:
		updateChatGPU< 64> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _rho, _nBus);
		break;
	case 32:
		updateChatGPU< 32> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _rho, _nBus);
		break;
	case 16:
		updateChatGPU< 16> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _rho, _nBus);
		break;
	case  8:
		updateChatGPU<  8> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _rho, _nBus);
		break;
	case  4:
		updateChatGPU<  4> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _rho, _nBus);
		break;
	case  2:
		updateChatGPU<  2> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _rho, _nBus);
		break;
	case  1:
		updateChatGPU<  1> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _rho, _nBus);
		break;
	}

	updateBpt2 << < _numBlocksB, _blockSize >> > (Bpt2._matrixGPU, Chat._matrixGPU, _nAgentByBus._matrixGPU, _nBus);


}

void OPFADMMGPU::CommunicationX()
{
	// X = { Pi, Qi, vi, li, pi, qi, vAi, (Pci, Qci, lci) for all child Ci }

	communicateX << <_nBus, _blockSize >> > (X._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nBus);

	//Y      (Pi, Qi, vi, li, pi, qi, vai, Pji, Qji, lji)
	updateQ << <_numBlocksM, _blockSize >> > (Q._matrixGPU, X._matrixGPU, Mu._matrixGPU, _rho, _sizeOPFTotal);
}



float OPFADMMGPU::updateRes(int indice)
{

	float resS = Y.max2(&Ypre);
	float resR = Y.max2(&X);
	float resV = 0;

	float oldrho = _rho;
	resF.set(0, indice, resR);
	resF.set(1, indice, oldrho * resS);
	resF.set(2, indice, resV);

	
	/*std::cout << " Y " << std::endl;
	Y.display(true);
	std::cout << " X " << std::endl;
	X.display(true);*/

	if (_tau > 1) {
		if (resR > _mu * resS) {
			_rho = _tau * _rho;
			Apt2.multiply(_tau);
			Hinv.divide(_tau);
		
			//std::cout << _iterGlobal << "rho augmente " << _rho << std::endl;
		}
		else if (resS > _mu * resR) {// rho = rho / tau_inc;
			_rho = _rho / _tau;
			Apt2.divide(_tau);

			Hinv.multiply(_tau);
			//std::cout << _tau << " " << _mu << std::endl;
			//std::cout << _iterGlobal << "rho diminue " << _rho << std::endl;
		}/**/
	}

	


	return MYMAX(MYMAX(resV, oldrho * resS), resR);
}

int OPFADMMGPU::feasiblePoint()
{
	bool mustTrans = false;
	if (X.getPos()) {
		X.transferCPU();
		_indiceBusBegin.transferCPU();
		mustTrans = true;
	}
	// X  (Pi, Qi, li, vi, pi, qi, vai, Pci ..., Qci... , lci...) !!!!!

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

MatrixCPU OPFADMMGPU::getPb(){
	bool transferToDo = false;
	if(Y.getPos()){
		Y.transferCPU();
		_indiceBusBegin.transferCPU();
		transferToDo = true;
	}
	MatrixCPU Pb(2*_nBus, 1);
	
	for (int i = 0; i <_nBus; i++)
	{
		Pb.set(i,0, Y.get(_indiceBusBegin.get(i, 0) + 4, 0));
		Pb.set(i + _nLine, 0, Y.get(_indiceBusBegin.get(i, 0) + 5, 0));
	}
	if(transferToDo){
		Y.transferGPU();
		_indiceBusBegin.transferGPU();
	}
	return Pb;
}
MatrixCPU OPFADMMGPU::getPhi(){
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
MatrixCPU OPFADMMGPU::getE(){
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



void OPFADMMGPU::display() {

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
	PnTmax.transferCPU();
	PnTmin.transferCPU();
	PnTilde.transferCPU();
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
			<< X.get(begining + 4, 0) << "|" << std::setw(16) << X.get(begining + 5, 0)
			<< "|" << std::setw(16) << Mu.get(begining + 3, 0) << "|" << std::setw(16)
			<< Mu.get(begining + 4, 0) << "|" << std::setw(16) << Mu.get(begining + 5, 0) << "|" << std::endl;

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
	std::cout << "  #  | Mag(pu) | MIN(pu) |  MYMAX(pu)|  P (pu) | Pmin (pu) | Pmax (pu) |  Q (pu)  | Qmin (pu) | Qmax (pu) |" << std::endl;
	std::cout << "-----|---------|---------|---------|---------|-----------|-----------|----------|-----------|-----------|" << std::endl;


	for (int b = 0; b < _nBus; b++) {
		int begining = _indiceBusBegin.get(b, 0);
		int nb = _nAgentByBus.get(b, 0);
		std::cout << std::setw(5) << b << "|" << std::setw(8) << sqrt(Y.get(begining + 3, 0)) << " |" << std::setw(9)
			<< VoltageLimitReal.get(0, b) << "|" << std::setw(9) << VoltageLimitReal.get(1, b)
			<< "|" << std::setw(9) << Y.get(begining + 4, 0) << "|" << std::setw(11)
			<< PnTmin.get(b, 0) * nb << "|" << std::setw(11) << PnTmax.get(b, 0) * nb << "|" << std::setw(10) << Y.get(begining + 5, 0)
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
__global__ void updateChatGPU(float* Chat, float* Y, float* MU, float* nChild, float* Ancestor, float* posChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, float _rho, int nBus) {

	int bus = blockIdx.x;
	int index = threadIdx.x;
	int step = blockDim.x;

	__shared__ float shArr[_blockSizeSmall]; // c'est grand pour pas grand chose...


	int indice = indiceBusBegin[bus];
	int indiceChild = (bus < (nBus - 1)) ? indiceChildBegin[bus] : 0;
	int nb = nChild[bus];
	int Ai = Ancestor[bus];
	int c = posChild[bus];
	float var = 0;

	if (index < 6) {
		//float Phat, Qhat, lhat, phat, qhat;
		var = Y[indice + index] / 2 - MU[indice + index] / (2 * _rho);
		if (bus > 0) {
			if (index < 3) {
				int nAi = nChild[Ai];
				int indiceAncBus = indiceBusBegin[Ai] + 7 + nAi * index + c;
				//var = indiceAncBus;
				var += Y[indiceAncBus] / 2 - MU[indiceAncBus] / (2 * _rho);
			}
		}
	}
	float vhat = 0;
	float muhat = 0;
	for (int i = index; i < nb; i += step) {
		int Bus2 = Childs[indiceChild + i];
		int indiceBusChild = indiceBusBegin[Bus2];
		muhat += MU[indiceBusChild + 6]; // pas du tout coalescent
		vhat += Y[indiceBusChild + 6]; // pas du tout coalescent
	}
	shArr[index] = vhat / (nb + 1) - muhat / (_rho * (nb + 1));
	__syncthreads();
	for (int size = _blockSizeSmall / 2; size > 0; size /= 2) { //uniform
		if (index < size) {
			shArr[index] += shArr[index + size];
		}
		__syncthreads();
	}

	if (index < 6) {
		if (index == 3) {
			var = shArr[0] + Y[indice + 3] / (nb + 1) - MU[indice + 3] / (_rho * (nb + 1)); //shArr[0];
		}
		Chat[index * nBus + bus] = var; // pas coalescent mais bon perdu pour perdu
		// pour p et q, on a un /2 en trop !
	}
}

__global__ void updateBpt2(float* Bpt2, float* Chat, float* nAgentByBus, int nBus) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;


	for (int b = index; b < nBus; b += step) {
		int nA = nAgentByBus[b];

		Bpt2[b] = nA > 0 ? 2 * Chat[b + 4 * nBus] / nA : 0; // �criture coalescente et lecture coalescente
		Bpt2[b + nBus] = nA > 0 ? 2 * Chat[b + 5 * nBus] / nA : 0;

	}

}


template <unsigned int _blockSizeSmall>
__global__ void updatePnPGPUSharedResidual(float* Pn, float* PnPre, float* PnMoy, float* PnTilde, float* MUL, float* nAgentByBus, float _rhol, float* Ap2, float* Cp, float* Pmin,
	float* Pmax, float* Apt1, float* Apt2, float* Bpt2, float* CoresSoloBusAgent, float* CoresBusAgent, float* CoresBusAgentBegin, float eps, int nIterLMax, int nAgent, int nBus) {


	//Definition de toutes les variables locales
	int i = blockIdx.x; // c'est aussi l'identifiant du bus !
	unsigned int thIdx = threadIdx.x;

	// ne change pas

	float Ap2local[2];
	float Ap12local[2];
	float Cplocal[2];
	float Pminlocal[2];
	float Pmaxlocal[2];

	// constant et commun � tous les thread d'un bloc
	__shared__ float Apt1Shared;
	__shared__ float Apt2Shared;
	__shared__ float Apt12Shared;
	__shared__ float Bpt2Shared[2];
	__shared__ int nAgentShared;
	__shared__ float at1Shared;

	// change
	float Pnlocal[2];
	float Pnprelocal[2]; // change

	float bpt, MULOCAL, moy, p;
	float m, r, ub, lb, t;
	// le changement doit �tre partag� par tous les threads du bloc

	__shared__ float MuShared[2];
	__shared__ float PnMoyShared[2];
	__shared__ float PnTildeShared[2];
	__shared__ bool mustContinue;



	__shared__ float shArrP[_blockSizeSmall];
	__shared__ float shArrQ[_blockSizeSmall];


	if (thIdx == 0) {
		Apt1Shared = Apt1[i]; // rho_l *Ni, m�me pour les 2
		Apt2Shared = Apt2[i]; // rho * Ni^2, m�me pour les 2
		Apt12Shared = Apt1Shared + Apt2Shared; // m�me pour les 2
		nAgentShared = nAgentByBus[i];
		at1Shared = _rhol;
		mustContinue = false;
	}


	if (thIdx < 2) {
		Bpt2Shared[thIdx] = Bpt2[i + nBus * thIdx];
		MuShared[thIdx] = MUL[i + nBus * thIdx];
		PnMoyShared[thIdx] = PnMoy[i + nBus * thIdx];
		PnTildeShared[thIdx] = PnTilde[i + nBus * thIdx];
	}
	__syncthreads();

	int iter = 0;
	if (nAgentShared > 0) { // sinon il n'y a rien � faire
		const int CoresAgentLinLocal = CoresBusAgentBegin[i];
		const int j = CoresAgentLinLocal + thIdx;
		//const int endLocal = CoresAgentLinLocal + nAgentShared;
		double res = 0;
		if (nAgentShared == 1) { // cas trivial s'il n'y a qu'un agent, la divergence est entre les blocs donc c'est ok
			if (thIdx == 0) {
				int agent = CoresSoloBusAgent[i];
				// Cplocal et Ap12local, Pmaxlocal, Pminlocal � definir
				Cplocal[0] = Cp[agent];
				Ap2local[0] = Ap2[agent];
				ub = Pmax[agent];
				lb = Pmin[agent];
				r = (Apt2Shared * Bpt2Shared[0] - Cplocal[0]) / (Apt2Shared + Ap2local[0]); //pn = (_rho * Bpt2.get(b, 0) - Cost2.get(n, 0)) / ((_rho + Cost1.get(n, 0)));
				t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
				Pnlocal[0] = t;
				Pnprelocal[0] = t;
				PnMoyShared[0] = t;
				PnTildeShared[0] = t;

				// Q 
				Cplocal[1] = Cp[agent + nAgent];
				Ap2local[1] = Ap2[agent + nAgent];
				ub = Pmax[agent + nAgent];
				lb = Pmin[agent + nAgent];
				r = (Apt2Shared * Bpt2Shared[1] - Cplocal[1]) / (Apt2Shared + Ap2local[1]); //pn = (_rho * Bpt2.get(b, 0) - Cost2.get(n, 0)) / ((_rho + Cost1.get(n, 0)));
				t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
				Pnlocal[1] = t;
				Pnprelocal[1] = t;
				PnMoyShared[1] = t;
				PnTildeShared[1] = t;
			}
		}
		else {
			Pnlocal[0] = 0;
			Pnlocal[1] = 0;
			if (thIdx < nAgentShared)
			{
				int agent = CoresBusAgent[j];
				// P & Q
				Ap2local[0] = Ap2[agent];
				Ap12local[0] = Ap2local[0] + _rhol;
				Cplocal[0] = Cp[agent];
				Pminlocal[0] = Pmin[agent];
				Pmaxlocal[0] = Pmax[agent];
				Pnlocal[0] = Pn[agent];

				Ap2local[1] = Ap2[agent + nAgent];
				Ap12local[1] = Ap2local[1] + _rhol;
				Cplocal[1] = Cp[agent + nAgent];
				Pminlocal[1] = Pmin[agent + nAgent];
				Pmaxlocal[1] = Pmax[agent + nAgent];
				Pnlocal[1] = Pn[agent + nAgent];
			}

			//Calcul des it�rations

			for (iter = 0; iter < nIterLMax; iter++) {
				__syncthreads();
				if (thIdx < nAgentShared) {
					// P
					MULOCAL = MuShared[0];
					moy = PnMoyShared[0];
					p = PnTildeShared[0];

					Pnprelocal[0] = Pnlocal[0];
					m = Pnlocal[0] - moy + p - MULOCAL; // Pn.get(n, 0) - PnMoy.get(bus, 0) + PnTilde.get(bus, 0) - MuL.get(bus, 0);
					r = (m * at1Shared - Cplocal[0]) / Ap12local[0]; // pn = (Bp1.get(n, 0) * _rhol - Cost2.get(n, 0)) / Ap12.get(n, 0);
					ub = Pmaxlocal[0];
					lb = Pminlocal[0];
					t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
					Pnlocal[0] = t;
					
					res = (double) t - Pnprelocal[0];
					res = res * res;
					if (res > eps) {
						mustContinue = true; // pas de race condition, car l'ordre n'importe pas,
					}


					// Q
					MULOCAL = MuShared[1];
					moy = PnMoyShared[1];
					p = PnTildeShared[1];

					Pnprelocal[1] = Pnlocal[1];
					m = Pnlocal[1] - moy + p - MULOCAL; // Pn.get(n + _nAgent, 0) - PnMoy.get(bus + _nBus, 0) + PnTilde.get(bus + _nBus, 0) - MuL.get(bus + _nBus, 0);
					r = (m * at1Shared - Cplocal[1]) / Ap12local[1]; // pn = (Bp1.get(n + _nAgent, 0) * _rhol - Cost2.get(n + _nAgent, 0)) / Ap12.get(n+ _nAgent, 0);
					ub = Pmaxlocal[1];
					lb = Pminlocal[1];
					t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
					Pnlocal[1] = t;
					res = (double)t - Pnprelocal[1];
					res = res * res;
					if (res > eps) {
						mustContinue = true; // pas de race condition, car l'ordre n'importe pas,
					}

				}

				shArrP[thIdx] = Pnlocal[0];
				shArrQ[thIdx] = Pnlocal[1];


				__syncthreads();
				if (_blockSizeSmall >= 512) {
					if (thIdx < 256) {
						shArrP[thIdx] += shArrP[thIdx + 256];
						shArrQ[thIdx] += shArrQ[thIdx + 256];
					}
					__syncthreads();
				}
				if (_blockSizeSmall >= 256) {
					if (thIdx < 128) {
						shArrP[thIdx] += shArrP[thIdx + 128];
						shArrQ[thIdx] += shArrQ[thIdx + 128];
					}
					__syncthreads();
				}
				if (_blockSizeSmall >= 128) {
					if (thIdx < 64) {
						shArrP[thIdx] += shArrP[thIdx + 64];
						shArrQ[thIdx] += shArrQ[thIdx + 64];
					}
					__syncthreads();
				}
				if (_blockSizeSmall >= 64) {
					if (thIdx < 32) {
						warpReduce<_blockSizeSmall>(shArrP, thIdx);
						warpReduce<_blockSizeSmall>(shArrQ, thIdx);
					}
				}else if (_blockSizeSmall >= 32) {
					warpReduce<_blockSizeSmall>(shArrP, thIdx);
					warpReduce<_blockSizeSmall>(shArrQ, thIdx);
				}
				__syncthreads();

				if (thIdx == 0) {
					// P
					moy = shArrP[0] / nAgentShared;
					PnMoyShared[0] = moy;
					bpt = moy + MuShared[0]; //Bpt1.set(b, 0, MuL.get(b, 0) + PnMoy.get(b, 0));
					p = (Apt1Shared * bpt + Apt2Shared * Bpt2Shared[0]) / Apt12Shared; //pn = (Bpt1.get(b, 0) * Apt1.get(b, 0) + Bpt2.get(b, 0) * Apt2.get(b, 0)) / Apt12.get(b, 0);
					PnTildeShared[0] = p;
					res = p - moy;
					res = res * res;
					if (res > eps) {
						mustContinue = true;
					}
					MuShared[0] = MuShared[0] + moy - p; // mu = MuL.get(b, 0) + PnMoy.get(b, 0) - PnTilde.get(b, 0);
					
					// Q
					moy = shArrQ[0] / nAgentShared;
					PnMoyShared[1] = moy;
					bpt = moy + MuShared[1]; //Bpt1.set(b, 0, MuL.get(b, 0) + PnMoy.get(b, 0));
					p = (Apt1Shared * bpt + Apt2Shared * Bpt2Shared[1]) / Apt12Shared; //pn = (Bpt1.get(b + _nBus, 0) * Apt1.get(b + _nBus, 0) + Bpt2.get(b + _nBus, 0) * Apt2.get(b + _nBus, 0)) / Apt12.get(b + _nBus, 0);

					PnTildeShared[1] = p;
					res = p - moy;
					res = res * res;
					if (res > eps) {
						mustContinue = true;
					}
					MuShared[1] = MuShared[1] + moy - p; // mu = MuL.get(b, 0) + PnMoy.get(b, 0) - PnTilde.get(b, 0);

				}
				__syncthreads();
				if (!mustContinue) {
					break;
				}
				else {
					__syncthreads();
					if (thIdx == 0) {
						mustContinue = false;
					}
				}
			}
		}
		//Ecriture des it�rations
		__syncthreads();

		if (thIdx < nAgentShared)
		{
			int agent = CoresBusAgent[j];

			Pn[agent] = Pnlocal[0];
			PnPre[agent] = Pnprelocal[0];

			Pn[agent + nAgent] = Pnlocal[1];
			PnPre[agent + nAgent] = Pnprelocal[1];

		}
		if (thIdx == 0) {
			PnMoy[blockIdx.x] = PnMoyShared[0];// TMoyShared;
			PnTilde[blockIdx.x] = PnTildeShared[0];// PShared;
			MUL[blockIdx.x] = MuShared[0];// MuShared;

			PnMoy[blockIdx.x + nBus] = PnMoyShared[1];// TMoyShared;
			PnTilde[blockIdx.x + nBus] = PnTildeShared[1];// PShared;
			MUL[blockIdx.x + nBus] = MuShared[1];// MuShared;
		}
	}
}

