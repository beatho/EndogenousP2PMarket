#include "../head/MarketEndoDirectGPU.cuh"
#define MAX(X, Y) X * (X >= Y) + Y * (Y > X)

// get loss
// comuunication loss
// feasible point
// init 




MarketEndoDirectGPU::MarketEndoDirectGPU() : Method()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " MarketEndoDirectGPU Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1ab , Fb2, Fb3, Fb4, Fb5,FB6, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu
}


MarketEndoDirectGPU::MarketEndoDirectGPU(float rho) : Method()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default MarketEndoDirectGPU Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb11, FB12, Fb2, Fb3, Fb4, Fb5, FB6, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu
}

MarketEndoDirectGPU::~MarketEndoDirectGPU()
{
	
}
void MarketEndoDirectGPU::setParam(float rho)
{
	_rho = rho;
}

bool MarketEndoDirectGPU::chekcase()
{
	if (_nBus != (_nLine + 1)) {
		std::cout << "wrong number of line " << _nLine << "against " << _nBus << std::endl;
		return false;
	}
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
		ZsRe.display();
		ZsIm.display();
		return false;
	}

	return true;
}

void MarketEndoDirectGPU::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
		CHECK_LAST_CUDA_ERROR();
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 0, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION
	}
	
	_iterG = sim.getIterG();
	_nIterL = sim.getIterL();
	_stepG = sim.getStepG();
	
	
	float epsG = sim.getEpsG();
	_epsL = sim.getEpsL() / 20;

	_epsL *= _epsL;
	
	float fc = 0;
	float resG = 2 * epsG;
	_iterGlobal = 0;
	
	
	while (((_iterGlobal < _iterG) && (resG>epsG)) || (_iterGlobal <= _stepG)) {
		/*std::cout << "---------------------------------" << std::endl;
		std::cout << " X avant" << std::endl;
		X.display(true);
		
		//LAMBDALin.display();
		//Bt1.display();
		//TradeLin.display();	
		
		std::cout << " Q "   << std::endl;
		Q.display(true);
		std::cout << " Y "   << std::endl;
		Y.display(true);
		std::cout << " Mu "   << std::endl;
		Mu.display(true);
		std::cout << " Chat "   << std::endl;
		Chat.display(true);
		std::cout << " Bp2 " << std::endl;
		Bp2.display(true);
		std::cout << " P " << std::endl;
		P.display(true);*/
		//Pn.saveCSVForce("TestPnGPU2.csv", 11, 1);
		//X.saveCSVForce("TestXGPU2.csv", 11, 1);
		//Y.saveCSVForce("TestYGPU2.csv", 11, 1);
		//Chat.saveCSVForce("TestChatGPU2.csv", 11, 1);
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		updatePMarket(); // puissance et trade (ind�pendament du bus, m�me si on aurait pu r�soudre bus par bus)
		//std::cout << " P " << std::endl;
		//P.display(true);
		
		
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 1, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

		//updateXWOCurrent(); // flux dans le r�seau, tension
		updateXWOCurrentCPU();
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 2, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		//Q.saveCSVForce("TestQGPU.csv", 11, 1);
		CommunicationX();
	
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 3, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		//Y.saveCSVForce("TestYGPU.csv", 11, 1);
		updateGlobalProb();
		
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 4, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		
		updateMu();
		
		updateChat();

#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 5, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		// FB 4
		if (!(_iterGlobal % _stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateRes(_iterGlobal / _stepG);
			//std::cout << _iterGlobal << " " << resF.get(0, _iterGlobal / _stepG) << " " << resF.get(1, _iterGlobal / _stepG) << " " << resF.get(2, _iterGlobal / _stepG) << std::endl;
			//resG = 1;
#ifdef INSTRUMENTATION
			cudaDeviceSynchronize();
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 6, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		}
		//std::cout << iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;
		_iterGlobal++;
	}
	//std::cout << _iterGlobal << " " << resF.get(0, (_iterGlobal - 1) / _stepG) << " " << resF.get(1, (_iterGlobal - 1) / _stepG) << " " << resF.get(2, (_iterGlobal - 1) / _stepG) << std::endl;


#ifdef INSTRUMENTATION	
	occurencePerBlock.increment(0, 1, _iterGlobal);
	occurencePerBlock.increment(0, 2, _iterGlobal);
	occurencePerBlock.increment(0, 3, _iterGlobal);
	occurencePerBlock.increment(0, 4, _iterGlobal);
	occurencePerBlock.increment(0, 5, _iterGlobal);
	occurencePerBlock.increment(0, 6, _iterGlobal / _stepG);

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	

	//setPnFromX << < _nBus, _blockSizeSmall >> > (Pn._matrixGPU, X._matrixGPU, _indiceBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgent);


	MatrixCPU tradeLinCPU;
	TradeLin.toMatCPU(tradeLinCPU);
	MatrixCPU LAMBDALinCPU;
	LAMBDALin.toMatCPU(LAMBDALinCPU);
	MatrixCPU PnCPU;
	Pn.toMatCPU(PnCPU);
	

	int indice = 0;
	for (int idAgent = 0; idAgent < _nAgentTrue; idAgent++) {
		MatrixCPU omega(cas.getVoisin(idAgent));
		int Nvoisinmax = nVoisinCPU.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			int idVoisin = omega.get(voisin, 0);
			trade.set(idAgent, idVoisin, tradeLinCPU.get(indice, 0));
			LAMBDA.set(idAgent, idVoisin, LAMBDALinCPU.get(indice, 0));
			indice = indice + 1;
		}
	}
	for (int idAgent = _nAgentTrue; idAgent < _nAgent; idAgent++) {
		for (int idVoisin = 0; idVoisin < _nAgentTrue; idVoisin++) {
			if (idVoisin != (idAgent - _nAgentTrue)) {
				trade.set(idAgent, idVoisin, tradeLinCPU.get(indice, 0));
				LAMBDA.set(idAgent, idVoisin, LAMBDALinCPU.get(indice, 0));
				indice = indice + 1;
			}

		}
	}
	
	
	/*PnCPU.set(0, 0, getPLoss());
	PnCPU.set(_nAgentTrue, 0, getQLoss());
	*/

	//Ct.display();
	//Tlocal.display();

	fc = calcFc(&Cost1, &Cost2, &Tlocal, &Pn, &Ct, &tempN1, &tempNN);
	// FB 5
	
	result->setResF(&resF);
	result->setIter(_iterGlobal);
	result->setPn(&PnCPU);
	result->setFc(fc);
	result->setLAMBDA(&LAMBDA);
	result->setTrade(&trade);

#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();  
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 7, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 7, 1);

	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
	tall = clock() - tall;
	timeMarketEndo = tall;
	//display();
	result->setTime((float)tall / CLOCKS_PER_SEC);
	
}

void MarketEndoDirectGPU::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    
	Pmin = cas.getPmin();
	Pmax = cas.getPmax();
	Cost2 = cas.getb();
	Cp = cas.getb();


	MatrixCPU Lb(cas.getLb());
	MatrixCPU Ub(cas.getUb());
	matLb.transferCPU();
	matUb.transferCPU();

	int indice = 0;

	for (int idAgent = 0; idAgent < _nAgent; idAgent++) {
		int Nvoisinmax = nVoisinCPU.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			matLb.set(indice, 0, Lb.get(idAgent, 0));
			matUb.set(indice, 0, Ub.get(idAgent, 0));
			indice = indice + 1;
		}
	}
	matLb.transferGPU();
	matUb.transferGPU();
	
	float Ploss = Pn.get(0, 0, false);
	float Qloss = Pn.get(_nAgentTrue, 0, false);
	// pour essayer que cela marche
	Pn.add(&Pmin, &Pmax);
	Pn.divide(2);

	// unleash powe
	Pmin.set(0, 0, -POWERLIMIT, true);
	Pmax.set(_nAgentTrue, 0, POWERLIMIT, true);
	Pmin.set(_nAgentTrue, 0, -POWERLIMIT, true);

	
	Pn.set(0, 0, Ploss, true);
	Pn.set(_nAgentTrue, 0, Qloss, true);

	Pb.set(0.0);
	Pbmax.set(0.0);
	Pbmin.set(0.0);
	ComputePFromAgentToBus();

	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);
	Cp.multiplyT(&nVoisin);

	initPQAgent << < _nBus, _blockSizeSmall >> > (X._matrixGPU, _indiceBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, Pn._matrixGPU, _nAgent);
	CHECK_LAST_CUDA_ERROR();
	
	

	initDFSPQ << <1, _nBus, _nBus* (8*sizeof(bool) + sizeof(int)) >> > (X._matrixGPU, Pb._matrixGPU, nChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nBus);
	CHECK_LAST_CUDA_ERROR();
	
	computeLoss();
	CHECK_LAST_CUDA_ERROR();
	
	CommunicationX();
	CHECK_LAST_CUDA_ERROR();
	//Y.set(&X);
	
	updateChat();
	CHECK_LAST_CUDA_ERROR();
	
	//Cp.display(true);

#ifdef INSTRUMENTATION
	cudaDeviceSynchronize(); 
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION

}

void MarketEndoDirectGPU::init(const Simparam& sim, const StudyCase& cas)
{
	 

	if (_CoresChatBegin.getPos()) {

		//_CoresChatBegin.transferCPU();
		//_indiceBusBegin.transferCPU();

		Ancestor.transferCPU();
		PosChild.transferCPU();

		_indiceChildBegin.transferCPU();
		Childs.transferCPU();

		ZsIm.transferCPU();
		ZsRe.transferCPU();
	}
	


	// intitilisation des matrixs et variables 
	
	clock_t t = clock();
	//std::cout << "init " << std::endl;
	_rho = sim.getRho();
	
	initMarket(sim, cas);
	 

	_iterG = sim.getIterG();
	_stepG = sim.getStepG();
	
	float epsG = sim.getEpsG();
	float epsGC = sim.getEpsGC();
	_ratioEps = epsG / epsGC;
	
	_nBus = cas.getNBus();
	_nBusWLoss = _nBus + 1;
	_nLine = cas.getNLine(true); // ne doit pas �tre r�duit ici !!!

	_debutloss =  3 * _nLine + 5 * _nBus + 2 * (_nAgentTrue - 1); // L = nChild.sum()
	_sizeEndoMarketTotal = _debutloss;
	_sizeChat = 4 * _nBus;
	
	if (losstype == LossType::CURRENT) {
		_sizeEndoMarketTotal += (_nBus + 2); // pertes et courants sauf premier bus ou + 2
	}
	else if (losstype == LossType::POWER) {
		_sizeEndoMarketTotal += _nAgent;
	}
	
	_numBlocksB = ceil((_nBus + _blockSize - 1) / _blockSize);
	_numBlocksH = ceil((_sizeEndoMarketTotal + _blockSize - 1) / _blockSize);
	
	
	_numLineByBlockY = _sizeEndoMarketTotal / 100 + 1; // on veut maimum de 100 blocks !!!


	//std::cout << _nAgentTrue << " " << _nBus << " " << _nLine << " " << _sizeEndoMarketTotal << std::endl;

	_CoresAgentBus = MatrixGPU(cas.getCoresAgentBusLin(), 1);
	_CoresAgentBusBegin = MatrixGPU(cas.getCoresAgentBusLinBegin(), 1);
	_nAgentByBus = MatrixGPU(cas.getNagentByBus(), 1);
	_nAgentByBusCPU = cas.getNagentByBus();
	PosAgent = MatrixGPU(_nAgentTrue, 1, 0, 1);

	 
	nChildCPU = MatrixCPU(_nBus, 1);

	CoresLineBusCPU = cas.getCoresLineBus(true);
	CoresLineBus = MatrixGPU(CoresLineBusCPU, 1);

	_CoresBusAgent = MatrixGPU(cas.getCoresBusAgentLin(), 1); // Cores[n] = b

	Ancestor = MatrixGPU(_nBus, 1, 0); // A_i = bus ant�c�dent de i
	PosChild = MatrixGPU(_nBus, 1, 0); // indice du bus i dans Child[Ai]
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
		int busFrom = CoresLineBusCPU.get(lold, 0);
		Ancestor.set(busTo, 0, busFrom);
		nChildCPU.increment(busFrom, 0, 1);
		ZsNorm.set(lold, 0, ZsRe.get(lold, 0) * ZsRe.get(lold, 0) + ZsIm.get(lold, 0) * ZsIm.get(lold, 0));
	}
	nChild = MatrixGPU(nChildCPU, 1);


	_rhoInv = 1 / _rho;
	resF = MatrixCPU(3, (_iterG / _stepG) + 1);

	
	
	
	//std::cout << " local resolution " << std::endl;
	// local resolution
	tempN2 = MatrixGPU(_nAgent, 1, 0, 1);
	tempB2 = MatrixGPU(2 * _nBus, 1, 0, 1);
	
	Pbmax = MatrixGPU(2 * _nBus, 1, 0, 1);
	Pbmin = MatrixGPU(2 * _nBus, 1, 0, 1);
	Pb = MatrixGPU(2 * _nBus, 1, 0, 1);
	CoresSoloBusAgent = MatrixGPU(_nBus, 1, -1, 1);
	
	_nAgentByBusCPU.increment(0, 0, -1);
	_nAgentOn0 = _nAgentByBusCPU.get(0, 0);
	
	removeLossAgent << <1, 1 >> > (_nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU);
	initPosAgent << <_nBus, _blockSizeSmall >> > (PosAgent._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU);

	//cudaDeviceSynchronize();
	//CHECK_LAST_CUDA_ERROR();
	//Pb.saveCSVForce("TestPbGPU2.csv", 11, 1);
	ComputePFromAgentToBus();
	//Pb.saveCSVForce("TestPbGPU2.csv", 11, 1);
	//cudaDeviceSynchronize();
	//CHECK_LAST_CUDA_ERROR();
	//_nAgentByBus.display();

	//std::cout << " creation " << std::endl;
	X = MatrixGPU(_sizeEndoMarketTotal, 1, 0, 1);
	Ypre = MatrixGPU(_sizeEndoMarketTotal, 1, 0, 1);
	Y = MatrixGPU(_sizeEndoMarketTotal, 1, 0, 1);
	Y.preallocateReduction();

	Mu = MatrixGPU(_sizeEndoMarketTotal, 1, 0, 1);
	
	tempN1 = MatrixGPU(_nAgent, 1, 0, 1);
	tempNN = MatrixGPU(_nTrade, 1, 0, 1);
	//tempM1 = new MatrixCPU[_nAgent];
	tempM = MatrixGPU(_sizeEndoMarketTotal, 1, 0, 1);
	
	sizeMarketEndoDirectGPU = MatrixGPU(_nBusWLoss, 1, 0, 1);
	sizeMarketEndoDirectGPU.preallocateReduction();
	sizeMarketEndoDirectGPUBig = MatrixGPU(_sizeEndoMarketTotal, 1, 0, 1);

	indiceBusBeginCPU = MatrixCPU(_nBusWLoss, 1);
	_indiceBusBeginBig = MatrixGPU(_sizeEndoMarketTotal, 1, 0, 1);
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
		debutChat += 4;
	}
	indiceBusBeginCPU.set(_nBus, 0, debut);
	CoresChatBeginCPU.set(_nBus, 0, debutChat);

	
	_CoresChatBegin = MatrixGPU(CoresChatBeginCPU, 1);
	_indiceBusBegin = MatrixGPU(indiceBusBeginCPU, 1);
	defineSizeBig << <_nBusWLoss, _blockSize >> > (sizeMarketEndoDirectGPUBig._matrixGPU, nChild._matrixGPU, _indiceBusBegin._matrixGPU, sizeMarketEndoDirectGPU._matrixGPU, _indiceBusBeginBig._matrixGPU, _nAgentByBus._matrixGPU, losstype, _nBus, _nAgentTrue);
	 
	
	_sizeEndoMarketMax = sizeMarketEndoDirectGPU.max2();
	Hinv = MatrixGPU(_sizeEndoMarketTotal, _sizeEndoMarketMax, 0, 1);
	Q = MatrixGPU(_sizeEndoMarketTotal, 1, 0, 1);
	
	Childs = MatrixGPU(_nLine, 1);
	PosChild = MatrixGPU(_nBus, 1, -1);

	Chat = MatrixGPU(_sizeChat, 1, 0, 1);


	MatrixGPU lowerBound(cas.getLowerBound(), 1); //voltage angle, voltage, line...
	MatrixGPU upperBound(cas.getUpperBound(), 1); //voltage angle, voltage, line...
	VoltageLimit = MatrixGPU( 2, _nBus, 0, 1); // min, max
	VoltageLimitReal = MatrixGPU( 2, _nBus, 0, 1); // min, max
	
	initVoltageBound << < _numBlocksB, _blockSize >> > (VoltageLimitReal._matrixGPU, VoltageLimit._matrixGPU, lowerBound._matrixGPU, upperBound._matrixGPU, nChild._matrixGPU, _nBus);
	VoltageLimit.toMatCPU(VoltageLimitCPU);

	MatrixCPU nChildTemp(_nBus, 1, 0);
	
	int debutChild = 0;
	
	_indiceChildBegin = MatrixGPU(_nBus, 1);
	
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
	
	//std::cout << " Hinv " << std::endl;
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
	float sumP = 0;
	float sumQ = 0;
	int sizeA = 0;
	MatrixGPU A;
	switch (losstype)
	{
	case LossType::POWER:
		sizeA = _nAgent;
		A = MatrixGPU(sizeA, sizeA);
		A.setEyes(-1);
		A.set(0, 0, 0);
		A.set(_nAgentTrue, _nAgentTrue, 0);
		for (int i = 1; i < _nAgentTrue; i++) {
			A.set(0, i, 1); // sum(p) + Ploss = 0
			A.set(_nAgentTrue, i + _nAgentTrue, 1); // Qloss + sum(q) = 0
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
	
	//Hinv.display(true);
	_indiceChildBegin.transferGPU();
	Childs.transferGPU();
	ZsIm.transferGPU();
	ZsRe.transferGPU();
	
	
	//std::cout << " init valeur " << std::endl;
	 
	initPQAgentV << < _nBus, _blockSizeSmall >> > (X._matrixGPU, _indiceBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, Pn._matrixGPU, _nAgentTrue);
	 
	
	
	initDFSPQ << <1, _nBus, _nBus* (8*sizeof(bool) + sizeof(int)) >> > (X._matrixGPU, Pb._matrixGPU, nChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nBus);
	 
	communicateX << <_nBusWLoss, _blockSize >> > (X._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nAgentByBus._matrixGPU, _CoresBusAgent._matrixGPU, PosAgent._matrixGPU, losstype, _nBus, _nAgentTrue);
	 
	computeLoss();
	 
	
	/*
	X[_nBus].set(0, 0, getPLoss());
	X[_nBus].set(1, 0, getQLoss());*/
	
	Y.set(&X);

	
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
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	//std::cout << "---------------------------------------------------------------------------------------" << std::endl;
}

void MarketEndoDirectGPU::initMarket(const Simparam& sim, const StudyCase& cas)
{
	if (CoresMatLin.getPos()) { // une copie en trop mais pour l'instant c'est ok...
		CoresMatLin.transferCPU();

		CoresLinAgent.transferCPU();
		CoresAgentLin.transferCPU();
		CoresLinVoisin.transferCPU();
		CoresLinTrans.transferCPU();

		Tlocal_pre.transferCPU();
		TradeLin.transferCPU();
		LAMBDALin.transferCPU();

		matLb.transferCPU();
		matUb.transferCPU();
		Ct.transferCPU();

		Pn.transferCPU();
		Pmin.transferCPU();
		Pmax.transferCPU();
	}


	_nAgentTrue = cas.getNagent();
	_nAgent = 2 * _nAgentTrue;
	
	if (_rhol == 0) {
		_rhol = _rho;
	}
	nVoisinCPU = cas.getNvoi();
	nVoisin = MatrixGPU(nVoisinCPU, 1);
	nVoisin.preallocateReduction();

	_nTrade = nVoisin.sum();
	_nTradeP = nVoisin.sum(0, _nAgentTrue);
	//std::cout << "nTrade " << _nTrade << " " << _nTradeP << std::endl;
	
	_nTradeQ = _nTrade - _nTradeP;
	if (_nTradeQ != (_nAgentTrue * (_nAgentTrue - 1))) {
		std::cout << "err MarketEndoDirectGPU : " << _nAgent << " " << _nAgentTrue << " " << _nTrade << " " << _nTradeP << " " << _nTradeQ << std::endl;
		throw std::invalid_argument("Agent must be fully conected for the Q echanges, WIP");
	}
	int nVoisinMax = nVoisin.max2();
	if (_blockSize * NMAXPEERPERTRHREAD < nVoisinMax) {
		std::cout << _blockSize << " " << NMAXPEERPERTRHREAD << " " << nVoisinMax << std::endl;
		throw std::invalid_argument("For this Method, an agent must not have more than 5120 peers");
	}
	_numBlocksN = ceil((_nAgent + _blockSize - 1) / _blockSize);
	_numBlocksM = ceil((_nTrade + _blockSize - 1) / _blockSize);
	
	
	
	
	if (initWithMarketClear) {
		ADMMMarketGPU market;
		Simparam res(sim);
		market.solve(&res, sim, cas);
		//res.display();
		LAMBDA = res.getLambda();
		trade = res.getTrade();
		Pn = res.getPn();
		//Pn.display();
		//std::cout << "****" << std::endl;
		
	}
	else {
		LAMBDA = sim.getLambda();
		trade = sim.getTrade();
		Pn = sim.getPn(); // somme des trades
	}
	Pmin = cas.getPmin();
	Pmax = cas.getPmax();
	 // unleash powe
	Pmin.set(0, 0, -POWERLIMIT);
	Pmax.set(_nAgentTrue, 0, POWERLIMIT);
	Pmin.set(_nAgentTrue, 0, -POWERLIMIT);
	
	Pn.transferGPU();
	Pn.preallocateReduction();
	Pmin.transferGPU();
	Pmax.transferGPU();

	if (Pn.max2() == 0) {
		Pn.add(&Pmin, &Pmax);
		Pn.divide(2);
		Pn.set(0, 0, 0, true);
	}
	/*Pn.display();
	LAMBDA.display();
	trade.display();*/
	 
	//CHECK_LAST_CUDA_ERROR();

	_at1 = _rho; // car that apparait 2 fois 
	_at2 = _rhol;

	MatrixCPU BETA(cas.getBeta());
	MatrixCPU Ub(cas.getUb());
	MatrixCPU Lb(cas.getLb());

	/*if (Ub.get(_nAgentTrue, 0) == 0) { // unleash power
		Ub.set(_nAgentTrue, 0, POWERLIMIT);
		Lb.set(_nAgentTrue, 0, -POWERLIMIT);
	}*/


	
	
	//std::cout << "mise sous forme lin�aire" << std::endl;
	


	CoresMatLin = MatrixGPU(_nAgent, _nAgentTrue, -1);
	CoresAgentLin = MatrixGPU(_nAgent + 1, 1);
	CoresLinAgent = MatrixGPU(_nTrade, 1);
	CoresLinVoisin = MatrixGPU(_nTrade, 1);
	CoresLinTrans = MatrixGPU(_nTrade, 1);

	Tlocal_pre = MatrixGPU(_nTrade, 1);
	TradeLin = MatrixGPU(_nTrade, 1);
	LAMBDALin = MatrixGPU(_nTrade, 1);

	matLb = MatrixGPU(_nTrade, 1);
	matUb = MatrixGPU(_nTrade, 1);
	Ct = MatrixGPU(_nTrade, 1);
	

	int indice = 0;
	//std::cout << " P " << std::endl;
	for (int idAgent = 0; idAgent < _nAgentTrue; idAgent++) { // P
		MatrixCPU omega(cas.getVoisin(idAgent));
		int Nvoisinmax = nVoisinCPU.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			int idVoisin = omega.get(voisin, 0);
			matLb.set(indice, 0, Lb.get(idAgent, 0));
			matUb.set(indice, 0, Ub.get(idAgent, 0));
			Ct.set(indice, 0, BETA.get(idAgent, idVoisin));
			TradeLin.set(indice, 0, trade.get(idAgent, idVoisin));
			Tlocal_pre.set(indice, 0, trade.get(idAgent, idVoisin));
			LAMBDALin.set(indice, 0, LAMBDA.get(idAgent, idVoisin));
			CoresLinAgent.set(indice, 0, idAgent);
			CoresLinVoisin.set(indice, 0, idVoisin);
			CoresMatLin.set(idAgent, idVoisin, indice);
			indice = indice + 1;
		}
		CoresAgentLin.set(idAgent + 1, 0, indice);
	}
	//std::cout << " Q " << std::endl;
	for (int idAgent = _nAgentTrue; idAgent < _nAgent; idAgent++) { // Q
		for (int idVoisin = 0; idVoisin < _nAgentTrue; idVoisin++) {
			if (idVoisin != (idAgent - _nAgentTrue)) {
				matLb.set(indice, 0, Lb.get(idAgent, 0));
				matUb.set(indice, 0, Ub.get(idAgent, 0));
				TradeLin.set(indice, 0, trade.get(idAgent, idVoisin));
				Tlocal_pre.set(indice, 0, trade.get(idAgent, idVoisin));
				LAMBDALin.set(indice, 0, LAMBDA.get(idAgent, idVoisin));
				CoresLinAgent.set(indice, 0, idAgent);
				CoresLinVoisin.set(indice, 0, idVoisin + _nAgentTrue);
				CoresMatLin.set(idAgent, idVoisin, indice);
				indice = indice + 1;
			}
		}

		CoresAgentLin.set(idAgent + 1, 0, indice);
	}
	for (int lin = 0; lin < _nTrade; lin++) {
		int i = CoresLinAgent.get(lin, 0);
		int j = CoresLinVoisin.get(lin, 0);
		if (lin >= _nTradeP) {
			i -= _nAgentTrue;
		}

		int k = CoresMatLin.get(j, i);
		CoresLinTrans.set(lin, 0, k);
	}
	// transfert des mises lineaire
	matUb.transferGPU();
	matLb.transferGPU();
	Ct.transferGPU();

	Tlocal_pre.transferGPU();
	TradeLin.transferGPU();
	LAMBDALin.transferGPU();

	CoresAgentLin.transferGPU();
	CoresLinAgent.transferGPU();
	CoresLinVoisin.transferGPU();
	CoresMatLin.transferGPU();
	CoresLinTrans.transferGPU();

	//std::cout << "autres donn�e sur CPU" << std::endl;
	tempNN = MatrixGPU(_nTrade, 1, 0, 1);
	tempNN.preallocateReduction();
	tempN1 = MatrixGPU(_nAgent, 1, 0, 1); // plut�t que de re-allouer de la m�moire � chaque utilisation
	//MatrixCPU temp1N(1, _nAgent, 0, 1);

	/**/

	

	P = Pn; // moyenne des trades, ici c'est juste pour qu'il ait la m�me taille sans avoir besoin de se poser de question
	P.divideT(&nVoisin);
	Tlocal = MatrixGPU(_nTrade, 1, 0, 1);
	Tlocal.preallocateReduction();

	Cost1 = MatrixGPU(cas.geta(), 1);
	Cost2 = MatrixGPU(cas.getb(), 1);
	Ap1 = nVoisin;
	Ap2 = nVoisin;
	Ap3 = Cost1;	
	Ap123 = MatrixGPU(_nAgent, 1, 0, 1);

	
	
	Cp = Cost2;
	Bp1 = MatrixGPU(_nAgent, 1, 0, 1);
	Bp2 = MatrixGPU(_nAgent, 1, 0, 1);
	Bt1 = MatrixGPU(_nTrade, 1, 0, 1);
	Bt2 = MatrixGPU(_nTrade, 1, 0, 1);

	MU = MatrixGPU(sim.getMU(), 1); // facteur reduit i.e lambda_l/_rho
	Tmoy = P;


	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);

	/*std::cout << "Power bound" << std::endl;
	Pmin.display();
	Pmax.display();*/


	Ap1.multiply(_rhol);
	Ap2.multiply( _rho);
	Ap2.multiplyT(&nVoisin);
	Cp.multiplyT(&nVoisin);
	Tmoy.divideT(&nVoisin);

	Ap3.multiplyT(&nVoisin);
	Ap3.multiplyT(&nVoisin);
	Ap123.add(&Ap1, &Ap2);
	Ap123.add(&Ap3);

	//CoresLinTrans.display();
	//CoresAgentLin.display();
	/*std::cout << _at1 << " " << _at2 << std::endl;
	
	Ct.display(true);
	Ap1.display(true);
	Ap2.display(true);
	Ap3.display(true);
	Ap123.display(true);
	Cp.display(true);

	Pmin.display(true);
	Pmax.display(true);
	matLb.display(true);
	matUb.display(true);
	std::cout << "fin init market" << std::endl;*/
	
	//CHECK_LAST_CUDA_ERROR();
}




void MarketEndoDirectGPU::updateGlobalProb() {
	
	
	Ypre.swap(&Y);
	int numBlock = _sizeEndoMarketTotal/ _numLineByBlockY;
	switch (_blockSize) {
	case 512:
		updateY<512> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeMarketEndoDirectGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeEndoMarketMax, _sizeEndoMarketTotal);
		break;
	case 256:
		updateY<256> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeMarketEndoDirectGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeEndoMarketMax, _sizeEndoMarketTotal);
		break;
	case 128:
		updateY<128> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeMarketEndoDirectGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeEndoMarketMax, _sizeEndoMarketTotal);
		break;
	case 64:
		updateY< 64> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeMarketEndoDirectGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeEndoMarketMax, _sizeEndoMarketTotal);
		break;
	case 32:
		updateY< 32> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeMarketEndoDirectGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeEndoMarketMax, _sizeEndoMarketTotal);
		break;
	case 16:
		updateY< 16> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeMarketEndoDirectGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeEndoMarketMax, _sizeEndoMarketTotal);
		break;
	case  8:
		updateY<  8> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeMarketEndoDirectGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeEndoMarketMax, _sizeEndoMarketTotal);
		break;
	case  4:
		updateY<  4> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeMarketEndoDirectGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeEndoMarketMax, _sizeEndoMarketTotal);
		break;
	case  2:
		updateY<  2> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeMarketEndoDirectGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeEndoMarketMax, _sizeEndoMarketTotal);
		break;
	case  1:
		updateY<  1> << <numBlock, _blockSize >> > (Y._matrixGPU, Hinv._matrixGPU, Q._matrixGPU, sizeMarketEndoDirectGPUBig._matrixGPU, _indiceBusBeginBig._matrixGPU, _sizeEndoMarketMax, _sizeEndoMarketTotal);
		break;
	}

	Y.set(3, 0, 1, 1);
	Y.set(4 + 2 * _nAgentOn0, 0, 1, 1);


	// communication of y, mu

}

void MarketEndoDirectGPU::updateX()
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

			gamma = k2 * x4 - (x1 * x1 + x2 * x2) / x3; // ce n'est pas vraiment gamma, doit �tre positif
			//std::cout << "x 1 : " << x1 << " " << x2 << " " << x3 * k2 << " " << x4 << " " << (x1 * x1 + x2 * x2) / x3  - k2 * x4 << std::endl;

			if (gamma >= 0) {
				// the solution is good !
				typeSol = 1;
				goodSol = true;
			}
			else {
				if (c1122 == 0) { // cas d�g�n�r�
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
		
	}
	*/

}

void MarketEndoDirectGPU::updateXWOCurrent()
{
	updateXEndoMarket <<<_numBlocksB, _blockSize >> > (X._matrixGPU, Chat._matrixGPU, VoltageLimit._matrixGPU, nChild._matrixGPU, _CoresChatBegin._matrixGPU, _indiceBusBegin._matrixGPU, _nBus);
}

void MarketEndoDirectGPU::updateXWOCurrentCPU()
{

	X.transferCPU();
	Chat.transferCPU();
	double x1, x2, x3, x4, c1, c2, c3, c4, lambdaLo, lambdaUp, x3min, x3max, gamma, k2;
	double c1122;
	int nSol = 0;
	
	double p = 0;
	int nRoot = 0;

	for (int i = 1; i < _nBus; i++) {
		int typeSol = 0;
		int BestRoot = 0;
		double bestGamma = -1;
		bool goodSol = false;

		int begining = indiceBusBeginCPU.get(i, 0);
		int nC = nChildCPU.get(i, 0);
		int beginChat = CoresChatBeginCPU.get(i, 0);

		k2 = sqrt(2.0 / (nC + 1));
		
		c1 = -2 * Chat.get(beginChat, 0);
		c2 = -2 * Chat.get(beginChat + 1, 0);
		c4 = -2 * Chat.get(beginChat + 2, 0);
		c3 = -2 * Chat.get(beginChat + 3, 0) / k2;
		c1122 = c1 * c1 + c2 * c2;


		x3min = VoltageLimitCPU.get(0, i);
		x3max = VoltageLimitCPU.get(1, i);

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
		
		if (!goodSol) { // cas d�g�n�r�
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
					break;
				}
				if (gamma > bestGamma && lambdaUp > bestGamma) {
					typeSol = 2;
					bestGamma = min(gamma, lambdaUp);
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
		X.set(begining, 0, x1);
		X.set(begining + 1, 0, x2);
		X.set(begining + 2, 0, x4);
		X.set(begining + 3, 0, x3* k2);

		//std::cout << "x F : " << x1 << " " << x2 << " " << x3*k2 << " " << x4 << " " << gamma << std::endl;

	}
	X.transferGPU();
	Chat.transferGPU();

}

void MarketEndoDirectGPU::updatePMarket()
{
	/*std::cout << "objective loss P " << Y.get(17, 0, false) << " "<< Y.get(18, 0, false) << std::endl;
	std::cout << Ap1.get(1, 0) << " " << Ap2.get(1, 0) << " " << Ap3.get(1, 0) << " "  << Bp2.get(1, 0) << " " << Cp.get(1, 0) << std::endl;
	std::cout << _at1 << " " << _at2 << std::endl;
	Ct.display();
	
	std::cout << "Bt1" << std::endl;
	Bt1.display(true);*/
	
	if (getQLoss() > 0) {
		matUb.setBloc(_nTradeP, _nTradeP + _nAgentTrue - 1, 0, 1, _Qloss);
		matLb.setBloc(_nTradeP, _nTradeP + _nAgentTrue - 1, 0, 1, 0.0);
	}
	else {
		matLb.setBloc(_nTradeP, _nTradeP + _nAgentTrue - 1, 0, 1, _Qloss);
		matUb.setBloc(_nTradeP, _nTradeP + _nAgentTrue - 1, 0, 1, 0.0);
	}
	
	/*std::cout << "limit" << std::endl;
	matLb.display(true);
	matUb.display(true);*/

	
	updateLocalProb();
	/*std::cout << "result" << std::endl;
	Tlocal.display(true);
	P.display(true);*/
	
	
	TradeLin.swap(&Tlocal);
	/*std::cout << Ap1.get(0, 0) << " " << Bp1.get(0, 0) << std::endl;
	std::cout << " P" << std::endl;
	P.display();
	Tmoy.display();
	TradeLin.display();*/
	
	updateXPn<<<_nBusWLoss, _blockSizeSmall>>>(X._matrixGPU, Pn._matrixGPU, P._matrixGPU, nVoisin._matrixGPU, _indiceBusBegin._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU,  _CoresAgentBus._matrixGPU, losstype, _nAgentTrue, _nBus);
	
	
	
}

void MarketEndoDirectGPU::updateMu()
{
	updateMUGPU << <_numBlocksH, _blockSize >> > (Mu._matrixGPU, Y._matrixGPU, X._matrixGPU, _rho, _sizeEndoMarketTotal);

}



float MarketEndoDirectGPU::getPLoss()
{
	_Ploss =  Y.get(_debutloss, 0, false);
	
	return _Ploss;
}

float MarketEndoDirectGPU::getQLoss()
{
	int indice = 1;
	if (losstype == LossType::POWER) {
		indice = _nAgentTrue;
	}

	_Qloss = Y.get(_debutloss + indice, 0, false);
	
	return _Qloss;
}

void MarketEndoDirectGPU::computeLoss()
{

	int numBlock = 1;
	switch (_blockSize) {
	case 512:
		ComputeLoss<512> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgentTrue, _nBus);
		break;
	case 256:
		ComputeLoss<256> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgentTrue, _nBus);
		break;
	case 128:
		ComputeLoss<128> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgentTrue, _nBus);
		break;
	case 64:
		ComputeLoss< 64> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgentTrue, _nBus);
		break;
	case 32:
		ComputeLoss< 32> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgentTrue, _nBus);
		break;
	case 16:
		ComputeLoss< 16> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgentTrue, _nBus); 
		break;
	case  8:
		ComputeLoss<  8> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgentTrue, _nBus);
		break;
	case  4:
		ComputeLoss<  4> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgentTrue, _nBus);
		break;
	case  2:
		ComputeLoss<  2> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgentTrue, _nBus);
		break;
	case  1:
		ComputeLoss<  1> << <numBlock, _blockSize >> > (X._matrixGPU, Pn._matrixGPU, _indiceBusBegin._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, losstype, _nAgentTrue, _nBus);
		break;
	}
}



void MarketEndoDirectGPU::updateChat()
{
	int numBlock = _nBus;
	switch (_blockSizeSmall) {
	case 512:
		updateChatGPU3<512> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case 256:
		updateChatGPU3<256> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case 128:
		updateChatGPU3<128> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case 64:
		updateChatGPU3< 64> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case 32:
		updateChatGPU3< 32> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case 16:
		updateChatGPU3< 16> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case  8:
		updateChatGPU3<  8> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case  4:
		updateChatGPU3<  4> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case  2:
		updateChatGPU3<  2> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	case  1:
		updateChatGPU3<  1> << <numBlock, _blockSizeSmall >> > (Chat._matrixGPU, Y._matrixGPU, Mu._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, PosChild._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _CoresChatBegin._matrixGPU, _nAgentByBus._matrixGPU, _rho, losstype, _nBus);
		break;
	}


	

	
	// pour puissance
	updateBp2();
	
	// pour �changes
	updateLAMBDABt1GPU << <_numBlocksM, _blockSize >> > (Bt1._matrixGPU, LAMBDALin._matrixGPU, TradeLin._matrixGPU, _rho, CoresLinTrans._matrixGPU, _nTrade);

	//Bt1.display();
}

void MarketEndoDirectGPU::CommunicationX()
{
/**/ // X = { Pi, Qi, li, vi,, (pn, qn), vAi (Pci, Qci, lci) for all child Ci }
	
	communicateX << <_nBusWLoss, _blockSize >> > (X._matrixGPU, nChild._matrixGPU, Ancestor._matrixGPU, Childs._matrixGPU, _indiceBusBegin._matrixGPU, _indiceChildBegin._matrixGPU, _nAgentByBus._matrixGPU, _CoresBusAgent._matrixGPU, PosAgent._matrixGPU, losstype, _nBus, _nAgentTrue);

	
	// Q udate in argmin 0.5yHy + Qy

	// Y = { Pi, Qi, vi, li, (pn ...), qn..., vAi,  Pci ... , Qci ... , lci ... for all child Ci }

	updateQ << <_numBlocksH, _blockSize >> > (Q._matrixGPU, X._matrixGPU, Mu._matrixGPU, _rho, _sizeEndoMarketTotal);


}


float MarketEndoDirectGPU::updateRes(int indice) 
{
	
	updateDiffGPU << <_numBlocksM, _blockSize >> > (tempNN._matrixGPU, Tlocal._matrixGPU, CoresLinTrans._matrixGPU, _nTrade);
	float resR = tempNN.max2();

	float resS = Tlocal.max2(&TradeLin); // nomalement * _rhog mais si _rhog est tres grand impossible que cela converge !!!




	float resSTemp = _rho *Y.max2(&Ypre); // 

	if (resSTemp > resS) {
		resS = resSTemp;
	}

	float resV = Y.max2(&X) * _ratioEps;
	
		

	resF.set(0, indice, resR);
	resF.set(1, indice, resS);
	resF.set(2, indice, resV);



	return MAX(MAX(resV, resS), resR);
}

float MarketEndoDirectGPU::updateResRhoFixe(int indice)
{
	updateDiffGPU << <_numBlocksM, _blockSize >> > (tempNN._matrixGPU, Tlocal._matrixGPU, CoresLinTrans._matrixGPU, _nTrade);
	float resR = tempNN.max2();

	float resS = Tlocal.max2(&TradeLin); // nomalement * _rhog mais si _rhog est tres grand impossible que cela converge !!!




	float resSTemp = _rho * Y.max2(&Ypre);

	if (resSTemp > resS) {
		resS = resSTemp;
	}

	float resV = Y.max2(&X);



	resF.set(0, indice, resR);
	resF.set(1, indice, resS);
	resF.set(2, indice, resV);


	return MAX(MAX(resV, resS), resR);
}

int MarketEndoDirectGPU::feasiblePoint()
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

void MarketEndoDirectGPU::ComputePFromAgentToBus()
{
	int numBlock = _nBus;
	switch (_blockSizeSmall) {
	case 512:
		ComputePFromAgentToBusGPU<512> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgentTrue, _nBus);
		break;
	case 256:
		ComputePFromAgentToBusGPU<256> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgentTrue, _nBus);
		break;
	case 128:
		ComputePFromAgentToBusGPU<128> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgentTrue, _nBus);
		break;
	case 64:
		ComputePFromAgentToBusGPU< 64> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgentTrue, _nBus);
		break;
	case 32:
		ComputePFromAgentToBusGPU< 32> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgentTrue, _nBus);
		break;
	case 16:
		ComputePFromAgentToBusGPU< 16> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgentTrue, _nBus);
		break;
	case  8:
		ComputePFromAgentToBusGPU<  8> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgentTrue, _nBus);
		break;
	case  4:
		ComputePFromAgentToBusGPU<  4> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgentTrue, _nBus);
		break;
	case  2:
		ComputePFromAgentToBusGPU<  2> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgentTrue, _nBus);
		break;
	case  1:
		ComputePFromAgentToBusGPU<  1> << <numBlock, _blockSizeSmall >> > (Pb._matrixGPU, Pbmin._matrixGPU, Pbmax._matrixGPU, CoresSoloBusAgent._matrixGPU, Pn._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, _CoresAgentBusBegin._matrixGPU, _nAgentTrue, _nBus);
		break;
	}
}


// Market !!!!

void MarketEndoDirectGPU::updateLocalProb() {
	// FB 1a
	int numBlocks = _nAgent;
	switch (_blockSize) {
	case 512:
		updateTradePGPUSharedResidual<512> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap2._matrixGPU, Ap123._matrixGPU, Bp2._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, _epsL, _nIterL);
		break;
	case 256:
		updateTradePGPUSharedResidual<256> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap2._matrixGPU, Ap123._matrixGPU, Bp2._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, _epsL, _nIterL);
		break;
	case 128:
		updateTradePGPUSharedResidual<128> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap2._matrixGPU, Ap123._matrixGPU, Bp2._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, _epsL, _nIterL);
		break;
	case 64:
		updateTradePGPUSharedResidual< 64> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap2._matrixGPU, Ap123._matrixGPU, Bp2._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, _epsL, _nIterL);
		break;
	case 32:
		updateTradePGPUSharedResidual< 32> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap2._matrixGPU, Ap123._matrixGPU, Bp2._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, _epsL, _nIterL);
		break;
	case 16:
		updateTradePGPUSharedResidual< 16> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap2._matrixGPU, Ap123._matrixGPU, Bp2._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, _epsL, _nIterL);
		break;
	case  8:
		updateTradePGPUSharedResidual<  8> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap2._matrixGPU, Ap123._matrixGPU, Bp2._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, _epsL, _nIterL);
		break;
	case  4:
		updateTradePGPUSharedResidual<  4> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap2._matrixGPU, Ap123._matrixGPU, Bp2._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, _epsL, _nIterL);
		break;
	case  2:
		updateTradePGPUSharedResidual<  2> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap2._matrixGPU, Ap123._matrixGPU, Bp2._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, _epsL, _nIterL);
		break;
	case  1:
		updateTradePGPUSharedResidual<  1> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap2._matrixGPU, Ap123._matrixGPU, Bp2._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, _epsL, _nIterL);
		break;
	}
}


void MarketEndoDirectGPU::updateBp2()
{
	
	updateBp2GPU << <_nBusWLoss, _blockSizeSmall >> > (Bp2._matrixGPU, Y._matrixGPU, Mu._matrixGPU, _indiceBusBegin._matrixGPU, _CoresAgentBusBegin._matrixGPU, _CoresAgentBus._matrixGPU, _nAgentByBus._matrixGPU, losstype, _rho, _nBus, _nAgentTrue);

	Bp2.divideT(&nVoisin);

}

// autre

void MarketEndoDirectGPU::display() {

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
	Pb.set(0.0);
	Pb.transferCPU();

	Cost1.transferCPU();
	Cost2.transferCPU();
	for (int i = 0; i < _nBus; i++) {
		int Nb = _nAgentByBus.get(i, 0);
		int begin = _CoresAgentBusBegin.get(i, 0);
		for (int In = 0; In < Nb; In++) {
			int n = _CoresAgentBus.get(In + begin, 0);
			Pb.set(i, 0, Pb.get(i, 0) + Pn.get(n, 0));
			Pb.set(i + _nBus, 0, Pb.get(i + _nBus, 0) + Pn.get(n + _nAgentTrue, 0));
		}
	}



	if (_iterGlobal == 0) {
		std::cout << "algorithm not launch" << std::endl;
	}
	else if (_iterGlobal < _iterG) {
		std::cout << "method " << _name << " converged in " << _iterGlobal << " iterations." << std::endl;
		std::cout << "Converged in " << (float) timeMarketEndo / CLOCKS_PER_SEC << " seconds" << std::endl;

	}
	else {
		std::cout << "method " << _name << " not converged in " << _iterGlobal << " iterations." << std::endl;
		std::cout << "time taken " << (float) timeMarketEndo / CLOCKS_PER_SEC << " seconds" << std::endl;
	}
	std::cout << "The power error of this state is (constraint) " << resF.get(0, _iterGlobal / _stepG) << " and convergence " << resF.get(1, _iterGlobal / _stepG) << std::endl;
	std::cout << "===============================================================|" << std::endl;
	std::cout << "      System Summary                                           |" << std::endl;
	std::cout << "===============================================================|" << std::endl;
	std::cout << "Buses            " << _nBus << std::endl;
	std::cout << "Branches         " << _nLine << std::endl;
	std::cout << "Agent            " << _nAgentTrue << std::endl;
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
	std::cout << " Line |    From     |    To      |                Upstream flow                      |    Constraint    |" << std::endl;
	std::cout << "  #   |    Bus      |    Bus     |   P (pu)   |   Q (pu)   |   l (pu)   |  Loss (pu) |      lmax        |" << std::endl;
	std::cout << "------|-------------|------------|------------|------------|------------|------------|------------------|" << std::endl;

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
__global__ void updateChatGPU3(float* Chat, float* Y, float* MU, float* nChild, float* Ancestor, float* posChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, float* CoresChatBegin, float* nAgentByBus, float _rho, int losstype, int nBus) {

	int bus = blockIdx.x;
	int index = threadIdx.x;
	int step = blockDim.x;
	
	int beginChat = CoresChatBegin[bus];

	__shared__ float shArr[_blockSizeSmall]; // c'est grand pour pas grand chose...


	int beginBus = indiceBusBegin[bus];
	int beginChild = (bus < (nBus - 1)) ? indiceChildBegin[bus] : 0;
	int childCount = nChild[bus];
	int Ai = Ancestor[bus];
	int nAgent = nAgentByBus[bus];
	int c = posChild[bus];
	float var = 0;
	int borne = 4;
	int divideVar = 2 + ((losstype == 1) && (index == 2)) * 1;
	int indiceLoss = indiceBusBegin[nBus];

	if (index < borne) {
		//float Phat, Qhat, lhat, vihat, pnhat..., qnhat...;
		var = Y[beginBus + index] / divideVar - MU[beginBus + index] / (divideVar * _rho);
		if (bus > 0) {
			if (index < 3) {
				int nAi = nChild[Ai];
				int nAgentAi = nAgentByBus[Ai];
				int indiceAncBus = indiceBusBegin[Ai] + 5 + 2 * nAgentAi + nAi * index + c;
				//var = indiceAncBus;
				var += Y[indiceAncBus] / divideVar - MU[indiceAncBus] / (divideVar * _rho);
				var += ((index == 2) && (losstype == 1)) ? (Y[indiceLoss + 2 + bus] / divideVar - MU[indiceLoss + 2 + bus] / (divideVar * _rho)) : 0;
			}
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


/*
float phat, qhat;
	int divideP = 1;
	if (losstype == LossType::POWER) { // POWER
		divideP += 1;
	}

	for (int i = 0; i < _nBus; i++) {
		int Nb = _nAgentByBus.get(i, 0);
		int begin = _CoresAgentBusBegin.get(i, 0);
		for (int In = 0; In < Nb; In++) {
			int n = _CoresAgentBus.get(In + begin, 0);

			//std::cout << "bus " << i << " agent " << n << " en pos " << In << " Y " << Y[i].get(5 + 2 * In, 0) << " " <<
			//	Y[_nBus].get(n, 0) << " " << Mu[i].get(5 + 2 * In, 0) << " " << Mu[_nBus].get(n, 0);

			phat = Y[i].get(5 + 2 * In, 0) - Mu[i].get(5 + 2 * In, 0) / _rho;
			qhat = Y[i].get(6 + 2 * In, 0) - Mu[i].get(6 + 2 * In, 0) / _rho;
			if (losstype == LossType::POWER) {
				phat += Y[_nBus].get(n, 0) - Mu[_nBus].get(n, 0) / (_rho);
				qhat += Y[_nBus].get(n + _nAgentTrue, 0) - Mu[_nBus].get(n + _nAgentTrue, 0) / (_rho);
			}
			//std::cout <<  " phat " << phat/divideP << " Bp2 " << phat / (divideP * nVoisin.get(n,0)) << std::endl;
			//phat = (Y[i].get(5 + 2 * In, 0) + Y[_nBus].get(n, 0)) / 2			    - (Mu[i].get(5 + 2 * In, 0) + Mu[_nBus].get(n, 0)) / (2 * _rho);
			//qhat = (Y[i].get(6 + 2 * In, 0) + Y[_nBus].get(n + _nAgentTrue, 0)) / 2 - (Mu[i].get(6 + 2 * In, 0) + Mu[_nBus].get(n + _nAgentTrue, 0)) / (2 * _rho);


			Bp2.set(n, 0, phat / divideP);
			Bp2.set(n + _nAgentTrue, 0, qhat / divideP);
		}
	}
	// Y     (Ploss, Pn, Qloss Qn)

	float phatLoss = 0; // Y[_nBus].get(0, 0) - Mu[_nBus].get(0, 0) / _rho;
	float qhatLoss = 0; // Y[_nBus].get(1, 0) - Mu[_nBus].get(1, 0) / _rho;
	switch (losstype)
	{
	case LossType::POWER:
		phatLoss = Y[_nBus].get(0, 0) - Mu[_nBus].get(0, 0) / _rho;
		qhatLoss = Y[_nBus].get(_nAgentTrue, 0) - Mu[_nBus].get(_nAgentTrue, 0) / _rho;


		break;
	case LossType::CURRENT:
		phatLoss = Y[_nBus].get(0, 0) - Mu[_nBus].get(0, 0) / _rho;
		qhatLoss = Y[_nBus].get(1, 0) - Mu[_nBus].get(1, 0) / _rho;
		break;
	}


Bp2.set(0, 0, phatLoss);
Bp2.set(_nAgentTrue, 0, qhatLoss);

Bp2.divideT(&nVoisin);

*/

__global__ void updateBp2GPU(float* Bp2, float* Y, float* MU, float* indiceBusBegin, float* indiceAgentBegin, float* CoresAgentBus, float* nAgentByBus, int losstype, float rho, int nBus, int nAgent) {
	
	int bus	    = blockIdx.x;
	int thIdx   = threadIdx.x;
	int step	= blockDim.x;

	int begin   = indiceBusBegin[bus];
	int divideP = 1 + 1 * (losstype == 0);
	
	float phat = 0;
	float qhat = 0;

	if(bus < nBus) // bus normaux 
	{
		int Nb      = nAgentByBus[bus];
		int beginAgent = indiceAgentBegin[bus];
		int beginLoss = indiceBusBegin[nBus];
		
		for (int In = thIdx; In < Nb; In += step) {
			int n = CoresAgentBus[In + beginAgent];

			phat = Y[begin + 4 + In] - MU[begin + 4 + In] / rho;
			qhat = Y[begin + 4 + Nb + In] - MU[begin + 4 + Nb + In] / rho;
			
			if (losstype == 0) {
				phat += Y[beginLoss + n] - MU[beginLoss + n] / rho;
				qhat += Y[beginLoss + n + Nb] - MU[beginLoss + n + Nb] / rho;
			}
			

			Bp2[n] = phat / divideP;
			Bp2[n + nAgent] = qhat / divideP;
		}
	}
	else { // bus des pertes
		if (thIdx == 0) {
			if (losstype == 0) {
				phat = Y[begin] - MU[begin] / rho;
				qhat = Y[begin + nAgent] - MU[begin + nAgent] / rho;
			}
			else {
				phat = Y[begin] - MU[begin] / rho;
				qhat = Y[begin + 1] - MU[begin + 1] / rho;
			}
			Bp2[0] = phat;
			Bp2[nAgent] = qhat;
		}
	}

}

