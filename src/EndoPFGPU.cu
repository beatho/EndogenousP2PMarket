#include "../head/EndoPFGPU.cuh"
#define MAX(X, Y) X * (X >= Y) + Y * (Y > X)

// On prend la transpos�e de G !!! (ie G(n,i) = G[n*Nvar + i] )


EndoPFGPU::EndoPFGPU() : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " EndoPFGPU Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu
}


EndoPFGPU::EndoPFGPU(float rho) : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default EndoPFGPU Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu
}

EndoPFGPU::~EndoPFGPU()
{
	DELETEB(pf);
}
void EndoPFGPU::setParam(float rho)
{
	_rho = rho;
}

void EndoPFGPU::setTau(float tau)
{
	if (tau < 1) {
		throw std::invalid_argument("tau must be greater than 1");
	}
	_tau = tau;
}



void EndoPFGPU::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
{
#ifdef DEBUG_SOLVE
	cas.display();
	sim.display(1);
#endif // DEBUG_SOLVE
	
	timeEndoPF = clock();
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
		timePerBlock.set(0, 0, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		occurencePerBlock.set(0, 0, 1);
#endif // INSTRUMENTATION
	}
	_rhog = sim.getRho();
	_at1 = _rhog;
	
	int iterL = sim.getIterL();
	int stepL = sim.getStepL()/20;

	float epsL = sim.getEpsL();
	float epsG = sim.getEpsG();
	
	float resG = 2 * epsG;
	float epsL2 = epsL * epsL;
	_iterGlobal = 0;
	//CHECK_LAST_CUDA_ERROR();
	//Pn.display(true);
	//std::cout << "*******" << std::endl;
	while ((_iterGlobal < _iterG) && (resG>epsG) || (_iterGlobal <= _stepG)) {
		
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		updateLocalProbGPU(epsL2, iterL);
		//CHECK_LAST_CUDA_ERROR();
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 1, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		
		tradeLin.swap(&Tlocal); // echange juste les pointeurs	


		updateGlobalProbGPU();
		//CHECK_LAST_CUDA_ERROR();
		if (!(_iterGlobal % _stepG)) {
#ifdef INSTRUMENTATION
			cudaDeviceSynchronize();
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateResBis(&resF, &Tlocal, _iterGlobal / _stepG, &tempNN);
			//CHECK_LAST_CUDA_ERROR();
#ifdef INSTRUMENTATION
			cudaDeviceSynchronize();
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		}
		_iterGlobal++;
	}
	//std::cout << _iterGlobal << " " << resG << " " << resF.get(0, (_iterGlobal - 1) / _stepG) << " " << resF.get(1, (_iterGlobal - 1) / _stepG) << " " << resF.get(2, (_iterGlobal - 1) / _stepG) << std::endl;
#ifdef INSTRUMENTATION
	occurencePerBlock.increment(0, 1, _iterGlobal);
	occurencePerBlock.increment(0, 3, _iterGlobal);
	occurencePerBlock.increment(0, 4, _iterGlobal);
	occurencePerBlock.increment(0, 5, _iterGlobal);
	occurencePerBlock.increment(0, 6, _iterGlobal / _stepG);

	cudaDeviceSynchronize();
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	float fc = calcFc(&a, &b, &tradeLin, &Pn, &Ct, &tempN1, &tempNN);
	//Method::calcFc(MatrixGPU* cost1, MatrixGPU* cost2, MatrixGPU* trade, MatrixGPU* Pn, MatrixGPU* BETA, MatrixGPU* tempN1, MatrixGPU* tempNN)
	MatrixCPU tradeLinCPU;
	tradeLin.toMatCPU(tradeLinCPU);
	MatrixCPU LAMBDALinCPU;
	LAMBDALin.toMatCPU(LAMBDALinCPU);
	MatrixCPU PnCPU;
	Pn.toMatCPU(PnCPU);
	MatrixCPU MUCPU;
	MU.toMatCPU(MUCPU);
	//CHECK_LAST_CUDA_ERROR();
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
	


	// FB 5
	result->setResF(&resF);
	result->setLAMBDA(&LAMBDA);
	result->setTrade(&trade);
	//result->setDelta1(&delta1);
	//result->setDelta2(&delta2);
	result->setIter(_iterGlobal);
	
	result->setPn(&PnCPU);
	result->setFc(fc);
	result->setMU(&MUCPU);

#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 7, 1);
	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
	timeEndoPF = clock() - timeEndoPF;
	result->setTime((float)timeEndoPF / CLOCKS_PER_SEC);
	
}

void EndoPFGPU::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	matLb.transferCPU();
	Pmin = MatrixGPU(cas.getPmin());
	Pmax = MatrixGPU(cas.getPmax());


	MatrixCPU Lb(cas.getLb());

	b = cas.getb();
	Cp1 = b;
	int indice = 0;

	for (int idAgent = 0; idAgent < _nAgentTrue; idAgent++) {
		int Nvoisinmax = nVoisinCPU.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			matLb.set(indice, 0, Lb.get(idAgent, 0));
			indice = indice + 1;
		}
	}
	for (int idAgent = _nAgentTrue; idAgent < _nAgent; idAgent++) {
		for (int voisin = 0; voisin < (_nAgentTrue - 1); voisin++) {
			matLb.set(indice, 0, Lb.get(idAgent, 0));
			indice = indice + 1;
		}
	}

	matLb.transferGPU();

	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);
	Cp1.multiplyT(&nVoisin);
	

#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION


}

void EndoPFGPU::init(const Simparam& sim, const StudyCase& cas)
{
	DELETEB(pf);
	if (CoresMatLin.getPos()) { 
		CoresMatLin.transferCPU();
		CoresLinAgent.transferCPU();
		CoresAgentLin.transferCPU();
		CoresLinVoisin.transferCPU();
		CoresLinTrans.transferCPU();

		Tlocal_pre.transferCPU();
		tradeLin.transferCPU();
		LAMBDALin.transferCPU();

		matLb.transferCPU();
		matUb.transferCPU();
		Ct.transferCPU();
	}

	// intitilisation des matrixs et variables 
	
	//std::cout << "init " << std::endl;
	if (!cas.isAC()) {
		throw std::invalid_argument("Wrong studyCase must be AC");
	}

	_rhog = sim.getRho();
	_rho1 = sim.getRho1();
	_iterG = sim.getIterG();
	_stepG = sim.getStepG();
	float epsG = sim.getEpsG();
	float epsGC = sim.getEpsGC();
	_ratioEps = epsG / epsGC;
	isRadial = cas.isRadial();
	_nAgentTrue = sim.getNAgent();
	_nAgent = _nAgentTrue + _nAgentTrue;

	_rhol = _rho; //*nAgent
	//std::cout << "rho " << _rho << std::endl;
	if (_rho == 0) {
		_rhol = _rhog;
	}
	
	nVoisinCPU = cas.getNvoi();
	nVoisin = MatrixGPU(nVoisinCPU, 1);
	nVoisin.preallocateReduction();
	int nVoisinMax = nVoisin.max2();
	if (_blockSize * NMAXPEERPERTRHREAD < nVoisinMax) {
		std::cout << _blockSize << " " << NMAXPEERPERTRHREAD << " " << nVoisinMax << std::endl;
		throw std::invalid_argument("For this Method, an agent must not have more than 5120 peers");
	}

	//CHECK_LAST_CUDA_ERROR();
	_nLine = cas.getNLine(true);
	//std::cout << "_nLine " << _nLine << std::endl;
	_nBus = cas.getNBus();
	_nVarPF = _nLine + 2 * _nBus;

	//std::cout << _nVarPF << std::endl;
	
	_nTrade = nVoisin.sum();
	_nTradeP = nVoisin.sum(0,_nAgentTrue);
	//std::cout << "nTrade " << _nTrade << " " << _nTradeP << std::endl;
	//CHECK_LAST_CUDA_ERROR();
	_nTradeQ = _nTrade - _nTradeP;
	//std::cout << "nTrade " << _nTradeQ << " " << nVoisin.sum(_nAgentTrue, _nAgent) << std::endl;
	if (_nTradeQ != (_nAgentTrue * (_nAgentTrue - 1))) {
		std::cout << "err EndoPFGPU : " << _nAgent << " " << _nAgentTrue << " " << _nTrade << " " << _nTradeP << " " << _nTradeQ << std::endl;
		throw std::invalid_argument("Agent must be fully conected for the Q echanges, WIP");
	}
	_numBlocksN = ceil((_nAgent + _blockSize - 1) / _blockSize);
	_numBlocksM = ceil((_nTrade + _blockSize - 1) / _blockSize);
	_numBlocksL = ceil((_nLine + _blockSize - 1) / _blockSize);
	_numBlocksBL = ceil((_nVarPF + _blockSize - 1) / _blockSize);
	
	
	//std::cout << _numBlocksN << " " << _numBlocksM << " " << _numBlocksL << " " << _numBlocksBL << std::endl;
	
	//std::cout <<  _blockSize << std::endl;
	if (initWithMarketClear) {
		ADMMMarketGPU market;
		Simparam res(sim);
		market.solve(&res, sim, cas);
		//res.display();
		LAMBDA = res.getLambda();
		trade = res.getTrade();
		Pnpre = MatrixGPU(res.getPn(), 1);

	}
	else {
		LAMBDA = sim.getLambda();
		trade = sim.getTrade();
		Pnpre = MatrixGPU(sim.getPn(), 1);
	}
	//Pnpre.display(true);
	Tmoy = Pnpre;
	Tmoy.divideT(&nVoisin);
	//std::cout << "*******" << std::endl;
	
	
	_at1 = _rhog; // represente en fait 2*a
	_at2 = _rhol;

	resF = MatrixCPU(3, (_iterG / _stepG) + 1);
	resX = MatrixCPU(4, (_iterG / _stepG) + 1);

	MatrixGPU BETA(cas.getBeta());
	
	MatrixGPU Ub(cas.getUb());
	MatrixGPU Lb(cas.getLb());
	
	 
	//std::cout << "mise sous forme lin�aire" << std::endl;
	


	CoresMatLin = MatrixGPU(_nAgent, _nAgentTrue, -1);
	CoresLinAgent = MatrixGPU(_nTrade, 1);
	CoresAgentLin = MatrixGPU(_nAgent + 1, 1);
	CoresLinVoisin = MatrixGPU(_nTrade, 1);
	CoresLinTrans = MatrixGPU(_nTrade, 1);

	Tlocal_pre = MatrixGPU(_nTrade, 1);
	tradeLin = MatrixGPU(_nTrade, 1);
	LAMBDALin = MatrixGPU(_nTrade, 1);

	matLb = MatrixGPU(_nTrade, 1);
	matUb = MatrixGPU(_nTrade, 1);
	Ct = MatrixGPU(_nTrade, 1);
	

	int indice = 0;
	//std::cout << " P " << std::endl;
	for (int idAgent = 0; idAgent < _nAgentTrue; idAgent++) { // P
		MatrixGPU omega(cas.getVoisin(idAgent));
		int Nvoisinmax = nVoisinCPU.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			int idVoisin = omega.get(voisin, 0);
			if(Lb.getNCol()== 1){
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
	//std::cout << " Q " << std::endl;
	for (int idAgent = _nAgentTrue; idAgent < _nAgent; idAgent++) { // Q
		for (int idVoisin = 0; idVoisin < _nAgentTrue; idVoisin++) {
			if (idVoisin != (idAgent - _nAgentTrue)) {
				matLb.set(indice, 0, Lb.get(idAgent, 0));
				matUb.set(indice, 0, Ub.get(idAgent, 0));
				//Ct.set(indice, 0, BETA.get(idAgent, idVoisin));
				tradeLin.set(indice, 0, trade.get(idAgent, idVoisin));
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
	tradeLin.transferGPU();
	LAMBDALin.transferGPU();

	CoresAgentLin.transferGPU();
	CoresLinAgent.transferGPU();
	CoresLinVoisin.transferGPU();
	CoresMatLin.transferGPU();
	CoresLinTrans.transferGPU();

	
	//CHECK_LAST_CUDA_ERROR();
	//std::cout << "donnees sur CPU pour le grid" << std::endl;
	delta1 = MatrixGPU(_nVarPF, 1, 0, 1);
	delta2 = MatrixGPU(_nVarPF, 1, 0, 1);
	Z1 = MatrixGPU(_nVarPF, 1, 0, 1);
	Z2 = MatrixGPU(_nVarPF, 1, 0, 1);
	Y = MatrixGPU(_nVarPF, 1, 0, 1);
	Ypre = MatrixGPU(_nVarPF, 1, 0, 1);
	dY = MatrixGPU(_nVarPF, 1, 0, 1);

	Ylimit = MatrixGPU(_nVarPF, 1, 0, 1); // angle, amplitude, flux
	YOffset = MatrixGPU(_nVarPF, 1, 0, 1); // angle, amplitude, flux
	G = MatrixGPU(_nAgent, _nVarPF, 0, 1);
	SensiBis = MatrixGPU(_nVarPF, 1, 0, 1);
	
	
	MatrixGPU LimitsUb(cas.getUpperBound(), 1); // angle, amplitude, flux
	MatrixGPU LimitsLb(cas.getLowerBound(), 1);
	
	initLimits << <_numBlocksBL, _blockSize >> > (Ylimit._matrixGPU, YOffset._matrixGPU, LimitsLb._matrixGPU, LimitsUb._matrixGPU, _nVarPF);
	//Ylimit.display(true);
	//YOffset.display(true);
	
	CHECK_LAST_CUDA_ERROR();

	//std::cout << " PF " << std::endl;
	if (isRadial) {
		pf = new GPUPFdistPQ;
	}
	else {
		pf = new GPUPF;
	}
	//Ylimit.display();
	//YOffset.display();

	pf->init(cas, &Pnpre);

	
	//CHECK_LAST_CUDA_ERROR();
	//std::cout << "autres donn�e sur CPU" << std::endl;
	tempNN = MatrixGPU(_nTrade, 1, 0, 1);
	tempNN.preallocateReduction();
	tempN1 = MatrixGPU(_nAgent, 1, 0, 1); // plut�t que de re-allouer de la m�moire � chaque utilisation
	tempL1 = MatrixGPU(_nVarPF, 1, 0, 1);
	tempL1.preallocateReduction();
	//MatrixGPU temp1N(1, _nAgent, 0, 1);

	Tlocal = MatrixGPU(_nTrade, 1, 0, 1);
	Tlocal.preallocateReduction();
	

	P = MatrixGPU(_nAgent, 1, 0, 1); // moyenne des trades
	Pn = MatrixGPU(sim.getPn(), 1);
	dP = MatrixGPU(_nAgent, 1, 0, 1);

	a = MatrixGPU(cas.geta(), 1);
	b = MatrixGPU(cas.getb(), 1);
	Ap2 = a;
	Ap1 = nVoisin;
	Ap12 = MatrixGPU(_nAgent, 1, 0, 1);

	Bt1 = MatrixGPU(_nTrade, 1, 0, 1);
	Bt2 = MatrixGPU(_nTrade, 1, 0, 1);
	Cp = MatrixGPU(_nAgent, 1, 0, 1);
	Cp2 = MatrixGPU(_nAgent, 1, 0, 1);
	Cp1 = b;
	Bp1 = MatrixGPU(_nAgent, 1, 0, 1);

	Pmin = MatrixGPU(cas.getPmin(), 1);
	Pmax = MatrixGPU(cas.getPmax(), 1);
	

	MU = MatrixGPU(sim.getMU(), 1); // facteur reduit i.e lambda_l/_rho
	

	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);
	
	


	Ap1.multiply(_rhol);
	Ap2.multiplyT(&nVoisin);
	Ap2.multiplyT(&nVoisin);
	Ap12.add(&Ap1, &Ap2);

	Cp1.multiplyT(&nVoisin);
	
	
	

	//std::cout << "update Global" << std::endl;
	updateGlobalProbGPU();
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	//std::cout << "fin init " << std::endl;
	CHECK_LAST_CUDA_ERROR();

}

void EndoPFGPU::updateGlobalProbGPU() {
	
	// FB 2a
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	//std::cout << "*************" << std::endl;
	
	Pn.swap(&Pnpre);
	updatePn(&Pn,&Tmoy,&nVoisin);
	//std::cout << " Tmoy " << std::endl;
	
	//std::cout << " Pn " << std::endl;
	//Pn.display(true);
//	tradeLin.display(true);
	
	//std::cout << "update PF" << std::endl;
	pf->updatePQ(&Pn);
	
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	timePerBlock.increment(0, 3, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 3b
	cudaDeviceSynchronize();
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION	
	 

	//std::cout << "solve PF" << std::endl;
	
	// FB 3
	pf->solve();
	 

	tempL1 = pf->getY();
	Ypre.swap(&Y);
	Y.subtract(&tempL1, &YOffset);
	
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 4, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 3c
	cudaDeviceSynchronize();
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	//std::cout << "update PF" << std::endl;
	//std::cout << "Y " << std::endl;
	//Y.display(true);
	
	float Ploss = - pf->getPloss() ;
	float Qloss = - pf->getQloss();	
	//std::cout << " Ploss " << Ploss << " Qloss " << Qloss << std::endl;
	
	Pmin.set(0, 0, Ploss / nVoisinCPU.get(0, 0), true);
	Pmax.set(0, 0, Ploss / nVoisinCPU.get(0, 0), true);
	
	Pmin.set(_nAgentTrue, 0, Qloss / nVoisinCPU.get(_nAgentTrue, 0), true);
	Pmax.set(_nAgentTrue, 0, Qloss / nVoisinCPU.get(_nAgentTrue, 0), true);
	
	if (Qloss > 0) {
		matUb.setBloc(_nTradeP, _nTradeP + _nAgentTrue - 1, 0, 1, Qloss);
	}
	else {
		matLb.setBloc(_nTradeP, _nTradeP + _nAgentTrue - 1, 0, 1, Qloss);
	}/**/
	
 
	


	updateZDeltaGPU << <_numBlocksBL, _blockSize >> > (Z1._matrixGPU, Z2._matrixGPU, Ylimit._matrixGPU, delta1._matrixGPU, delta2._matrixGPU, Y._matrixGPU, _nVarPF);

	updateSensi();
	
	updateCp2();
	Cp.add(&Cp1, &Cp2);
/*std::cout << " Z " << std::endl;
	Z1.display(true);
	Z2.display(true);

	std::cout << " Delta " << std::endl;
	delta1.display(true);
	delta2.display(true);
std::cout << " Sensi " << std::endl;
	SensiBis.display(true);
	G.display(true);	
	std::cout << " Cp " << std::endl;
	Cp2.display(true);
	//Cp.display(true);
*/
	////CHECK_LAST_CUDA_ERROR();

	

	// FB 3c
	updateLAMBDABt1GPU << <_numBlocksM, _blockSize >> > (Bt1._matrixGPU, LAMBDALin._matrixGPU, tradeLin._matrixGPU, _rhog, CoresLinTrans._matrixGPU, _nTrade);
	 
	
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 5, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
	
}

void EndoPFGPU::updateLocalProbGPU(float epsL, int nIterL) {
	// FB 1a
	int numBlocks = _nAgent;
	/*std::cout << "probl�me local" << std::endl;
	std::cout << _at1 << " " << _at2 << std::endl;
	Bt1.display(true);
	Ct.display(true);
	matLb.display(true);
	matUb.display(true);
	Ap1.display(true);
	Ap12.display(true);
	Bp1.display(true);
	Cp.display(true);
	Pmin.display(true);
	Pmax.display(true);*/


	switch (_blockSize) {
	case 512:
		updateTradePGPUSharedResidual<512> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case 256:
		updateTradePGPUSharedResidual<256> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case 128:
		updateTradePGPUSharedResidual<128> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case 64:
		updateTradePGPUSharedResidual< 64> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case 32:
		updateTradePGPUSharedResidual< 32> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case 16:
		updateTradePGPUSharedResidual< 16> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case  8:
		updateTradePGPUSharedResidual<  8> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case  4:
		updateTradePGPUSharedResidual<  4> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case  2:
		updateTradePGPUSharedResidual<  2> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case  1:
		updateTradePGPUSharedResidual<  1> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	}

	

}


void EndoPFGPU::updateCp2()
{
	int numBlocks = _nAgent;
	switch (_blockSize) {
	case 512:
		updateCp2GPU<512> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, SensiBis._matrixGPU, G._matrixGPU, nVoisin._matrixGPU, _rho1, _nVarPF);
		break;
	case 256:
		updateCp2GPU<256> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, SensiBis._matrixGPU, G._matrixGPU, nVoisin._matrixGPU, _rho1, _nVarPF);
		break;
	case 128:
		updateCp2GPU<128> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, SensiBis._matrixGPU, G._matrixGPU, nVoisin._matrixGPU, _rho1, _nVarPF);
		break;
	case 64:
		updateCp2GPU< 64> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, SensiBis._matrixGPU, G._matrixGPU, nVoisin._matrixGPU, _rho1, _nVarPF);
		break;
	case 32:
		updateCp2GPU< 32> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, SensiBis._matrixGPU, G._matrixGPU, nVoisin._matrixGPU, _rho1, _nVarPF);
		break;
	case 16:
		updateCp2GPU< 16> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, SensiBis._matrixGPU, G._matrixGPU, nVoisin._matrixGPU, _rho1, _nVarPF);
		break;
	case  8:
		updateCp2GPU<  8> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, SensiBis._matrixGPU, G._matrixGPU, nVoisin._matrixGPU, _rho1, _nVarPF);
		break;
	case  4:
		updateCp2GPU<  4> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, SensiBis._matrixGPU, G._matrixGPU, nVoisin._matrixGPU, _rho1, _nVarPF);
		break;
	case  2:
		updateCp2GPU<  2> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, SensiBis._matrixGPU, G._matrixGPU, nVoisin._matrixGPU, _rho1, _nVarPF);
		break;
	case  1:
		updateCp2GPU<  1> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, SensiBis._matrixGPU, G._matrixGPU, nVoisin._matrixGPU, _rho1, _nVarPF);
		break;
	}



}



void EndoPFGPU::updateSensi()
{

	updateSensiBis<<<_numBlocksBL,_blockSize>>> (SensiBis._matrixGPU, Y._matrixGPU, Z1._matrixGPU, Z2._matrixGPU, delta1._matrixGPU, delta2._matrixGPU, _nVarPF);
	dY.subtract(&Y, &Ypre);
	dP.subtract(&Pn, &Pnpre);
	//dY.display(true);
	//dP.display(true);
	updateSensiGPU << <_nAgent, _blockSize >> > (G._matrixGPU, dY._matrixGPU, dP._matrixGPU, _nVarPF);
	//G.display();
}




float EndoPFGPU::updateResBis(MatrixCPU* res, MatrixGPU* Tlocal, int iter, MatrixGPU* tempNN)
{
	
	float resS = Tlocal->max2(&tradeLin);

	updateDiffGPU << <_numBlocksM, _blockSize >> > (tempNN->_matrixGPU, tradeLin._matrixGPU, CoresLinTrans._matrixGPU, _nTrade);
	float resR = tempNN->max2();
	// Residus reseau
	
	tempL1.subtract(&Ylimit, &Y);
	tempL1.projectNeg();

	float resXf = _ratioEps * tempL1.max2();
	res->set(0, iter, resR);
	res->set(1, iter, resS);
	res->set(2, iter, resXf);
	return MAX(MAX(resXf, resS), resR);
}

void EndoPFGPU::display() {

	a.transferCPU();
	b.transferCPU();
	Pn.transferCPU();
	Pmin.transferCPU();
	Pmax.transferCPU();

	if (_iterGlobal == 0) {
		std::cout << "algorithm not launch" << std::endl;
	}
	else if (_iterGlobal < _iterG) {
		std::cout << "method " << _name << " converged in " << _iterGlobal << " iterations." << std::endl;
		std::cout << "Converged in " << (float) timeEndoPF / CLOCKS_PER_SEC << " seconds" << std::endl;

	}
	else {
		std::cout << "method " << _name << " not converged in " << _iterGlobal << " iterations." << std::endl;
		std::cout << "time taken " << (float) timeEndoPF / CLOCKS_PER_SEC << " seconds" << std::endl;
	}
	std::cout << "The power error of this state is (constraint) " << resF.get(0, _iterGlobal / _stepG) << " and convergence " << resF.get(1, _iterGlobal / _stepG) << std::endl;
	std::cout << "===============================================================|" << std::endl;
	std::cout << "      System Summary                                           |" << std::endl;
	std::cout << "===============================================================|" << std::endl;
	std::cout << "Agent            " << _nAgentTrue << std::endl;



	std::cout << std::endl << std::endl;
	std::cout << std::endl << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "      Agent Data                                                                                        |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Agent  |  Cost    |  Cost    |          Power Injection           |           Power Injection          |" << std::endl;
	std::cout << "  #     |   a (pu) |   b (pu) |  P (pu)  | Pmin (pu)  | Pmax (pu)  |  Q (pu)   | Qmin (pu)  | Qmax (pu) |" << std::endl;
	std::cout << "--------|----------|----------|----------|------------|------------|-----------|------------|-----------|" << std::endl;

	for (int n = 0; n < _nAgentTrue; n++) {

		std::cout << std::setw(8) << n << "|" << std::setw(9) << a.get(n, 0) << " |" << std::setw(10)
			<< b.get(n, 0) << "|" << std::setw(10) << Pn.get(n, 0) << "|" << std::setw(12)
			<< Pmin.get(n, 0) * nVoisinCPU.get(n, 0) << "|" << std::setw(12) << Pmax.get(n, 0) * nVoisinCPU.get(n, 0)
			<< "|" << std::setw(11) << Pn.get(n + _nAgentTrue, 0) << "|" << std::setw(12) << Pmin.get(n + _nAgentTrue, 0) * nVoisinCPU.get(n + _nAgentTrue, 0)
			<< "|" << std::setw(11) << Pmax.get(n + _nAgentTrue, 0) * nVoisinCPU.get(n + _nAgentTrue, 0) << "|" << std::endl;
	}


	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "                      END PRINT                                                                         |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;


	a.transferGPU();
	b.transferGPU();
	Pn.transferGPU();
	Pmin.transferGPU();
	Pmax.transferGPU();
}





__global__ void initLimits(float* Ylimit, float* Yoffset, float* limitsLb, float* limitsUb, int nVarPF) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;

	for (int i = index; i < nVarPF; i += step) {
		float ub = limitsUb[i];
		float lb = limitsLb[i];
		float mid = (ub + lb) / 2;
		float lim = ub - mid;
		Ylimit[i] = lim;
		Yoffset[i] = mid;
	}
}



/*
	Z1.add(&Ylimit, &delta1);
	Z1.subtract(&Y);
	Z1.projectPos();

	Z2.add(&Ylimit, &delta2);
	Z2.add(&Y);
	Z2.projectPos();
*/

__global__ void updateZGPU(float* Z1, float* Z2, float* Ylimit, float* delta1, float* delta2, float* Y, int nVarPF) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;

	for (int i = index; i < nVarPF; i += step) {

		float z1 = Ylimit[i] + delta1[i] - Y[i];
		float z2 = Ylimit[i] + delta2[i] + Y[i];
		Z1[i] = (z1 > 0) * z1;
		Z2[i] = (z2 > 0) * z2;

	}
}


/*
	delta1.add(&Ylimit);
	delta1.subtract(&Z1);
	delta1.subtract(&Y);

	delta2.add(&Ylimit);
	delta2.subtract(&Z2);
	delta2.add(&Y);
*/
__global__ void updateDeltaGPU(float* delta1, float* delta2, float* Z1, float* Z2, float* Ylimit, float* Y, int nVarPF) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;

	for (int i = index; i < nVarPF; i += step) {
		float d1 = Ylimit[i] + delta1[i] - Y[i];
		float d2 = Ylimit[i] + delta2[i] + Y[i];

		delta1[i] = (d1 < 0) * d1;//delta1[i] + Ylimit[i] - Z1[i] - Y[i];
		delta2[i] = (d2 < 0) * d2;//delta2[i] + Ylimit[i] - Z2[i] + Y[i];
	}
}


__global__ void updateZDeltaGPU(float* Z1, float* Z2, float* Ylimit, float* delta1, float* delta2, float* Y, int nVarPF) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;

	for (int i = index; i < nVarPF; i += step) {

		float z1 = Ylimit[i] + delta1[i] - Y[i];
		float z2 = Ylimit[i] + delta2[i] + Y[i];
		Z1[i] = (z1 > 0) * z1;
		Z2[i] = (z2 > 0) * z2;
		delta1[i] = (z1 < 0) * z1;
		delta2[i] = (z2 < 0) * z2;
	}

}



/*
for (int n = 0; n < _nAgent; n++) {
		float sum = 0;
		for (int i = 0; i < _nVarPF; i++) {
			sum += SensiBis.get(i, 0) * G.get(i, n);
		}
		Cp2.set(n, 0, sum * _rho1 * nVoisin.get(n, 0));
	}

*/

template <unsigned int blockSize>
__global__ void updateCp2GPU(float* Cp2, float* SensiBis, float* G, float* nVoisin, float rho1, int nVarPF) {
	int index = threadIdx.x;
	int step = blockDim.x;
	int agent = blockIdx.x; // un bloc par agent
	__shared__ float shArr[blockSize];

	float sum = 0;
	for (int i = index; i < nVarPF; i += step) {
		sum += SensiBis[i] * G[agent * nVarPF + i];
	}
	shArr[index] = sum;
	__syncthreads();

	if (blockSize >= 512) { if (index < 256) { shArr[index] += shArr[index + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (index < 128) { shArr[index] += shArr[index + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (index < 64) { shArr[index] += shArr[index + 64]; } __syncthreads(); }
	if (index < 32) {
		warpReduce<blockSize>(shArr, index);
	}
	if (index == 0) {
		Cp2[agent] = shArr[0] * rho1 * nVoisin[agent];
	}

}


/*
	SensiBis.set(&Y);
	SensiBis.multiply(2);
	SensiBis.add(&Z1);
	SensiBis.subtract(&Z2);
	SensiBis.add(&delta2);
	SensiBis.subtract(&delta1);
*/

__global__ void updateSensiBis(float* sensiBis, float* Y, float* Z1, float* Z2, float* delta1, float* delta2, int nVarPF) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;

	for (int i = index; i < nVarPF; i += step) {
		sensiBis[i] = 2 * Y[i] + Z1[i] - Z2[i] + delta2[i] - delta1[i];
	}

}



/*

	for (int i = 0; i < _nVarPF; i++) {
		for (int n = 1; n < _nAgent; n++) {
			if (abs(dP.get(n,0)) > 0.01) {
				G.set(i, n, dY.get(i, 0) / dP.get(n, 0));
			}
		}
	}

*/

__global__ void updateSensiGPU(float* G, float* dY, float* dP, int nVar) {
	int index = threadIdx.x;
	int step = blockDim.x;
	int agent = blockIdx.x; // un bloc par agent
	float dp = dP[agent];

	if (dp > 0.01 || dp < -0.01) {
		for (int i = index; i < nVar; i += step) {
			G[agent * nVar + i] = dY[i] / dp;
		}
	}



}







/*
template <unsigned int blockSize>
__global__ void updateTradePGPUSharedResidual(float* Tlocal, float* Tlocal_pre, float* Tmoy, float* P, float* MU, float* nVoisin, float at1, float at2, float* Bt1, float* Ct,
	float* matLb, float* matUb, float* Ap1, float* Ap12, float* Cp, float* Pmin, float* Pmax, float* CoresAgentLin, float eps, int nStepL) {

	//Definition de toutes les variables locales
	int i = blockIdx.x; // c'est aussi l'identifiant de l'agent !
	unsigned int thIdx = threadIdx.x;
	const int step = blockSize;
	// ne change pas


	float Bt1local[NMAXPEERPERTRHREAD];
	float Ctlocal[NMAXPEERPERTRHREAD];
	float matUblocal[NMAXPEERPERTRHREAD];
	float matLblocal[NMAXPEERPERTRHREAD];

	float Tlocallocal[NMAXPEERPERTRHREAD]; // change
	float Tlocalprelocal[NMAXPEERPERTRHREAD]; // change
	float sum;
	float bp, MULOCAL, moy, p;
	float m, r, ub, lb, t;
	// le changement doit �tre partag� par tous les threads du bloc

	__shared__ float MuShared;
	__shared__ float TMoyShared;
	__shared__ float PShared;


	// constant et commun � tous les thread d'un bloc
	__shared__ float Ap1Shared;
	__shared__ float CpShared;
	__shared__ float Ap12Shared;
	__shared__ float PmaxShared;
	__shared__ float PminShared;
	__shared__ float nVoisinShared;
	__shared__ float at1Shared;
	__shared__ float at2Shared;
	__shared__ float at12Shared;
	__shared__ bool mustContinue;


	if (thIdx == 0) {
		Ap1Shared = Ap1[i];
		CpShared = Cp[i];
		Ap12Shared = Ap12[i];
		PmaxShared = Pmax[i];
		PminShared = Pmin[i];
		nVoisinShared = nVoisin[i];
		at1Shared = at1;
		at2Shared = at2;
		at12Shared = at1 + at2;
		MuShared = MU[i];
		TMoyShared = Tmoy[i];
		PShared = P[i];
		mustContinue = false;
	}
	int k = 0;
	__syncthreads();
	const int CoresAgentLinLocal = CoresAgentLin[i];
	const int beginLocal = CoresAgentLinLocal + thIdx;
	const int endLocal = CoresAgentLinLocal + nVoisinShared;
	float res;
	for (int j = beginLocal; j < endLocal; j += step) {
		Bt1local[k] = Bt1[j];
		Ctlocal[k] = Ct[j];
		matUblocal[k] = matUb[j];
		matLblocal[k] = matLb[j];
		//Tlocalprelocal[k] = Tlocal_pre[j];
		Tlocallocal[k] = Tlocal_pre[j];
		k = k + 1;
	}

	__shared__ float shArr[blockSize];

	//Calcul des it�rations

	for (int iter = 0; iter < nStepL; iter++) {

		MULOCAL = MuShared; // tous lisent le m�me : broadcast !
		moy = TMoyShared;
		p = PShared;
		sum = 0;
		k = 0;
		for (int j = beginLocal; j < endLocal; j += step) {
			Tlocalprelocal[k] = Tlocallocal[k];
			m = Tlocallocal[k] - moy + p - MULOCAL;
			r = (Bt1local[k] * at1Shared + m * at2Shared - Ctlocal[k]) / (at12Shared);
			ub = matUblocal[k];
			lb = matLblocal[k];
			t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
			Tlocallocal[k] = t;
			sum += t;
			res = (t - Tlocalprelocal[k]);
			res = (double) res*res;
			if (res > eps) {
				mustContinue = true; // pas de race condition, car l'ordre n'importe pas,
				//mais est ce que cela ne va pas physiquement bloquer ?
			}
			k = k + 1;
		}

		shArr[thIdx] = sum;
		__syncthreads();
		if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
		if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
		if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
		if (thIdx < 32) {
			warpReduce<blockSize>(shArr, thIdx);
		}
		__syncthreads();

		if (thIdx == 0) {
			moy = shArr[0] / nVoisinShared;
			TMoyShared = moy;
			bp = moy + MuShared;
			p = (Ap1Shared * bp - CpShared) / Ap12Shared;
			p = (PmaxShared - p) * (p > PmaxShared) + (PminShared - p) * (p < PminShared) + p;
			PShared = p;
			res = p - moy;
			res = (double) res* res;
			if (res > eps) {
				mustContinue = true;
			}
			MuShared = MULOCAL + moy - p;
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
	//Ecriture des it�rations
	__syncthreads();
	k = 0;
	for (int j = beginLocal; j < endLocal; j += step) {
		Tlocal[j] = Tlocallocal[k];
		Tlocal_pre[j] = Tlocalprelocal[k];
		k = k + 1;
	}
	if (thIdx == 0) {
		Tmoy[blockIdx.x] = TMoyShared;// TMoyShared;
		P[blockIdx.x] = PShared;// PShared;
		MU[blockIdx.x] = MuShared;// MuShared;
	}

}











*/









