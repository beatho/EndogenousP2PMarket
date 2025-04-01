#include "../head/MarEndoConsGPU.cuh"




MarEndoConsGPU::MarEndoConsGPU() : MethodP2PGPU()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " MarEndoConsGPU Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu
}


MarEndoConsGPU::MarEndoConsGPU(float rho) : MethodP2PGPU()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default MarEndoConsGPU Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu
}

MarEndoConsGPU::~MarEndoConsGPU()
{
	DELETEB(OPF);
	//DELETEB(OPFCPU);
}


void MarEndoConsGPU::setParam(float rho)
{
	_rho = rho;
}

void MarEndoConsGPU::setTau(float tau)
{
	if (tau < 1) {
		throw std::invalid_argument("tau must be greater than 1");
	}
	_tau = tau;
}



void MarEndoConsGPU::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
		timePerBlock.set(0, 0, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		occurencePerBlock.set(0, 0, 1);
#endif // INSTRUMENTATION
	}
	
	_rhog = sim.getRho();
	_at1 = _rhog;
	_iterG = sim.getIterG();
	int iterL = sim.getIterL();
	int stepL = sim.getStepL();

	
	float epsL = sim.getEpsL() / 20;
	float epsG = sim.getEpsG();
	

	float fc = 0;

	
	_resG = 2 * epsG;
	float epsL2 = epsL * epsL;

	_iterGlobal = 0;
	while ((_iterGlobal < _iterG) && (_resG > epsG)  || (_iterGlobal <=_stepG)) { // || (_iterGlobal <= _stepG)
		//std::cout << "*";
		/*std::cout << "---------------------------------" << std::endl;
		std::cout << " Pn " << std::endl;
		Pn.display(true);
		std::cout << " Pso " << std::endl;
		PSO.display(true);
		std::cout << " Bp3 " << std::endl;
		Bp3.display(true);*/

		/*std::cout << " Tlocal " << std::endl;
		Tlocal.display(true);
		std::cout << " Bt1 " << std::endl;
		Bt1.display(true);
		std::cout << " matlb " << std::endl;
		matLb.display(true);
		std::cout << " matUb " << std::endl;
		matUb.display(true); */
		

#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		//CHECK_LAST_CUDA_ERROR();
		updateLocalProbGPU(epsL2, iterL);
		//CHECK_LAST_CUDA_ERROR();
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 1, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

#endif // INSTRUMENTATION
		//std::cout << _iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, _iterGlobal / _stepG) << " " << resF.get(1, _iterGlobal / _stepG) << std::endl;

		tradeLin.swap(&Tlocal);
		updateGlobalProb();
		//CHECK_LAST_CUDA_ERROR();
		// FB 4
		if (!(_iterGlobal % _stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			_resG = updateResBis( _iterGlobal / _stepG );
			//std::cout << _iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, _iterGlobal / _stepG)
				//<< " " << resF.get(1, _iterGlobal / _stepG) << " " << resF.get(2, _iterGlobal / _stepG) << std::endl;
#ifdef INSTRUMENTATION
			cudaDeviceSynchronize();
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION		
		}
		//CHECK_LAST_CUDA_ERROR();

		_iterGlobal++;
	}
	
	//std::cout << _iterGlobal << " " << iterLocal  << " " << resF.get(0, (_iterGlobal - 1) / _stepG) << " " << resF.get(1, (_iterGlobal - 1) / _stepG) << " " << resF.get(2, (_iterGlobal - 1) / _stepG) << std::endl;
#ifdef INSTRUMENTATION
	occurencePerBlock.increment(0, 1, _iterGlobal);
	occurencePerBlock.increment(0, 3, _iterGlobal);
	occurencePerBlock.increment(0, 4, _iterGlobal);
	occurencePerBlock.increment(0, 5, _iterGlobal);
	occurencePerBlock.increment(0, 6, _iterGlobal / _stepG);

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	
	fc = calcFc();
	
	MatrixCPU tradeLinCPU;
	tradeLin.toMatCPU(tradeLinCPU);
	MatrixCPU LAMBDALinCPU;
	LAMBDALin.toMatCPU(LAMBDALinCPU);
	MatrixCPU PnCPU;
	Pn.toMatCPU(PnCPU);
	MatrixCPU MUCPU;
	MU.toMatCPU(MUCPU);


	int indice = 0;
	for (int idAgent = 0; idAgent < _nAgentTrue; idAgent++) {
		MatrixGPU omega(cas.getVoisin(idAgent));
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
	//trade.display();
	/*std::cout << "PSO, Pn" << std::endl;
	PSO.display();
	Pn.display();*/

	// FB 5
	
	MatrixCPU Pb(OPF->getPb());
	MatrixCPU Phi(OPF->getPhi());
	MatrixCPU E(OPF->getE());
	
	result->setE(&E);
	result->setPhi(&Phi);
	result->setPb(&Pb);

	result->setResF(&resF);
	
	result->setLAMBDA(&LAMBDA);
	
	result->setTrade(&trade);
	result->setIter(_iterGlobal);
	
	result->setPn(&PnCPU);
	result->setFc(fc);
	result->setMU(&MUCPU);

#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 7, 1);
	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
	timeMarketEndo = clock() - tall;
	result->setTime((float) timeMarketEndo / CLOCKS_PER_SEC);


	//std::cout << "****" << std::endl;
	//OPF->display();
	//display();
}

void MarEndoConsGPU::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	
	matLb.transferCPU();
	matUb.transferCPU();

	Pmin = cas.getPmin();
	Pmax = cas.getPmax();


	MatrixCPU Lb(cas.getLb());
	MatrixCPU Ub(cas.getUb());

	b = cas.getb();
	Cp = b;
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

	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);
	Cp.multiplyT(&nVoisin);
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION


}

void MarEndoConsGPU::init(const Simparam& sim, const StudyCase& cas)
{
	DELETEB(OPF);
	DELETEB(OPFCPU);
	


	if (CoresMatLin.getPos()) { // une copie en trop mais pour l'instant c'est ok...
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
	//std::cout << "init " << std::endl;
	//cudaDeviceSynchronize();
	//CHECK_LAST_CUDA_ERROR();
	// intitilisation des matrixs et variables 
	
	if (!cas.isAC()) {
		throw std::invalid_argument("Wrong studyCase must be AC");
	}
	if (!cas.isRadial()) {
		throw std::invalid_argument("Wrong studyCase must be radial, dont have OPF on non-radial cases");
	}
	////CHECK_LAST_CUDA_ERROR();
	_rhog = sim.getRho();
	_rhoSO = _rhog;
	//_rhoSO = sim.getRho1();
	_iterG = sim.getIterG();
	_stepG = sim.getStepG();
	float epsG = sim.getEpsG();
	float epsGC = sim.getEpsGC();
	_epsLim = sim.getEpsIntern();
	_stepL = sim.getStepL();
	_ratioEps = epsG / epsGC;
	_nAgentTrue = sim.getNAgent();
	if (_nAgentTrue != cas.getNagent()) {
		throw std::invalid_argument("nAgent different on Simparam and study case");
	}
	_nAgent = 2 * _nAgentTrue;

	paramOPF = sim;
	paramOPF.setItG(sim.getIterIntern());
	paramOPF.setEpsG(sim.getEpsIntern());
	_stepIntern = sim.getStepIntern();


	_rhol = _rho; //*nAgent
	
	//std::cout << "rho " << _rho << std::endl;
	if (_rho == 0) {
		_rhol = _rhog;
	}
	
	CHECK_LAST_CUDA_ERROR();
	nVoisinCPU = cas.getNvoi();
	nVoisin = MatrixGPU(nVoisinCPU, 1);
	nVoisin.preallocateReduction();
	int nVoisinMax = nVoisin.max2();
	if (_blockSize * NMAXPEERPERTRHREAD < nVoisinMax) {
		std::cout << _blockSize << " " << NMAXPEERPERTRHREAD << " " << nVoisinMax << std::endl;
		throw std::invalid_argument("For this Method, an agent must not have more than 5120 peers");
	}
	//std::cout << "Trades" << std::endl;
	_nTrade = nVoisin.sum();
	_nTradeP = nVoisin.sum(0, _nAgentTrue);
	
	_nTradeQ = _nTrade - _nTradeP;
	if (_nTradeQ != (_nAgentTrue * (_nAgentTrue - 1))) {
		std::cout << "err MarEndoGPU : " << _nAgent << " " << _nAgentTrue << " " << _nTrade << " " << _nTradeP << " " << _nTradeQ << std::endl;
		throw std::invalid_argument("Agent must be fully conected for the Q echanges, WIP");
	}
	_numBlocksN = ceil((_nAgent + _blockSize - 1) / _blockSize);
	_numBlocksM = ceil((_nTrade + _blockSize - 1) / _blockSize);

	/*cudaDeviceSynchronize();
	//CHECK_LAST_CUDA_ERROR();*/
	//std::cout << "Market" << std::endl;
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	
	MatrixCPU PnCPU(2 * _nAgent, 1);
	if (initWithMarketClear) {
		ADMMMarketGPU market;
		Simparam res(sim);
		market.solve(&res, sim, cas);
		//res.display();
		LAMBDA = res.getLambda();
		trade = res.getTrade();
		PnCPU = res.getPn();
	}
	else {
		LAMBDA = sim.getLambda();
		trade = sim.getTrade();
		PnCPU = sim.getPn();
	}
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	
	//std::cout << "time : " <<  (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / 1000000000 << std::endl;
	//PnCPU.display();
	Pn = MatrixGPU(PnCPU, 1);
	//Pn.display(true);
	paramOPF.setPn(&PnCPU);
	
	_at1 = _rhog; 
	_at2 = _rhol;

	resF = MatrixCPU(3, (_iterG / _stepG) + 1);
	
	/*cudaDeviceSynchronize();
	//CHECK_LAST_CUDA_ERROR();*/

	MatrixGPU BETA(cas.getBeta());
	
	MatrixGPU Ub(cas.getUb());
	MatrixGPU Lb(cas.getLb());
	
	
	//std::cout << "mise sous forme lin�aire" << std::endl;
	
	CoresMatLin = MatrixGPU(_nAgent, _nAgentTrue, -1);
	CoresAgentLin = MatrixGPU(_nAgent + 1, 1);
	CoresLinAgent = MatrixGPU(_nTrade, 1);
	CoresLinVoisin = MatrixGPU(_nTrade, 1);
	CoresLinTrans = MatrixGPU(_nTrade, 1);

	Tlocal_pre = MatrixGPU(_nTrade, 1);
	tradeLin = MatrixGPU(_nTrade, 1);
	LAMBDALin = MatrixGPU(_nTrade, 1);

	matLb = MatrixGPU(_nTrade, 1);
	matUb = MatrixGPU(_nTrade, 1);
	Ct = MatrixGPU(_nTrade, 1);
	//Bt2 = MatrixGPU(_nTrade, 1);
	/*cudaDeviceSynchronize();
	//CHECK_LAST_CUDA_ERROR();*/
	int indice = 0;

	for (int idAgent = 0; idAgent < _nAgentTrue; idAgent++) { // P
		MatrixGPU omega(cas.getVoisin(idAgent));
		int Nvoisinmax = nVoisinCPU.get(idAgent, 0);
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
		CoresAgentLin.set(idAgent + 1, 0, indice);
	}
	//std::cout << " Q " << std::endl;
	for (int idAgent = _nAgentTrue; idAgent < _nAgent; idAgent++) { // Q
		for (int idVoisin = 0; idVoisin < _nAgentTrue; idVoisin++) {
			if (idVoisin != (idAgent - _nAgentTrue)) {
				if(Lb.getNCol()==1){
					matLb.set(indice, 0, Lb.get(idAgent, 0));
					matUb.set(indice, 0, Ub.get(idAgent, 0));
				} else {
					matLb.set(indice, 0, Lb.get(idAgent, idVoisin));
					matUb.set(indice, 0, Ub.get(idAgent, idVoisin));
				}
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
	CHECK_LAST_CUDA_ERROR();
	
	//CHECK_LAST_CUDA_ERROR();
	//std::cout << "donnees sur CPU pour le grid" << std::endl;

	Ap3 = nVoisin;
	Ap3.multiplyT(&nVoisin);
	Ap3.multiply(_rhoSO);
	PSO = MatrixGPU(_nAgent, 1, 0, 1);
	PSO.preallocateReduction();
	etaSO = MatrixGPU(_nAgent, 1, 0, 1);
	Bp3 = MatrixGPU(_nAgent, 1, 0, 1);
	_radial = cas.isRadial();

	 
	if (OPFonCPU) {
		OPFCPU = new OPFADMMCons;
		OPFCPU->initConsensus(paramOPF, cas, _rhoSO);
	}
	else {
		OPF = new OPFADMMConsGPU;
		OPF->initConsensus(paramOPF, cas, _rhoSO);
	}
	CHECK_LAST_CUDA_ERROR();

	
	 

	
	//CHECK_LAST_CUDA_ERROR();
	//std::cout << "autres donn�e sur CPU" << std::endl;
	tempNN = MatrixGPU(_nTrade, 1, 0, 1);
	tempNN.preallocateReduction();
	tempN1 = MatrixGPU(_nAgent, 1, 0, 1); // plut�t que de re-allouer de la m�moire � chaque utilisation
	

	Tlocal = MatrixGPU(_nTrade, 1, 0, 1);
	Tlocal.preallocateReduction();
	P = Pn; // moyenne des trades
	P.divideT(&nVoisin);

	a = MatrixGPU(cas.geta(), 1);
	b = MatrixGPU(cas.getb(), 1);

	// on enn veut pas que l'agent des pertes consomme plus que n�cessaire !!!
	//a.set(0, 0, 1);
	//a.set(_nAgentTrue, 0, 1);

	//CHECK_LAST_CUDA_ERROR();
	Ap2 = a;
	Ap1 = nVoisin;
	Ap123 = Ap3;

	Bt1 = MatrixGPU(_nTrade, 1, 0, 1);
	Cp = b;
	
	Bp1 = MatrixGPU(_nAgent, 1, 0, 1);

	Pmin = MatrixGPU(cas.getPmin(), 1);
	
	Pmin.set(0, 0, -100000, true); // unleash power !!!
	Pmin.set(_nAgentTrue, 0, -100000, true); // unleash power !!!	
	Pmax = MatrixGPU(cas.getPmax(), 1);
	Pmax.set(_nAgentTrue, 0, 100000, true); // unleash power !!!
	
	MU = MatrixGPU(sim.getMU(), 1); // facteur reduit i.e lambda_l/_rho
	Tmoy = P;

	
	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);
	Ap1.multiply(_rhol);
	Cp.multiplyT(&nVoisin);
	Tmoy.divideT(&nVoisin);
	
	Ap2.multiplyT(&nVoisin);
	Ap2.multiplyT(&nVoisin);
	Ap123.add(&Ap1);
	Ap123.add(&Ap2);
	//CHECK_LAST_CUDA_ERROR();
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
	matUb.display(true);*/

	/*PSO.display();
	Pn.display();*/
	//std::cout << "******" << std::endl;
	//std::cout << "updateGlobal" << std::endl;
	CHECK_LAST_CUDA_ERROR();
	updateGlobalProb();
	CHECK_LAST_CUDA_ERROR();
	
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	//std::cout << "fin init " << std::endl;
}

void MarEndoConsGPU::updateGlobalProb() {

	// FB 3a
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	std::chrono::high_resolution_clock::time_point t2;
#endif // INSTRUMENTATION

	
	float eps = MYMIN(_resG * _delta, _epsLim);
	
	
	//std::cout << "SolveOPF" << std::endl;
	if (_iterGlobal % _stepIntern == 0) {
		if (OPFonCPU) {
			PSO.toMatCPU(PSOCPU);
			OPFCPU->solveConsensus(eps, &PSOCPU);
			PSO = PSOCPU;
		}
		else {
			OPF->solveConsensus(eps, &PSO);
		}
		CHECK_LAST_CUDA_ERROR();
	



#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 4, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	// FB 3b
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	//std::cout << "update OPF" << std::endl;

	updatePn();
	CHECK_LAST_CUDA_ERROR();
	if (OPFonCPU) {
		Pn.toMatCPU(PnCPU);
		OPFCPU->updateConsensus(&PnCPU);
	}
	else {
		OPF->updateConsensus(&Pn);
	}
	CHECK_LAST_CUDA_ERROR();
	//MatrixCPU PnCPU;
	

	
	//
	//Pn.display();
	//PSO.display();
	
	//
	
	
	

#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 3, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	//Agent des pertes
	//std::cout << "update Market" << std::endl;
	/*Pmin.set(0, 0, Ploss / nVoisin.get(0,0));
	Pmax.set(0, 0, Ploss / nVoisin.get(0, 0));
	
	Pmin.set(_nAgentTrue, 0, Qloss / nVoisin.get(_nAgentTrue, 0));
	Pmax.set(_nAgentTrue, 0, Qloss / nVoisin.get(_nAgentTrue, 0));*/
	

	float Qloss = PSO.get(_nAgentTrue, 0, false);
	CHECK_LAST_CUDA_ERROR();
	if (Qloss > 0) {
		matUb.setBloc(_nTradeP, _nTradeP + _nAgentTrue - 1, 0, 1, Qloss);
		matLb.setBloc(_nTradeP, _nTradeP + _nAgentTrue - 1, 0, 1, 0.0);
	}
	else {
		matLb.setBloc(_nTradeP, _nTradeP + _nAgentTrue - 1, 0, 1, Qloss);
		matUb.setBloc(_nTradeP, _nTradeP + _nAgentTrue - 1, 0, 1, 0.0);
	}
	/**/
	CHECK_LAST_CUDA_ERROR();
	// FB 3c
	
	
		updateEtaSO();
		updateBp3();
	}
	CHECK_LAST_CUDA_ERROR();
	updateLAMBDABt1GPU << <_numBlocksM, _blockSize >> > (Bt1._matrixGPU, LAMBDALin._matrixGPU, tradeLin._matrixGPU, _rhog, CoresLinTrans._matrixGPU, _nTrade);
	//Bp3.display();
	CHECK_LAST_CUDA_ERROR();
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 5, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
	
}

void MarEndoConsGPU::updateLocalProbGPU(float epsL, int nIterL)
{
	int numBlocks = _nAgent;
	switch (_blockSize) {
	case 512:
		updateTradePGPUSharedResidual<512> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU,  Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case 256:
		updateTradePGPUSharedResidual<256> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case 128:
		updateTradePGPUSharedResidual<128> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case 64:
		updateTradePGPUSharedResidual< 64> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case 32:
		updateTradePGPUSharedResidual< 32> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case 16:
		updateTradePGPUSharedResidual< 16> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case  8:
		updateTradePGPUSharedResidual<  8> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case  4:
		updateTradePGPUSharedResidual<  4> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case  2:
		updateTradePGPUSharedResidual<  2> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case  1:
		updateTradePGPUSharedResidual<  1> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	}


}


void MarEndoConsGPU::updateEtaSO()
{
	tempN1.subtract(&Pn, &PSO);
	tempN1.multiply(0.5);
	etaSO.add(&tempN1);

	/*for (int n = 0; n < _nAgent; n++) {
		float eta = 0.5 * (Pn.get(n, 0) - PSO.get(n, 0));
		etaSO.set(n, 0, etaSO.get(n, 0) + eta);
	}*/
}



void MarEndoConsGPU::updateBp3()
{
	Bp3.add(&PSO, &Pn);
	Bp3.multiply(0.5);
	Bp3.subtract(&etaSO);
	Bp3.divideT(&nVoisin);
}



float MarEndoConsGPU::updateResBis(int iter)
{
	
	//tradeLin.display();
	float resS = Tlocal.max2(&tradeLin);

	updateDiffGPU << <_numBlocksM, _blockSize >> > (tempNN._matrixGPU, Tlocal._matrixGPU, CoresLinTrans._matrixGPU, _nTrade);
	float resR = tempNN.max2();


	 
	float resXf = PSO.max2(&Pn);
	/*for (int i = 1; i < _nAgentTrue; i++) {
		resXf = MYMAX(abs(PSO.get(i,0) - Pn.get(i,0)), resXf);
		resXf = MYMAX(abs(PSO.get(i + _nAgentTrue, 0) - Pn.get(i + _nAgentTrue, 0)), resXf);
	}*/

	
	
	resF.set(0, iter, resR);
	resF.set(1, iter, resS);
	resF.set(2, iter, resXf);
	return MYMAX(MYMAX(resXf * _ratioEps, resS), resR);
}







void MarEndoConsGPU::display() {

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


// updateConsensusGPU << <_numBlocksN, _blockSize >> > (Cost2._matrixGPU, etaSO._matrixGPU, Pn._matrixGPU, Pmarket->_matrixGPU, _rhoSO, _nAgent);
