#include "../head/ADMMGPUConst5.cuh"

ADMMGPUConst5::ADMMGPUConst5() : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "Constructeur ADMMGPUConst5" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilisé pendant la simu

}


ADMMGPUConst5::ADMMGPUConst5(float rho) : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "Constructeur ADMMGPUConst5 defaut" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
	timePerBlock = MatrixCPU(1, 8, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 8, 0); //nb de fois utilisé pendant la simu

}

ADMMGPUConst5::~ADMMGPUConst5()
{

}

void ADMMGPUConst5::setParam(float rho)
{
	_rho = rho;
}

void ADMMGPUConst5::setTau(float tau)
{
	if (tau < 1) {
		throw std::invalid_argument("tau must be greater than 1");
	}
	_tau = tau;
}

void ADMMGPUConst5::init(const Simparam& sim, const StudyCase& cas)
{
	// intitilisation des matrixs et variables 
	_rhog = sim.getRho();
	_rho1 = sim.getRho1();
	//std::cout << "rho initial " << _rhog << std::endl;
	_nAgent = sim.getNAgent();
	
	_rhol = _rho;
	if (_rho == 0) {
		_rhol = _rhog;
	}
	//std::cout << _rhog << " " << _rho1 << " " << _rhol << std::endl;
	const int iterG = sim.getIterG();
	const int stepG = sim.getStepG();
	float epsG = sim.getEpsG();
	float epsGC = sim.getEpsGC();
	_ratioEps = epsG / epsGC;
	
	nVoisinCPU = cas.getNvoi();
	nVoisin = MatrixGPU(nVoisinCPU, 1);
	nVoisin.preallocateReduction();

	int nVoisinMax = nVoisin.max2();
	if (_blockSize * NMAXPEERPERTRHREAD < nVoisinMax) {
		std::cout << _blockSize << " " << NMAXPEERPERTRHREAD << " " << nVoisinMax << std::endl;
		throw std::invalid_argument("For this Method, an agent must not have more than 5120 peers");
	}

	_nLine = cas.getNLine();
	
	_nBus = cas.getNBus();

	_nTrade = nVoisin.sum();
	_numBlocksN = ceil((_nAgent + _blockSize - 1) / _blockSize);
	_numBlocksM = ceil((_nTrade + _blockSize - 1) / _blockSize);
	_numBlocksL = ceil((_nLine + _blockSize - 1) / _blockSize);
	_numBlocksNL = ceil((_nAgent * _nLine + _blockSize - 1) / _blockSize);
	_at1 = _rhog; // represente en fait 2*a
	_at2 = _rhol;

	resF = MatrixCPU(3, (iterG / stepG) + 1);
	resX = MatrixCPU(4, (iterG / stepG) + 1);

	MatrixCPU BETA(cas.getBeta());
	MatrixGPU Ub(cas.getUb());
	MatrixGPU Lb(cas.getLb());
	LAMBDA = sim.getLambda();
	trade = sim.getTrade();
	
	//std::cout << "mise sous forme linéaire" << std::endl;
	// Rem : si matrice déjà existante, elles sont déjà sur GPU donc bug pour les get
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




	CoresMatLin = MatrixGPU(_nAgent, _nAgent, -1);
	
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

	for (int idAgent = 0; idAgent < _nAgent; idAgent++) {
		MatrixCPU omega(cas.getVoisin(idAgent));
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
	for (int lin = 0; lin < _nTrade; lin++) {
		int i = CoresLinAgent.get(lin, 0);
		int j = CoresLinVoisin.get(lin, 0);
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
	//std::cout << "donnees sur GPU pour le grid" << std::endl;
	Kappa1 = MatrixGPU(_nLine, 1, 0, 1);
	Kappa2 = MatrixGPU(_nLine, 1, 0, 1);
	Kappa1_pre = MatrixGPU(_nLine, 1, 0, 1);
	Kappa2_pre = MatrixGPU(_nLine, 1, 0, 1);
	Qpart = MatrixGPU(_nAgent, _nLine, 0, 1);
	Qtot = MatrixGPU(_nLine, 1, 0, 1);
	alpha = MatrixGPU(_nAgent, _nLine, 0, 1);

	G = MatrixGPU(cas.getPowerSensi());
	lLimit = MatrixGPU(cas.getLineLimit(), 1);
	
	GTrans = MatrixGPU(_nAgent, _nLine);
	if (GTrans.getPos()) {
		GTrans.transferCPU();
	}
	GTrans.setTrans(&G);
	//G.transferGPU();
	GTrans.transferGPU();
	G2 = GTrans;
	G2.multiplyT(&GTrans);

	//std::cout << "autres donnée sur GPU" << std::endl;
	tempNN = MatrixGPU(_nTrade, 1, 0, 1);
	tempN1 = MatrixGPU(_nAgent, 1, 0, 1); // plutôt que de re-allouer de la mémoire à chaque utilisation
	tempL1 = MatrixGPU(_nLine, 1, 0, 1);
	tempL2 = MatrixGPU(_nLine, 1, 0, 1);
	//MatrixGPU temp1N(1, _nAgent, 0, 1);

	Tlocal = MatrixGPU(_nTrade, 1, 0, 1);
	P = MatrixGPU(_nAgent, 1, 0, 1); // moyenne des trades
	Pn = MatrixGPU(sim.getPn(), 1); // somme des trades


	a = MatrixGPU(cas.geta(), 1);
	b = MatrixGPU(cas.getb(), 1);
	Ap2 = a;
	Ap1 = nVoisin;
	Ap12 = MatrixGPU(_nAgent, 1, 0, 1);
	Ap2a = a;
	Ap2b = MatrixGPU(_nAgent, 1, 0, 1);

	Bt1 = MatrixGPU(_nTrade, 1, 0, 1);
	Cp = MatrixGPU(_nAgent, 1, 0, 1);
	Cp2 = MatrixGPU(_nAgent, 1, 0, 1);
	Cp1 = b;

	Pmin = MatrixGPU(cas.getPmin(), 1);
	Pmax = MatrixGPU(cas.getPmax(), 1);
	MU = MatrixGPU(sim.getMU(), 1); // facteur reduit i.e lambda_l/_rho
	Tmoy = MatrixGPU(sim.getPn(), 1);

	tempNN.preallocateReduction();
	Tlocal.preallocateReduction();
	tempL1.preallocateReduction();

	P.preallocateReduction();


	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);
	Ap1.multiply(_rhol);
	Cp1.multiplyT(&nVoisin);
	Tmoy.divideT(&nVoisin);


	tempN1.sum(&G2);
	tempN1.multiply(2 * _rho1);
	Ap2.add(&tempN1);
	Ap2.multiplyT(&nVoisin);
	Ap2.multiplyT(&nVoisin);
	Ap12.add(&Ap1, &Ap2);
	
	
	
	updateGlobalProbGPU();

}



void ADMMGPUConst5::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	matLb.transferCPU();
	//std::cout << "update" << std::endl;
	Pmin = MatrixGPU(cas.getPmin());
	Pmax = MatrixGPU(cas.getPmax());


	MatrixGPU Lb(cas.getLb());

	b = cas.getb();
	Cp1 = cas.getb();
	int indice = 0;

	for (int idAgent = 0; idAgent < _nAgent; idAgent++) {
		int Nvoisinmax = nVoisinCPU.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
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
	timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION

	//std::cout << "fin update temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
}


void ADMMGPUConst5::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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

	if (_id == 0) {
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
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
	//std::cout << _numBlocks2 << " " <<  _blockSize << std::endl;
	
	_at1 = _rhog; // represente en fait 2*a
	
	float epsG = sim.getEpsG();
	float epsL = sim.getEpsL();
	const int stepL = sim.getStepL();
	const int stepG = sim.getStepG();
	const int iterG = sim.getIterG();
	const int iterL = sim.getIterL();
	

	float resG = 2 * epsG;
	float epsL2 = epsL * epsL;
	int iterGlobal = 0;
	int iterLocal = 0;
	int realOccurence = 0;
	
	//std::cout << iterG << " " << iterL << " " << epsL << " " << epsG << std::endl;
	while ((iterGlobal < iterG) && (resG > epsG)) {
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

		updateLocalProbGPU(epsL2, iterL);
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 1, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION

		tradeLin.swap(&Tlocal); // echange juste les pointeurs	


		updateGlobalProbGPU();
		if (!(iterGlobal % stepG)) {
#ifdef INSTRUMENTATION
			cudaDeviceSynchronize();
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateRes(&resF, &Tlocal, iterGlobal / stepG, &tempNN);
#ifdef INSTRUMENTATION
			cudaDeviceSynchronize();
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 6, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		}
		//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;
		iterGlobal++;
	}
#ifdef INSTRUMENTATION
	occurencePerBlock.increment(0, 1, iterGlobal);
	occurencePerBlock.increment(0, 3, iterGlobal);
	occurencePerBlock.increment(0, 4, iterGlobal);
	occurencePerBlock.increment(0, 5, iterGlobal);
	occurencePerBlock.increment(0, 6, iterGlobal / stepG);

	cudaDeviceSynchronize();
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	
	Kappa1.projectNeg(); //delta1
	Kappa2.projectNeg(); // delta2

	float fc = calcFc(&a, &b, &tradeLin, &Pn, &Ct, &tempN1, &tempNN);
	//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resG << std::endl;
	MatrixCPU tradeLinCPU;
	tradeLin.toMatCPU(tradeLinCPU);
	MatrixCPU LAMBDALinCPU;
	LAMBDALin.toMatCPU(LAMBDALinCPU);
	MatrixCPU PnCPU;
	Pn.toMatCPU(PnCPU);
	MatrixCPU MUCPU;
	MU.toMatCPU(MUCPU);
	MatrixCPU delta1CPU;
	Kappa1.toMatCPU(delta1CPU);
	MatrixCPU delta2CPU;
	Kappa2.toMatCPU(delta2CPU);

	int indice = 0;
	for (int idAgent = 0;idAgent < _nAgent; idAgent++) {
		MatrixCPU omega(cas.getVoisin(idAgent));
		int Nvoisinmax = nVoisinCPU.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			int idVoisin = omega.get(voisin, 0);
			trade.set(idAgent, idVoisin, tradeLinCPU.get(indice, 0));
			LAMBDA.set(idAgent, idVoisin, LAMBDALinCPU.get(indice, 0));
			indice = indice + 1;
		}
	}
	result->setResF(&resF);
	result->setLAMBDA(&LAMBDA);
	result->setTrade(&trade);
	result->setDelta1(&delta1CPU);
	result->setDelta2(&delta2CPU);
	result->setIter(iterGlobal);
	result->setPn(&PnCPU);
	result->setFc(fc);
	result->setMU(&MUCPU);
	result->setRho(_rhog);
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 7, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 7, 1);
	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION

	tall = clock() - tall;
	result->setTime((float)tall / CLOCKS_PER_SEC);
}

void ADMMGPUConst5::updateLocalProbGPU(float epsL, int nIterL) {
	int numBlocks = _nAgent;
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
	//cudaStreamSynchronize(streamCalculation);
}



void ADMMGPUConst5::updateGlobalProbGPU()
{
	//Rem : tout calcul qui est de taille N ou M peut être fait par les agents
		// Si le calcul est de taile L, soit c'est calculé par un/des superviseurs, soit tous les agents le calcul (un peu absurde)

#ifdef INSTRUMENTATION
// FB 3a
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION


	updatePnGPU << <_numBlocksN, _blockSize >> > (Pn._matrixGPU, Tmoy._matrixGPU, nVoisin._matrixGPU, _nAgent);
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	timePerBlock.increment(0, 3, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 3b
	cudaDeviceSynchronize();
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION	
	updateAlphaTrans << < _numBlocksNL, _blockSize >> > (alpha._matrixGPU, GTrans._matrixGPU, Pn._matrixGPU, _nLine, _nAgent);
	updateQpartTrans << < _nLine, _blockSize, _nAgent * sizeof(float) >> > (Qpart._matrixGPU, alpha._matrixGPU, _nAgent, _nLine);
	updateQtotTrans << <_numBlocksL, _blockSize >> > (Qtot._matrixGPU, Qpart._matrixGPU, alpha._matrixGPU, _nLine);


	
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 4, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 3c
	cudaDeviceSynchronize();
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	Kappa1_pre.set(&Kappa1);
	Kappa2_pre.set(&Kappa2);
	updateKappaGPU << <_numBlocksL, _blockSize >> > (Kappa1._matrixGPU, Kappa2._matrixGPU, lLimit._matrixGPU, Qtot._matrixGPU, _nLine);
	diffKappa << <_numBlocksL, _blockSize >> > (tempL1._matrixGPU, Kappa1._matrixGPU, Kappa2._matrixGPU, _nLine);
	int numBlocks = _nAgent;
	switch (_blockSize) {
	case 512:
		updateCp2GPUTrans<512> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		break;
	case 256:
		updateCp2GPUTrans<256> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		break;
	case 128:
		updateCp2GPUTrans<128> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		break;
	case 64:
		updateCp2GPUTrans<64> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		break;
	case 32:
		updateCp2GPUTrans<32> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		break;
	case 16:
		updateCp2GPUTrans<16> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		break;
	case  8:
		updateCp2GPUTrans<8> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		break;
	case  4:
		updateCp2GPUTrans<4> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		break;
	case  2:
		updateCp2GPUTrans<2> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		break;
	case  1:
		updateCp2GPUTrans<1> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		break;
	}



	updateLAMBDABt1GPU << <_numBlocksM, _blockSize >> > (Bt1._matrixGPU, LAMBDALin._matrixGPU, tradeLin._matrixGPU, _rhog, CoresLinTrans._matrixGPU, _nTrade);
	updateCp << <_numBlocksN, _blockSize >> > (Cp._matrixGPU, Cp1._matrixGPU, Cp2._matrixGPU, _nAgent);

#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 5, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
}



float ADMMGPUConst5::updateRes(MatrixCPU* res, MatrixGPU* Tlocal, int iter, MatrixGPU* tempNN)
{

	float resS = Tlocal->max2(&tradeLin);

	updateDiffGPU <<<_numBlocksM, _blockSize >> > (tempNN->_matrixGPU, Tlocal->_matrixGPU, CoresLinTrans._matrixGPU, _nAgent);
	float resR = tempNN->max2();

	// version de l'article
	/*tempL1.set(&Kappa1);
	tempL2.set(&Kappa2);
	Kappa1_pre.projectNeg();
	Kappa2_pre.projectNeg();
	tempL1.projectNeg();
	tempL2.projectNeg();
	tempL1.subtract(&Kappa1_pre);
	tempL2.subtract(&Kappa2_pre);
	tempL1.multiplyT(&tempL1);
	tempL2.multiplyT(&tempL2);
	tempL1.add(&tempL2);*/
	updateResX << <_numBlocksL, _blockSize >> > (tempL1._matrixGPU, Kappa1._matrixGPU, Kappa2._matrixGPU, Kappa1_pre._matrixGPU, Kappa2_pre._matrixGPU, _nLine);

	float resXf = _ratioEps * sqrt(tempL1.max2());
	res->set(0, iter, resR);
	res->set(1, iter, resS);
	res->set(2, iter, resXf);
	return MAX(MAX(resXf, resS), resR);

}





void ADMMGPUConst5::display() {

	std::cout << _name << std::endl;
}

