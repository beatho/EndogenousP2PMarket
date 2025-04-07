#include "../head/ADMMGPUConstCons3.cuh"

ADMMGPUConstCons3::ADMMGPUConstCons3() : MethodP2PGPU()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "Constructeur ADMMGPUConstCons3" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu

}


ADMMGPUConstCons3::ADMMGPUConstCons3(float rho) : MethodP2PGPU()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "Constructeur ADMMGPUConstCons3 defaut" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu

}

ADMMGPUConstCons3::~ADMMGPUConstCons3()
{
	if (alpha != nullptr) {
		//std::cout << "delete alpha" << std::endl;
		cudaFree(alpha);
		alpha = nullptr;
	}
}

void ADMMGPUConstCons3::setParam(float rho)
{
	_rho = rho;
}

void ADMMGPUConstCons3::setTau(float tau)
{
	if (tau < 1) {
		throw std::invalid_argument("tau must be greater than 1");
	}
	_tau = tau;
}

void ADMMGPUConstCons3::init(const Simparam& sim, const StudyCase& cas)
{
	// intitilisation des matrixs et variables 
	isAC = false;
	initSize(cas);
	initSimParam(sim);

	int nVoisinMax = nVoisin.max2();
	if (_blockSize * NMAXPEERPERTRHREAD < nVoisinMax) {
		std::cout << _blockSize << " " << NMAXPEERPERTRHREAD << " " << nVoisinMax << std::endl;
		throw std::invalid_argument("For this Method, an agent must not have more than 5120 peers");
	}

	
	L2 = 2 * _nLine;
	_Msize = _nAgent + L2 + 1;
	_Asize = L2 * _nAgent;
	//std::cout << _nAgent << " " << _nLine << " " << _Msize << std::endl;
	initCaseParam(sim, cas);
	
	//std::cout << "mise sous forme lineaire" << std::endl;
	initLinForm(cas);
		
	//std::cout << "donnees sur GPU pour le grid" << std::endl;
	if (_nLine) {
		cudaMalloc((void**)&alpha, sizeof(float));
		tempL21 = MatrixGPU(L2 + 1, 1, 0, 1);

		H = MatrixGPU(_nAgent, _nAgent, 0, 1);
		H.setEyes(_rho1);

		q = MatrixGPU(_nAgent, 1, 0, 1); // 0.5x^THx + q^T*x
		diffPso = MatrixGPU(_nAgent, 1, 0, 1);

		c = MatrixGPU(L2 + 1, 1, 0, 1); // contrainte Ax+b>0 ou = 0 pour egalit�
		Ai = MatrixGPU(L2 + 1, _nAgent, 0, 1);
		MatrixGPU ones(1, _nAgent, 1, 1);
		MatrixGPU temp(cas.getPowerSensi(), 1);
		Ai.setBloc(0, _nLine, 0, _nAgent, &temp, -1);
		Ai.setBloc(_nLine, L2, 0, _nAgent, &temp);
		Ai.setBloc(L2, L2 + 1, 0, _nAgent, &ones);
		bi = MatrixGPU(L2 + 1, 1, 0, 1);
		lLimit = MatrixGPU(cas.getLineLimit(), 1);
		bi.setBloc(0, _nLine, 0, 1, &lLimit);
		bi.setBloc(_nLine, L2, 0, 1, &lLimit);

		M = MatrixGPU(_Msize, _Msize, 0, 1); // M*pas = R
		Minv = MatrixGPU(_Msize, _Msize, 0, 1); // M*pas = R
		pas = MatrixGPU(_Msize, 1, 0, 1);
		R = MatrixGPU(_Msize, 1, 0, 1);


		ZA = MatrixGPU(L2 + 1, _nAgent, 0, 1); // M = (H -Atrans ZA W)
		Z = MatrixGPU(L2 + 1, L2 + 1, 0, 1);
		Zvect = MatrixGPU(L2 + 1, 1, 0, 1);
		W = MatrixGPU(L2 + 1, L2 + 1, 0, 1);
		Wvect = MatrixGPU(L2 + 1, 1, 0, 1);
		Atrans = MatrixGPU(_nAgent, L2 + 1, 0, 1);
		Atrans.setTrans(&Ai);

		M.setBloc(0, _nAgent, 0, _nAgent, &H);
		M.setBloc(0, _nAgent, _nAgent, _Msize, &Atrans, -1);

		Rx1 = MatrixGPU(_nAgent, 1, 0, 1); // Hx+q
		Rx2 = MatrixGPU(_nAgent, 1, 0, 1); // -Ai^T*U

		Ru = MatrixGPU(L2 + 1, 1, 0, 1); // Ru = W*(U-PI)

		U = MatrixGPU(L2 + 1, 1, 0, 1);
		PI = MatrixGPU(L2 + 1, 1, 0, 1);
		
		Apas = MatrixGPU(L2 + 1, 1, 0, 1);
	}
	Pso = MatrixGPU(_nAgent, 1, 0, 1); // = Pn ? risque de non respect des contraintes
	etaP = MatrixGPU(_nAgent, 1, 0, 1); 
		
	R.preallocateReduction();
	Pso.preallocateReduction();

	//std::cout << "autres donn�e sur GPU" << std::endl;
	
	initDCEndoMarket();

	Ap3 = nVoisin;
	Bp3 = MatrixGPU(_nAgent, 1, 0, 1); // 1/Mn * (Pso + P)/2 - eta/rho1
	Ap3.multiplyT(&nVoisin);
	Ap3.multiply(_rho1);
	Ap123.add(&Ap12, &Ap3);
	
	
	updateGlobalProbGPU();

	//std::cout << " end init " << std::endl;
}

void ADMMGPUConstCons3::solveOPF()
{
	// update q
    diffPso.set(&Pso);
	q.set(&etaP);
	diffPso.add(&Pn);
	diffPso.multiply(-_rho1/2);
	q.add(&diffPso);
	
	//init
	int k = 0;
		
	float err = 2 * _epsIntern;
	mu = 10;
	float valMin = 0.0000001;
	//boucle
	while (k< _iterIntern && err > _epsIntern) {
	// update c
		
		c.linearOperation(&Ai, &Pso, &bi);
		
	// update PI
		/*c.transferCPU();
		PI.transferCPU();
		for (int l = 0; l < L2; l++) {
			if (c.get(l, 0) < valMin) {
				PI.set(l, 0, mu / valMin); // eviter division par O
			}
			else {
				PI.set(l, 0, mu / c.get(l, 0)); // eviter division par O
			}
		}
		PI.set(L2, 0, -c.get(L2, 0) / mu);
		PI.transferGPU();
		c.transferGPU();*/
		updatePI << <_numBlocksL, _blockSize >> > (PI._matrixGPU, c._matrixGPU, mu, valMin, L2);
	
	// update M
		// update Zvect
		Zvect.set(&U); 
		
		Zvect.set(L2, 0, 1, true); // egalite
		// update Z
		Z.setEyes(&Zvect);
		// update ZA
		ZA.multiplyMat(&Z, &Ai);
		// update W
		Wvect.set(&c);
		Wvect.set(L2, 0, mu, true);
		W.setEyes(&Wvect);

		M.setBloc(_nAgent, _Msize, 0, _nAgent, &ZA);
		M.setBloc(_nAgent, _Msize, _nAgent, _Msize, &W);
		try
		{
			Minv.invertGaussJordan(&M);
		}
		catch (const std::exception& e)
		{
			std::cout << e.what() << std::endl;
			float alphaCPU;
			cudaMemcpy(&alphaCPU, alpha, sizeof(float), cudaMemcpyDeviceToHost);
			std::cout << "k = " << k << " err= " << err << " alpha = " << alphaCPU  << " mu=" << mu << std::endl;
			c.display(true);
			Pn.display(true);
			Pso.display(true);
			std::cout << "---------------------------------" << std::endl;
			Pso.set(&Pn);
			return;
		}
		
	
	//update R
		// Rx
		Rx1.linearOperation(&H, &Pso, &q);
		Rx2.multiply(&Atrans, &U);
		Rx2.subtract(&Rx1);
		// Ru
		tempL21.subtract(&PI, &U);
		Ru.multiply(&W, &tempL21);

		R.setBloc(0, _nAgent, 0, 1, &Rx2);
		R.setBloc(_nAgent, _Msize, 0, 1, &Ru);
		//update pas
		pas.multiply(&Minv, &R);
		// find alpha
		findalpha();
		/*float alphaCPU;
		cudaMemcpy(&alphaCPU, alpha, sizeof(float), cudaMemcpyDeviceToHost);
		
		if (alphaCPU < valMin) { // trop proche de la froncti�re
			std::cout << "alpha " << alphaCPU << std::endl;
			c.display(true);
			Pso.display(true);
			Pn.display(true);
			break;
		}*/
		
		// update P, U
		updatePso << <_numBlocksN, _blockSize >> > (Pso._matrixGPU, pas._matrixGPU, alpha, _nAgent);
		updateU << <_numBlocksL, _blockSize >> > (U._matrixGPU, pas._matrixGPU, alpha, _nAgent, L2 + 1);

		

		// update mu
		mu *= 0.8;
		mu = MYMAX(mu, valMin);
		//R.transferCPU();
		err = R.distance2();
		//R.transferGPU();
		k++;
	}
	//std::cout << k << " " << err << std::endl;
}



void ADMMGPUConstCons3::findalpha()
{
	

	int numBlock = L2;
	switch (_blockSize) {
	case 512:
		updateAPas<512> << <numBlock, _blockSize >> > (Apas._matrixGPU, Ai._matrixGPU, pas._matrixGPU, _nAgent);
		updateAlpha<512> << <1, _blockSize >> > (alpha, U._matrixGPU, pas._matrixGPU, c._matrixGPU, Apas._matrixGPU, _nAgent, L2);
		break;
	case 256:
		updateAPas<256> << <numBlock, _blockSize >> > (Apas._matrixGPU, Ai._matrixGPU, pas._matrixGPU, _nAgent);
		updateAlpha<256> << <1, _blockSize >> > (alpha, U._matrixGPU, pas._matrixGPU, c._matrixGPU, Apas._matrixGPU, _nAgent, L2);
		break;
	case 128:
		updateAPas<128> << <numBlock, _blockSize >> > (Apas._matrixGPU, Ai._matrixGPU, pas._matrixGPU, _nAgent);
		updateAlpha<128> << <1, _blockSize >> > (alpha, U._matrixGPU, pas._matrixGPU, c._matrixGPU, Apas._matrixGPU, _nAgent, L2);
		break;
	case 64:
		updateAPas< 64> << <numBlock, _blockSize >> > (Apas._matrixGPU, Ai._matrixGPU, pas._matrixGPU, _nAgent);
		updateAlpha< 64> << <1, _blockSize >> > (alpha, U._matrixGPU, pas._matrixGPU, c._matrixGPU, Apas._matrixGPU, _nAgent, L2);
		break;
	case 32:
		updateAPas< 32> << <numBlock, _blockSize >> > (Apas._matrixGPU, Ai._matrixGPU, pas._matrixGPU, _nAgent);
		updateAlpha< 32> << <1, _blockSize >> > (alpha, U._matrixGPU, pas._matrixGPU, c._matrixGPU, Apas._matrixGPU, _nAgent, L2);
		break;
	case 16:
		updateAPas< 16> << <numBlock, _blockSize >> > (Apas._matrixGPU, Ai._matrixGPU, pas._matrixGPU, _nAgent);
		updateAlpha< 16> << <1, _blockSize >> > (alpha, U._matrixGPU, pas._matrixGPU, c._matrixGPU, Apas._matrixGPU, _nAgent, L2);
		break;
	case  8:
		updateAPas<  8> << <numBlock, _blockSize >> > (Apas._matrixGPU, Ai._matrixGPU, pas._matrixGPU, _nAgent);
		updateAlpha< 8> << <1, _blockSize >> > (alpha, U._matrixGPU, pas._matrixGPU, c._matrixGPU, Apas._matrixGPU, _nAgent, L2);
		break;
	case  4:
		updateAPas<  4> << <numBlock, _blockSize >> > (Apas._matrixGPU, Ai._matrixGPU, pas._matrixGPU, _nAgent);
		updateAlpha<  4> << <1, _blockSize >> > (alpha, U._matrixGPU, pas._matrixGPU, c._matrixGPU, Apas._matrixGPU, _nAgent, L2);
		break;
	case  2:
		updateAPas<  2> << <numBlock, _blockSize >> > (Apas._matrixGPU, Ai._matrixGPU, pas._matrixGPU, _nAgent);
		updateAlpha<  2> << <1, _blockSize >> > (alpha, U._matrixGPU, pas._matrixGPU, c._matrixGPU, Apas._matrixGPU, _nAgent, L2);
		break;
	case  1:
		updateAPas<  1> << <numBlock, _blockSize >> > (Apas._matrixGPU, Ai._matrixGPU, pas._matrixGPU, _nAgent);
		updateAlpha<  1> << <1, _blockSize >> > (alpha, U._matrixGPU, pas._matrixGPU, c._matrixGPU, Apas._matrixGPU, _nAgent, L2);
		break;
	}
	

	// truc avec alpha et beta pour optimiser sa valeur


}




void ADMMGPUConstCons3::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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

	if (_id == 0) {
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
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
	//std::cout << _numBlocks2 << " " <<  _blockSize << std::endl;
	_rhog = sim.getRho();
	_at1 = _rhog; // represente en fait 2*a
	
	
	float resG = 2 * _epsG;
	float epsL2 = _epsL * _epsL;
	_iterGlobal = 0;
	
	//std::cout << iterG << " " << iterL << " " << epsL << " " << epsG << std::endl;
	while ((_iterGlobal < _iterG) && (resG > _epsG)) {
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

		updateLocalProbGPU(epsL2, _iterL);
#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 1, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION

		tradeLin.swap(&Tlocal); // echange juste les pointeurs	
		//std::cout << "-";
		
		updateGlobalProbGPU();
		
		if (!(_iterGlobal % _stepG)) {
#ifdef INSTRUMENTATION
			cudaDeviceSynchronize();
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateResEndo(_iterGlobal / _stepG);
#ifdef INSTRUMENTATION
			cudaDeviceSynchronize();
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		}
		//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;
		_iterGlobal++;
	}
#ifdef INSTRUMENTATION
	occurencePerBlock.increment(0, 1, iterGlobal);
	occurencePerBlock.increment(0, 5, iterGlobal);
	occurencePerBlock.increment(0, 6, iterGlobal);
	occurencePerBlock.increment(0, 7, iterGlobal);
	occurencePerBlock.increment(0, 8, iterGlobal / stepG);

	cudaDeviceSynchronize();
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	
	

	//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resG << std::endl;
	std::cout << "valeur finale des contraintes de l'opf : " << std::endl;
	c.display(true);
	
	setResult(result, cas.isAC());


#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 9, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 9, 1);
	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
}

void ADMMGPUConstCons3::updateLocalProbGPU(float epsL, int nIterL) {
	int numBlocks = _nAgent;
	switch (_blockSize) {
	case 512:
		updateTradePGPUSharedResidualCons<512> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case 256:
		updateTradePGPUSharedResidualCons<256> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case 128:
		updateTradePGPUSharedResidualCons<128> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case 64:
		updateTradePGPUSharedResidualCons< 64> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case 32:
		updateTradePGPUSharedResidualCons< 32> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case 16:
		updateTradePGPUSharedResidualCons< 16> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case  8:
		updateTradePGPUSharedResidualCons<  8> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case  4:
		updateTradePGPUSharedResidualCons<  4> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case  2:
		updateTradePGPUSharedResidualCons<  2> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	case  1:
		updateTradePGPUSharedResidualCons<  1> << <numBlocks, _blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap3._matrixGPU, Ap123._matrixGPU, Bp3._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, epsL, nIterL);
		break;
	}
	//cudaStreamSynchronize(streamCalculation);
}



void ADMMGPUConstCons3::updateGlobalProbGPU()
{
	//Rem : tout calcul qui est de taille N ou M peut �tre fait par les agents
		// Si le calcul est de taile L, soit c'est calcul� par un/des superviseurs, soit tous les agents le calcul (un peu absurde)

#ifdef INSTRUMENTATION
// FB 3a
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION


	updatePnGPU << <_numBlocksN, _blockSize >> > (Pn._matrixGPU, Tmoy._matrixGPU, nVoisin._matrixGPU, _nAgent);
	
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	timePerBlock.increment(0, 5, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 3b
	cudaDeviceSynchronize();
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION	

	// Resolution de l'OPF
	if (_nLine) {
		//std::cout << " Pn :" << std::endl;
		//Pn.display(true);
		//std::cout << " etaP :" << std::endl;
		//etaP.display(true);

		solveOPF();
		//std::cout << " Pso :" << std::endl;
		//Pso.display(true);
	}
	else {
		Pso = Pn;
		/*tempN1.set(&etaP);
		tempN1.divide(_rho1);
		Pso.add(&Pn);
		Pso.divide(2);
		Pso.subtract(&tempN1);*/
	}
	
	
	

#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 3c
	cudaDeviceSynchronize();
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	// update Bp3
	updateEtaPBp3 << <_numBlocksN, _blockSize >> > (Bp3._matrixGPU, etaP._matrixGPU, nVoisin._matrixGPU, Pso._matrixGPU, Pn._matrixGPU, _rho1, _nAgent);
	updateLAMBDABt1GPU << <_numBlocksM, _blockSize >> > (Bt1._matrixGPU, LAMBDALin._matrixGPU, tradeLin._matrixGPU, _rhog, CoresLinTrans._matrixGPU, _nTrade);
	

#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
}



float ADMMGPUConstCons3::updateResEndo(int iter)
{

	float resS = Tlocal.max2(&tradeLin);

	updateDiffGPU <<<_numBlocksM, _blockSize >> > (tempNN._matrixGPU, Tlocal._matrixGPU, CoresLinTrans._matrixGPU, _nAgent);
	float resR = tempNN.max2();

	float resXf = _ratioEps * Pso.max2(&Pn);
	resF.set(0, iter, resR);
	resF.set(1, iter, resS);
	resF.set(2, iter, resXf);
	return MYMAX(MYMAX(resXf, resS), resR);

}



