#include "../head/ADMMGPUConst3.cuh"

ADMMGPUConst3::ADMMGPUConst3() : MethodP2PGPU()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "Constructeur ADMMGPUConst3" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu

}


ADMMGPUConst3::ADMMGPUConst3(float rho) : MethodP2PGPU()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "Constructeur ADMMGPUConst3 defaut" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu

}

ADMMGPUConst3::~ADMMGPUConst3()
{
}

void ADMMGPUConst3::setParam(float rho)
{
	_rho = rho;
}

void ADMMGPUConst3::setTau(float tau)
{
	if (tau < 1) {
		throw std::invalid_argument("tau must be greater than 1");
	}
	_tau = tau;
}

void ADMMGPUConst3::init(const Simparam& sim, const StudyCase& cas)
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
	initCaseParam(sim, cas);
	//std::cout << "mise sous forme lineaire" << std::endl;
	initLinForm(cas);

	//std::cout << "donnees sur GPU pour le grid" << std::endl;
	initDCEndoGrid(cas);
	
	//std::cout << "autres donnee sur GPU" << std::endl;
	

	initDCEndoMarket();

	
	updateGlobalProbGPU();	
}


void ADMMGPUConst3::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
	/*_mu = _muInit;
	_mu1 = _muInit;
	_rhog = sim.getRho();
	float oldRho1 = _rho1;
	_rho1 = sim.getRho1();
	float tau = _rho1 / oldRho1;
	Ap2b.multiply(_tau);
	Ap2.add(&Ap2a, &Ap2b);
	Ap12.add(&Ap1, &Ap2);*/
	_rhog = sim.getRho();
	_at1 = _rhog; // represente en fait 2*a
	
	
	float resG = 2 * _epsG;
	float resL = 2 * _epsL;
	_iterGlobal = 0;
	int iterLocal = 0;

	
	//std::cout << iterG << " " << iterL << " " << epsL << " " << epsG << std::endl;
	while ((_iterGlobal < _iterG) && (resG > _epsG)) {
		resL = 2 * _epsL;
		iterLocal = 0;
		while (iterLocal < _iterL && resL > _epsL) {
#ifdef INSTRUMENTATION
			cudaDeviceSynchronize();
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			
			updateLocalProbGPU(&Tlocal, &P);
#ifdef INSTRUMENTATION
			cudaDeviceSynchronize();
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 1, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION

			if (!(iterLocal % _stepL)) {
#ifdef INSTRUMENTATION
				cudaDeviceSynchronize();
				t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
				resL = calcRes();
#ifdef INSTRUMENTATION
				cudaDeviceSynchronize();
				t2 = std::chrono::high_resolution_clock::now();
				timePerBlock.increment(0, 4, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
				
			}
			//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resG << std::endl;
			Tlocal.swap(&Tlocal_pre); 
			iterLocal++;
		}
		if (iterLocal == _iterL) {
			std::cout << _iterGlobal << " " << iterLocal << " " << resL << " " << resG << std::endl;
		}
#ifdef INSTRUMENTATION
		occurencePerBlock.increment(0, 1, iterLocal);
		occurencePerBlock.increment(0, 4, iterLocal / stepL);
#endif // INSTRUMENTATION


		Tlocal.swap(&Tlocal_pre); // on �viter d'echanger lorsque qu'il ne faut pas
		tradeLin.swap(&Tlocal); // echange juste les pointeurs	
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
	occurencePerBlock.increment(0, 5, iterGlobal);
	occurencePerBlock.increment(0, 6, iterGlobal);
	occurencePerBlock.increment(0, 7, iterGlobal);
	occurencePerBlock.increment(0, 8, iterGlobal / stepG);

	//std::cout << "fin simu temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	cudaDeviceSynchronize();
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION


	setResult(result, cas.isAC());


#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 9, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 9, 1);
	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
	

	

}

void ADMMGPUConst3::updateLocalProbGPU( MatrixGPU* Tlocal, MatrixGPU* P) {
	int numBlocks = _nAgent;
	switch (_blockSize) {
	case 512:
		updateTradePGPUShared<512> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU);
		break;
	case 256:
		updateTradePGPUShared<256> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU);
		break;
	case 128:
		updateTradePGPUShared<128> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU);
		break;
	case 64:
		updateTradePGPUShared< 64> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU);
		break;
	case 32:
		updateTradePGPUShared< 32> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU);
		break;
	case 16:
		updateTradePGPUShared< 16> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU);
		break;
	case  8:
		updateTradePGPUShared<  8> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU);
		break;
	case  4:
		updateTradePGPUShared<  4> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU);
		break;
	case  2:
		updateTradePGPUShared<  2> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU);
		break;
	case  1:
		updateTradePGPUShared<  1> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU);
		break;
	}
}



void ADMMGPUConst3::updateGlobalProbGPU()
{
	//Rem : tout calcul qui est de taille N ou M peut �tre fait par les agents
		// Si le calcul est de taile L, soit c'est calcul� par un/des superviseurs, soit tous les agents le calcul (un peu absurde)


#ifdef INSTRUMENTATION
		// FB 3a
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

#endif // INSTRUMENTATION

	updatePnGPU << <_numBlocksN, _blockSize >> > (Pn._matrixGPU, Tmoy._matrixGPU, nVoisin._matrixGPU, _nAgent);
	updateAlphaTrans << < _numBlocksNL, _blockSize >> > (alpha._matrixGPU, GTrans._matrixGPU, Pn._matrixGPU, _nLine, _nAgent);
	updateQpartTrans << < _nLine, _blockSize, _nAgent * sizeof(float) >> > (Qpart._matrixGPU, alpha._matrixGPU, _nAgent, _nLine);
	updateQtotTrans << <_numBlocksL, _blockSize >> > (Qtot._matrixGPU, Qpart._matrixGPU, alpha._matrixGPU, _nLine);
	
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	timePerBlock.increment(0, 5, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 3b
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
		//updateCp2aTrans<512> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, _nLine, _nAgent);
		//updateCp2bTrans<512> << <numBlocks, _blockSize >> > (tempN1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, _nLine, _nAgent);

		break;
	case 256:
		updateCp2GPUTrans<256> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		//updateCp2aTrans<256> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, _nLine, _nAgent);
		//updateCp2bTrans<256> << <numBlocks, _blockSize >> > (tempN1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, _nLine, _nAgent);

		break;
	case 128:
		updateCp2GPUTrans<128> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		//updateCp2aTrans<128> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, _nLine, _nAgent);
		//updateCp2bTrans<128> << <numBlocks, _blockSize >> > (tempN1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, _nLine, _nAgent);

		break;
	case 64:
		updateCp2GPUTrans<64> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		//updateCp2aTrans<64> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, _nLine, _nAgent);
		//updateCp2bTrans<64> << <numBlocks, _blockSize >> > (tempN1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, _nLine, _nAgent);

		break;
	case 32:
		updateCp2GPUTrans<32> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		//updateCp2aTrans<32> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, _nLine, _nAgent);
		//updateCp2bTrans<32> << <numBlocks, _blockSize >> > (tempN1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, _nLine, _nAgent);

		break;
	case 16:
		updateCp2GPUTrans<16> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		//updateCp2aTrans<16> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, _nLine, _nAgent);
		//updateCp2bTrans<16> << <numBlocks, _blockSize >> > (tempN1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, _nLine, _nAgent);

		break;
	case  8:
		updateCp2GPUTrans<8> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		//updateCp2aTrans<8> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, _nLine, _nAgent);
		//updateCp2bTrans<8> << <numBlocks, _blockSize >> > (tempN1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, _nLine, _nAgent);

		break;
	case  4:
		updateCp2GPUTrans<4> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		//updateCp2aTrans<4> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, _nLine, _nAgent);
		//updateCp2bTrans<4> << <numBlocks, _blockSize >> > (tempN1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, _nLine, _nAgent);

		break;
	case  2:
		updateCp2GPUTrans<2> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		//updateCp2aTrans<2> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, _nLine, _nAgent);
		//updateCp2bTrans<2> << <numBlocks, _blockSize >> > (tempN1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, _nLine, _nAgent);

		break;
	case  1:
		updateCp2GPUTrans<1> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, nVoisin._matrixGPU, _rho1, _nLine, _nAgent);
		//updateCp2aTrans<1> << <numBlocks, _blockSize >> > (Cp2._matrixGPU, tempL1._matrixGPU, GTrans._matrixGPU, _nLine, _nAgent);
		//updateCp2bTrans<1> << <numBlocks, _blockSize >> > (tempN1._matrixGPU, GTrans._matrixGPU, Qpart._matrixGPU, _nLine, _nAgent);
		break;
	}

#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 3c
	cudaDeviceSynchronize();
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION



	
	updateLAMBDABt1GPU << <_numBlocksM, _blockSize >> > (Bt1._matrixGPU, LAMBDALin._matrixGPU, tradeLin._matrixGPU, _rhog, CoresLinTrans._matrixGPU, _nTrade);
	updateCp << <_numBlocksN, _blockSize >> > (Cp._matrixGPU, Cp1._matrixGPU, Cp2._matrixGPU, _nAgent);

#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
	

}



