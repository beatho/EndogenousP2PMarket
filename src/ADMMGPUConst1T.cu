#include "../head/ADMMGPUConst1T.cuh"




ADMMGPUConst1T::ADMMGPUConst1T() : MethodP2PGPU()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "Constructeur ADMMGPUConst1T" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu

}


ADMMGPUConst1T::ADMMGPUConst1T(float rho) : MethodP2PGPU()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "Constructeur ADMMGPUConst1T defaut" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu

}

ADMMGPUConst1T::~ADMMGPUConst1T()
{
}

void ADMMGPUConst1T::setParam(float rho)
{
	_rho = rho;
}

void ADMMGPUConst1T::setTau(float tau)
{
	if (tau < 1) {
		throw std::invalid_argument("tau must be greater than 1");
	}
	_tau = tau;
}

void ADMMGPUConst1T::init(const Simparam& sim, const StudyCase& cas)
{
	// intitilisation des matrixs et variables 

	
	//std::cout << "init " << std::endl;
	// reste sur CPU 
	isAC = false;
	initSize(cas);
	initSimParam(sim);


	initLinForm(sim, cas);


	//std::cout << "donnees sur GPU pour le grid" << std::endl;
	
	initDCEndoGrid(cas);
	

	//std::cout << "autres donn�e sur GPU" << std::endl;
	initCaseParam(sim, cas);

	initDCEndoMarket();

	updateGlobalProbGPU();
	//std::cout << "fin init" << std::endl;
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
}


void ADMMGPUConst1T::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
	float epsG = sim.getEpsG();
	float epsL = sim.getEpsL();
	
	const int stepL = sim.getStepL();
	const int stepG = sim.getStepG();
	const int iterG = sim.getIterG();
	const int iterL = sim.getIterL();
	

	float resG = 2 * epsG;
	float resL = 2 * epsL;
	int iterGlobal = 0;
	int iterLocal = 0;

	
	while ((iterGlobal < iterG) && (resG>epsG)) {
		resL = 2 * epsL;
		iterLocal = 0;
		while (iterLocal< iterL && resL>epsL) {
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

			if (!(iterLocal % stepL)) {
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
		if (iterLocal == iterL) {
			std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resG << std::endl;
		}
#ifdef INSTRUMENTATION
		occurencePerBlock.increment(0, 1, iterLocal);
		occurencePerBlock.increment(0, 4, iterLocal / stepL);
#endif // INSTRUMENTATION
		

		Tlocal.swap(&Tlocal_pre); // on �viter d'echanger lorsque qu'il ne faut pas
		tradeLin.swap(&Tlocal); // echange juste les pointeurs
		updateGlobalProbGPU();
		
		if (!(iterGlobal % stepG)) {
#ifdef INSTRUMENTATION
			cudaDeviceSynchronize();
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			
			resG = updateResEndo(iterGlobal / stepG);
#ifdef INSTRUMENTATION
			cudaDeviceSynchronize();
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
			
			//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG)
			//<< " " << resF.get(1, iterGlobal / stepG) << " " << resF.get(2, iterGlobal / stepG) << std::endl;
			/*if (iterLocal == iterL) {
				float d1 = tradeLin.max2(&Tlocal_pre);
				float d2 = P.max2(&Tmoy);
				std::cout << " local : resS " << d1 << " resR " << d2 << std::endl;
			}*/
		}
		
		iterGlobal++;
	}
	//std::cout << "fin simu temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
#ifdef INSTRUMENTATION
	occurencePerBlock.increment(0, 5, iterGlobal);
	occurencePerBlock.increment(0, 6, iterGlobal);
	occurencePerBlock.increment(0, 7, iterGlobal);
	occurencePerBlock.increment(0, 8, iterGlobal / stepG);
	cudaDeviceSynchronize();
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	
	Kappa1.projectNeg(); //delta1
	Kappa2.projectNeg(); // delta2


	float fc = calcFc();
	//float d1 = tradeLin.max2(&Tlocal_pre);
	//float d2 = P.max2(&Tmoy);
	//std::cout << " local : resS " << d1 << " resR " << d2 << std::endl;
	//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, (iterGlobal-1) / stepG) << " " << resF.get(1, (iterGlobal - 1) / stepG) << " " << resF.get(2, (iterGlobal - 1) / stepG) << std::endl;
	//std::cout << "Desequilibre " << Pn.sum() << std::endl;
	
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
	for (int idAgent = 0; idAgent < _nAgent; idAgent++) {
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
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 9, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 9, 1);
	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
	

	tall = clock() - tall;
	result->setTime((float)tall / CLOCKS_PER_SEC);

	// � changer � la main pour l'instant
	//std::string name = cas.getName();
	//std::string fileName = "AllResidual"+ name +".csv";
	//std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	//resF.saveCSV(fileName, mode);
	//resX.saveCSV(fileName, mode);


}

void ADMMGPUConst1T::updateLocalProbGPU( MatrixGPU* Tlocal, MatrixGPU* P) {
	int numBlocks = _nAgent;
	
	switch (_blockSize) {
	case 512:
		updateTradePGPU<512> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, CoresLinAgent._matrixGPU, _nAgent);
		break;
	case 256:
		updateTradePGPU<256> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, CoresLinAgent._matrixGPU, _nAgent);
		break;
	case 128:
		updateTradePGPU<128> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, CoresLinAgent._matrixGPU, _nAgent);
		break;
	case 64:
		updateTradePGPU< 64> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, CoresLinAgent._matrixGPU, _nAgent);
		break;
	case 32:
		updateTradePGPU< 32> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, CoresLinAgent._matrixGPU, _nAgent);
		break;
	case 16:
		updateTradePGPU< 16> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, CoresLinAgent._matrixGPU, _nAgent);
		break;
	case  8:
		updateTradePGPU<  8> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, CoresLinAgent._matrixGPU, _nAgent);
		break;
	case  4:
		updateTradePGPU<  4> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, CoresLinAgent._matrixGPU, _nAgent);
		break;
	case  2:
		updateTradePGPU<  2> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, CoresLinAgent._matrixGPU, _nAgent);
		break;
	case  1:
		updateTradePGPU<  1> << <numBlocks, _blockSize >> > (Tlocal->_matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P->_matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, _at1, _at2, Bt1._matrixGPU, Ct._matrixGPU,
			matLb._matrixGPU, matUb._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, CoresAgentLin._matrixGPU, CoresLinAgent._matrixGPU, _nAgent);
		break;
	}
}



void ADMMGPUConst1T::updateGlobalProbGPU()
{
	//Rem : tout calcul qui est de taille N ou M peut �tre fait par les agents
	// Si le calcul est de taile L, soit c'est calcul� par un/des superviseurs, soit tous les agents le calcul (un peu absurde)
	
	// FB 3a
#ifdef INSTRUMENTATION
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
	//updateCpOld << <_numBlocksN, _blockSize >> > (Cp._matrixGPU, Cp1._matrixGPU, Cp2._matrixGPU, tempN1._matrixGPU, nVoisin._matrixGPU, _rho1, _nAgent);

#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
}



