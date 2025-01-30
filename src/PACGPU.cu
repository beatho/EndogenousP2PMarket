	#include "../head/PACGPU.cuh"
#define MAX(X, Y) X * (X >= Y) + Y * (Y > X)
// betaLin, CoresAgentLinBig, indiceNu

PACGPU::PACGPU() : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " PACGPU Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
}


PACGPU::PACGPU(float rho) : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default PACGPU Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
}

PACGPU::~PACGPU()
{
	/*DELETEA(tempM1);
	DELETEA(tempM);
	DELETEA(X);
	DELETEA(Xpre);
	DELETEA(Xhat);
	DELETEA(Mu);
	DELETEA(Muhat);
	DELETEA(Nu);
	DELETEA(Nuhat);
	DELETEA(Hinv);
	DELETEA(H);
	DELETEA(matLb);
	DELETEA(Q);
	DELETEA(Qinit);
	DELETEA(matUb);*/
}
void PACGPU::setParam(float rho)
{
	_rho = rho;
}

void PACGPU::setGamma(float gamma)
{
	_gamma = gamma;
	if (_gammahat > _gamma) {
		std::cout << "Warning : gammahat should be smaller than gamma";
	}
}

void PACGPU::setGammahat(float gammahat)
{
	_gammahat = gammahat;
	if (_gammahat > _gamma) {
		std::cout << "Warning : gammahat should be smaller than gamma";
	}
}

void PACGPU::setInitCoef(float alpha, float phi, float theta)
{
	if (alpha <= 0 || phi < 0 || theta < 0 || alpha >= 1 || phi >= 1 || theta >= 1) {
		throw std::invalid_argument("coefficient must be positive and <1");
	}
	_alpha = alpha;
	_phi = phi;
	_theta = theta;
	
}

void PACGPU::setBestRhoGamma(float lambdaMax, float lambdaMin, const StudyCase& cas)
{
	throw std::runtime_error("setBestRhoGamma : WIP not implemented");
	/*MatrixGPU a = cas.geta();
	float alpha = 1 * a.min2();
	float L = 1 * a.max2();
	_gamma = (2 * alpha * L) / (2 * lambdaMax + lambdaMin);
	_gammahat = _gamma;
	if (augmente && _alpha>1) {
		_rho = 1 / (sqrt(_gamma * _alpha * lambdaMax));
	}
	else {
		_rho = 1 / (sqrt(_gamma * lambdaMax));
	}*/
	
	std::cout << "best gamma " << _gamma << " and rho " << _rho << std::endl;

}

void PACGPU::setBestRhoGammaHeuristic(const StudyCase& cas)
{
	int N = cas.getNagent();
	int Y = cas.getNCons();
	float lambdaMax = 1;
	float lambdaMin = 1;
	if (cas.isAC()) {
		lambdaMax = 0.9995 * N + 2.0455;
		lambdaMin = (1.9927 - 3.3289 / N) / N;
	}
	else {
		lambdaMax = abs(Y - 0.5 * N) + 0.4961 * N + 3.3561;
		if (N > 10) {
			lambdaMin = (-0.00675 * N + 0.9675) / (N * N) * (Y - 0.5 * N) * (Y - 0.5 * N) + (log(831 * N - 5755)) / (3 * N);
		}
		else {
			lambdaMin = 0.35 + 0.4676 / N;
		}
	}
	if (lambdaMin < 0) {
		lambdaMin = 0.000001;
	}

	MatrixCPU a = cas.geta();
	float alpha = 1 * a.min2Nnull(0.0001);
	float L = 1 * a.max2();
	//std::cout << "alpha " << alpha << " and L " << L << std::endl;
	_gamma = (2 * alpha * L) / (2 * lambdaMax + lambdaMin);
	_gammahat = _gamma;
	if (augmente && _alpha > 1) {
		_rho = 1 / (sqrt(_gamma * _alpha * lambdaMax));
	}
	else {
		_rho = 1 / (sqrt(_gamma * lambdaMax));
	}

	//std::cout << "best gamma " << _gamma << " and rho " << _rho << std::endl;
}

void PACGPU::updateCoef()
{
	throw std::runtime_error("updateCoef : WIP not implemented");
	/*if (augmente) {

		_rho = MAX(0.99 * _rho, 0.1);
		_rhoInv = 1 / _rho;
		for (int i = 0; i < _nAgent; i++) {
			H[i].set(0, 0, Cost1.get(i, 0) + _rhoInv);
			if (augmente) {
				H[i].increment(0, 0, _rho * _gamma);
			}
			int M = nVoisin.get(i, 0);
			for (int m = 0; m < M; m++) {

				H[i].set(m + 1, m + 1, _rhoInv); // diag tnm
				H[i].set(M + m + 1, M + m + 1, _rhoInv); // diag anm
				if (augmente) {
					H[i].set(m + 1, 0, -_rho * _gamma); // first column pn <-> tnm
					H[i].set(0, m + 1, -_rho * _gamma); // fisrt row    pn <-> tnm

					H[i].increment(m + 1, m + 1, 2 * _rho * _gamma); // diag tnm
					H[i].increment(M + m + 1, M + m + 1, 2 * _rho * _gamma); // diag anm
				}
			}
			for (int k = 1; k < M + 1; k++) {
				for (int j = 1; j < M + 1; j++) {
					if (k != j) {
						H[i].set(k, j, _rho * _gamma);
					}
				}
			}
			for (int k = 1; k < M + 1; k++) {
				for (int j = M + 1; j < 2 * M + 1; j++) {
					H[i].set(k, j, _rho * _gamma);
					H[i].set(j, k, _rho * _gamma);
				}
			}
			Hinv[i].invertEigen(&H[i]);
		}
	}
	*/
	
	/**/
	/*_alpha = 0.9 * _alpha;
	_theta = 0.9 * _theta;
	_phi = 0.9 * _phi;*/
}



void PACGPU::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
	_rhog = sim.getRho();
	
	int iterG = sim.getIterG();
	int stepG = sim.getStepG();
	
	float epsG = sim.getEpsG()/2;
	
	
	float resG = 2 * epsG;
	int iterGlobal = 0;
	while ((iterGlobal < iterG) && (resG>epsG)) {

		updateLocalProb();
		
#ifdef INSTRUMENTATION
		occurencePerBlock.increment(0, 1, 1);
		//occurencePerBlock.increment(0, 2, 0);
		//occurencePerBlock.increment(0, 3, 0);

		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

		updateGlobalProb();

#ifdef INSTRUMENTATION
		cudaDeviceSynchronize();
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 5, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		// FB 4
		if (!((iterGlobal -1) % stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateRes(iterGlobal / stepG);
			//resG = 1;
#ifdef INSTRUMENTATION
			cudaDeviceSynchronize();
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 6, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		}
		//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;

		iterGlobal++;
	}
	//std::cout << iterGlobal  << " " << resF.get(0, (iterGlobal - 1) / stepG) << " " << resF.get(1, (iterGlobal - 1) / stepG) << " " << resF.get(2, (iterGlobal - 1) / stepG) << std::endl;


	
	
#ifdef INSTRUMENTATION	
	occurencePerBlock.increment(0, 5, iterGlobal);
	occurencePerBlock.increment(0, 6, iterGlobal / stepG);
	

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	MatrixCPU tradeLinCPU;
	X.toMatCPU(tradeLinCPU);
	MatrixCPU CoresAgentLinCPU;
	CoresAgentLin.toMatCPU(CoresAgentLinCPU);
	MatrixCPU CoresLinVoisinCPU;
	CoresLinVoisin.toMatCPU(CoresLinVoisinCPU);

	//MatrixCPU MUCPU;
	//MU.toMatCPU(MUCPU);
	
	
	for (int idAgent = 0; idAgent < _nAgentTrue; idAgent++) {
		int begin = CoresAgentLinCPU.get(idAgent, 0);
		int Nvoisinmax = nVoisinCPU.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			int idVoisin = CoresLinVoisinCPU.get(begin + voisin + 1,0);
			trade.set(idAgent, idVoisin, tradeLinCPU.get(begin + voisin + 1, 0));
		}
		Pn.set(idAgent, 0, tradeLinCPU.get(begin, 0));
	}
	int indice = 0;
	for (int idAgent = _nAgentTrue; idAgent < _nAgent; idAgent++) {
		indice = 0;
		int begin = CoresAgentLinCPU.get(idAgent, 0);
		for (int idVoisin = 0; idVoisin < _nAgentTrue; idVoisin++) {
			if (idVoisin != (idAgent - _nAgentTrue)) {
				trade.set(idAgent, idVoisin, tradeLinCPU.get(begin + indice + 1, 0));
				indice = indice + 1;
			}
			Pn.set(idAgent, 0, tradeLinCPU.get(begin, 0));
		}
	}
	//std::cout << "Pn :" << std::endl;
	//Pn.display();
	
	
	calculFcPAC << <_nAgent, _blockSize >> > (tempN1._matrixGPU, tempNN._matrixGPU, Cost1._matrixGPU, Cost2._matrixGPU, Qinit._matrixGPU, X._matrixGPU, CoresAgentLin._matrixGPU, nVoisin._matrixGPU);
	
	float fc = tempN1.sum();

	fc += tempNN.sum();
	
	 //fc = calcFc(&Cost1, &Cost2, &trade, &Pn, &BETA, &tempN1, &tempNN);
	// FB 5
	
	//std::cout << "set end" << std::endl;


	result->setResF(&resF);
	result->setTrade(&trade); 
	result->setIter(iterGlobal);
	result->setPn(&Pn);
	result->setFc(fc);
	
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

void PACGPU::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	

	MatrixGPU Pmin(cas.getPmin(), 1);
	MatrixGPU Pmax(cas.getPmax(), 1);
	MatrixGPU Cost2(cas.getb(), 1);
	MatrixGPU Lb(cas.getLb(), 1);
	MatrixGPU Ub(cas.getUb(), 1);

	updateP0PAC << <_nAgent, _blockSize >> > (matLb._matrixGPU, matUb._matrixGPU, Q._matrixGPU, Qinit._matrixGPU, Pmin._matrixGPU, Pmax._matrixGPU, Cost2._matrixGPU, Lb._matrixGPU, Ub._matrixGPU, CoresAgentLin._matrixGPU, nVoisin._matrixGPU);
	
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION

}

void PACGPU::init(const Simparam& sim, const StudyCase& cas)
{
	// initialisation des matrixs et variables 
	
	clock_t t = clock();
	//std::cout << "init " << std::endl;
	if (_rho == 0) {
		_rho = sim.getRho();
	}

	_rhog = _rho;
	const int iterG = sim.getIterG();
	const int stepG = sim.getStepG();
	float epsG = sim.getEpsG();
	isAC = cas.isAC();
	
	_nAgentTrue = sim.getNAgent();
	_nAgent = _nAgentTrue + isAC * _nAgentTrue;
	_numBlocksN = ceil((_nAgent + _blockSize - 1) / _blockSize);
	nVoisin = MatrixGPU(cas.getNvoi(), 1);
	nVoisin.preallocateReduction();
	nVoisinCPU = cas.getNvoi();


	_nTrade = nVoisin.sum();
	_nTradeP = 0;
	if (isAC) {
		for (int n = 0; n < _nAgentTrue; n++) {
			_nTradeP += nVoisinCPU.get(n, 0);
		}
		_nTradeQ = _nTrade - _nTradeP;
		if (_nTradeQ != (_nAgentTrue * (_nAgentTrue - 1))) {
			std::cout << "err PACGPU : " << _nAgent << " " << _nAgentTrue << " " << _nTrade << " " << _nTradeP << " " << _nTradeQ << std::endl;

			throw std::invalid_argument("Agent must be fully connected for the Q echanges, WIP");
		}
	}
	else {
		_nTradeP = _nTrade;
	}
	//_sizePACGPU = _nAgent + 2 * _nTrade;


	_rhoInv = 1 / _rhog;
	resF = MatrixCPU(3, (iterG / stepG) + 1);

	

	_sizePACX = _nAgent + 2 * _nTrade; //_nAgent * (1 + 2*Mn)
	_sizePACMu = _nAgent + _nTrade; // _nAgent * (1 + Mn)
	_sizePACNu = _nTrade; // _nAgent * Mn
	_sizeHinv = 1 + 2 * nVoisin.max2();


	if (CoresMatLin.getPos()) { // si sur GPU il faut refaire les transferts
		CoresIndiceNu.transferCPU();
		CoresAgentLin.transferCPU();
		CoresAgentLinBig.transferCPU();
		_sizeQ.transferCPU();
		CoresLinAgent.transferCPU();
		CoresLinVoisin.transferCPU();
		CoresLinTrans.transferCPU();
		CoresMatLin.transferCPU();
		CoresLinTransLocal.transferCPU();

		X.transferCPU();
		matLb.transferCPU();
		matUb.transferCPU();
		Hinv.transferCPU();
		Qinit.transferCPU();
		Cost1.transferCPU();
		Cost2.transferCPU();
	}
	BETA = MatrixGPU(cas.getBeta());
	MatrixCPU Ub(cas.getUb());
	MatrixCPU Lb(cas.getLb());
	MatrixCPU Pmin(cas.getPmin());
	MatrixCPU Pmax(cas.getPmax());
	Cost1 = cas.geta();
	Cost2 = cas.getb();

	trade = sim.getTrade();
	Pn = sim.getPn();

	CoresMatLin = MatrixGPU(_nAgent, _nAgentTrue, -1);
	CoresLinAgent = MatrixGPU(_sizePACX, 1, -1);
	CoresAgentLin = MatrixGPU(_nAgent, 1, -1);
	CoresAgentLinBig = MatrixGPU(_sizePACX, 1, -1);
	CoresLinVoisin = MatrixGPU(_sizePACX, 1, -1);
	CoresLinTrans = MatrixGPU(_sizePACX, 1, -1);
	CoresLinTransLocal = MatrixGPU(_sizePACX, 1, -1);
	CoresIndiceNu = MatrixGPU(_sizePACX, 1, -1);
	_sizeQ = MatrixGPU(_sizePACX, 1, -1);
	
	// Who is the peer ? 
	int indice = 0;
	int debutNu = 0;
	
	//std::cout << " P " << std::endl;
	for (int idAgent = 0; idAgent < _nAgentTrue; idAgent++) { // P
		MatrixGPU omega(cas.getVoisin(idAgent));
		int Nvoisinmax = nVoisinCPU.get(idAgent, 0);
		int begin = indice;
		
		
		CoresIndiceNu.set(indice, 0, debutNu); // le debut de Mu[agent]
		CoresAgentLin.set(idAgent, 0, indice);

		// Pn
		CoresAgentLinBig.set(indice, 0, begin);
		_sizeQ.set(indice, 0, 1 + 2 * Nvoisinmax);
		indice += 1; 
		debutNu += Nvoisinmax;
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			int idVoisin = omega.get(voisin, 0);
			CoresAgentLinBig.set(indice, 0, begin);
			CoresLinAgent.set(indice, 0, idAgent);
			CoresLinVoisin.set(indice, 0, idVoisin); // 
			CoresAgentLinBig.set(indice + Nvoisinmax, 0, begin);
			CoresLinAgent.set(indice + Nvoisinmax, 0, idVoisin);
			CoresLinVoisin.set(indice + Nvoisinmax, 0, idAgent); // 
			_sizeQ.set(indice, 0, 1 + 2 * Nvoisinmax);
			_sizeQ.set(indice + Nvoisinmax, 0, 1 + 2 * Nvoisinmax);
			CoresMatLin.set(idAgent, idVoisin, indice);
			indice = indice + 1;
		}
		indice += Nvoisinmax;
	}
	//std::cout << " Q " << std::endl;
	for (int idAgent = _nAgentTrue; idAgent < _nAgent; idAgent++) { // Q
		int begin = indice;
		int Nvoisinmax = (_nAgentTrue - 1);
		
		CoresIndiceNu.set(indice, 0, debutNu); // le debut de Mu[agent]
		CoresAgentLin.set(idAgent, 0, indice);

		CoresAgentLinBig.set(indice, 0, begin);
		_sizeQ.set(indice, 0, 1 + 2 * Nvoisinmax);
		indice += 1;
		debutNu += Nvoisinmax;
		for (int idVoisin = 0; idVoisin < _nAgentTrue; idVoisin++) {
			if (idVoisin != (idAgent - _nAgentTrue)) {
				CoresAgentLinBig.set(indice, 0, begin);
				CoresAgentLinBig.set(indice + Nvoisinmax, 0, begin);
				CoresLinAgent.set(indice, 0, idAgent);
				CoresLinVoisin.set(indice, 0, idVoisin + _nAgentTrue);
				CoresLinAgent.set(indice + Nvoisinmax, 0, idVoisin + _nAgentTrue);
				CoresLinVoisin.set(indice + Nvoisinmax, 0, idAgent);
				CoresMatLin.set(idAgent, idVoisin, indice);
				_sizeQ.set(indice, 0, 1 + 2 * Nvoisinmax);
				_sizeQ.set(indice + Nvoisinmax, 0, 1 + 2 * Nvoisinmax);
				indice = indice + 1;
			}
		}
		indice += Nvoisinmax;
	}
	
	/*CoresLinAgent.display();
	CoresLinVoisin.display();*/

	for (int idAgent = 0; idAgent < _nAgent; idAgent++) {
		int debut = CoresAgentLin.get(idAgent, 0, 0);
		int Nvoisinmax = nVoisinCPU.get(idAgent, 0);

		for (int m = 0; m < Nvoisinmax; m++) {
			int lin = debut + m + 1;
			int voisin = CoresLinVoisin.get(lin, 0);
			int k = 0;
			if (idAgent < _nAgentTrue) // P
				k = CoresMatLin.get(voisin, idAgent);
			else {
				k = CoresMatLin.get(voisin, idAgent - _nAgentTrue);
			}

			//std::cout << " trade num " << lin << " entre " << idAgent << " et " << voisin;
			//std::cout << " le trade symetrique le trade num " << k << std::endl;
			CoresLinTrans.set(lin, 0, k);
		}
	}

	for (int i = 0; i < _nAgent; i++) {
		int M = nVoisinCPU.get(i, 0);
		for (int m = 0; m < M; m++) {
			int lin = CoresAgentLin.get(i, 0) + m + 1; // indice global de tim = tip (m est le num�ro du voisin, p est le num�ro de l'agent)
			int p = CoresLinVoisin.get(lin, 0); // valeur de p
			int lin2 = CoresLinTrans.get(lin, 0); // indice global de tpi
			int linLoc = lin2 - CoresAgentLin.get(p, 0) - 1; //indice local de tpi
			CoresLinTransLocal.set(lin, 0, linLoc);

			int debutNuTrans = CoresIndiceNu.get(CoresAgentLin.get(p, 0), 0);
			CoresIndiceNu.set(lin, 0, debutNuTrans + linLoc);

			//std::cout << lin << ": trade entre " << i << " et " << p << " num local " << m ;
			//std::cout << " le symetrique est " << lin2 <<  " ou " << linLoc << std::endl;
		}
	}
	//CoresAgentLin.display();
	//CoresLinTrans.display();
	//CoresIndiceNu.display();
	

	//std::cout << " creation" << std::endl;
	
	X = MatrixGPU(_sizePACX, 1);
	Xpre = MatrixGPU(_sizePACX, 1, 0, 1);
	Xhat = MatrixGPU(_sizePACX, 1, 0, 1);
	matLb = MatrixGPU(_sizePACX, 1);
	matUb = MatrixGPU(_sizePACX, 1);

	Mu = MatrixGPU(_sizePACMu, 1, 0, 1);
	Muhat = MatrixGPU(_sizePACMu, 1, 0, 1);
	Nu = MatrixGPU(_sizePACNu, 1, 0, 1);
	Nuhat = MatrixGPU(_sizePACNu, 1, 0, 1);

	tempN1 = MatrixGPU(_nAgent, 1, 0, 1);
	tempNN = MatrixGPU(_nAgent, _nAgentTrue, 0, 1);
	tempM1 = MatrixGPU(_sizePACMu, 1, 0, 1);
	tempM = MatrixGPU(_sizePACNu, 1, 0, 1);


	Hinv = MatrixGPU(_sizePACX, _sizeHinv);
	Q = MatrixGPU(_sizePACX, 1, 0, 1);
	Qinit = MatrixGPU(_sizePACX, 1);
	
	//std::cout << _sizePACX << " " << _sizeHinv << std::endl;

	//std::cout << "problem setup " << std::endl;
	for (int i = 0; i < _nAgent; i++) {
		//std::cout << "********* Agent " << i<< "**********"<< std::endl;
		// def
		int M = nVoisinCPU.get(i, 0);
		int begin = CoresAgentLin.get(i, 0);
		MatrixCPU Htemp(1 + 2 * M, 1 + 2 * M);
		H = MatrixCPU(1 + 2 * M, 1 + 2 * M);
		
		
		// init
		X.set(begin, 0, Pn.get(i, 0));
		matLb.set(begin, 0, Pmin.get(i, 0));
		matUb.set(begin, 0, Pmax.get(i, 0));
		H.set(0, 0, Cost1.get(i, 0) + _rhoInv);
		if (augmente) {
			H.increment(0, 0, _rho * _gamma);
		}
		
		Qinit.set(begin, 0, Cost2.get(i, 0));
		
		for (int m = 0; m < M; m++) {
			int indice = begin + m + 1;
			int voisin = CoresLinVoisin.get(indice, 0);
			_sizeQ.set(indice, 0, 1 + 2 * M);
			_sizeQ.set(indice + M, 0, 1 + 2 * M);
			//std::cout << "voisin num " << m << " is " << voisin << std::endl;

			X.set(indice, 0, trade.get(i, voisin % _nAgentTrue)); // tnm
			X.set(indice + M , 0, trade.get(voisin, i % _nAgentTrue)); //amn
			matLb.set(indice, 0, Lb.get(i, 0));
			matLb.set(indice + M, 0, -Ub.get(i, 0)); // est ce que cela g�ne la convergence ou est ce que cela l'aide ?
			matUb.set(indice, 0, Ub.get(i, 0));
			matUb.set(indice + M, 0, -Lb.get(i, 0));
			
			H.set(m + 1, m + 1, _rhoInv); // diag tnm
			H.set(M + m + 1, M + m + 1, _rhoInv); // diag anm
			if (augmente) {
				H.set(m + 1, 0, -_rho * _gamma); // first column pn <-> tnm
				H.set(0, m + 1, -_rho * _gamma); // fisrt row    pn <-> tnm

				H.increment(m + 1, m + 1, 2 * _rho * _gamma); // diag tnm
				H.increment(M + m + 1, M + m + 1, 2 * _rho * _gamma); // diag anm
			}
			if (i < _nAgentTrue) {
				Qinit.set(indice, 0, BETA.get(i, voisin));
			}
		}
		
		if (augmente) {
			for (int k = 1; k < M + 1; k++) {
				for (int j = 1; j < M + 1; j++) {
					if (k != j) {
						H.set(k, j, _rho * _gamma); // sum(sum(tnk*tnj)) -> ||Ga||
					}
				}
			}
			for (int k = 1; k < M + 1; k++) {
				int j = M + k;
				H.set(k, j, _rho * _gamma);
				H.set(j, k, _rho * _gamma);
			}
		}	
		//H.display(); 
		
		Htemp.invertGaussJordan(&H);
		
		
		
		Hinv.setBloc(begin, begin + 2 * M + 1, 0, 1 + 2 * M, &Htemp);
		
		
		
	}
	//Hinv.display();
	//std::cout << std::endl;
	//CoresAgentLinBig.display();
	//std::cout << " Transfert GPU " <<std::endl;
	
	//CoresAgentLin.display();
	//CoresLinTrans.display();



	CoresMatLin.transferGPU();
	CoresLinAgent.transferGPU();
	CoresAgentLin.transferGPU();
	CoresAgentLinBig.transferGPU();
	CoresLinVoisin.transferGPU();
	CoresLinTrans.transferGPU();
	CoresLinTransLocal.transferGPU();
	CoresIndiceNu.transferGPU();
	
	_sizeQ.transferGPU();
	X.transferGPU();

	X.preallocateReduction();
	Muhat.preallocateReduction();
	Nuhat.preallocateReduction();
	tempM1.preallocateReduction();
	tempM.preallocateReduction();
	tempN1.preallocateReduction();
	tempNN.preallocateReduction();

	matLb.transferGPU();
	matUb.transferGPU();
	Hinv.transferGPU();
	
	Qinit.transferGPU();

	Xpre.set(&X);

	Cost1.transferGPU();
	Cost2.transferGPU();
	


	//matLb.display(true);
	//matUb.display(true);

	//std::cout << "Global update" << std::endl;
	
	
	updateGlobalProb();
	


	//std::cout << "rho " << _rhog << " _alpha " << _alpha << " _phi " << _phi << std::endl;
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
}

void PACGPU::setBestParam(const StudyCase& cas)
{
	setBestRhoGammaHeuristic(cas);
}


void PACGPU::updateGlobalProb() {

	// communication of xhat
	updateNu();
	//Nu.display(true);
	//Nuhat.display(true);
	// communication of nuhat
	updateQ();
	if (augmente) {
		//updateCoef();
	}
	//Q.display(true);
	//std::cout << "******" << std::endl;
}

void PACGPU::updateLocalProb() {
	// FB 1a
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		
	
	Xpre.swap(&X);
	int numBlock = _sizePACX;
	switch (_blockSize) {
	case 512:
		updateLocalProblPAC<512> <<<numBlock, _blockSize >> > (X._matrixGPU, Q._matrixGPU, Hinv._matrixGPU, matLb._matrixGPU, matUb._matrixGPU, CoresAgentLinBig._matrixGPU, _sizeQ._matrixGPU, _sizeHinv);
		break;
	case 256:
		updateLocalProblPAC<256> << <numBlock, _blockSize >> > (X._matrixGPU, Q._matrixGPU, Hinv._matrixGPU, matLb._matrixGPU, matUb._matrixGPU, CoresAgentLinBig._matrixGPU, _sizeQ._matrixGPU, _sizeHinv);
		break;
	case 128:
		updateLocalProblPAC<128> << <numBlock, _blockSize >> > (X._matrixGPU, Q._matrixGPU, Hinv._matrixGPU, matLb._matrixGPU, matUb._matrixGPU, CoresAgentLinBig._matrixGPU, _sizeQ._matrixGPU, _sizeHinv);
		break;
	case 64:
		updateLocalProblPAC< 64> << <numBlock, _blockSize >> > (X._matrixGPU, Q._matrixGPU, Hinv._matrixGPU, matLb._matrixGPU, matUb._matrixGPU, CoresAgentLinBig._matrixGPU, _sizeQ._matrixGPU, _sizeHinv);
		break;
	case 32:
		updateLocalProblPAC< 32> << <numBlock, _blockSize >> > (X._matrixGPU, Q._matrixGPU, Hinv._matrixGPU, matLb._matrixGPU, matUb._matrixGPU, CoresAgentLinBig._matrixGPU, _sizeQ._matrixGPU, _sizeHinv);
		break;
	case 16:
		updateLocalProblPAC< 16> << <numBlock, _blockSize >> > (X._matrixGPU, Q._matrixGPU, Hinv._matrixGPU, matLb._matrixGPU, matUb._matrixGPU, CoresAgentLinBig._matrixGPU, _sizeQ._matrixGPU, _sizeHinv);
		break;
	case  8:
		updateLocalProblPAC<  8> << <numBlock, _blockSize >> > (X._matrixGPU, Q._matrixGPU, Hinv._matrixGPU, matLb._matrixGPU, matUb._matrixGPU, CoresAgentLinBig._matrixGPU, _sizeQ._matrixGPU, _sizeHinv);
		break;
	case  4:
		updateLocalProblPAC<  4> << <numBlock, _blockSize >> > (X._matrixGPU, Q._matrixGPU, Hinv._matrixGPU, matLb._matrixGPU, matUb._matrixGPU, CoresAgentLinBig._matrixGPU, _sizeQ._matrixGPU, _sizeHinv);
		break;
	case  2:
		updateLocalProblPAC<  2> << <numBlock, _blockSize >> > (X._matrixGPU, Q._matrixGPU, Hinv._matrixGPU, matLb._matrixGPU, matUb._matrixGPU, CoresAgentLinBig._matrixGPU, _sizeQ._matrixGPU, _sizeHinv);
		break;
	case  1:
		updateLocalProblPAC<  1> << <numBlock, _blockSize >> > (X._matrixGPU, Q._matrixGPU, Hinv._matrixGPU, matLb._matrixGPU, matUb._matrixGPU, CoresAgentLinBig._matrixGPU, _sizeQ._matrixGPU, _sizeHinv);
		break;
	}
	
	//X.display(true);
	updateXhat();
	//Xhat.display(true);
	updateMu();
	//Mu.display(true);
	//Muhat.display(true);
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 1, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION	

	

}

void PACGPU::updateXhat()
{
	
	if (augmente) {
		//Xhat.subtract(&X, &Xhat);
		Xhat.subtract(&X, &Xpre);
		Xhat.multiply(_alpha);
		Xhat.add(&X);
	}
	else {
		Xhat.set(&X);
	}

}

void PACGPU::updateMu()
{
	int numBlock = _nAgent;

	if (augmente) {
		switch (_blockSize) {
		case 512:
			updatePACMuAug<512> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _phi);
			break;
		case 256:
			updatePACMuAug<256> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _phi);
			break;
		case 128:
			updatePACMuAug<128> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _phi);
			break;
		case 64:
			updatePACMuAug< 64> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _phi);
			break;
		case 32:
			updatePACMuAug< 32> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _phi);
			break;
		case 16:
			updatePACMuAug< 16> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _phi);
			break;
		case  8:
			updatePACMuAug<  8> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _phi);
			break;
		case  4:
			updatePACMuAug<  4> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _phi);
			break;
		case  2:
			updatePACMuAug<  2> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _phi);
			break;
		case  1:
			updatePACMuAug<  1> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _phi);
			break;
		}
	}
	else {
		switch (_blockSize) {
		case 512:
			updatePACMu<512> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _gammahat);
			break;
		case 256:
			updatePACMu<256> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _gammahat);
			break;
		case 128:
			updatePACMu<128> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _gammahat);
			break;
		case 64:
			updatePACMu< 64> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _gammahat);
			break;
		case 32:
			updatePACMu< 32> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _gammahat);
			break;
		case 16:
			updatePACMu< 16> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _gammahat);
			break;
		case  8:
			updatePACMu<  8> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _gammahat);
			break;
		case  4:
			updatePACMu<  4> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _gammahat);
			break;
		case  2:
			updatePACMu<  2> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _gammahat);
			break;
		case  1:
			updatePACMu<  1> << <numBlock, _blockSize >> > (Mu._matrixGPU, Muhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _gammahat);
			break;
		}
	}

	
}



void PACGPU::updateNu()
{
	if (augmente) {
		updateNuAug <<<_nAgent, _blockSize >> > (Nu._matrixGPU, Nuhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresLinTrans._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _theta);
	}
	else {
		updateNuGPU << <_nAgent, _blockSize >> > (Nu._matrixGPU, Nuhat._matrixGPU, Xhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresLinTrans._matrixGPU, CoresIndiceNu._matrixGPU, _rho, _gamma, _gammahat);
	}
}


void PACGPU::updateQ()
{
	updateQAug << <_nAgent, _blockSize >> > (Q._matrixGPU, Qinit._matrixGPU, Xhat._matrixGPU, Muhat._matrixGPU, Nuhat._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresLinTrans._matrixGPU, CoresIndiceNu._matrixGPU, _rhoInv, _rho, _gamma, augmente);
}

float PACGPU::updateRes(int indice)
{
	float resS = 0;
	float resR = 0;
	float resV = 0;
	
	float resTempS = X.max2(&Xhat);
	float resTempR = Muhat.max2(&Mu);
	float resTempV = Nuhat.max2(&Nu);

	if (resTempS > resS) {
		resS = resTempS;
	}
	if (resTempR > resR) {
		resR = resTempR;
	}
	if (resTempV > resV) {
		resV = resTempV;
	}
	//std::cout << "iter : " << indice << " " << resTempS << " " << resTempR << " " << resTempV << " | ";

	int numBlock = _nAgent;
	switch (_blockSize) {
	case 512:
		calcConstraintPAC<512> << <numBlock, _blockSize >> > (tempM._matrixGPU, tempM1._matrixGPU, X._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresLinTrans._matrixGPU, CoresIndiceNu._matrixGPU);
		break;
	case 256:
		calcConstraintPAC<256> << <numBlock, _blockSize >> > (tempM._matrixGPU, tempM1._matrixGPU, X._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresLinTrans._matrixGPU, CoresIndiceNu._matrixGPU);
		break;
	case 128:
		calcConstraintPAC<128> << <numBlock, _blockSize >> > (tempM._matrixGPU, tempM1._matrixGPU, X._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresLinTrans._matrixGPU, CoresIndiceNu._matrixGPU);
		break;
	case 64:
		calcConstraintPAC< 64> << <numBlock, _blockSize >> > (tempM._matrixGPU, tempM1._matrixGPU, X._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresLinTrans._matrixGPU, CoresIndiceNu._matrixGPU);
		break;
	case 32:
		calcConstraintPAC< 32> << <numBlock, _blockSize >> > (tempM._matrixGPU, tempM1._matrixGPU, X._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresLinTrans._matrixGPU, CoresIndiceNu._matrixGPU);
		break;
	case 16:
		calcConstraintPAC< 16> << <numBlock, _blockSize >> > (tempM._matrixGPU, tempM1._matrixGPU, X._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresLinTrans._matrixGPU, CoresIndiceNu._matrixGPU);
		break;
	case  8:
		calcConstraintPAC<  8> << <numBlock, _blockSize >> > (tempM._matrixGPU, tempM1._matrixGPU, X._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresLinTrans._matrixGPU, CoresIndiceNu._matrixGPU);
		break;
	case  4:
		calcConstraintPAC<  4> << <numBlock, _blockSize >> > (tempM._matrixGPU, tempM1._matrixGPU, X._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresLinTrans._matrixGPU, CoresIndiceNu._matrixGPU);
		break;
	case  2:
		calcConstraintPAC<  2> << <numBlock, _blockSize >> > (tempM._matrixGPU, tempM1._matrixGPU, X._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresLinTrans._matrixGPU, CoresIndiceNu._matrixGPU);
		break;
	case  1:
		calcConstraintPAC<  1> << <numBlock, _blockSize >> > (tempM._matrixGPU, tempM1._matrixGPU, X._matrixGPU, nVoisin._matrixGPU, CoresAgentLin._matrixGPU, CoresLinTrans._matrixGPU, CoresIndiceNu._matrixGPU);
		break;
	}/**/
	resTempR = tempM1.max2();
	resTempV = tempM.max2();
	
	//std::cout << X.max2(&Xpre) << std::endl;

	resTempS = X.max2(&Xpre);
	
	//std::cout << " | " << resTempS << " " << resTempR << " " << resTempV << std::endl;

	
	if (resTempS > resS) {
		resS = resTempS;
	}
	if (resTempR > resR) {
		resR = resTempR;
	}
	if (resTempV > resV) {
		resV = resTempV;
	}

	
	resF.set(0, indice, resR);
	resF.set(1, indice, resS);
	resF.set(2, indice, resV);
	return MAX(MAX(resV, resS), resR);
}



void PACGPU::display() {

	std::cout << _name << std::endl;
}


__global__ void updateNuAug(float* Nu, float* Nuhat, float* Xhat, float* nVoisin, float* CoresAgentLin, float* CoresLinTrans, float* CoresindiceNu, float rho, float gamma, float theta) {

	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int agent = blockIdx.x;
	int begin = CoresAgentLin[agent];
	int Mn = nVoisin[agent];
	int indiceNu = CoresindiceNu[begin];
	

	for (int lin = thIdx; lin < Mn; lin+= step) {

		float nuOld = Nu[indiceNu + lin];
		int linPeer = CoresLinTrans[begin + lin + 1];
		float nu = Nuhat[indiceNu + lin] + rho * gamma * (Xhat[begin + lin + 1 + Mn] - Xhat[linPeer]);
		
		
		Nu[indiceNu + lin] = nu;
		Nuhat[indiceNu + lin] = nu + theta * (nu - nuOld);

	}
}


__global__ void updateNuGPU(float* Nu, float* Nuhat, float* Xhat, float* nVoisin, float* CoresAgentLin, float* CoresLinTrans, float* CoresindiceNu, float rho, float gamma, float gammahat) {

	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int agent = blockIdx.x;
	int begin = CoresAgentLin[agent];
	int Mn = nVoisin[agent];
	int indiceNu = CoresindiceNu[begin];

	for (int lin = thIdx; lin < Mn; lin += step) {

		float nuOld = Nu[indiceNu + lin];
		int linPeer = CoresLinTrans[begin + lin + 1];
		float dX = (Xhat[begin + lin + 1 + Mn] - Xhat[linPeer]);
		float nu = rho * gamma * dX + nuOld;
		Nu[indiceNu + lin] =  nu;
		Nuhat[indiceNu + lin] = nu + rho * gammahat * dX;
	}
}


__global__ void updateQAug(float* Q, float* Qinit, float* Xhat, float* Muhat, float* Nuhat, float* nVoisin, float* CoresAgentLin, float* CoresLinTrans, float* CoresIndiceNu, float rhoInv, float rho, float gamma, bool augmente) {
	int agent = blockIdx.x;
	int thIdx = threadIdx.x;
	int step = blockDim.x;

	int begin = CoresAgentLin[agent];
	int Mn = nVoisin[agent];
	int indiceNu = CoresIndiceNu[begin];
	int indiceMu = indiceNu + agent;
	float MuLocal = Muhat[indiceMu];
	
	if (thIdx == 0) {
		Q[begin] = Qinit[begin] + MuLocal - rhoInv * Xhat[begin];
	}

	for (int i = thIdx; i < Mn; i += step) {
		int lin = begin + i + 1;
		int indiceNuTrans = CoresIndiceNu[lin];
		float q  = Qinit[lin]      + Muhat[indiceMu + thIdx + 1] - MuLocal - Nuhat[indiceNuTrans] - rhoInv * Xhat[lin];
		float q2 = Qinit[lin + Mn] + Muhat[indiceMu + thIdx + 1] + Nuhat[indiceNu + i] - rhoInv * Xhat[lin + Mn];

		if (augmente) {	
			int linTrans = CoresLinTrans[lin];
			float tpi = Xhat[linTrans];
			q2 -= rho*gamma* tpi;
		}
		Q[lin] = q;
		Q[lin + Mn] = q2;

	}	
}



template <unsigned int blockSize>
__global__ void updateLocalProblPAC(float* X, float* Q, float* Hinv, float* matLb, float* matUb, float* CoresAgentLinBig, float* sizeQ, int sizeOPFmax) {

	// un bloc par ligne
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int l = blockIdx.x;
	float sum = 0;
	__shared__ float shArr[blockSize];
	int N = sizeQ[l];
	int indiceBegin = CoresAgentLinBig[l];
	for (int i = thIdx; i < N; i += step) {
		sum += Hinv[l * sizeOPFmax + i] * Q[indiceBegin + i];
	}

	shArr[thIdx] = sum;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduce<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		float res = -shArr[0];
		
		float ub = matUb[l];
		float lb = matLb[l];
		X[l] = res + (ub - res) * (res > ub) + (lb - res) * (res < lb);
	}
}


__global__ void calculFcPAC(float* tempN1, float* tempNN, float* a, float* b, float* betaLin, float* X, float* CoresAgentLin, float* nVoisin) {
	int agent = blockIdx.x;
	int nAgent = gridDim.x;
	int thI = threadIdx.x;
	int step = blockDim.x;
	int begin = CoresAgentLin[agent];
	int Mn = nVoisin[agent];

	for (int i = thI; i < Mn; i += step) {
		tempNN[agent * nAgent + thI] = betaLin[begin + i + 1] * X[begin + i + 1];
	}
	
	if (thI == 0) {
		tempN1[agent] = X[begin] * (0.5*a[agent] * X[begin] + b[agent]);
	}

}






__global__ void updateP0PAC(float* matLb, float* matUb, float* Q, float* Qinit, float* Pmin, float* Pmax, float* Cost2, float* Lb, float* Ub, float* CoresAgentLin, float* nVoisin) {
	int agent = blockIdx.x;
	int thI = threadIdx.x;
	int step = blockDim.x;
	int begin = CoresAgentLin[agent];
	int Mn = nVoisin[agent];

	for (int lin = thI; lin < Mn; lin++) {
		int i = lin + begin;
		matLb[i] = Lb[agent];
		matLb[i + Mn] = -Ub[agent];
		matUb[i] = Ub[agent];
		matUb[i + Mn] = -Lb[agent];
	}

	if (thI == 0) {
		float q = Q[begin] - Qinit[begin];
		float a = Cost2[agent];
		Qinit[begin] = a;
		Q[begin] = q + a;
	}

}

