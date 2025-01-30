#include "../head/PAC.h"
#define MAX(X, Y) X * (X >= Y) + Y * (Y > X)


PAC::PAC() : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " PAC Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
}


PAC::PAC(float rho) : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default PAC Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
}

PAC::~PAC()
{
	DELETEA(tempM1);
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
	DELETEA(matUb);
}
void PAC::setParam(float rho)
{
	_rho = rho;
}

void PAC::setGamma(float gamma)
{
	_gamma = gamma;
	if (_gammahat > _gamma) {
		std::cout << "Warning : gammahat should be smaller than gamma";
	}
}

void PAC::setGammahat(float gammahat)
{
	_gammahat = gammahat;
	if (_gammahat > _gamma) {
		std::cout << "Warning : gammahat should be smaller than gamma";
	}
}

void PAC::setInitCoef(float alpha, float phi, float theta)
{
	if (alpha <= 0 || phi < 0 || theta < 0 || alpha >= 1 || phi >= 1 || theta >= 1) {
		throw std::invalid_argument("coefficient must be positive and <1");
	}
	_alpha = alpha;
	_phi = phi;
	_theta = theta;
	
}

void PAC::setBestRhoGamma(float lambdaMax, float lambdaMin, const StudyCase& cas)
{
	MatrixCPU a = cas.geta();
	float eps = 0;
	float alpha = (1 - eps) * a.min2Nnull();
	float L = (1 + eps) * a.max2();
	std::cout << "alpha " << alpha << " and L " << L << std::endl;
	_gamma = (2 * alpha * L) / (2 * lambdaMax + lambdaMin);
	_gammahat = _gamma;
	if (augmente && _alpha>1) {
		_rho = 1 / (sqrt(_gamma * _alpha * lambdaMax));
	}
	else {
		_rho = 1 / (sqrt(_gamma * lambdaMax));
	}
	
	std::cout << "best gamma " << _gamma << " and rho " << _rho << std::endl;

}

void PAC::setBestRhoGammaHeuristic(const StudyCase& cas)
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
	if (lambdaMin <= 0) {
		lambdaMin = 0.000001;
	}
	//std::cout << "lambdaMax " << lambdaMax << " and lambdaMin " << lambdaMin << std::endl;

	MatrixCPU a = cas.geta();
	//a.display();
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


void PAC::setBestRhoGammaHeuristic(const StudyCase& cas, float factor)
{
	if (factor > 1 || factor <= 0) {
		throw std::invalid_argument("setBestRhoGammaHeuristic : factor must be between 0 and 1 exclude");
	}
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
	if (lambdaMin <= 0) {
		lambdaMin = 0.000001;
	}
	//std::cout << "lambdaMax " << lambdaMax << " and lambdaMin " << lambdaMin << std::endl;

	MatrixCPU a = cas.geta();
	//a.display();
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

void PAC::updateCoef()
{
	if (augmente) {

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
			Hinv[i].invertGaussJordan(&H[i]);
		}
	}
	
	
	/**/
	/*_alpha = 0.9 * _alpha;
	_theta = 0.9 * _theta;
	_phi = 0.9 * _phi;*/
}



void PAC::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 0, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION
	}
	_rhog = sim.getRho();
	
	int iterG = sim.getIterG();
	int stepG = sim.getStepG();
	
	float epsG = sim.getEpsG()/2;
	
	float fc = 0;
	float resG = 2 * epsG;
	int iterGlobal = 0;
	//std::cout << "iterG " << iterG <<std::endl;
	while ((iterGlobal < iterG) && (resG>epsG)) {

		
		updateLocalProb();
		
#ifdef INSTRUMENTATION
		occurencePerBlock.increment(0, 1, 1);
		occurencePerBlock.increment(0, 2, 1);
		occurencePerBlock.increment(0, 3, 1);

		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

		updateGlobalProb();

#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 5, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		// FB 4
		if (!(iterGlobal % stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateRes(iterGlobal / stepG);
			//resG = 1;
#ifdef INSTRUMENTATION
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 6, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		}
		//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;

		iterGlobal++;
	}
	//std::cout << iterGlobal  << " " << resF.get(0, (iterGlobal - 1) / stepG) << " " << resF.get(1, (iterGlobal - 1) / stepG) << " " << resF.get(2, (iterGlobal - 1) / stepG) << std::endl;

	/*for (int i = 0; i < _nAgent; i++) {
		std::cout << "agent " << i << std::endl;
		X[i].display();
	}*/



#ifdef INSTRUMENTATION	
	occurencePerBlock.increment(0, 5, iterGlobal);
	occurencePerBlock.increment(0, 6, iterGlobal / stepG);
	

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	
	for (int idAgent = 0; idAgent < _nAgentTrue; idAgent++) {
		MatrixCPU omega(cas.getVoisin(idAgent));
		int Nvoisinmax = nVoisin.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			int idVoisin = omega.get(voisin, 0);
			trade.set(idAgent, idVoisin, X[idAgent].get(voisin + 1, 0));
		}
		Pn.set(idAgent, 0, X[idAgent].get(0, 0));
	}
	int indice = 0;
	for (int idAgent = _nAgentTrue; idAgent < _nAgent; idAgent++) {
		indice = 0;
		for (int idVoisin = 0; idVoisin < _nAgentTrue; idVoisin++) {
			if (idVoisin != (idAgent - _nAgentTrue)) {
				trade.set(idAgent, idVoisin, X[idAgent].get(indice + 1, 0));
				indice = indice + 1;
			}
			Pn.set(idAgent, 0, X[idAgent].get(0, 0));
		}
	}
	//std::cout << "Pn :" << std::endl;
	//Pn.display();
	fc = calcFc(&Cost1, &Cost2, &trade, &Pn, &BETA, &tempN1, &tempNN);
	// FB 5
	
	//std::cout << "set end" << std::endl;


	result->setResF(&resF);
	result->setTrade(&trade); 
	result->setIter(iterGlobal);
	result->setPn(&Pn);
	result->setFc(fc);
	
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 7, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 7, 1);

	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
	tall = clock() - tall;

	result->setTime((float)tall / CLOCKS_PER_SEC);
	
}

void PAC::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	MatrixCPU Pmin(cas.getPmin());
	MatrixCPU Pmax(cas.getPmax());
	MatrixCPU Cost2(cas.getb());

	MatrixCPU Lb(cas.getLb());
	MatrixCPU Ub(cas.getUb());

	for (int i = 0; i < _nAgent; i++) {
		int M = nVoisin.get(i, 0);
		matLb[i].set(0, 0, Pmin.get(i, 0));
		matUb[i].set(0, 0, Pmax.get(i, 0));
		
		float q = Qinit[i].get(0, 0);
		float a = Cost2.get(i, 0);
		Q[i].set(0, 0, Q[i].get(0, 0) - q + a);
		Qinit[i].set(0, 0, a);
		
		for (int m = 0; m < M; m++) {
			matLb[i].set(m + 1, 0, Lb.get(i, 0));
			matLb[i].set(M + m + 1, 0, -Ub.get(i, 0));	
			matUb[i].set(m + 1, 0, Ub.get(i, 0));
			matUb[i].set(M + m + 1, 0, -Lb.get(i, 0));
		}
	}

#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION

}

void PAC::init(const Simparam& sim, const StudyCase& cas)
{
	DELETEA(tempM1);
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
	DELETEA(matUb);


	// intitilisation des matrixs et variables 
	
	
	clock_t t = clock();
	//std::cout << "init " << std::endl;
	if (_rho == 0) {
		_rho = sim.getRho();
	}
	
	
	_rhog = _rho;
	const int iterG = sim.getIterG();
	const int stepG = sim.getStepG();
	float epsG = sim.getEpsG();
	//std::cout << "iterG " << iterG << std::endl;
	isAC = cas.isAC();
	
	_nAgentTrue = sim.getNAgent();
	_nAgent = _nAgentTrue + isAC * _nAgentTrue;
	nVoisin = cas.getNvoi();
	_nTrade = nVoisin.sum();
	_nTradeP = 0;
	if (isAC) {
		for (int n = 0; n < _nAgentTrue; n++) {
			_nTradeP += nVoisin.get(n, 0);
		}
		_nTradeQ = _nTrade - _nTradeP;
		if (_nTradeQ != (_nAgentTrue * (_nAgentTrue - 1))) {
			std::cout << "err PAC : " << _nAgent << " " << _nAgentTrue << " " << _nTrade << " " << _nTradeP << " " << _nTradeQ << std::endl;

			throw std::invalid_argument("Agent must be fully conected for the Q echanges, WIP");
		}
	}
	else {
		_nTradeP = _nTrade;
	}
	//_sizePAC = _nAgent + 2 * _nTrade;


	_rhoInv = 1 / _rhog;
	resF = MatrixCPU(3, (iterG / stepG) + 1);

	BETA = MatrixCPU(cas.getBeta());
	MatrixCPU Ub(cas.getUb());
	MatrixCPU Lb(cas.getLb());
	MatrixCPU Pmin(cas.getPmin());
	MatrixCPU Pmax(cas.getPmax());
	Cost1 = cas.geta();
	Cost2 = cas.getb();

	trade = sim.getTrade();
	Pn = sim.getPn();

	CoresMatLin = MatrixCPU(_nAgent, _nAgentTrue, -1);
	CoresLinAgent = MatrixCPU(_nTrade, 1);
	CoresAgentLin = MatrixCPU(_nAgent + 1, 1);
	CoresLinVoisin = MatrixCPU(_nTrade, 1);
	CoresLinTrans = MatrixCPU(_nTrade, 1);
	CoresLinTransLocal = MatrixCPU(_nTrade, 1);

	// Who is the peer ? 
	int indice = 0;
	//std::cout << " P " << std::endl;
	for (int idAgent = 0; idAgent < _nAgentTrue; idAgent++) { // P
		MatrixCPU omega(cas.getVoisin(idAgent));
		int Nvoisinmax = nVoisin.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			int idVoisin = omega.get(voisin, 0);
			CoresLinAgent.set(indice, 0, idAgent);
			CoresLinVoisin.set(indice, 0, idVoisin); // 
			CoresMatLin.set(idAgent, idVoisin, indice);
			indice = indice + 1;
		}
		CoresAgentLin.set(idAgent + 1, 0, indice);
	}
	//std::cout << " Q " << std::endl;
	for (int idAgent = _nAgentTrue; idAgent < _nAgent; idAgent++) { // Q
		for (int idVoisin = 0; idVoisin < _nAgentTrue; idVoisin++) {
			if (idVoisin != (idAgent - _nAgentTrue)) {
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
		int k = 0;
		if(i < _nAgentTrue) // P
			k = CoresMatLin.get(j, i);
		else {
			k = CoresMatLin.get(j, i % _nAgentTrue);
		}
		
		//std::cout << " trade num " << lin << " entre " << i % _nAgentTrue << " et " << j;
		//std::cout << " le trade symetrique le trade num " << k << std::endl;
		CoresLinTrans.set(lin, 0, k);
	}

	for (int i = 0; i < _nAgent; i++) {
		int M = nVoisin.get(i, 0);
		for (int m = 0; m < M; m++) {
			int lin = CoresAgentLin.get(i, 0) + m; // indice global de tim = tip (m est le num�ro du voisin, p est le num�ro de l'agent)
			int p = CoresLinVoisin.get(lin, 0); // valeur de p
			int lin2 = CoresLinTrans.get(lin, 0); // indice global de tpi
			int linLoc = lin2 - CoresAgentLin.get(p, 0); //indice local de tpi
			CoresLinTransLocal.set(lin, 0, linLoc);
			//std::cout << lin << ": trade entre " << i << " et " << p << " num local " << m ;
			//std::cout << " le symetrique est " << lin2 <<  " ou " << linLoc << std::endl;
		}
	}
	/*CoresMatLin.display();
	CoresAgentLin.display();
	CoresLinAgent.display();
	CoresLinVoisin.display();
	CoresLinTrans.display();
	CoresLinTransLocal.display();*/

	//std::cout << "creation " << std::endl;
	X = new MatrixCPU[_nAgent];
	Xpre = new MatrixCPU[_nAgent];
	Xhat = new MatrixCPU[_nAgent];
	Mu = new MatrixCPU[_nAgent];
	Muhat = new MatrixCPU[_nAgent];
	Nu = new MatrixCPU[_nAgent];
	Nuhat = new MatrixCPU[_nAgent];


	tempM1 = new MatrixCPU[_nAgent];
	tempM = new MatrixCPU[_nAgent];


	Hinv = new MatrixCPU[_nAgent];
	H = new MatrixCPU[_nAgent];
	Q = new MatrixCPU[_nAgent];
	Qinit = new MatrixCPU[_nAgent];

	matLb = new MatrixCPU[_nAgent];
	matUb = new MatrixCPU[_nAgent];
	
	tempN1 = MatrixCPU(_nAgent, 1);
	tempNN = MatrixCPU(_nAgent, _nAgentTrue);
	//std::cout << "problem setup " << std::endl;
	for (int i = 0; i < _nAgent; i++) {
		//std::cout << "********* Agent " << i << "**********" << std::endl;
		// def
		int M = nVoisin.get(i, 0);
		int begin = CoresAgentLin.get(i, 0);
		X[i] = MatrixCPU(1 + 2 * M, 1);
		Xpre[i] = MatrixCPU(1 + 2 * M, 1);
		Xhat[i] = MatrixCPU(1 + 2 * M, 1);
		matLb[i] = MatrixCPU(1 + 2 * M, 1);
		matUb[i] = MatrixCPU(1 + 2 * M, 1);

		Mu[i] = MatrixCPU(1 + M, 1);
		Muhat[i] = MatrixCPU(1 + M, 1);
		Nu[i] = MatrixCPU(M, 1);
		Nuhat[i] = MatrixCPU(M, 1);	

		tempM1[i] = MatrixCPU(1 + M, 1);
		tempM[i] = MatrixCPU(M, 1);
		
		H[i] = MatrixCPU(1 + 2 * M, 1 + 2 * M);
		Hinv[i] = MatrixCPU(1 + 2 * M, 1 + 2 * M);
		Q[i] = MatrixCPU(1 + 2 * M, 1, 0);
		Qinit[i] = MatrixCPU(1 + 2 * M, 1, 0);

		// init
		X[i].set(0, 0, Pn.get(i, 0));
		matLb[i].set(0, 0, Pmin.get(i, 0));
		matUb[i].set(0, 0, Pmax.get(i, 0));
		H[i].set(0, 0, Cost1.get(i, 0) + _rhoInv);
		if (augmente) {
			H[i].increment(0, 0, _rho * _gamma);
		}
				
		Qinit[i].set(0, 0, Cost2.get(i, 0));
		
		for (int m = 0; m < M; m++) {
			int voisin = CoresLinVoisin.get(begin + m, 0);
			//std::cout << "voisin num " << m << " is " << voisin << std::endl;

			X[i].set(m + 1, 0, trade.get(i, voisin % _nAgentTrue)); // tnm
			X[i].set(M + m + 1, 0, trade.get(voisin, i % _nAgentTrue)); //amn
			matLb[i].set(m + 1, 0, Lb.get(i, 0));
			matLb[i].set(M + m + 1, 0, -Ub.get(i, 0)); // est ce que cela g�ne la convergence ou est ce que cela l'aide ?
			matUb[i].set(m + 1, 0, Ub.get(i, 0)); 
			matUb[i].set(M + m + 1, 0, -Lb.get(i, 0));
			
			H[i].set(m + 1, m + 1, _rhoInv); // diag tnm
			H[i].set(M + m + 1, M + m + 1, _rhoInv); // diag anm
			if (augmente) {
				H[i].set(m + 1, 0, -_rho * _gamma); // first column pn <-> tnm
				H[i].set(0, m + 1, -_rho * _gamma); // fisrt row    pn <-> tnm

				H[i].increment(m + 1, m + 1, 2 * _rho * _gamma); // diag tnm
				H[i].increment(M + m + 1, M + m + 1, 2 * _rho * _gamma); // diag anm
			}
			if (i < _nAgentTrue) {
				Qinit[i].set(m + 1, 0, BETA.get(i, voisin));
			}
		}
		
		if (augmente) {
			for (int k = 1; k < M + 1; k++) {
				for (int j = 1; j < M + 1; j++) {
					if (k != j) {
						H[i].set(k, j, _rho * _gamma); // sum(sum(tnk*tnj)) -> ||Ga||
					}
				}
			}
			for (int k = 1; k < M + 1; k++) {
				int j = M + k;
				H[i].set(k, j, _rho * _gamma);
				H[i].set(j, k, _rho * _gamma);
			}
		}	
		Xpre[i].set(&X[i]);
		//H[i].display();
		Hinv[i].invertGaussJordan(&H[i]);
		//Hinv[i].display();
		//std::cout << "****" << std::endl;
	}
	/*for (int i = 0; i < _nAgent; i++) {
		std::cout << "agent " << i << std::endl;
		Qinit[i].display();
	}*/
	
	//std::cout << "Global update" << std::endl;

	updateGlobalProb();

	

	
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
}

float PAC::calcFc(MatrixCPU* cost1, MatrixCPU* cost2, MatrixCPU* trade, MatrixCPU* Pn, MatrixCPU* BETA, MatrixCPU* tempN1, MatrixCPU* tempNN)
{
	float fc = 0;
	tempN1->set(cost1);
	tempN1->multiply(0.5);
	tempN1->multiplyT(Pn);
	tempN1->add(cost2);
	tempN1->multiplyT(Pn);
	

	fc = fc + tempN1->sum();
	for (int i = 0; i < _nAgentTrue; i++) {
		for (int j = 0; j < _nAgentTrue; j++) {
			tempNN->set(i, j, trade->get(i, j) * BETA->get(i, j));
		}
	}
	for (int i = _nAgentTrue; i < _nAgent; i++) {
		for (int j = 0; j < _nAgentTrue; j++) {
			tempNN->set(i, j, 0);
		}
	}

	fc = fc + tempNN->sum();

	return fc;
}

void PAC::setBestParam(const StudyCase& cas)
{
	setBestRhoGammaHeuristic(cas);
}

void PAC::updateGlobalProb() {

	
	// communication of xhat
	updateNu();
	/*for (int i = 0; i < _nAgent; i++) {
		int m = nVoisin.get(i, 0);
		for (int j = 0; j < m ; j++) {
			std::cout << Nu[i].get(j, 0) << " ";
		}
	}
	std::cout << std::endl;
	for (int i = 0; i < _nAgent; i++) {
		int m = nVoisin.get(i, 0);
		for (int j = 0; j < m; j++) {
			std::cout << Nuhat[i].get(j, 0) << " ";
		}
	}
	std::cout << std::endl;*/
	// communication of nuhat
	updateQ();
	if (augmente) {
		//updateCoef();
	}
	/*for (int i = 0; i < _nAgent; i++) {
		int m = nVoisin.get(i, 0);
		for (int j = 0; j < 2 * m + 1; j++) {
			std::cout << Q[i].get(j, 0) << " ";
		}
	}
	std::cout << std::endl;
	std::cout << "******" << std::endl;*/
}

void PAC::updateLocalProb() {
	// FB 1a
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	

	for (int i = 0; i < _nAgent; i++) { 
		Xpre[i].swap(&X[i]);
		X[i].MultiplyMatVec(&Hinv[i], &Q[i]); // solve system by using the inverse
		X[i].multiply(-1);
	}
	

#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 1, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 1b
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION	
	for (int i = 0; i < _nAgent; i++) {
		X[i].project(&matLb[i], &matUb[i]);
	}
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 2, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION	
	
	/*for (int i = 0; i < _nAgent; i++) {
		int m = nVoisin.get(i, 0);
		for (int j = 0; j < 2*m + 1; j++) {
			std::cout << X[i].get(j, 0) << " ";
		}
	}
	std::cout << std::endl;*/
	updateXhat();
	/*for (int i = 0; i < _nAgent; i++) {
		int m = nVoisin.get(i, 0);
		for (int j = 0; j < 2 * m + 1; j++) {
			std::cout << Xhat[i].get(j, 0) << " ";
		}
	}
	std::cout << std::endl;*/
	updateMu();
	/*for (int i = 0; i < _nAgent; i++) {
		int m = nVoisin.get(i, 0);
		for (int j = 0; j < m + 1; j++) {
			std::cout << Mu[i].get(j, 0) << " ";
		}
	}
	std::cout << std::endl;
	for (int i = 0; i < _nAgent; i++) {
		int m = nVoisin.get(i, 0);
		for (int j = 0; j < m + 1; j++) {
			std::cout << Muhat[i].get(j, 0) << " ";
		}
	}
	std::cout << std::endl;*/

#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 3, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
}

void PAC::updateXhat()
{
	for (int i = 0; i < _nAgent; i++) {
		if (augmente) {
			//Xhat[i].subtract(&X[i], &Xhat[i]);
			Xhat[i].subtract(&X[i], &Xpre[i]);
			Xhat[i].multiply(_alpha);
			Xhat[i].add(&X[i]);
		}
		else {
			Xhat[i].set(&X[i]);
		}
	}
}

void PAC::updateMu()
{
	for (int i = 0; i < _nAgent; i++) {
		tempM1[i].set(&Mu[i]); // mu^k
		double pSum = 0;
		int M = nVoisin.get(i, 0);
		for (int m = 0; m < M; m++) {// Gn*xhat
			pSum += Xhat[i].get(1 + m, 0);
			//Mu[i].set(m + 1, 0, (Xhat[i].get(1 + m, 0) - Xhat[i].get(1 + m + M, 0)) / 2 - Xpre[i].get(1 + m, 0)); ne marche pas ?
			Mu[i].set(m + 1, 0, Xhat[i].get(1 + m, 0) + Xhat[i].get(1 + m + M, 0));
		}
		Mu[i].set(0, 0, Xhat[i].get(0, 0) - pSum);
		

		if (augmente) {
			Mu[i].multiply(_rho * _gamma);
			Mu[i].add(&Muhat[i]);

			Muhat[i].subtract(&Mu[i], &tempM1[i]);
			Muhat[i].multiply(_phi);
			Muhat[i].add(&Mu[i]);
		}
		else {
			Muhat[i].set(&Mu[i]);
			Mu[i].multiply(_rho * _gamma);
			Muhat[i].multiply(_rho * _gammahat);
			Mu[i].add(&tempM1[i]);
			Muhat[i].add(&Mu[i]);
		}
	}
}



void PAC::updateNu()
{
	for (int i = 0; i < _nAgent; i++) {
		
		tempM[i].set(&Nu[i]); // nu^k
		// Gn*xhat
		
		int M = nVoisin.get(i, 0);
		for (int m = 0; m < M; m++) {
			
			// find t^[p]mi ???
			int lin = CoresAgentLin.get(i, 0) + m; // indice global de tim = tip (m est le num�ro du voisin, p est le num�ro de l'agent)
			int p = CoresLinVoisin.get(lin, 0); // valeur de p

			int linLoc = CoresLinTransLocal.get(lin, 0);
			//int lin2 = CoresLinTrans.get(lin, 0); // indice global de tpi
			//int linLoc = lin2 - CoresAgentLin.get(p, 0); //indice local de tpi

			float tpi = Xhat[p].get(linLoc + 1, 0);
			// Bj aj
			Nu[i].set(m, 0, Xhat[i].get(1 + m + M, 0) - tpi); // 
		}
		
		if (augmente) {
			Nu[i].multiply(_rho * _gamma);
			Nu[i].add(&Nuhat[i]);
			Nuhat[i].subtract(&Nu[i], &tempM[i]);
			Nuhat[i].multiply(_theta);
			Nuhat[i].add(&Nu[i]);
		}
		else {
			Nuhat[i].set(&Nu[i]);
			Nu[i].multiply(_rho * _gamma);
			Nuhat[i].multiply(_rho * _gammahat);

			Nu[i].add(&tempM[i]);
			Nuhat[i].add(&Nu[i]);
		}/**/
	}

}


void PAC::updateQ()
{
	for (int i = 0; i < _nAgent; i++) {
		Q[i].set(&Qinit[i]);
		Q[i].increment(0, 0, Muhat[i].get(0, 0) - _rhoInv * Xhat[i].get(0, 0));

		int M = nVoisin.get(i, 0);
		for (int m = 0; m < M; m++) {
			int lin = CoresAgentLin.get(i, 0) + m; // indice global de tim = tip (m est le num�ro du voisin, p est le num�ro de l'agent)
			int p = CoresLinVoisin.get(lin, 0); // valeur de p

			int linLoc = CoresLinTransLocal.get(lin, 0);

			Q[i].increment(m + 1    , 0, Muhat[i].get(m + 1, 0) - Muhat[i].get(0, 0) - Nuhat[p].get(linLoc, 0) - _rhoInv * Xhat[i].get(m + 1, 0));
			Q[i].increment(M + m + 1, 0, Muhat[i].get(m + 1, 0) + Nuhat[i].get(m, 0) - _rhoInv * Xhat[i].get(M + m + 1, 0));
			//Q[i].increment(m + 1, 0, Muhat[i].get(m + 1, 0) - Muhat[i].get(0, 0)  - _rhoInv * X[i].get(m + 1, 0));
			//Q[i].increment(M + m + 1, 0, Muhat[i].get(m + 1, 0) + Nuhat[i].get(m, 0) - _rhoInv * X[i].get(M + m + 1, 0));

			if (augmente) {	
				float tpi = Xhat[p].get(linLoc + 1, 0);
				Q[i].increment(M + m + 1, 0, -_rho * _gamma * tpi); 
			}
		}
	}
}

float PAC::updateRes(int indice)
{
	float resS = 0;
	float resR = 0;
	float resV = 0;
	double pSum = 0;
	for (int i = 0; i < _nAgent; i++) {
		float resTempS = Xhat[i].max2(&X[i]);
		float resTempR = Muhat[i].max2(&Mu[i]);
		float resTempV = Nuhat[i].max2(&Nu[i]);

		if (resTempS > resS) {
			resS = resTempS;
		}
		if (resTempR > resR) {
			resR = resTempR;
		}
		if (resTempV > resV) {
			resV = resTempV;
		}/**/
		//std::cout << "iter : " << indice << " agent " << i << " " << resTempS << " " << resTempR << " " << resTempV << " | ";
		
		pSum = 0;
		int M = nVoisin.get(i, 0);
		for (int m = 0; m < M; m++) {
			pSum += X[i].get(1 + m, 0);
			tempM1[i].set(m + 1, 0, X[i].get(1 + m, 0) + X[i].get(1 + m + M, 0));

			int lin = CoresAgentLin.get(i, 0) + m; // indice global de tim = tip (m est le num�ro du voisin, p est le num�ro de l'agent)
			int p = CoresLinVoisin.get(lin, 0); // valeur de p
			int linLoc = CoresLinTransLocal.get(lin, 0);
			float tpi = X[p].get(linLoc + 1, 0);
			tempM[i].set(m, 0, X[i].get(1 + m + M, 0) - tpi); // 
		}
		tempM1[i].set(0, 0, X[i].get(0, 0) - pSum);
		resTempR = tempM1[i].max2();
		resTempV = tempM[i].max2();
		

		resTempS = X[i].max2(&Xpre[i]);
		

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
	}
	resF.set(0, indice, resR);
	resF.set(1, indice, resS);
	resF.set(2, indice, resV);

	return MAX(MAX(resV, resS), resR);
}




void PAC::display() {

	std::cout << _name << std::endl;
}
