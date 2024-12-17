#include "../head/PACConst.h"
#define MAX(X, Y) X * (X >= Y) + Y * (Y > X)


PACConst::PACConst() : Method()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " PACConst Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
}


PACConst::PACConst(float rho) : Method()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default PACConst Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
}

PACConst::~PACConst()
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


void PACConst::setParam(float rho)
{
	_rho = rho;
}

void PACConst::setGamma(float gamma)
{
	_gamma = gamma;
	if (_gammahat > _gamma) {
		std::cout << "Warning : gammahat should be smaller than gamma";
	}
}

void PACConst::setGammahat(float gammahat)
{
	_gammahat = gammahat;
	if (_gammahat > _gamma) {
		std::cout << "Warning : gammahat should be smaller than gamma";
	}
}

void PACConst::setInitCoef(float alpha, float phi, float theta)
{
	if (alpha <= 0 || phi < 0 || theta < 0 || alpha >= 1 || phi >= 1 || theta >= 1) {
		throw std::invalid_argument("coefficient must be positive and <1");
	}
	_alpha = alpha;
	_phi = phi;
	_theta = theta;
	
}

void PACConst::setBestRhoGamma(float lambdaMax, float lambdaMin, const StudyCase& cas)
{
	MatrixCPU a = cas.geta();
	float alpha = 1 * a.min2();
	float L = 1 * a.max2();
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

void PACConst::updateCoef()
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
			Hinv[i].invertEigen(&H[i]);
		}
	}
	
	
	/**/
	/*_alpha = 0.9 * _alpha;
	_theta = 0.9 * _theta;
	_phi = 0.9 * _phi;*/
}



void PACConst::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
	/*std::cout << iterGlobal << " Copies " << resF.get(0, (iterGlobal - 1) / stepG) << " Convergence " << resF.get(1, (iterGlobal - 1) / stepG) << " Contraintes " << resF.get(2, (iterGlobal - 1) / stepG) << std::endl;

	std::cout << "X :" << std::endl;
	for (int i = 0; i < _nAgent + 1; i++) {
		X[i].display();
	}*/



#ifdef INSTRUMENTATION	
	occurencePerBlock.increment(0, 5, iterGlobal);
	occurencePerBlock.increment(0, 6, iterGlobal / stepG);
	

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	int indice = 0;
	for (int idAgent = 0; idAgent < _nAgent; idAgent++) {
		MatrixCPU omega(cas.getVoisin(idAgent));
		int Nvoisinmax = nVoisin.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			int idVoisin = omega.get(voisin, 0);
			trade.set(idAgent, idVoisin, X[idAgent].get(voisin + 1, 0));
		}
		Pn.set(idAgent, 0, X[idAgent].get(0, 0));
	}
	//std::cout << "Pn :" << std::endl;
	//Pn.display();
	fc = calcFc(&Cost1, &Cost2, &trade, &Pn, &BETA, &tempN1, &tempNN);
	// FB 5
	
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

void PACConst::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	MatrixCPU Pmin(cas.getPmin());
	MatrixCPU Pmax(cas.getPmax());
	MatrixCPU Cost2(cas.getb());

	MatrixCPU Lb(cas.getLb());

	for (int i = 0; i < _nAgent; i++) {
		int M = nVoisin.get(i, 0);
		matLb[i].set(0, 0, Pmin.get(i, 0));
		matUb[i].set(0, 0, Pmax.get(i, 0));
		
		Q[i].set(0, 0, Cost2.get(i, 0) + Mu[i].get(0, 0) - _rhoInv * X[i].get(0, 0));
		Qinit[i].set(0, 0, Cost2.get(i, 0));
		
		for (int m = 0; m < M; m++) {
			matLb[i].set(m + 1, 0, Lb.get(i, 0));
			matLb[i].set(M + m + 1, 0, Lb.get(i, 0));	
		}
	}

	
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);

}

void PACConst::init(const Simparam& sim, const StudyCase& cas)
{
	// intitilisation des matrixs et variables 
	
	clock_t t = clock();
	std::cout << "init " << std::endl;
	_rhog = sim.getRho();
	
	const int iterG = sim.getIterG();
	const int stepG = sim.getStepG();
	float epsG = sim.getEpsG();
	
	_nAgent = sim.getNAgent();
	nVoisin = cas.getNvoi();
	_nLine = cas.getNLine();
	_nTrade = nVoisin.sum();
	_sizePACConst = _nAgent + 2 * _nTrade;


	_rhoInv = 1 / _rho;
	resF = MatrixCPU(3, (iterG / stepG) + 1);

	BETA = MatrixCPU(cas.getBeta());
	MatrixCPU Ub(cas.getUb());
	MatrixCPU Lb(cas.getLb());
	MatrixCPU Pmin(cas.getPmin());
	MatrixCPU Pmax(cas.getPmax());
	Cost1 = MatrixCPU(cas.geta());
	Cost2 = MatrixCPU(cas.getb());

	trade = sim.getTrade();
	Pn = sim.getPn();

	X = new MatrixCPU[_nAgent + 1];
	Xpre = new MatrixCPU[_nAgent + 1];
	Xhat = new MatrixCPU[_nAgent + 1];
	Mu = new MatrixCPU[_nAgent];
	Muhat = new MatrixCPU[_nAgent];
	Nu = new MatrixCPU[_nAgent + 1];
	Nuhat = new MatrixCPU[_nAgent + 1];
	delta = MatrixCPU(2*_nLine, 1);
	deltahat = MatrixCPU(2*_nLine, 1);

	tempN1 = MatrixCPU(_nAgent, 1);
	tempNN = MatrixCPU(_nAgent, _nAgent);
	tempL = MatrixCPU(2*_nLine, 1);
	tempM1 = new MatrixCPU[_nAgent];
	tempM = new MatrixCPU[_nAgent];
	

	Hinv = new MatrixCPU[_nAgent + 1];
	H = new MatrixCPU[_nAgent + 1];
	Q = new MatrixCPU[_nAgent + 1];
	Qinit = new MatrixCPU[_nAgent + 1];


	CoresMatLin = MatrixCPU(_nAgent, _nAgent, -1);
	CoresLinAgent = MatrixCPU(_nTrade, 1);
	CoresAgentLin = MatrixCPU(_nAgent + 1, 1);
	CoresLinVoisin = MatrixCPU(_nTrade, 1);
	CoresLinTrans = MatrixCPU(_nTrade, 1);
	CoresLinTransLocal = MatrixCPU(_nTrade, 1);
	G = cas.getPowerSensi();
	Phi = MatrixCPU(_nLine, 1);
	lLimit = cas.getLineLimit();
	
	matLb = new MatrixCPU[_nAgent];
	matUb = new MatrixCPU[_nAgent];

	int indice = 0;
	for (int i = 0; i < _nAgent; i++) {
		// def
		int M = nVoisin.get(i, 0);
		MatrixCPU omega(cas.getVoisin(i));
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
			int j = omega.get(m, 0);
			X[i].set(m + 1, 0, trade.get(i, j)); // tnm
			X[i].set(M + m + 1, 0, trade.get(j, i)); //amn
			matLb[i].set(m + 1, 0, Lb.get(i, 0));
			matLb[i].set(M + m + 1, 0, -Ub.get(i, 0)); // est ce que cela gêne la convergence ou est ce que cela l'aide ?
			matUb[i].set(m + 1, 0, Ub.get(i, 0)); 
			matUb[i].set(M + m + 1, 0, -Lb.get(i, 0));
			
			H[i].set(m + 1, m + 1, _rhoInv); // diag tnm
			H[i].set(M + m + 1, M + m + 1, _rhoInv); // diag anm
			if (augmente) {
				H[i].set(m + 1, 0, -_rho * _gamma); // first column pn <-> tnm ||Ga||
				H[i].set(0, m + 1, -_rho * _gamma); // fisrt row    pn <-> tnm ||Ga||

				H[i].increment(m + 1, m + 1, 2 * _rho * _gamma); // diag tnm
				H[i].increment(M + m + 1, M + m + 1, 2 * _rho * _gamma); // diag anm
			}
			Qinit[i].set(m + 1, 0, BETA.get(i, j));

			CoresLinAgent.set(indice, 0, i);
			CoresLinVoisin.set(indice, 0, j);
			CoresMatLin.set(i, j, indice);
			indice = indice + 1;
		}
		CoresAgentLin.set(i + 1, 0, indice);
		if (augmente) {
			for (int k = 1; k < M + 1; k++) {
				for (int j = 1; j < M + 1; j++) {
					if (k != j) {
						H[i].set(k, j, _rho * _gamma);
					}
				}
			}
			for (int k = 1; k < M + 1; k++) {
				int j = M + k;
				H[i].set(k, j, _rho * _gamma); // sum(amn*tnm)
				H[i].set(j, k, _rho * _gamma);
			}
		}		
		Xpre[i].set(&X[i]);
		Hinv[i].invertEigen(&H[i]);
	}

	X[_nAgent] = MatrixCPU(_nAgent, 1);
	Xpre[_nAgent] = MatrixCPU(_nAgent, 1);
	Xhat[_nAgent] = MatrixCPU(_nAgent, 1);
	Nu[_nAgent] = MatrixCPU(_nAgent, 1);
	Nuhat[_nAgent] = MatrixCPU(_nAgent, 1);
	
	Hinv[_nAgent] = MatrixCPU(_nAgent, _nAgent);
	H[_nAgent] = MatrixCPU(_nAgent, _nAgent);
	Q[_nAgent] = MatrixCPU(_nAgent, 1);
	Qinit[_nAgent] = MatrixCPU(_nAgent, 1);
	if (augmente) {
		H[_nAgent].multiplyTrans(&G, &G, 0);
		H[_nAgent].multiply(2 * _rho * _gamma);
	}

	for (int i = 0; i < _nAgent; i++) {
		H[_nAgent].increment(i, i, _rhoInv);
		X[_nAgent].set(i, 0, Pn.get(i, 0));
		Xpre[_nAgent].set(i, 0, Pn.get(i, 0));
		if (augmente) {
			H[_nAgent].increment(i, i, _rho * _gamma);
		}
	}
	
	
	Hinv[_nAgent].invertEigen(&H[_nAgent]);
	//std::cout << "mise sous forme linéaire" << std::endl;
	
	for (int lin = 0; lin < _nTrade; lin++) {
		int i = CoresLinAgent.get(lin, 0);
		int j = CoresLinVoisin.get(lin, 0);
		int k = CoresMatLin.get(j, i);
		CoresLinTrans.set(lin, 0, k);
	}

	for (int i = 0; i < _nAgent; i++) {
		int M = nVoisin.get(i, 0);
		for (int m = 0; m < M; m++) {
			int lin = CoresAgentLin.get(i, 0) + m; // indice global de tim = tip (m est le numéro du voisin, p est le numéro de l'agent)
			int p = CoresLinVoisin.get(lin, 0); // valeur de p

			int lin2 = CoresLinTrans.get(lin, 0); // indice global de tpi
			int linLoc = lin2 - CoresAgentLin.get(p, 0); //indice local de tpi

			CoresLinTransLocal.set(lin, 0, linLoc);
		}
	}
	updateGlobalProb();
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
}

void PACConst::updateGlobalProb()  {
	//std::cout << "Xhat" << std::endl;
	updateXhat();
	//std::cout << "Mu" << std::endl;
	updateMu();
	//std::cout << "Delta" << std::endl;
	updateDelta();
	// communication of xhat
	//std::cout << "Nu" << std::endl;
	updateNu();
	// communication of nuhat
	//std::cout << "Q" << std::endl;
	updateQ();
	if (augmente) {
		//updateCoef();
	}

}

void PACConst::updateLocalProb() {
	// FB 1a
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	for (int i = 0; i < _nAgent + 1; i++) { 
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
	
	// FB 1c
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION	
 
	// nothing 

#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 3, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
}

void PACConst::updateXhat()
{
	for (int i = 0; i < _nAgent + 1; i++) {
		if (augmente) {
			Xhat[i].subtract(&X[i], &Xhat[i]);
			//Xhat[i].subtract(&X[i], &Xpre[i]);
			Xhat[i].multiply(_alpha);
			Xhat[i].add(&X[i]);
		}
		else {
			Xhat[i].set(&X[i]);
		}
	}
}

void PACConst::updateMu()
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


void PACConst::updateNu()
{

	tempN1.set(&Nu[_nAgent]); // nu^k
	for (int i = 0; i < _nAgent; i++) {
		tempM[i].set(&Nu[i]); // nu^k
		int M = nVoisin.get(i, 0);
		if (i < _nAgent) {
			for (int m = 0; m < M; m++) {
				// find t^[p]mi ???
				int lin = CoresAgentLin.get(i, 0) + m; // indice global de tim = tip (m est le numéro du voisin, p est le numéro de l'agent)
				int p = CoresLinVoisin.get(lin, 0); // valeur de p

				int lin2 = CoresLinTrans.get(lin, 0); // indice global de tpi
				int linLoc = lin2 - CoresAgentLin.get(p, 0); //indice local de tpi

				float tpi = Xhat[p].get(linLoc + 1, 0);
				// Bj aj
				Nu[i].set(m, 0, Xhat[i].get(1 + m + M, 0) - tpi); // 
			}
		}
		Nu[_nAgent].set(i, 0, Xhat[_nAgent].get(i, 0) - Xhat[i].get(0, 0)); // 
	}

	for (int i = 0; i < _nAgent + 1; i++) {
		if (augmente) {
			Nu[i].multiply(_rho * _gamma);
			Nu[i].add(&Nuhat[i]);
			if (i == _nAgent) {
				Nuhat[i].subtract(&Nu[i], &tempN1);
			}
			else {
				Nuhat[i].subtract(&Nu[i], &tempM[i]);
			}
			
			Nuhat[i].multiply(_theta);
			Nuhat[i].add(&Nu[i]);
		}
		else {
			Nuhat[i].set(&Nu[i]);
			Nu[i].multiply(_rho * _gamma);
			Nuhat[i].multiply(_rho * _gammahat);
			if (i == _nAgent) {
				Nu[i].add(&tempN1);
			}
			else {
				Nu[i].add(&tempM[i]);
			}
			
			Nuhat[i].add(&Nu[i]);
		}
	}
}


void PACConst::updateQ()
{
	for (int i = 0; i < _nAgent; i++) {
		Q[i].set(&Qinit[i]);
		Q[i].increment(0, 0, Muhat[i].get(0, 0) - _rhoInv * Xhat[i].get(0, 0) - Nuhat[_nAgent].get(i,0));

		int M = nVoisin.get(i, 0);
		for (int m = 0; m < M; m++) {
			int lin = CoresAgentLin.get(i, 0) + m; // indice global de tim = tip (m est le numéro du voisin, p est le numéro de l'agent)
			int p = CoresLinVoisin.get(lin, 0); // valeur de p

			int lin2 = CoresLinTrans.get(lin, 0); // indice global de tpi
			int linLoc = lin2 - CoresAgentLin.get(p, 0); //indice local de tpi

			float tpi = Xhat[p].get(linLoc + 1, 0);

			Q[i].increment(m + 1, 0, Muhat[i].get(m + 1, 0) - Muhat[i].get(0, 0) - Nuhat[p].get(linLoc,0) - _rhoInv * Xhat[i].get(m + 1, 0));
			Q[i].increment(M + m + 1, 0, Muhat[i].get(m + 1, 0) + Nuhat[i].get(m, 0) - _rhoInv * Xhat[i].get(M + m + 1, 0));
			

			if (augmente) {
				Q[i].increment(M + m + 1, 0, -_rho * _gamma * tpi);
			}
		}

		double sumdeltaG = 0;

		for (int l = 0; l < _nLine; l++) {
			sumdeltaG += (delta.get(l, 0) - delta.get(_nLine + l, 0)) * G.get(l, i);
		}

		Q[_nAgent].set(i, 0, sumdeltaG + Nuhat[_nAgent].get(i, 0) - _rhoInv * Xhat[_nAgent].get(i, 0));
		if (augmente) {
			Q[_nAgent].increment(i, 0, -_rho * _gamma * Xhat[i].get(0,0) );
		}


	}
	


}

void PACConst::updateDelta()
{
	Phi.MultiplyMatVec(&G, &Xhat[_nAgent]);
	tempL.set(&delta);

	for (int l = 0; l < _nLine; l++) {
		delta.set(l, 0, max((Phi.get(l, 0) - lLimit.get(l, 0)), 0));
		delta.set(_nLine + l, 0, max((-Phi.get(l, 0) - lLimit.get(l, 0)), 0));
	}

	if (augmente) {
		delta.multiply(_rho * _gamma);
		delta.add(&deltahat);
		deltahat.subtract(&delta, &tempL);
		deltahat.multiply(_theta);
		deltahat.add(&delta);
	}
	else {
		deltahat.set(&delta);
		delta.multiply(_rho * _gamma);
		deltahat.multiply(_rho * _gammahat);

		delta.add(&tempL);
		deltahat.add(&delta);
	}


}

float PACConst::updateRes(int indice)
{
	float resS = 0;
	float resR = 0;
	float resV = 0;
	for (int i = 0; i < _nAgent; i++) {
		float resTempS = 0;// Xhat[i].max2(&X[i]);
		float resTempR = 0;// Muhat[i].max2(&Mu[i]);
		float resTempV = 0;//Nuhat[i].max2(&Nu[i]);

		if (resTempS > resS) {
			resS = resTempS;
		}
		if (resTempR > resR) {
			resR = resTempR;
		}
		if (resTempV > resV) {
			resV = resTempV;
		}
		//std::cout << resTempS << " " << resTempR << " " << resTempV << " ";
		// est égale à la puissance échangé, petit problème quelque part...
		double pSum = 0;
		int M = nVoisin.get(i, 0);
		for (int m = 0; m < M; m++) {
			pSum += Xhat[i].get(1 + m, 0);
			tempM1[i].set(m + 1, 0, Xhat[i].get(1 + m, 0) + Xhat[i].get(1 + m + M, 0));

			int lin = CoresAgentLin.get(i, 0) + m; // indice global de tim = tip (m est le numéro du voisin, p est le numéro de l'agent)
			int p = CoresLinVoisin.get(lin, 0); // valeur de p
			int lin2 = CoresLinTrans.get(lin, 0); // indice global de tpi
			int linLoc = lin2 - CoresAgentLin.get(p, 0); //indice local de tpi
			float tpi = Xhat[p].get(linLoc + 1, 0);
			tempM[i].set(m, 0, Xhat[i].get(1 + m + M, 0) - tpi); // 
		}
		tempM1[i].set(0, 0, Xhat[i].get(0, 0) - pSum);
		resTempR = MAX(tempM1[i].max2(), abs(Xhat[_nAgent].get(i, 0) - Xhat[i].get(0, 0))); // Bx = 0
		resTempV = tempM[i].max2(); // Gx = 0
		resTempS = X[i].max2(&Xpre[i]); // x^{k+1} = x^k
		//std::cout << resTempS << " " << resTempR << " " << resTempV;

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
	resS = MAX(resS, Xhat[_nAgent].max2(&X[_nAgent]));
	resS = MAX(resS, Xpre[_nAgent].max2(&X[_nAgent]));
	resR = MAX(resR, deltahat.max2(&delta));
	resV = MAX(resV, Nuhat[_nAgent].max2(&Nu[_nAgent]));
	//std::cout << resS << " " << resR << " " << resV << std::endl;

	resF.set(0, indice, resR);
	resF.set(1, indice, resS);
	resF.set(2, indice, resV);
	return MAX(MAX(resV, resS), resR);
}




void PACConst::display() {

	std::cout << _name << std::endl;
}
