#ifdef OSQP
#include "../head/OSQP.h"
 

OSQP::OSQP() : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " OSQP Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_alpha = 1;
	_name = NAME;
}


OSQP::OSQP(float alpha) : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default OSQP Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	
	_alpha = alpha;
	_name = NAME;
}

OSQP::~OSQP()
{
	for (int agent = 0; agent < _nAgent; agent++) {
		if (work[agent]) {
			osqp_cleanup(work[agent]);
			work[agent] = nullptr;
		}
#ifdef OSQPGPU
		DELETEB(Acsv[agent]);
		DELETEB(Pcsv[agent]);
		DELETEB(Pu[agent]);
#endif // OSQPGPU

		DELETEA(l[agent]);
		DELETEA(u[agent]);
		DELETEA(Q[agent]);
		//DELETEA(xResult[agent]); delete by work
		DELETEA(Pdata[agent]);
		DELETEA(Pidx[agent]);
		DELETEA(Pptr[agent]);
		DELETEA(Adata[agent]);
		DELETEA(Aidx[agent]);
		DELETEA(Aptr[agent]);
	}
	DELETEA(work);
#ifdef OSQPGPU
	DELETEA(Acsv);
	DELETEA(Pcsv);
	DELETEA(Pu);
#endif // OSQPGPU
	DELETEB(settings);
}
void OSQP::setParam(float alpha)
{
	if (alpha > 0 && alpha < 2) {
		_alpha = alpha;
	}
	std::cout << " Warning alpha not changed, must be between 0 and 2 and alpha = " << alpha << std::endl;
	
}

void OSQP::setTau(float tau)
{
	throw std::domain_error("tau is not define for this method");
}



void OSQP::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 0, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION
	}
	
	int iterG = sim.getIterG();
	int iterL = sim.getIterL();
	int stepG = sim.getStepG();

	float epsG = sim.getEpsG();
	float epsL = sim.getEpsL();
		

	float resG = 2 * epsG;
	float resL = 2 * epsL;
	float fc = 0;

	
	int iterGlobal = 0;
	int iterLocal = 0;
	c_float* xResult = nullptr;
	MatrixCPU resF(sim.getRes());
	//MatrixCPU timeAgent(_nAgent, iterG);
	MatrixCPU Pn(sim.getPn()); // trades' sum for each agent
	MatrixCPU tempN1(_nAgent, 1);
	MatrixCPU tempNN(_nAgent, _nAgentTrue);
	MatrixCPU tempM(_nTrade, 1);
	
	while ((iterGlobal < iterG) && (resG>epsG)) {
		resL = 2 * epsL;
		iterLocal = 0;
		
		for (int agent = 0; agent < _nAgent; agent++) {
			// FB 3b
#ifdef INSTRUMENTATION
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			updateQ(Q[agent], agent);
			osqp_update_lin_cost(work[agent], Q[agent]);;
#ifdef INSTRUMENTATION
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 5, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			// FB 1
			osqp_solve(work[agent]);
			xResult = work[agent]->solution->x;
			iterLocal += work[agent]->info->iter;
			updateTrade(xResult, agent);
#ifdef INSTRUMENTATION
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 1, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		}
		
		trade.swap(&Tlocal);
		// FB 3
#ifdef INSTRUMENTATION
		occurencePerBlock.increment(0, 1, 1);
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		updateLAMBDA();
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 5, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		// FB 4
		if (!(iterGlobal % stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateResBis(&resF, (iterGlobal / stepG), &tempM);
#ifdef INSTRUMENTATION
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		}
		
		iterGlobal++;
	}
#ifdef INSTRUMENTATION	
	occurencePerBlock.increment(0, 5, iterGlobal);
	occurencePerBlock.increment(0, 6, iterGlobal / stepG);


	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	// FB 5
	result->setResF(&resF);
	result->setLAMBDA(&LAMBDA);
	result->setTrade(&trade);
	result->setIter(iterGlobal);
	

	updatePn();
	result->setPn(&Pn);
	fc = calcFc();
	result->setFc(fc);
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 7, 1);

	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION

	tall = clock() - tall;
	result->setTime((float)tall / CLOCKS_PER_SEC);
}

float OSQP::updateResBis(MatrixCPU* res, int iter, MatrixCPU* tempM)
{
	int idTrade = 0;
	for (int agent = 0; agent < _nAgent; agent++) {
		int Mn = nVoisin.get(agent, 0);
		int begin = CoresAgentLin.get(agent, 0);
		for (int m = 0; m < Mn; m++) {
			int voisin = CoresLinVoisin.get(begin + m, 0);
			int agent2 = agent - _nAgentTrue * (agent >= _nAgentTrue);
			int voisin2 = voisin - _nAgentTrue * (voisin >=_nAgentTrue);
			
			tempM->set(idTrade, 0, trade.get(agent, voisin2) + trade.get(voisin, agent2));
			idTrade++;
		}

	}
	float resR = tempM->max2();


	float resS = _rho * Tlocal.max2(&trade);
	

	res->set(0, iter, resR);
	res->set(1, iter, resS);

	return MYMAX(resS, resR);
}
void OSQP::updateTrade(c_float* xResult, int agent) {

	int Mn = nVoisin.get(agent, 0);
	int begin = CoresAgentLin.get(agent, 0);

	for (int m = 0; m < Mn; m++) {
		int voisin = CoresLinVoisin.get(begin + m, 0);
		int voisin2 = voisin - _nAgentTrue * (voisin >= _nAgentTrue);

		Tlocal.set(agent, voisin2, xResult[m + 1]);
	}
	
}

void OSQP::init(const Simparam& sim, const StudyCase& cas)
{
	// d�truire avant de cr�er pour �viter les fuites m�moires
	if (settings) {
		for (int agent = 0; agent < _nAgent; agent++) {
			if (work[agent]) {
				osqp_cleanup(work[agent]);
				work[agent] = nullptr;
			}
#ifdef OSQPGPU
			DELETEB(Acsv[agent]);
			DELETEB(Pcsv[agent]);
			DELETEB(Pu[agent]);
#endif // OSQPGPU

			DELETEA(l[agent]);
			DELETEA(u[agent]);
			DELETEA(Q[agent]);
			//DELETEA(xResult[agent]); deleted by work
			DELETEA(Pdata[agent]);
			DELETEA(Pidx[agent]);
			DELETEA(Pptr[agent]);
			DELETEA(Adata[agent]);
			DELETEA(Aidx[agent]);
			DELETEA(Aptr[agent]);
		}
		DELETEA(work);
#ifdef OSQPGPU
		DELETEA(Acsv);
		DELETEA(Pcsv);
		DELETEA(Pu);
#endif // OSQPGPU
		DELETEB(settings);
	}



	//std::cout << "init" << std::endl;
	_rho = sim.getRho();
	isAC = cas.isAC();
	int iterL = sim.getIterL();
	float epsL = sim.getEpsL();
	
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
			std::cout << "err OSQP : " << _nAgent << " " << _nAgentTrue << " " << _nTrade << " " << _nTradeP << " " << _nTradeQ << std::endl;
			throw std::invalid_argument("Agent must be fully conected for the Q echanges, WIP");
		}
	}
	else {
		_nTradeP = _nTrade;
	}


	GAMMA = cas.getBeta();
	
	a = cas.geta();
	b = cas.getb();
	Ub = cas.getUb();
	MatrixCPU Lb(cas.getLb());
	MatrixCPU Pmin(cas.getPmin());
	MatrixCPU Pmax(cas.getPmax());
	

	LAMBDA = sim.getLambda();
	trade = sim.getTrade();
	Tlocal = sim.getTrade();
	MatrixCPU Pn(sim.getPn());


	CoresAgentLin = MatrixCPU(_nAgent + 1, 1);
	CoresLinVoisin = MatrixCPU(_nTrade, 1);
	

	// Who is the peer ? 
	int indice = 0;
	//std::cout << " P " << std::endl;
	for (int idAgent = 0; idAgent < _nAgentTrue; idAgent++) { // P
		MatrixCPU omega(cas.getVoisin(idAgent));
		int Nvoisinmax = nVoisin.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			int idVoisin = omega.get(voisin, 0);
			CoresLinVoisin.set(indice, 0, idVoisin); // 
			indice = indice + 1;
		}
		CoresAgentLin.set(idAgent + 1, 0, indice);
	}
	//std::cout << " Q " << std::endl;
	for (int idAgent = _nAgentTrue; idAgent < _nAgent; idAgent++) { // Q
		for (int idVoisin = 0; idVoisin < _nAgentTrue; idVoisin++) {
			if (idVoisin != (idAgent - _nAgentTrue)) {
				CoresLinVoisin.set(indice, 0, idVoisin + _nAgentTrue);
				indice = indice + 1;
			}
		}
		CoresAgentLin.set(idAgent + 1, 0, indice);
	}
	
	// OSQP object

	c_int exitflag = 0;// Exitflag

	
	// Workspace structures
	settings = new OSQPSettings;
#ifdef OSQPGPU
	work = new OSQPSolver * [_nAgent];
	Pcsv = new csc * [_nAgent];
	Acsv = new csc * [_nAgent];
	Pu = new csc * [_nAgent];
#else
	data = new OSQPData * [_nAgent];
	work = new OSQPWorkspace * [_nAgent];
#endif
	if (settings) {
		osqp_set_default_settings(settings);
		settings->alpha = _alpha;
		settings->eps_abs = epsL;
		settings->eps_rel = 0;
		settings->verbose = 0;
		settings->max_iter = iterL;
		//settings->linsys_solver = QDLDL_SOLVER;
	}
	
	
	Q = new c_float * [_nAgent];
	u = new c_float * [_nAgent];
	l = new c_float * [_nAgent];
	xResult = new c_float * [_nAgent];

	Pdata = new c_float * [_nAgent];
	Pidx = new c_int * [_nAgent];
	Pptr = new c_int * [_nAgent];

	Adata = new c_float * [_nAgent];
	Aidx = new c_int * [_nAgent];
	Aptr = new c_int * [_nAgent];

	for (int agent = 0; agent < _nAgent; agent++) {
		
		int n = nVoisin.get(agent, 0);
		float aAgent = a.get(agent, 0);
		int begin = CoresAgentLin.get(agent, 0);
#ifdef OSQPGPU
		Pcsv[agent] = new csc();
		Acsv[agent] = new csc();
		Pu[agent] = new csc();
#endif // OSQPGPU

		Q[agent] = new c_float[n + 1];
		u[agent] = new c_float[n + 2];
		l[agent] = new c_float[n + 2];
		xResult[agent] = new c_float[n + 1];

		l[agent][0] = 0; // sum(t) - p > 0
		u[agent][0] = 0; // sum(t) - p < 0
		l[agent][1] = Pmin.get(agent, 0); // p>pmin
		u[agent][1] = Pmax.get(agent, 0); // p<pmax 
		Q[agent][0] = (c_float) b.get(agent, 0);
		xResult[agent][0] = (c_float) Pn.get(agent, 0);

	
		for (int indV = 0; indV < n; indV++) {
			
			int voisin = CoresLinVoisin.get(begin + indV, 0);
			int agent2 = agent - _nAgentTrue * (agent >= _nAgentTrue);
			int voisin2 = voisin - _nAgentTrue * (voisin >= _nAgentTrue);
			
			c_float q = (c_float) LAMBDA.get(agent, voisin2) - _rho * (trade.get(agent, voisin2) - trade.get(voisin, agent2)) / 2.0;
			
			if (agent < _nAgentTrue) {
				q += (c_float) GAMMA.get(agent, voisin);
			}
			
			Q[agent][indV + 1] = q;
			l[agent][indV + 2] = Lb.get(agent, 0);
			u[agent][indV + 2] = Ub.get(agent, 0);
			xResult[agent][indV + 1] = trade.get(agent, voisin2);
		}


		MatrixCPU A(n + 2, n + 1);
		MatrixCPU P(n + 1, n + 1);
		
		P.setEyes(_rho);
		P.set(0, 0, aAgent);
		
		MatrixCPU diag(n + 1, n + 1);
		MatrixCPU ones(1, n + 1, 1);
		diag.setEyes(1);
		A.setBloc(0, 1, 0, n + 1, &ones);
		A.setBloc(1, n + 2, 0, n + 1, &diag);

		A.set(0, 0, -1);

		c_int P_nnz = P.getNNullHalf();
		int Pcol = P.getNCol();
		Pdata[agent] = new c_float[P_nnz];
		Pidx[agent] = new c_int[P_nnz];
		Pptr[agent] = new c_int[Pcol + 1];
		P.toCSCHalf(Pdata[agent], Pidx[agent], Pptr[agent]);


		c_int A_nnz = A.getNNull();
		int Acol = A.getNCol();
		Adata[agent] = new c_float[A_nnz];
		Aidx[agent] = new c_int[A_nnz];
		Aptr[agent] = new c_int[Acol + 1];
		A.toCSC(Adata[agent], Aidx[agent], Aptr[agent]);

#ifdef OSQPGPU
		/* Populate matrices */
		csc_set_data(Acsv[agent], n + 2, n + 1, A_nnz, Adata[agent], Aidx[agent], Aptr[agent]);
		csc_set_data(Pcsv[agent], n + 1, n + 1, P_nnz, Pdata[agent], Pidx[agent], Pptr[agent]);

		Pu[agent] = csc_to_triu(Pcsv[agent]);
		exitflag = osqp_setup(&work[agent], Pu[agent], Q[agent], Acsv[agent], l[agent], u[agent], n + 2, n + 1, settings);
		osqp_warm_start(work[agent], xResult[agent], NULL);
#else
		/* old OSQP*/
		
		data[agent] = new OSQPData;

		if (data[agent]) {
			data[agent]->n = n + 1;
			data[agent]->m = n + 2;
			data[agent]->P = csc_matrix(data[agent]->n, data[agent]->n, P_nnz, Pdata[agent], Pidx[agent], Pptr[agent]);
			data[agent]->q = Q[agent];
			data[agent]->A = csc_matrix(data[agent]->m, data[agent]->n, A_nnz, Adata[agent], Aidx[agent], Aptr[agent]);
			data[agent]->l = l[agent];
			data[agent]->u = u[agent];
		}
		exitflag = osqp_setup(&work[agent], data[agent], settings);

		osqp_warm_start_x(work[agent], xResult[agent]);
#endif // OSQPGPU
	}
	//std::cout << "fin init" << std::endl;
	//updateLAMBDA()
}

void OSQP::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	std::cout << "pas � jour si puissance r�active" << std::endl;
	// change  : ub, Pmin, Pmax, b pour les consomateurs
	c_int exitflag = 0;
	MatrixCPU Pmin(cas.getPmin());
	MatrixCPU Pmax(cas.getPmax());

	int nCons = cas.getNCons();
	MatrixGPU Lb(cas.getLb());

	b = cas.getb();

	for (int agent = 0; agent < nCons; agent++) {
		int n = nVoisin.get(agent, 0);
		
		l[agent][0] = 0;
		u[agent][0] = 0;
		l[agent][1] = Pmin.get(agent, 0);
		u[agent][1] = Pmax.get(agent, 0);
		MatrixCPU omega = cas.getVoisin(agent);

		for (int m = 0; m < n; m++) {
			// Q automatiquement mis � jour dans solve
			// int voisin = omega.get(m, 0);
			// c_float q = (c_float)(b.get(agent, 0) + GAMMA.get(agent, voisin) + LAMBDA.get(agent, voisin) - rho * (trade.get(agent, voisin) - trade.get(voisin, agent)) / 2);
			// Q[m] = q;
			l[agent][m + 2] = Lb.get(agent, 0);
			u[agent][m + 2] = Ub.get(agent, 0);
		}

		exitflag = osqp_update_bounds(work[agent], l[agent], u[agent]);

	}


#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION

}


void OSQP::updateQ(c_float* Q, int agent) {
	
	Q[0] = (c_float) b.get(agent, 0);
	int n = nVoisin.get(agent, 0);
	int begin = CoresAgentLin.get(agent, 0);
	for (int indV = 0; indV < n; indV++) {
		int voisin = CoresLinVoisin.get(begin + indV, 0);
		int agent2 = agent - _nAgentTrue * (agent >= _nAgentTrue);
		int voisin2 = voisin - _nAgentTrue * (voisin >= _nAgentTrue);
		c_float q = (c_float)(LAMBDA.get(agent, voisin2) - _rho * (trade.get(agent, voisin2) - trade.get(voisin, agent2)) / 2);

		if (agent < _nAgentTrue) {
			q += (c_float)GAMMA.get(agent, voisin);
		}

		Q[indV + 1] = q;
	}
}

void OSQP::updateLAMBDA()
{
	for (int agent = 0; agent < _nAgent; agent++) {
		int Mn = nVoisin.get(agent, 0);
		int begin = CoresAgentLin.get(agent, 0);
		for (int m = 0; m < Mn; m++) {
			int voisin = CoresLinVoisin.get(begin + m, 0);
			int agent2 = agent - _nAgentTrue * (agent >= _nAgentTrue);
			int voisin2 = voisin - _nAgentTrue * (voisin >= _nAgentTrue);
			float lambda = LAMBDA.get(agent, voisin2) + 0.5 * _rho * (trade.get(agent, voisin2) + trade.get(voisin, agent2));
			LAMBDA.set(agent, voisin2, lambda);
		}

	}
}




void OSQP::display() {

	std::cout << _name << std::endl;
}
#endif