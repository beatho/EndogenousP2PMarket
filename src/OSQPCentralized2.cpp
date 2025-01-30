
#ifdef OSQP
#include "../head/OSQPCentralized2.h"


OSQPCentralized2::OSQPCentralized2() : Method()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " OSQPCentralized2 Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_alpha = 1;
	_name = NAME;
	//def comme P2P m�me si cela n'en est pas
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1abcd, Fb2, Fb3, Fb5, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu
}


OSQPCentralized2::OSQPCentralized2(float alpha) : Method()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default OSQPCentralized2 Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	
	_alpha = alpha;
	_name = NAME;
}

OSQPCentralized2::~OSQPCentralized2()
{

	
	//if (settings) c_free(settings);
	DELETEB(settings);
	if (work) {
		osqp_cleanup(work);
		work = nullptr;
	}
	DELETEA(l);
	DELETEA(u);
	DELETEA(Q);
#ifdef OSQPGPU
	DELETEB(Acsv);
	DELETEB(Pcsv);
	DELETEB(Pu);
#else
	if (data) {
		if (data->A) c_free(data->A);
		if (data->P) c_free(data->P);
		DELETEB(data);
	}
#endif // OSQPGPU


	DELETEA(Pdata);
	DELETEA(Pidx);
	DELETEA(Pptr);

	DELETEA(Adata);
	DELETEA(Aidx);
	DELETEA(Aptr);
	//DELETEA(xResult); deja supprime par work
}
void OSQPCentralized2::setParam(float alpha)
{
	_alpha = alpha;
}

void OSQPCentralized2::setTau(float tau)
{
	throw std::domain_error("tau is not define for this method");
}



void OSQPCentralized2::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
		timePerBlock.increment(0, 0, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION
	}
	MatrixCPU Pn(sim.getPn()); // trades' sum for each agent
	MatrixCPU tempN1(_nAgent, 1);
	MatrixCPU tempNN(_nAgent, _nAgent);
	
	float fc = 0;
#ifdef INSTRUMENTATION
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	osqp_solve(work);
	
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 5, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 5, 1);
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION


	

	// FB 5
	xResult = work->solution->x;
	updateTrade();
	int iterG = work->info->iter;

	result->setTrade(&trade);
	result->setIter(iterG);
	
	updatePn(&Pn, &trade);
	result->setPn(&Pn);
	fc = calcFc(&a, &b, &trade, &Pn, &GAMMA,&tempN1,&tempNN);
	result->setFc(fc);

	tall = clock() - tall;
	result->setTime((float)tall / CLOCKS_PER_SEC);
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 7, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 7, 1);

	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
}

void OSQPCentralized2::updateTrade() {

	for (int lin = 0; lin < _nTrade; lin++) {
		int i = CoresLinAgent.get(lin, 0);
		int j = CoresLinVoisin.get(lin, 0);
		trade.set(i, j, xResult[lin]);
	}
	
}

void OSQPCentralized2::init(const Simparam& sim, const StudyCase& cas)
{

	
	float epsG = sim.getEpsG();
	int iterG = sim.getIterG();
	isAC = cas.isAC();
	_nAgentTrue = sim.getNAgent();
	_nAgent = _nAgentTrue + isAC * _nAgentTrue;

	//std::cout << _nAgent << std::endl;
	GAMMA = cas.getBeta();
	a = cas.geta();
	b = cas.getb();
	Ub = cas.getUb();
	MatrixCPU Lb(cas.getLb());
	MatrixCPU Pmin(cas.getPmin());
	MatrixCPU Pmax(cas.getPmax());
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
	nConstraint = _nTrade / 2 + _nTrade + _nAgent; // asym�trie, limite des trades limites des puissances
	GAMMALin = MatrixCPU(_nTradeP, 1);

	trade = sim.getTrade();
	//Lb.display();
	//Ub.display();
	// OSQPCentralized2 object


	// Exitflag
	c_int exitflag = 0;

	// Workspace structures
	
	settings = new OSQPSettings; //(OSQPSettings*)c_malloc(sizeof(OSQPSettings));
#ifdef OSQPGPU
	work = (OSQPSolver*)c_malloc(sizeof(OSQPSolver));
#else
	work = new OSQPWorkspace;// (OSQPWorkspace*)c_malloc(sizeof(OSQPWorkspace));
	data = new OSQPData;// (OSQPData*)c_malloc(sizeof(OSQPData));
#endif // OSQPGPU

	if (settings) {
		osqp_set_default_settings(settings);
		settings->alpha = _alpha;
		settings->eps_abs = epsG;
		settings->eps_rel = 0;
		settings->verbose = 0;
		settings->max_iter = iterG;
		//settings->linsys_solver = QDLDL_SOLVER;
		// Cette m�thode ne marche pas settings->linsys_solver = MKL_PARDISO_SOLVER; // MKL_PARDISO_SOLVER, QDLDL_SOLVER, CUDA_PCG_SOLVER
	}

	
	c_int P_nnz = 0;
	c_int Pcol = _nTrade;
		
	Q = new c_float[_nTrade];
	l = new c_float[nConstraint];
	u = new c_float[nConstraint];
	xResult = new c_float[_nTrade];

	CoresMatLin = MatrixGPU(_nAgent, _nAgent, -1);
	CoresLinAgent = MatrixGPU(_nTrade, 1);
	CoresAgentLin = MatrixGPU(_nAgent + 1, 1);
	CoresLinVoisin = MatrixGPU(_nTrade, 1);
	CoresLinTrans = MatrixGPU(_nTrade, 1);


	int indice = 0;
	
	for (int idAgent = 0; idAgent < _nAgentTrue; idAgent++) {
		MatrixCPU omega(cas.getVoisin(idAgent));
		int Nvoisinmax = nVoisin.get(idAgent, 0);
		P_nnz += Nvoisinmax * Nvoisinmax + Nvoisinmax;
		//A.setBloc(nTradeMax + idAgent, nTradeMax + idAgent + 1, indice, indice + Nvoisinmax, &ones); // pmin<P<Pmax
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			// l & u & q & warmstart
			int idVoisin = omega.get(voisin, 0);
			GAMMALin.set(indice, 0, GAMMA.get(idAgent, idVoisin));
			l[indice] = (c_float)Lb.get(idAgent, 0);
			u[indice] = (c_float)Ub.get(idAgent, 0);
			Q[indice] = (c_float)GAMMALin.get(indice, 0) + b.get(idAgent, 0);
			xResult[indice] = (c_float)trade.get(idAgent, idVoisin);
			CoresLinAgent.set(indice, 0, idAgent);
			CoresLinVoisin.set(indice, 0, idVoisin);
			CoresMatLin.set(idAgent, idVoisin, indice);
			indice = indice + 1;
		}
		l[_nTradeP + idAgent] = (c_float)Pmin.get(idAgent, 0);
		u[_nTradeP + idAgent] = (c_float)Pmax.get(idAgent, 0);
		CoresAgentLin.set(idAgent + 1, 0, indice);
	}
	indice += _nAgentTrue; // pmin < p < pma
	//std::cout << " Q " << std::endl;
	for (int idAgent = _nAgentTrue; idAgent < _nAgent; idAgent++) { // Q
		P_nnz += (_nAgentTrue - 1) * (_nAgentTrue - 1) + (_nAgentTrue - 1);
		for (int idVoisin = 0; idVoisin < _nAgentTrue; idVoisin++) {
			if (idVoisin != (idAgent - _nAgentTrue)) {
				l[indice] = (c_float)Lb.get(idAgent, 0);
				u[indice] = (c_float)Ub.get(idAgent, 0);
				Q[indice] = (c_float) b.get(idAgent, 0);
				xResult[indice] = (c_float)trade.get(idAgent, idVoisin);
				CoresLinAgent.set(indice, 0, idAgent);
				CoresLinVoisin.set(indice, 0, idVoisin + _nAgentTrue);
				CoresMatLin.set(idAgent, idVoisin, indice);
				indice = indice + 1;
			}
		}
		CoresAgentLin.set(idAgent + 1, 0, indice);
	}
	if (isAC) {
		indice = indice + _nAgent;// qmin < q < qma
	}
	


	P_nnz = P_nnz / 2;
	
	Pdata = new c_float[P_nnz];
	Pidx = new c_int[P_nnz];
	Pptr = new c_int[Pcol + 1];

	c_int A_nnz = 3 * _nTrade; 
	c_int Acol = _nTrade;

	Adata = new c_float[A_nnz];
	Aidx = new c_int[A_nnz];
	Aptr = new c_int[Acol + 1];

	int indData = 0;
	int PindData = 0;
	int offset = 0;
	int previousAgent = 0;
	int Nvoisinmax = 0;
	int indice2 = indice;

	int* indTranspo = new int[_nTrade];
	for (int lin = 0; lin < _nTrade; lin++) {
		int i = CoresLinAgent.get(lin, 0);
		int n = i;
		int j = CoresLinVoisin.get(lin, 0);
		if (lin >= _nTradeP) {
			i -= _nAgentTrue;
		}

		int k = CoresMatLin.get(j, i);
		CoresLinTrans.set(lin, 0, k);
		
		if (n != previousAgent) {
			offset += Nvoisinmax;
			previousAgent = n;
		}
		
		// A : parcourir les rows
			// lin = column
		// tmin<t<tmax
		Adata[indData] = 1;
		Aidx[indData] = lin;
		Aptr[lin] = indData;
		indData++;

		// pmin<p<pmax
		
		Adata[indData] = 1;
		Aidx[indData] = _nTrade + n;
		indData++;
		
		
		if (k > lin) { // antisym�trie
			//std::cout << indice << std::endl;
			l[indice] = (c_float) 0;
			u[indice] = (c_float) 0;
			//std::cout << l[indice] << " " << u[indice] << std::endl;
			Adata[indData] = 1;
			Aidx[indData] = indice;
			indTranspo[k] = indice;
			indData++;
			indice++;
			indice2++;
		}
		else {
			Adata[indData] = 1;
			Aidx[indData] = indTranspo[lin];
			//Aidx[indData] = indice2-nTradeMax/2; n'est vrai que si pas de prosumers
			indData++;
			indice2++;
		}
		// P
		Nvoisinmax = nVoisin.get(n, 0);
		float an = a.get(n, 0);
		Pptr[lin] = PindData;
		//std::cout << lin << " " << n << " " << Nvoisinmax << " " << offset << std::endl;
		for (int i = 0; i < Nvoisinmax; i++) {
			if (i  <= (lin - offset)) {
				Pdata[PindData] = an;
				Pidx[PindData] = i + offset;
				PindData++;
			}
		}
	}
	
	DELETEA(indTranspo);
	Aptr[Acol] = A_nnz;
	Pptr[Pcol] = P_nnz;
	/*for (int i = 0; i < nTradeMax; i++) {
		std::cout << xResult[i] << " ";
	}
	std::cout << std::endl;
	for (int i = 0; i < Acol + 1; i++) {
		std::cout << Aptr[i] << " ";
	}
	std::cout << std::endl;
	for (int i = 0; i < Pcol + 1; i++) {
		std::cout << Pptr[i] << " ";
	}
	std::cout << std::endl;
	
	for (int i = 0; i < A_nnz; i++) {
		std::cout << Adata[i] << " ";
	}
	std::cout << std::endl;
	for (int i = 0; i < A_nnz; i++) {
		std::cout << Aidx[i] << " ";
	}
	std::cout << std::endl;
	for (int i = 0; i < nConstraint; i++) {
		std::cout << l[i] << " ";
	}
	std::cout << std::endl;
	for (int i = 0; i < nConstraint; i++) {
		std::cout << u[i] << " ";
	}
	std::cout << std::endl;
	std::cout << PindData << " " << indData << std::endl;
	std::cout << P_nnz << " " << A_nnz << std::endl;*/
#ifdef OSQPGPU
	
	csc_set_data(Pcsv, nTradeMax, nTradeMax, P_nnz, Pdata, Pidx, Pptr);
	Pu = csc_to_triu(Pcsv);
	csc_set_data(Acsv, nConstraint, nTradeMax, A_nnz, Adata, Aidx, Aptr);
	exitflag = osqp_setup(&work, Pu, Q, Acsv, l, u, nConstraint, nTradeMax, settings);
	osqp_warm_start(work, xResult, NULL);
#else
	/*for (int i = 0; i < nConstraint; i++) {
		std::cout << l[i] << " ";
	}
	std::cout << std::endl;
	for (int i = 0; i < nConstraint; i++) {
		std::cout << u[i] << " ";
	}*/
	if (data) {
		data->n = _nTrade;
		data->m = nConstraint;
		data->P = csc_matrix(data->n, data->n, P_nnz, Pdata, Pidx, Pptr);
		data->q = Q;
		data->A = csc_matrix(data->m, data->n, A_nnz, Adata, Aidx, Aptr);
		data->l = l;
		data->u = u;
	}
	
	std::cout << std::endl;
	exitflag = osqp_setup(&work, data, settings);

	osqp_warm_start_x(work, xResult);
#endif // OSQPGPU
	
	
}

void OSQPCentralized2::updatePn(MatrixCPU* Pn, MatrixCPU* trade)
{
	Pn->set(0.0);
	for (int n = 0; n < _nAgent; n++) {
		for (int m = 0; m < _nAgentTrue; m++) {
			Pn->increment(n, 0, trade->get(n, m));
		}
	}

}

float OSQPCentralized2::calcFc(MatrixCPU* cost1, MatrixCPU* cost2, MatrixCPU* trade, MatrixCPU* Pn, MatrixCPU* BETA, MatrixCPU* tempN1, MatrixCPU* tempNN)
{
	float fc = 0;
	tempN1->set(cost1);
	tempN1->multiply(0.5);
	tempN1->multiplyT(Pn);
	tempN1->add(cost2);
	tempN1->multiplyT(Pn);

	fc = fc + tempN1->sum();

	for (int n = 0; n < _nAgentTrue; n++) {
		for (int m = 0; m < _nAgentTrue; m++) {
			fc += trade->get(n, m) * BETA->get(n, m);
		}
	}


	return fc;
}

void OSQPCentralized2::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	// change  : Lb, Pmin, Pmax, b pour les consomateurs
	c_int exitflag = 0;
	MatrixCPU Pmin(cas.getPmin());
	MatrixCPU Pmax(cas.getPmax());

	int nCons = cas.getNCons();
	MatrixGPU Lb(cas.getLb());

	MatrixCPU b(cas.getb());

	int indice = 0;
	for (int agent = 0; agent < nCons; agent++) { // mise � jour des consommateurs
		int m = nVoisin.get(agent, 0);
		for (int peer = 0; peer < m; peer++) {
			Q[indice] = b.get(agent, 0) + GAMMALin.get(indice, 0);
			l[indice] = Lb.get(agent, 0);
			indice = indice + 1;
		}
		l[_nTradeP + agent] = Pmin.get(agent, 0);
		u[_nTradeP + agent] = Pmax.get(agent, 0);
	}
	indice += _nAgentTrue;
	for (int agent = _nAgentTrue; agent < _nAgent; agent++) { // mise � jour des puissances r�actives
		int m = _nAgentTrue - 1;
		for (int peer = 0; peer < m; peer++) {
			Q[indice] = b.get(agent, 0);
			l[indice] = Lb.get(agent, 0);
			u[indice] = Ub.get(agent, 0);
			indice = indice + 1;
		}
		l[_nTrade + _nAgentTrue + agent] = Pmin.get(agent, 0);
		u[_nTrade + _nAgentTrue + agent] = Pmax.get(agent, 0);
	}

	exitflag = osqp_update_bounds(work, l, u);
	osqp_update_lin_cost(work, Q);

	


#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION

}


void OSQPCentralized2::display() {

	std::cout << _name << std::endl;
}


#endif