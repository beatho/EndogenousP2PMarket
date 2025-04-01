#ifdef OSQP
#include "../head/OSQPConst.h"
 

OSQPConst::OSQPConst() : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " OSQPConst Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_alpha = 1;
	_name = NAME;
}


OSQPConst::OSQPConst(float alpha) : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default OSQPConst Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	
	_alpha = alpha;
	_name = NAME;
}

OSQPConst::~OSQPConst()
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
		DELETEA(IdVoisin[agent]);
	}
	DELETEA(work);
#ifdef OSQPGPU
	DELETEA(Acsv);
	DELETEA(Pcsv);
	DELETEA(Pu);
#endif // OSQPGPU
	DELETEB(settings);
}
void OSQPConst::setParam(float alpha)
{
	_alpha = alpha;
}

void OSQPConst::setTau(float tau)
{
	throw std::domain_error("tau is not define for this method");
}



void OSQPConst::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
	int nAgent = sim.getNAgent();
	

	float resG = 2 * epsG;
	float resL = 2 * epsL;
	float fc = 0;

	
	int iterGlobal = 0;
	int iterLocal = 0;
	MatrixCPU resF(sim.getRes());
	//MatrixCPU timeAgent(_nAgent, iterG);
	
	MatrixCPU tempN1(_nAgent, 1);
	MatrixCPU tempNN(_nAgent, _nAgent);
	
	while ((iterGlobal < iterG) && (resG>epsG)) {
		resL = 2 * epsL;
		// FB 1
		iterLocal = 0;
		
		for (int agent = 0; agent < nAgent; agent++) {
			// FB 1a
#ifdef INSTRUMENTATION
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION		
			osqp_update_lin_cost(work[agent], Q[agent]);;
#ifdef INSTRUMENTATION
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 1, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			// FB 1b
			osqp_solve(work[agent]);
			xResult[agent] = work[agent]->solution->x;
			int M = nVoisin.get(agent, 0);
			for (int m = 0; m < M; m++) {
				int voisin = IdVoisin[agent][m];
				Tlocal.set(agent, voisin, xResult[agent][m+1]);
			}
			iterLocal += work[agent]->info->iter;
			//updateTrade(xResult, agent, &omega);
#ifdef INSTRUMENTATION
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 2, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		}
		// FB 3
#ifdef INSTRUMENTATION
		occurencePerBlock.increment(0, 1, iterLocal);
		occurencePerBlock.increment(0, 2, iterLocal);
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		updateLAMBDA();
		updatePn();
		updatePhi();
		Kappa1_pre.set(&Kappa1);
		Kappa2_pre.set(&Kappa2);
		updateKappa(&Kappa1, &Kappa2, &lLimit, &Phi);
		updateQ();
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 5, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		// FB 4
		if (!(iterGlobal % stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateResBis(&resF, (iterGlobal / stepG));
#ifdef INSTRUMENTATION
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		}
		trade.swap(&Tlocal);
		iterGlobal++;
	}
#ifdef INSTRUMENTATION	
	occurencePerBlock.increment(0, 5, iterGlobal);
	occurencePerBlock.increment(0, 6, iterGlobal / stepG);


	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	// FB 5
	
	for (int n = 0; n < _nAgent; n++) {
		int M = nVoisin.get(n, 0);
		for (int m = 0; m < M; m++) {
			int voisin = IdVoisin[n][m];
			trade.set(n, voisin, xResult[n][m + 1]);
		}
	}


	Kappa1.projectNeg(); //delta1
	Kappa2.projectNeg(); // delta2
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

void OSQPConst::updateTrade(c_float* xResult, int agent, MatrixCPU* omega) {

	int indice = 1;
	int M = nVoisin.get(agent, 0);
	for (int m = 0; m < M;m++) {
		Tlocal.set(agent, omega->get(m, 0), xResult[indice]);
		indice++;
	}
	
}

void OSQPConst::init(const Simparam& sim, const StudyCase& cas)
{
	
	_rho = sim.getRho();
	_rho1 = sim.getRho1();
	
	int iterL = sim.getIterL();
	float epsL = sim.getEpsL();
	float epsGC = sim.getEpsGC();
	float epsG = sim.getEpsG();
	_ratioEps = epsG / epsGC;
	_nAgent = sim.getNAgent();
	_nLine = cas.getNLine();
	_nBus = cas.getNBus();



	Beta = cas.getBeta();
	Ct = Beta;
	
	MatrixCPU connect(cas.getC());
	a = cas.geta();
	b = cas.getb();
	Ub = cas.getUb();
	MatrixCPU Lb(cas.getLb());
	MatrixCPU Pmin(cas.getPmin());
	MatrixCPU Pmax(cas.getPmax());
	nVoisin = cas.getNvoi();

	LAMBDA = sim.getLambda();
	trade = sim.getTrade();
	Tlocal = sim.getTrade();
	Pn = sim.getPn();


	Kappa1 = MatrixCPU(_nLine, 1, 0);
	Kappa2 = MatrixCPU(_nLine, 1, 0);
	Kappa1_pre = MatrixCPU(_nLine, 1, 0);
	Kappa2_pre = MatrixCPU(_nLine, 1, 0);
	Phipart = MatrixCPU(_nLine, _nAgent, 0);
	Phi = MatrixCPU(_nLine, 1);
	Ggrid = cas.getPowerSensi();
	Ggrid2 = cas.getPowerSensi();
	Ggrid2.multiplyT(&Ggrid);
	lLimit = MatrixCPU(cas.getLineLimit());
	
	updatePhi();
	updateKappa(&Kappa1, &Kappa2, &lLimit, &Phi);
	// OSQPConst object


	// Exitflag
	c_int exitflag = 0;

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

	IdVoisin = new int* [_nAgent];
	PosVoisin = MatrixCPU(_nAgent, _nAgent, -1);

	for (int agent = 0; agent < _nAgent; agent++) {
		int n = nVoisin.get(agent, 0);
		float aAgent = a.get(agent, 0);
#ifdef OSQPGPU
		Pcsv[agent] = new csc();
		Acsv[agent] = new csc();
		Pu[agent] = new csc();
#endif // OSQPGPU

		Q[agent] = new c_float[n + 1];
		u[agent] = new c_float[n + 2];
		l[agent] = new c_float[n + 2];
		xResult[agent] = new c_float[n + 1];
		IdVoisin[agent] = new int[n];

		l[agent][0] = 0; // sum(t) - p > 0
		u[agent][0] = 0; // sum(t) - p < 0
		l[agent][1] = Pmin.get(agent, 0); // p>pmin
		u[agent][1] = Pmax.get(agent, 0); // p<pmax 
		Q[agent][0] = (c_float) b.get(agent, 0);
		xResult[agent][0] = (c_float) Pn.get(agent, 0);

		MatrixCPU omega = cas.getVoisin(agent);

		for (int l = 0; l < _nLine; l++) {
			Q[agent][0] += _rho1 * (abs(Kappa1.get(l, 0)) - abs(Kappa2.get(l, 0)) + 2 * Phipart.get(l, agent)) * Ggrid.get(l, agent);
		}
	

		for (int m = 0; m < n; m++) {
			int voisin = omega.get(m, 0);
			IdVoisin[agent][m] = voisin;
			c_float q = (c_float)(Beta.get(agent, voisin) + LAMBDA.get(agent, voisin) - _rho * (trade.get(agent, voisin) - trade.get(voisin, agent)) / 2);
			Q[agent][m + 1] = q;
			l[agent][m + 2] = Lb.get(agent, 0);
			u[agent][m + 2] = Ub.get(agent, 0);
			xResult[agent][m + 1] = trade.get(agent, voisin);
			PosVoisin.set(agent, voisin, m + 1);
			
		}
		
		MatrixCPU A(n + 2, n + 1);
		MatrixCPU P(n + 1, n + 1);
		
		P.setEyes(_rho);
		float sumG = 0;
		for (int l = 0; l < _nLine; l++) {
			sumG += Ggrid2.get(l, agent);
		}
		P.set(0, 0, aAgent + 2 * _rho1 * sumG);
		
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
	
	/*for (int agent = 0; agent < _nAgent; agent++) {
		int n = nVoisin.get(agent, 0);
		for (int m = 0; m < n; m++) {
			std::cout << IdVoisin[agent][m] << " ";
		}
		std::cout << std::endl;
	}
	PosVoisin.display();*/
	updateLAMBDA();
}

void OSQPConst::updatePn()
{
	for (int n = 0; n < _nAgent; n++) {
		Pn.set(n, 0, xResult[n][0]);
	}
}

void OSQPConst::updateLAMBDA()
{
	for (int n = 0; n < _nAgent; n++) {
		int M = nVoisin.get(n, 0);
		for (int indV = 0; indV < M; indV++) {
			int voisin = IdVoisin[n][indV];
			int indA = PosVoisin.get(voisin, n);
			float lamb = 0.5 * _rho * (xResult[n][indV + 1] + xResult[voisin][indA]);
			LAMBDA.increment(n, voisin, lamb);
		}
	}
}

float OSQPConst::updateResBis(MatrixCPU* res, int iter)
{
	float resR = 0;
	for (int n = 0; n < _nAgent; n++) {
		int M = nVoisin.get(n, 0);
		for (int indV = 0; indV < M; indV++) {
			int voisin = IdVoisin[n][indV];
			int indA = PosVoisin.get(voisin, n);
			float resTemp = abs((xResult[n][indV + 1] + xResult[voisin][indA]));
			if (resTemp > resR) {
				resR = resTemp;
			}
		}
	}
	


	float resS = trade.max2(&Tlocal);

	// version de l'article
	MatrixCPU tempL(Kappa1);
	MatrixCPU tempL2(Kappa2);
	Kappa1_pre.projectNeg();
	Kappa2_pre.projectNeg();
	tempL.projectNeg();
	tempL2.projectNeg();
	tempL.subtract(&Kappa1_pre);
	tempL2.subtract(&Kappa2_pre);
	tempL.multiplyT(&tempL);
	tempL2.multiplyT(&tempL2);
	tempL.add(&tempL2);

	float resXf = _ratioEps * sqrt(tempL.max2());
	res->set(0, iter, resR);
	res->set(1, iter, resS);
	res->set(2, iter, resXf);
	return MYMAX(MYMAX(resXf, resS), resR);
}

void OSQPConst::updatePhi()
{
	for (int l = 0; l < _nLine; l++) {
		float qt = 0;
		for (int n = _nAgent - 1; n >= 0; n--) {
			qt += Ggrid.get(l, n) * Pn.get(n, 0);
			if (n > 0) {
				Phipart.set(l, n - 1, qt);
			}
		}
		Phi.set(l, 0, qt);
	}
}

void OSQPConst::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	// change  : ub, Pmin, Pmax, b pour les consomateurs
	c_int exitflag = 0;
	MatrixCPU Pmin(cas.getPmin());
	MatrixCPU Pmax(cas.getPmax());

	int nCons = cas.getNCons();
	MatrixGPU Lb(cas.getLb());

	b = cas.getb();
	updateQ();
	for (int agent = 0; agent < nCons; agent++) {
		osqp_update_lin_cost(work[agent], Q[agent]);;
		int n = nVoisin.get(agent, 0);
		
		l[agent][0] = 0;
		u[agent][0] = 0;
		l[agent][1] = Pmin.get(agent, 0);
		u[agent][1] = Pmax.get(agent, 0);
		MatrixCPU omega = cas.getVoisin(agent);

		for (int m = 0; m < n; m++) {
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


void OSQPConst::updateQ() {
	
	for (int n = 0; n < _nAgent; n++) {
		Q[n][0] = (c_float) b.get(n, 0);
		int M = nVoisin.get(n, 0);
		for (int l = 0; l < _nLine; l++) {
			Q[n][0] += _rho1 * (abs(Kappa1.get(l, 0)) - abs(Kappa2.get(l, 0)) + 2 * Phipart.get(l, n)) * Ggrid.get(l, n);
		}
		for (int indV = 0; indV < M; indV++) {
			int voisin = IdVoisin[n][indV];
			int indA = PosVoisin.get(voisin, n);
			c_float q = (Beta.get(n, voisin) + LAMBDA.get(n, voisin) - _rho * (xResult[n][indV + 1] - xResult[voisin][indA]) / 2);
			Q[n][indV + 1] = q;
		}
	}
	
}



void OSQPConst::display() {

	std::cout << _name << std::endl;
}


#endif