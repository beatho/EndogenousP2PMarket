#include "../head/DCOPFOSQP.h"


DCOPFOSQP::DCOPFOSQP() : Method()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "Constructeur DCOPFOSQP" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
}


DCOPFOSQP::~DCOPFOSQP()
{
#ifndef OSQPGPU
	
	if (work) {
		osqp_cleanup(work);
	}
if (data) {
		if (data->A) c_free(data->A);
		if (data->P) c_free(data->P);
		c_free(data);
	}

#endif
	
	
	if (settings) c_free(settings);
	
	DELETEA(L);
	
	DELETEA(U);
	DELETEA(Q);
	//DELETEA(xResult); deja supprime par work 
}



void DCOPFOSQP::init(const Simparam& sim, const StudyCase& cas)
{
	// intitilisation des matrixs et variables 
	_rhog = sim.getRho();
	
	
	//std::cout << "rho initial " << _rhog << std::endl;
	_nAgent = sim.getNAgent();
	_nLine = cas.getNLine();
	_nBus = cas.getNBus();
	a = cas.geta();
	b = cas.getb();
	Pmin = cas.getPmin();
	Pmax = cas.getPmax();
	Pn = sim.getPn();
	lLimit = cas.getLineLimit();
	tempN1 = MatrixCPU(_nAgent, 1);
	nConst = _nLine + _nAgent + 1;

	settings = new OSQPSettings;
#ifndef OSQPGPU
	work = new OSQPWorkspace;
	data = new OSQPData;
#endif
	if (settings) {
		osqp_set_default_settings(settings);
		settings->alpha = 1;
		settings->eps_abs = sim.getEpsG();
		settings->verbose = 0;
		settings->max_iter = sim.getIterG();
		settings->adaptive_rho_interval = 0;
	}
	
	
	Q = new c_float[_nAgent];
		
	c_float* U = new c_float[nConst];
	c_float* L = new c_float[nConst];
	xResult = new c_float[_nAgent];
	for (int l = 0; l < _nLine; l++) {

		L[l] = -lLimit.get(l, 0);
		U[l] =  lLimit.get(l, 0);
			
	}
	for (int n = 0; n < _nAgent; n++) {
		L[n + _nLine] = Pmin.get(n, 0);
		U[n + _nLine] = Pmax.get(n, 0);
		Q[n] = b.get(n, 0);
		xResult[n] = Pn.get(n, 0);
	}
	L[nConst - 1] = 0;
	U[nConst - 1] = 0;
	

	Aosqp = MatrixCPU(nConst, _nAgent);
	G = cas.getPowerSensi();
	Aosqp.setBloc(0, _nLine, 0, _nAgent, &G);
	MatrixCPU iden(_nAgent, _nAgent);
	iden.setEyes(1);
	MatrixCPU ones(1, _nAgent, 1);
	Aosqp.setBloc(_nLine, _nLine + _nAgent, 0, _nAgent, &iden);
	Aosqp.setBloc(nConst-1, nConst, 0, _nAgent, &ones);


	Hosqp = MatrixCPU(_nAgent, _nAgent);
	Hosqp.setEyes(&a);
			
	c_int H_nnz = Hosqp.getNNullHalf();
	int Hcol = Hosqp.getNCol();
			
	c_float* Hdata = new c_float[H_nnz];
	c_int* Hidx = new c_int[H_nnz];
	c_int* Hptr = new c_int[Hcol + 1];
	Hosqp.toCSCHalf(Hdata, Hidx, Hptr);

	c_int A_nnz = Aosqp.getNNull();
	int Acol = Aosqp.getNCol();
	c_float* Adata = new c_float[A_nnz];
	c_int* Aidx = new c_int[A_nnz];
	c_int* Aptr = new c_int[Acol + 1];
	Aosqp.toCSC(Adata, Aidx, Aptr);
	
#ifndef OSQPGPU



	data = new OSQPData;
	if (data) {
		data->n = _nAgent;
		data->m = nConst;
		data->P = csc_matrix(data->n, data->n, H_nnz, Hdata, Hidx, Hptr);
		data->q = Q;
		data->A = csc_matrix(data->m, data->n, A_nnz, Adata, Aidx, Aptr);
		data->l = L;
		data->u = U;
	}
	c_int exitflag = osqp_setup(&work, data, settings);

	osqp_warm_start_x(work, xResult);
	

	DELETEA(Hdata);
	DELETEA(Hidx);
	DELETEA(Hptr);
	DELETEA(Adata);
	DELETEA(Aidx);
	DELETEA(Aptr);
#endif // !OSQPGPU
	
}



void DCOPFOSQP::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	
	
	Pmin = cas.getPmin();
	Pmax = cas.getPmax();

	b = cas.getb();
	
	for (int n = 0; n < _nAgent; n++) {
		Q[n] = b.get(n, 0);
		L[_nLine + n] = Pmin.get(n, 0);
		U[_nLine + n] = Pmax.get(n, 0);
	}

#ifndef OSQPGPU
	osqp_update_lin_cost(work, Q);
	osqp_update_bounds(work, L, U);
#endif
	
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 10, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 10, 1);
#endif // INSTRUMENTATION

	//std::cout << "fin update temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
}



void DCOPFOSQP::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
	
#ifdef INSTRUMENTATION
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	int iterGlobal = 0;
#ifndef OSQPGPU
	osqp_solve(work);
	xResult = work->solution->x;
	iterGlobal = work->info->iter;
#endif
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 1, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 1, iterGlobal);
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	for (int n = 0; n < _nAgent; n++) {
		Pn.set(n, 0, xResult[n]);
	}

	float fc = calcFc(&a, &b, &Pn, &tempN1);
	
	
	result->setIter(iterGlobal);
	result->setPn(&Pn);
	result->setFc(fc);
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 9, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 9, 1);
	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION

	tall = clock() - tall;
	result->setTime((float)tall / CLOCKS_PER_SEC);
}



void DCOPFOSQP::display() {

	std::cout << _name << std::endl;
}


