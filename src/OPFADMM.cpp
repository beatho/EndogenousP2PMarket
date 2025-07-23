#include "../head/OPFADMM.h"
 
// old ? (Pi,Qi,vi,li,pi,qi,vai,pji,qji,lji)


OPFADMM::OPFADMM() : MethodOPF()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " OPFADMM Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
}

OPFADMM::OPFADMM(float rho) : MethodOPF(rho)
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default OPFADMM Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
}

void OPFADMM::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
{
#ifdef DEBUG_SOLVE
	cas.display();
	sim.display(1);
#endif // DEBUG_SOLVE
	
	timeOPF =clock();
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
	
		
	float fc = 0;
	float resG = 2 * _epsG;
	float resL = 2 * _epsL;
	_iterGlobal = 0;
	
	/*CoresSoloBusAgent.display();

	Apt1.display();
	Apt2.display();
	Bpt2.display();
	Cost1.display();
	Cost2.display();
	Pmin.display();
	Pmax.display();

	PnPre.display();
	PnMoy.display();
	PnTilde.display();
	MuL.display();
	_nAgentByBus.display();
	std::cout << _rhol << " " << epsL << " " << iterL << " " << _nAgent << " " << _nBus;
	std::cout << "------" << std::endl;*/

	while ((_iterGlobal < _iterG) && (resG>_epsG)) {
		/*std::cout << "--------" << std::endl;
		Pn.display();
		std::cout << " PnTilde " << std::endl;
		PnTilde.display();
		for (int i = 0; i < 2; i++) {
			std::cout << " X " << i << std::endl;
			X[i].display();
			std::cout << " Y " << i << std::endl;
			Y[i].display();
			std::cout << " Q " << i << std::endl;
			//Q[i].display();
			std::cout << " Mu " << i << std::endl;
			Mu[i].display();
		}
		Chat.display();
		Bpt2.display();*/
		
		_iterLocal = 0;
		resL = 2 * _epsL;
		while (_iterLocal< _iterL && resL> _epsL) {
			updateLocalProb();
			// FB 3
			if (!(_iterLocal % _stepL)) {
#ifdef INSTRUMENTATION
				t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
				resL = calcRes();
#ifdef INSTRUMENTATION
				t2 = std::chrono::high_resolution_clock::now();
				timePerBlock.increment(0, 4, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION			
			}
			_iterLocal++;
		}
		/*std::cout << " Pn " << std::endl;
		Pn.display();
		std::cout << " PnMoy " << std::endl;
		PnMoy.display();
		std::cout << " PnTilde " << std::endl;
		PnTilde.display();
		std::cout << " MuL " << std::endl;
		MuL.display();*/
		/*if (_iterLocal == iterL) {
			std::cout << _iterGlobal << " " << _iterLocal << " " << resL << " " << resG << std::endl;
		}*/
#ifdef INSTRUMENTATION
		occurencePerBlock.increment(0, 1, _iterLocal);
		occurencePerBlock.increment(0, 2, _iterLocal);
		occurencePerBlock.increment(0, 3, _iterLocal);
		occurencePerBlock.increment(0, 4, _iterLocal / stepL);
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

		updateXWOCurrent();
		updateXPQ();

#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 5, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		//std::cout << "---------" << std::endl;
		CommunicationX();
		
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		updateGlobalProb();
		updateMu();
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		updateChat();
		
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		// FB 4
		if (!(_iterGlobal % _stepG)) {
#ifdef INSTRUMENTATION
			t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
			resG = updateRes(_iterGlobal / _stepG);
			//std::cout << _iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, _iterGlobal / _stepG) << " " << resF.get(1, _iterGlobal / _stepG) << std::endl;
#ifdef INSTRUMENTATION
			t2 = std::chrono::high_resolution_clock::now();
			timePerBlock.increment(0, 9, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		}
		//std::cout << iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;

		_iterGlobal++;
	}
	//std::cout << iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, (iterGlobal - 1) / stepG) << " " << resF.get(1, (iterGlobal - 1) / stepG) << " " << resF.get(2, (iterGlobal - 1) / stepG) << std::endl;


#ifdef INSTRUMENTATION	
	occurencePerBlock.increment(0, 5, _iterGlobal);
	occurencePerBlock.increment(0, 6, _iterGlobal);
	occurencePerBlock.increment(0, 7, _iterGlobal);
	occurencePerBlock.increment(0, 8, _iterGlobal);
	occurencePerBlock.increment(0, 9, _iterGlobal / _stepG);

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	
	
	fc = calcFc();
	// FB 5
	
	Pn.set(0, 0, getPLoss());
	Pn.set(_nAgent, 0, getQLoss());
	result->setResF(&resF);
	

	/*std::cout << " Pn " << std::endl;
	Pn.display();
	std::cout << " PnTilde " << std::endl;
	PnTilde.display();*/
	
	MatrixCPU Pb(getPb());
	MatrixCPU Phi(getPhi());
	
	//Phi.display();
	MatrixCPU E(getE());
	
	result->setE(&E);
	result->setPhi(&Phi);
	result->setPb(&Pb);

	result->setIter(_iterGlobal);
	result->setPn(&Pn);
	result->setFc(fc);

#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 10, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 10, 1);

	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
	timeOPF = clock() - timeOPF;
	

	result->setTime((float) timeOPF / CLOCKS_PER_SEC);
	
}

void OPFADMM::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	Pmin = cas.getPmin();
	Pmax = cas.getPmax();
	Cost2 = cas.getb();

	// pour essayer que cela marche
	Pn.add(&Pmin, &Pmax);
	Pn.divide(2);
	PnTilde.set(0.0);

	for (int n = 1; n < _nAgent; n++) {
		int b = (int) _CoresBusAgent.get(n, 0);
		PnTilde.increment(b, 0, Pn.get(n, 0));
		PnTilde.increment(b + _nBus, 0, Pn.get(n + _nAgent, 0));
	}

	for (int b = 0; b < _nBus; b++) {
		int nb = (int) _nAgentByBus.get(b, 0);
		if (nb > 0) {
			PnTilde.set(b, 0, PnTilde.get(b, 0) / nb);
			PnTilde.set(b + _nBus, 0, PnTilde.get(b + _nBus, 0) / nb);
		}
		X[b].set(4, 0, nb * PnTilde.get(b, 0)); // pi = sum(pn)
		X[b].set(5, 0, nb * PnTilde.get(b + _nBus, 0)); // qi = sum(pn)
	}
	PnMoy.set(&PnTilde);
	DFSP(0); // Pi
	X[0].set(0, 0, 0);
	DFSQ(0); // Qi
	X[0].set(1, 0, 0);
	for (int i = 0; i < _nBus; i++) {
		//X[i].set(2, 0, 1);
		float Si = X[i].get(0, 0) * X[i].get(0, 0) + X[i].get(1, 0) * X[i].get(1, 0);
		X[i].set(3, 0, Si / X[i].get(2, 0)); // li = Si^2/vi
	}
	for (int i = 0; i < _nBus; i++) {
		int m = (int) nChild.get(i, 0);
		for (int j = 0; j < m; j++) {// (Pci, Qci, lci) for all child Ci
			int c = (int) Childs[i].get(j, 0);
			X[i].set(7 + 3 * j, 0, X[c].get(0, 0));
			X[i].set(8 + 3 * j, 0, X[c].get(1, 0));
			X[i].set(9 + 3 * j, 0, X[c].get(3, 0));
		}
		Y[i].set(&X[i]);
		//Mu[i].set(0.0);
	}/**/
	updateChat();
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 11, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 11, 1);
#endif // INSTRUMENTATION

}

void OPFADMM::init(const Simparam& sim, const StudyCase& cas)
{
	// intitilisation des matrixs et variables 
	//cas.display();
	//sim.display(1);
	clock_t t = clock();

	initSize(cas);
	
	initSimParam(sim);
	PnPre = sim.getPn();
	
	initLocalProb(cas);
	/*std::cout << "Apt1 " << std::endl;
	Apt1.display();
	std::cout << "Apt2 " << std::endl;
	Apt2.display();
	std::cout << " Pmin " << std::endl;
	Pmin.display();
	std::cout << " Pmax " << std::endl;
	Pmax.display();

	std::cout << " Pn  limits " << std::endl;
	PnTmin.display();
	PnTmax.display();*/
	
	
	initGrid(cas);
	
	
	sizeOPFADMM = MatrixCPU(_nBus, 1);
	for (int i = 0; i < _nBus; i++) {
		int sizeOPF =  3 * (int) nChild.get(i, 0) + 7;
		sizeOPFADMM.set(i, 0, sizeOPF);
		
	}
	
	allocateTab();
	
	for (int i = 0; i < _nBus; i++) {
		Chat[i] = MatrixCPU(4, 1);
		// vi = 1; vai = 1
		X[i].set(indvi, 0, 1); 
		X[i].set(indvai, 0, 1);

		// pi & qi
		int Nb = (int) _nAgentByBus.get(i, 0);
		int begin = (int) _CoresAgentBusBegin.get(i, 0);
		X[i].set(indpi, 0, Nb * PnTilde.get(i,0)); // pi = sum(pn)
		X[i].set(indqi, 0, Nb * PnTilde.get(i + _nBus, 0)); // qi = sum(pn)
		
		for (int In = 0; In < Nb; In++) {
			int n = (int) _CoresAgentBus.get(In + begin, 0);
			PosAgent.set(n, 0, In);
			Pb.increment(i, 0, Pn.get(n, 0));
			Pb.increment(i + _nBus, 0, Pn.get(n + _nAgent, 0));
			
			Pbmax.increment(i, 0, Pmax.get(n, 0));
			Pbmin.increment(i, 0, Pmin.get(n, 0));

			Pbmax.increment(i + _nBus, 0, Pmax.get(n + _nAgent, 0));
			Pbmin.increment(i + _nBus, 0, Pmin.get(n + _nAgent, 0));
		}
	
	}
	DFSP(0); // Pi
	X[0].set(indPi, 0, 0);
	DFSQ(0); // Qi
	X[0].set(indQi, 0, 0); 
	
	for (int i = 0; i < _nBus; i++) {
		float Si = X[i].get(indPi, 0) * X[i].get(indPi, 0) + X[i].get(indQi, 0) * X[i].get(indQi, 0);
		X[i].set(indli, 0, Si / X[i].get(indvi, 0)); // li = Si^2/vi
	}
	for (int i = 0; i < _nBus; i++) {
		int m = (int) nChild.get(i, 0);
		for (int j = 0; j < m; j++) { // (Pci, Qci, lci) for all child Ci
			int c = (int) Childs[i].get(j, 0);
			X[i].set(indPc + 3 * j, 0, X[c].get(indPi, 0));
			X[i].set(indQc + 3 * j, 0, X[c].get(indQi, 0));
			X[i].set(indlc + 3 * j, 0, X[c].get(indli, 0));
		}
		Y[i].set(&X[i]);
	}

	for (int i = 0; i < _nBus; i++) {
		if (i > 0) {
			A[i].set(2, indPi, 2 * ZsRe.get(i - 1, 0));
			A[i].set(2, indQi, 2 * ZsIm.get(i - 1, 0));
			A[i].set(2, indli, -ZsNorm.get(i - 1, 0));
			A[i].set(2, indvi, -1);
			A[i].set(2, indvai, 1);

			A[i].set(0, indPi, -1);
			A[i].set(1, indQi, -1);
		}
		A[i].set(0, indpi, 1);
		A[i].set(1, indqi, 1);
		for (int j = 0; j < nChild.get(i, 0); j++) {
			int c = (int) Childs[i].get(j, 0);
			A[i].set(0, indPc + 3 * j, 1); // Pci
			A[i].set(1, indQc + 3 * j, 1); // Qci
			A[i].set(0, indlc + 3 * j, -ZsRe.get(c - 1, 0)); // -R l
			A[i].set(1, indlc + 3 * j, -ZsIm.get(c - 1, 0)); // -X l
		}
		//A[i].display();
		MatrixCPU temp33(2 + 1 * (i > 0), 2 + 1 * (i > 0));
		MatrixCPU temp3M(2 + 1 * (i > 0), (int) sizeOPFADMM.get(i, 0));
		MatrixCPU tempMM((int)sizeOPFADMM.get(i, 0), (int) sizeOPFADMM.get(i, 0));
		temp33.multiplyTrans(&A[i], &A[i]); // (3*o_b) * (o_b*3) -> 9 * o_b
		temp33.invertGaussJordan(&temp33); // 3^3 = 27 (fixe !!!)
		temp3M.MultiplyMatMat(&temp33, &A[i]); // (3*3) * (3*o_b) -> 27 *o_b
		Hinv[i].multiplyTrans(&A[i], &temp3M, 0); // (o_b*3) * (3*o_b) -> 9 * o_b

		tempMM.setEyes(1);
		Hinv[i].subtract(&tempMM);
		Hinv[i].divide(_rho);
		//Hinv[i].display();
	}

	/*std::cout << "--------" << std::endl;
	std::cout << " Pn " << std::endl;
	Pn.display();
	PnMoy.display();
	std::cout << " PnTilde " << std::endl;
	PnTilde.display();
	for (int i = 0; i < 10; i++) {
		std::cout << " X " << i << std::endl;
		X[i].display();
		std::cout << " Y " << i << std::endl;
		Y[i].display();
		std::cout << " Q " << i << std::endl;
		Q[i].display();
		std::cout << " Mu " << i << std::endl;
		Mu[i].display();
	}*/
	std::cout << " updateChat " << std::endl;	
	updateChat();
	/*std::cout << " Chat " << std::endl;
	Chat.display();
	std::cout << " Bpt2 " << std::endl;
	Bpt2.display();
	std::cout << " Cp " << std::endl;
	Cost2.display();
	std::cout << " Ap2 " << std::endl;
	Cost1.display();
	std::cout << " Pmin " << std::endl;
	Pmin.display();
	std::cout << " Pma " << std::endl;
	Pmax.display();*/
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	//std::cout << "---------------------------------------------------------------------------------------" << std::endl;
}

void OPFADMM::initLocalProb(const StudyCase& cas){
	initAgent(cas);
	//std::cout << " local resolution " << std::endl;
	PnTmin = MatrixCPU(2 * _nBus, 1);
	PnTmax = MatrixCPU(2 * _nBus, 1);
	PnMoy = MatrixCPU(2 * _nBus, 1);
	MuL = MatrixCPU(2 * _nBus, 1);
	PnTilde = MatrixCPU(2 * _nBus, 1);
	Ap12 = Cost1;
	Ap12.add(_rhol);
	Bp1 = MatrixCPU(2 * _nAgent, 1);
	Bpt1 = MatrixCPU(2 * _nBus, 1);
	Bpt2 = MatrixCPU(2 * _nBus, 1);
	Apt1 = MatrixCPU(2 * _nBus, 1);
	Apt2 = MatrixCPU(2 * _nBus, 1);
	Apt12 = MatrixCPU(2 * _nBus, 1);
	for (int n = 1; n < _nAgent; n++) {
		int b = (int) _CoresBusAgent.get(n, 0);
		PnTilde.increment(b, 0, Pn.get(n, 0));
		PnTilde.increment(b + _nBus, 0, Pn.get(n + _nAgent, 0));
		
		PnTmax.increment(b, 0, Pmax.get(n, 0));
		PnTmax.increment(b + _nBus, 0, Pmax.get(n + _nAgent, 0));
		PnTmin.increment(b, 0, Pmin.get(n, 0));
		PnTmin.increment(b + _nBus, 0, Pmin.get(n + _nAgent, 0));

		if (_nAgentByBus.get(b, 0) == 1) {
			CoresSoloBusAgent.set(b, 0, n);
		}
	}
	
	for (int b = 0; b < _nBus; b++) {
		int nb = (int) _nAgentByBus.get(b, 0);
		if (nb > 0) {
			Apt1.set(b, 0, nb * _rhol);
			Apt1.set(b + _nBus, 0, nb * _rhol);

			Apt2.set(b, 0, nb * nb * _rho);
			Apt2.set(b + _nBus, 0, nb * nb * _rho);

			PnTilde.set(b, 0, PnTilde.get(b, 0) / nb);
			PnTilde.set(b + _nBus, 0, PnTilde.get(b + _nBus, 0) / nb);

			PnTmax.set(b, 0, PnTmax.get(b, 0) / nb);
			PnTmax.set(b + _nBus, 0, PnTmax.get(b + _nBus, 0) / nb);

			PnTmin.set(b, 0, PnTmin.get(b, 0) / nb);
			PnTmin.set(b + _nBus, 0, PnTmin.get(b + _nBus, 0) / nb);
		}
		
	}
	PnMoy.set(&PnTilde);
	
	Apt12.add(&Apt1, &Apt2);
}


void OPFADMM::updateLocalProb() {
	// FB 11a
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	updateBp1();
	updatepl();

#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 1, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 11b
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	PnMoy.set(0.0);
	for (int n = 0; n < _nAgent; n++) {
		int bus = (int) _CoresBusAgent.get(n, 0);
		PnMoy.increment(bus, 0, Pn.get(n, 0));
		PnMoy.increment(bus + _nBus, 0, Pn.get(n + _nAgent, 0));

	}
	for (int b = 0; b < _nBus; b++) {
		int nb = (int) _nAgentByBus.get(b, 0);
		if (nb > 0) {
			PnMoy.set(b, 0, PnMoy.get(b, 0) / nb );
			PnMoy.set(b + _nBus, 0, PnMoy.get(b + _nBus, 0) / nb);
		}
		
	}

#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 2, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());

	// FB 11c
	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	for (int b = 0; b < _nBus; b++) {
		int nb = (int) _nAgentByBus.get(b, 0);
		if (nb > 0) {
			Bpt1.set(b, 0, MuL.get(b, 0) + PnMoy.get(b, 0));
			Bpt1.set(b +_nBus, 0, MuL.get(b + _nBus, 0) + PnMoy.get(b + _nBus, 0));
		}
	}
	
	updatePTilde();
	float mu = 0;
	for (int b = 0; b < _nBus; b++) {
		int nb = (int) _nAgentByBus.get(b, 0);
		if (nb > 0) {
			mu = MuL.get(b, 0) + PnMoy.get(b, 0) - PnTilde.get(b, 0);
			MuL.set(b, 0, mu);
			mu = MuL.get(b + _nBus, 0) + PnMoy.get(b + _nBus, 0) - PnTilde.get(b + _nBus, 0);
			MuL.set(b + _nBus, 0, mu);
		}
	}
	
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 3, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
	
}

void OPFADMM::updateXPQ()
{	
	for (int i=0; i<_nBus; i++){
		X[i].set(indpi, 0, PnTilde.get(i, 0) * _nAgentByBus.get(i, 0));
		X[i].set(indqi, 0, PnTilde.get(i + _nBus, 0) * _nAgentByBus.get(i, 0));
	}
}

void OPFADMM::updateBp1()
{
	for (int i = 0; i < _nAgent; i++) {
		int bus = (int) _CoresBusAgent.get(i, 0);
		int Nb = (int) _nAgentByBus.get(bus, 0);
		if (Nb > 1) {
			float m = Pn.get(i, 0) - PnMoy.get(bus, 0) + PnTilde.get(bus, 0) - MuL.get(bus, 0);
			Bp1.set(i, 0, m);
			m = Pn.get(i + _nAgent, 0) - PnMoy.get(bus + _nBus, 0) + PnTilde.get(bus + _nBus, 0) - MuL.get(bus + _nBus, 0);
			Bp1.set(i + _nAgent, 0, m);
		}
	}
}

void OPFADMM::updatepl()
{
	float pn = 0;
	PnPre.swap(&Pn);
	for (int n = 0; n < _nAgent; n++) {
		int b = (int) _CoresBusAgent.get(n, 0);
		int Nb = (int) _nAgentByBus.get(b, 0);
		if (Nb == 1) { // � ne faire qu'� la premier it�ration...
			if (_iterLocal == 0) {
				float ub = Pmax.get(n, 0);
				float lb = Pmin.get(n, 0);
				pn = (_rho * Bpt2.get(b, 0) - Cost2.get(n, 0)) / ((_rho + Cost1.get(n, 0)));
				pn = (ub - pn) * (pn > ub) + (lb - pn) * (pn < lb) + pn;
				Pn.set(n, 0, pn);
				PnPre.set(n, 0, pn);

				ub = Pmax.get(n + _nAgent, 0);
				lb = Pmin.get(n + _nAgent, 0);
				pn = (_rho * Bpt2.get(b + _nBus, 0) - Cost2.get(n + _nAgent, 0)) / ((_rho + Cost1.get(n + _nAgent, 0)));
				pn = (ub - pn) * (pn > ub) + (lb - pn) * (pn < lb) + pn;
				//std::cout << _rho << " "<< Bpt2.get(b + _nBus, 0)<< " "<< Cost2.get(n + _nAgent, 0) << " "<< Cost1.get(n + _nAgent, 0)<<" "<< lb << " "<< ub << " " <<   pn << std::endl;
				Pn.set(n + _nAgent, 0, pn);
				PnPre.set(n + _nAgent, 0, pn);
			}
		}
		else {
			float ub = Pmax.get(n, 0);
			float lb = Pmin.get(n, 0);
			pn = (Bp1.get(n, 0) * _rhol - Cost2.get(n, 0)) / Ap12.get(n, 0);
			pn = (ub - pn) * (pn > ub) + (lb - pn) * (pn < lb) + pn;
			Pn.set(n, 0, pn);

			ub = Pmax.get(n + _nAgent, 0);
			lb = Pmin.get(n + _nAgent, 0);
			pn = (Bp1.get(n + _nAgent, 0) * _rhol - Cost2.get(n + _nAgent, 0)) / Ap12.get(n + _nAgent, 0);
			pn = (ub - pn) * (pn > ub) + (lb - pn) * (pn < lb) + pn;
			Pn.set(n + _nAgent, 0, pn);
		}
	}
	/*
	Pn.set(&Bp1);
	Pn.multiply(_rhol);

	Pn.subtract(&Cost2);
	Pn.divideT(&Ap12);
	Pn.project(&Pmin, &Pmax);*/
}

void OPFADMM::updatePTilde()
{
	float pn = 0;
	for (int b = 0; b < _nBus; b++) {
		int n = (int) CoresSoloBusAgent.get(b, 0);
		if (n != -1) {
			if (_iterLocal == 0) {
				PnTilde.set(b, 0, PnMoy.get(b, 0));
				PnTilde.set(b + _nBus, 0, PnMoy.get(b + _nBus, 0));
			}
		}
		else {
			int nb = (int) _nAgentByBus.get(b, 0);
			if (nb > 0) {
				pn = (Bpt1.get(b, 0) * Apt1.get(b, 0) + Bpt2.get(b, 0) * Apt2.get(b, 0)) / Apt12.get(b, 0);
				PnTilde.set(b, 0, pn);
				pn = (Bpt1.get(b + _nBus, 0) * Apt1.get(b + _nBus, 0) + Bpt2.get(b + _nBus, 0) * Apt2.get(b + _nBus, 0)) / Apt12.get(b + _nBus, 0);
				PnTilde.set(b + _nBus, 0, pn);
			}
		}
	}

	/*PnTilde.set(&Bpt1);
	PnTilde.multiplyT(&Apt1);
	tempB2.set(&Bpt2);
	tempB2.multiplyT(&Apt2);
	PnTilde.add(&tempB2);

	PnTilde.divideT(&Apt12);

	PnTilde.project(&PnTmin, &PnTmax);*/
}

float OPFADMM::calcRes()
{
	float d1 = Pn.max2(&PnPre);
	float d2 = PnMoy.max2(&PnTilde);

	return d1 * (d1 > d2) + d2 * (d2 >= d1);
}

void OPFADMM::updateChat()
{
	//Ytrans (Pi, Qi, vi, li, pi, qi, Pai, Qai, lai, vij) !!!
	//Y      (Pi, Qi, vi, li, pi, qi, vai, Pji, Qji, lji)
	for (int i = 0; i < _nBus; i++) {
		
		float Phat, Qhat, lhat, phat, qhat;
		float vhat = 0;
		float muhat = 0;
		float m = nChild.get(i, 0);
		// est ce que c'est y - MU/rho ou y + MU/rho (comme dans le papier mais le signe change je ne sais pas pourquoi)
		
		Phat = Y[i].get(indPi, 0) / 2 - Mu[i].get(indPi, 0) / (2 * _rho);
		Qhat = Y[i].get(indQi, 0) / 2 - Mu[i].get(indQi, 0) / (2 * _rho);
		lhat = Y[i].get(indli, 0) / 2 - Mu[i].get(indli, 0) / (2 * _rho);
		if (i > 0) {
			int Ai = (int) Ancestor.get(i, 0);
			int c  = (int) PosChild.get(i, 0);
			Phat += Y[Ai].get(indPc + 3 * c, 0) / 2 - Mu[Ai].get(indPc + 3 * c, 0) / (2 * _rho); // or Y[Ai].get(7 + 3 * c, 0)
			Qhat += Y[Ai].get(indQc + 3 * c, 0) / 2 - Mu[Ai].get(indQc + 3 * c, 0) / (2 * _rho);
			lhat += Y[Ai].get(indlc + 3 * c, 0) / 2 - Mu[Ai].get(indlc + 3 * c, 0) / (2 * _rho);
		}
		phat = Y[i].get(indpi, 0) - Mu[i].get(indpi, 0) / _rho;
		qhat = Y[i].get(indqi, 0) - Mu[i].get(indqi, 0) / _rho;
		for (int j = 0; j < m; j++) { // Vai of all childs
			int c = (int) Childs[i].get(j, 0);
			muhat += Mu[c].get(indvai, 0);
			vhat += Y[c].get(indvai, 0);
		}
		vhat = (Y[i].get(indvi, 0) + vhat) / (m + 1) - (Mu[i].get(indvi, 0) + muhat) / (_rho * (m + 1));
	
				

		Chat[i].set(indPi, 0, Phat);
		Chat[i].set(indQi, 0, Qhat);
		Chat[i].set(indvi, 0, vhat);
		Chat[i].set(indli, 0, lhat);
		
		int nA = (int) _nAgentByBus.get(i, 0);
		if (nA > 0) {
			Bpt2.set(i, 0, phat / nA);
			Bpt2.set(i + _nBus, 0, qhat / nA);
		}	
	}
}

void OPFADMM::CommunicationX()
{
/**/ // X = { Pi, Qi, vi, li, pi, qi, vAi, (Pci, Qci, lci) for all child Ci }
	
	for (int i = 0; i < _nBus; i++) {
		if (i > 0) {
			int Ai = (int) Ancestor.get(i, 0);
			X[i].set(indvai, 0, X[Ai].get(indvi, 0));
		}
		int m = (int) nChild.get(i, 0);
		for (int j = 0; j < m; j++) {
			int c = (int) Childs[i].get(j, 0);
			X[i].set(indPc + 3 * j, 0, X[c].get(indPi, 0));
			X[i].set(indQc + 3 * j, 0, X[c].get(indQi, 0));
			X[i].set(indlc + 3 * j, 0, X[c].get(indli, 0));
		}
	}
	//Ytrans (Pi, Qi, vi, li, pi, qi, pai, qai, lai, vij) !!!
	//Y      (Pi, Qi, vi, li, pi, qi, vai, Pji, Qji, lji)
	// Q udate in argmin 0.5yHy + Qy

	for (int i = 0; i < _nBus; i++) {
		for (int j = 0; j < sizeOPFADMM.get(i, 0); j++) {
			Q[i].set(j, 0, -(Mu[i].get(j, 0) + _rho * X[i].get(j, 0)));
		}
	}

}


float OPFADMM::updateRes(int indice) 
{
	float resS = 0;
	float resR = 0;
	float resV = 0;
	for (int i = 0; i < _nBus; i++) {
		
		float resTempS = Y[i].max2(&Ypre[i]);
		float resTempR = Y[i].max2(&X[i]);

		if (resTempS > resS) {
			resS = resTempS;
		}
		if (resTempR > resR) {
			resR = resTempR;
		}

	}
	float oldrho = _rho;
	resF.set(0, indice, resR);
	resF.set(1, indice, oldrho * resS);
	resF.set(2, indice, resV);

	if (_tau > 1) {
		if (resR > _mu * resS) {
			_rho = _tau * _rho;
			Apt2.multiply(_tau);
			for (int i = 0; i < _nBus; i++) {
				Hinv[i].divide(_tau);
			}
			Apt12.add(&Apt1, &Apt2);
			//std::cout << _iterGlobal << "rho augmente " << _rho << std::endl;
		}
		else if (resS > _mu * resR) {// rho = rho / tau_inc;
			_rho = _rho / _tau;
			Apt2.divide(_tau);
			for (int i = 0; i < _nBus; i++) {
				Hinv[i].multiply(_tau);
			}
			Apt12.add(&Apt1, &Apt2);
			//std::cout << _iterGlobal << "rho diminue " << _rho << std::endl;
		}/**/
	}
	


	return MYMAX(MYMAX(resV, oldrho *resS),  resR);
}





