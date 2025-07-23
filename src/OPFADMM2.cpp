#include "../head/OPFADMM2.h"
 


OPFADMM2::OPFADMM2() : MethodOPF()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " OPFADMM2 Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
}

OPFADMM2::OPFADMM2(float rho) : MethodOPF(rho)
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default OPFADMM2 Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
}

void OPFADMM2::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
		timePerBlock.increment(0, 0, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION
	}
	
	_iterG = sim.getIterG();
	int iterL = sim.getIterL();
	_stepG = sim.getStepG();
	int stepL = sim.getStepL();
	
	float epsG = sim.getEpsG();
	float epsL = sim.getEpsL();
	
	
	float fc = 0;
	float resG = 2 * epsG;
	_iterGlobal = 0;
	

	//Chat.display();
	//Bpt2.display();

	while ((_iterGlobal < _iterG) && (resG>epsG)) {
		/*std::cout << "---------------------------------" << std::endl;
		for (int i = 0; i < 2; i++) {
			std::cout << " X " << i << std::endl;
			X[i].display();
			std::cout << " Q " << i << std::endl;
			Q[i].display();
			std::cout << " Y " << i << std::endl;
			Y[i].display();
			std::cout << " Mu " << i << std::endl;
			Mu[i].display();
			std::cout << " Chat " << i << std::endl;
			Chat[i].display();
		}*/
#ifdef INSTRUMENTATION
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		updateXWOCurrent();
		updateXPQ();
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 1, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		
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
			//resG = 1;
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
	occurencePerBlock.increment(0, 1, _iterGlobal);
	occurencePerBlock.increment(0, 6, _iterGlobal);
	occurencePerBlock.increment(0, 7, _iterGlobal);
	occurencePerBlock.increment(0, 8, _iterGlobal);
	occurencePerBlock.increment(0, 9, _iterGlobal / _stepG);

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	for (int b = 0; b < _nBus; b++) {
		// pi & qi
		int Nb = (int) _nAgentByBus.get(b, 0);
		int begin = (int) _CoresAgentBusBegin.get(b, 0);
		for (int In = 0; In < Nb; In++) {
			int n = (int) _CoresAgentBus.get(In + begin, 0);
			// pn & qn
			Pn.set(n, 0, X[b].get(5 + 2 * In, 0));
			Pn.set(n + _nAgent, 0, X[b].get(6 + 2 * In, 0));
		}
	}

	fc = calcFc();

	Pn.set(0, 0, getPLoss());
	Pn.set(_nAgent, 0, getQLoss());
	// FB 5
	
	result->setResF(&resF);
	/*std::cout << "--------" << std::endl;
	std::cout << " Pn " << std::endl;
	Pn.display();
	for (int i = 0; i < 2; i++) {
		std::cout << " X " << i << std::endl;
		X[i].display();
	}
	*/
	MatrixCPU Pb(getPb());
	MatrixCPU Phi(getPhi());
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
	tall = clock() - tall;
	timeOPF = tall;

	result->setTime((float)tall / CLOCKS_PER_SEC);
	
}

void OPFADMM2::updateP0(const StudyCase& cas)
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

	// remove loss agent
	Pn.set(0, 0, 0);
	Pmin.set(0, 0, 0);
	Pmax.set(0, 0, 0);
	Pn.set(_nAgent, 0, 0);
	Pmin.set(_nAgent, 0, 0);
	Pmax.set(_nAgent, 0, 0);
	
	Pb.set(0.0);
	Pbmax.set(0.0);
	Pbmin.set(0.0);

	for (int b = 0; b < _nBus; b++) {
		int Nb = (int) _nAgentByBus.get(b, 0);
		int begin = (int) _CoresAgentBusBegin.get(b, 0);
		for (int In = 0; In < Nb; In++) {
			int n = (int) _CoresAgentBus.get(In + begin, 0);
			Pb.increment(b, 0, Pn.get(n, 0));
			Pb.increment(b + _nBus, 0, Pn.get(n + _nAgent, 0));
			Pbmax.increment(b, 0, Pmax.get(n, 0));
			Pbmin.increment(b, 0, Pmin.get(n, 0));

			Pbmax.increment(b + _nBus, 0, Pmax.get(n + _nAgent, 0));
			Pbmin.increment(b + _nBus, 0, Pmin.get(n + _nAgent, 0));
	
			// pn & qn
			X[b].set(5 + 2 * In, 0, Pn.get(n, 0));
			X[b].set(6 + 2 * In, 0, Pn.get(n + _nAgent, 0));
		}
	}
	
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
		int Nb = (int) _nAgentByBus.get(i, 0);
		for (int j = 0; j < m; j++) {// (Pci, Qci, lci) for all child Ci
			int c = (int) Childs[i].get(j, 0);
			X[i].set(5 + 2 * Nb + 3 * j, 0, X[c].get(0, 0));
			X[i].set(6 + 2 * Nb + 3 * j, 0, X[c].get(1, 0));
			X[i].set(7 + 2 * Nb + 3 * j, 0, X[c].get(3, 0));
		}
		Y[i].set(&X[i]);
		//Mu[i].set(0.0);
	}/**/

#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 11, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 11, 1);
	t1 = std::chrono::high_resolution_clock::now();
#endif

	updateChat();
#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);

#endif
}

void OPFADMM2::init(const Simparam& sim, const StudyCase& cas)
{
	// intitilisation des matrixs et variables 
	
	clock_t t = clock();
	//std::cout << "init " << std::endl;
	initSize(cas);
	

	initSimParam(sim);
	
	
	initAgent(cas);
	
	initGrid(cas);
	
	sizeOPFADMM = MatrixCPU(_nBus, 1);
	for (int i = 0; i < _nBus; i++) {
		int Nb = (int) _nAgentByBus.get(i, 0);
		int sizeOPF = 3 * (int) nChild.get(i, 0) + 5 + 2 * Nb;
		sizeOPFADMM.set(i, 0, sizeOPF);
	}

	allocateTab();
	
		
	
	
	//std::cout << " init valeur " << std::endl;
	for (int i = 0; i < _nBus; i++) {
		int Nb = (int) _nAgentByBus.get(i, 0);
		int begin = (int) _CoresAgentBusBegin.get(i, 0);
		Chat[i] = MatrixCPU(4 + 2 * Nb, 1);
		// vi = 1; vai = 1
		X[i].set(indvi, 0, 1); 
		X[i].set(indvai, 0, 1);

		// pi & qi
		
		for (int In = 0; In < Nb; In++) {
			int n = (int) _CoresAgentBus.get(In + begin, 0);
			PosAgent.set(n, 0, In);
			Pb.increment(i, 0, Pn.get(n, 0));
			Pb.increment(i + _nBus, 0, Pn.get(n + _nAgent, 0));
			//X[i].increment(4, 0, Pn.get(n, 0));
			//X[i].increment(5, 0, Pn.get(n + _nAgent, 0));
			
			Pbmax.increment(i, 0, Pmax.get(n, 0));
			Pbmin.increment(i, 0, Pmin.get(n, 0));

			Pbmax.increment(i + _nBus, 0, Pmax.get(n + _nAgent, 0));
			Pbmin.increment(i + _nBus, 0, Pmin.get(n + _nAgent, 0));

			if (Nb == 1) {
				CoresSoloBusAgent.set(i, 0, n);
			}
			// pn & qn
			X[i].set(indpi + 2 * In, 0, Pn.get(n, 0));
			X[i].set(indqi + 2 * In, 0, Pn.get(n + _nAgent, 0));
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
		int Ni = (int) _nAgentByBus.get(i, 0);
		for (int j = 0; j < m; j++) { // (Pci, Qci, lci) for all child Ci
			int c = (int) Childs[i].get(j, 0);
			X[i].set(indPc + 2 * Ni + 3 * j, 0, X[c].get(indPi, 0));
			X[i].set(indQc + 2 * Ni + 3 * j, 0, X[c].get(indQi, 0));
			X[i].set(indlc + 2 * Ni + 3 * j, 0, X[c].get(indli, 0));		
		}
		Y[i].set(&X[i]);
	}


	//std::cout << " Hinv " << std::endl;
	for (int i = 0; i < _nBus; i++) {
		int Ni = (int) _nAgentByBus.get(i, 0);
		if (i > 0) {
			A[i].set(2, indPi, 2 * ZsRe.get(i - 1, 0));
			A[i].set(2, indQi, 2 * ZsIm.get(i - 1, 0));
			A[i].set(2, indli, -ZsNorm.get(i - 1, 0));
			A[i].set(2, indvi, -1);
			A[i].set(2, indvai, 1);
			A[i].set(0, indPi, -1);
			A[i].set(1, indQi, -1);
		}
		
		// pi = sum(pn) & qi = sum(qn)
		for (int In = 0; In < Ni; In++) {
			A[i].set(0, indpi + 2 * In, 1);
			A[i].set(1, indqi + 2 * In, 1);
		}
		
		for (int j = 0; j < nChild.get(i, 0); j++) {
			int c = (int) Childs[i].get(j, 0);
			A[i].set(0, indPc + 2 * Ni + 3 * j, 1); // Pci
			A[i].set(1, indQc + 2 * Ni + 3 * j, 1); // Qci
			A[i].set(0, indlc + 2 * Ni + 3 * j, -ZsRe.get(c - 1, 0)); // -R l
			A[i].set(1, indlc + 2 * Ni + 3 * j, -ZsIm.get(c - 1, 0)); // -X l
		}
		//A[i].display();
		MatrixCPU temp33(2 + 1 * (i > 0), 2 + 1 * (i > 0));
		MatrixCPU temp3M(2 + 1 * (i > 0), (int) sizeOPFADMM.get(i, 0));
		MatrixCPU tempMM( (int) sizeOPFADMM.get(i, 0), (int) sizeOPFADMM.get(i, 0));
		temp33.multiplyTrans(&A[i], &A[i]); // (3*o_b) * (o_b*3) -> 9 * o_b^2
		temp33.invertGaussJordan(&temp33); // 3^3 = 27 (fixe !!!)
		temp3M.MultiplyMatMat(&temp33, &A[i]); // (3*3) * (3*o_b) -> 27 *o_b
		Hinv[i].multiplyTrans(&A[i], &temp3M, 0); // (o_b*3) * (3*o_b) -> 9 * o_b

		tempMM.setEyes(1);
		Hinv[i].subtract(&tempMM);
		Hinv[i].divide(_rho);
		//Hinv[i].display();
	}
	
	//std::cout << " updateChat " << std::endl;
	updateChat();
	/*std::cout << "--------" << std::endl;
	std::cout << " Pn " << std::endl;
	Pn.display();*/
	/*for (int i = 0; i < _nBus; i++) {
		std::cout << " X " << i << std::endl;
		X[i].display();
		std::cout << " Y " << i << std::endl;
		Y[i].display();
		std::cout << " Q " << i << std::endl;
		Q[i].display();
		std::cout << " Mu " << i << std::endl;
		Mu[i].display();
		std::cout << " Chat " << i << std::endl;
		Chat[i].display();
	}*/
	

	
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	//std::cout << "---------------------------------------------------------------------------------------" << std::endl;
}

void OPFADMM2::updateXPQ(){
	for (int i = 0; i < _nBus; i++) {
		int Nb = (int) _nAgentByBus.get(i, 0);
		int begin = (int) _CoresAgentBusBegin.get(i, 0);
		for (int In = 0; In < Nb; In++) {
			int n = (int) _CoresAgentBus.get(In + begin, 0);
			float ub = Pmax.get(n, 0);
			float lb = Pmin.get(n, 0);
			float pn = (_rho * Chat[i].get(indChatpi + 2 * In, 0) - Cost2.get(n, 0)) / (Cost1.get(n,0) + _rho);
			pn = (ub - pn) * (pn > ub) + (lb - pn) * (pn < lb) + pn;
		

			ub = Pmax.get(n + _nAgent, 0);
			lb = Pmin.get(n + _nAgent, 0);
			float qn = (_rho * Chat[i].get(indChatqi + 2 * In, 0) - Cost2.get(n + _nAgent, 0)) / (Cost1.get(n + _nAgent, 0) + _rho);
			qn = (ub - qn) * (qn > ub) + (lb - qn) * (qn < lb) + qn;
			
			// pn & qn
			X[i].set(indpi + 2 * In, 0, pn);
			X[i].set(indqi + 2 * In, 0, qn);
		}

	}
}


void OPFADMM2::updateChat()
{
	//Y      (Pi, Qi, vi, li, pi, qi, vai,(pn, qn), (Pji, Qji, lji))
	for (int i = 0; i < _nBus; i++) {
		//std::cout << "***** Bus " << i << std::endl;
		int Nb = (int) _nAgentByBus.get(i, 0);
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
			int c =  (int) PosChild.get(i, 0);
			int NAi = (int) _nAgentByBus.get(Ai, 0);
			Phat += Y[Ai].get(indPc + 2 * NAi + 3 * c, 0) / 2 - Mu[Ai].get(indPc + 2 * NAi + 3 * c, 0) / (2 * _rho);
			Qhat += Y[Ai].get(indQc + 2 * NAi + 3 * c, 0) / 2 - Mu[Ai].get(indQc + 2 * NAi + 3 * c, 0) / (2 * _rho);
			lhat += Y[Ai].get(indlc + 2 * NAi + 3 * c, 0) / 2 - Mu[Ai].get(indlc + 2 * NAi + 3 * c, 0) / (2 * _rho);
		}
		for (int j = 0; j < m; j++) { // Vai of all childs
			//std::cout << "Child " << j << std::endl;;
			int c = (int) Childs[i].get(j, 0);
			muhat += Mu[c].get(indvai, 0);
			vhat += Y[c].get(indvai, 0);
		}
		vhat = (Y[i].get(indvi, 0) + vhat) / (m + 1) - (Mu[i].get(indvi, 0) + muhat) / (_rho * (m + 1));
		
		Chat[i].set(indPi, 0, Phat);
		Chat[i].set(indQi, 0, Qhat);
		Chat[i].set(indvi, 0, vhat);
		Chat[i].set(indli, 0, lhat);
		
		for (int In = 0; In < Nb; In++) {
			//std::cout << "Agent " << In << std::endl;;
			
			phat = Y[i].get(indpi + 2 * In, 0) - Mu[i].get(indpi + 2 * In, 0) / _rho;
			qhat = Y[i].get(indqi + 2 * In, 0) - Mu[i].get(indqi + 2 * In, 0) / _rho;
			
			
			Chat[i].set(indChatpi + 2 * In, 0, phat);
			Chat[i].set(indChatqi + 2 * In, 0, qhat);
		}
		
	}


}

void OPFADMM2::CommunicationX()
{
/**/ // X = { Pi, Qi, vi, li, vAi, (pn, qn) (Pci, Qci, lci) for all child Ci }
	
	for (int i = 0; i < _nBus; i++) {
		int Ni = (int) _nAgentByBus.get(i, 0);
		if (i > 0) {
			int Ai = (int) Ancestor.get(i, 0);
			X[i].set(indvai, 0, X[Ai].get(indvi, 0));
		}

		int m = (int) nChild.get(i, 0);
		for (int j = 0; j < m; j++) {
			int c = (int) Childs[i].get(j, 0);
			X[i].set(indPc + 2 * Ni + 3 * j, 0, X[c].get(indPi, 0));
			X[i].set(indQc + 2 * Ni + 3 * j, 0, X[c].get(indQi, 0));
			X[i].set(indlc + 2 * Ni + 3 * j, 0, X[c].get(indli, 0));
		}
	}
	
	//Y      (Pi, Qi, vi, li, pi, qi, vai, (pn, qn) Pji, Qji, lji)
	// Q udate in argmin 0.5yHy + Qy

	for (int i = 0; i < _nBus; i++) {
		for (int j = 0; j < sizeOPFADMM.get(i, 0); j++) {
			Q[i].set(j, 0, -(Mu[i].get(j, 0) + _rho * X[i].get(j, 0)));
		}
	}   

}

// (Pi,Qi,vi,li,pi,qi,vai,pji,qji,lji)


