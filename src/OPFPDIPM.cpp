#include "../head/OPFPDIPM.h"
#define MAX(X, Y) X * (X >= Y) + Y * (Y > X)


OPFPDIPM::OPFPDIPM() : MethodOPF()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " OPFPDIPM Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	timePerBlock = MatrixCPU(1, 12, 0); // Fb0, Fb11abcd, FB12, Fb2, Fb3, Fb4, Fb5,FB6, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 12, 0); //nb de fois utilis� pendant la simu
}


OPFPDIPM::OPFPDIPM(float sigma) : MethodOPF()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default OPFPDIPM Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_sigma = sigma;
	timePerBlock = MatrixCPU(1, 12, 0); // Fb0, Fb11, FB12, Fb2, Fb3, Fb4, Fb5, FB6, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 12, 0); //nb de fois utilis� pendant la simu
}

OPFPDIPM::~OPFPDIPM()
{
		
}
void OPFPDIPM::setParam(float sigma)
{
	_sigma = sigma;
}


void OPFPDIPM::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
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
	
	
	float epsG = sim.getEpsG();
	
	
	CgapBest = Cgap;
	Xbest = X;
	float fc = 0;
	float resG = 2 * epsG;
	
	_iterGlobal = 0;
	


	while ((_iterGlobal < _iterG) && (resG>epsG)) {
		
		
		
		/*std::cout << " l , u " << std::endl;
		l.display();
		u.display();
		std::cout << " z , w " << std::endl;
		z.display();
		w.display();
				std::cout << "X" << std::endl;
		X.display();
		std::cout << "Y" << std::endl;
		Y.display();
		std::cout << " z , w " << std::endl;
		z.display();
		w.display();
		std::cout << "g Jac" << std::endl;
		gJac.display();
		std::cout << "X Jac" << std::endl;
		XJac.display();
		std::cout << "h Jac" << std::endl;
		hJac.display();
		std::cout << "g Jac" << std::endl;
		gJac.display();
			std::cout << "gHess" << std::endl;
		gHess.display();
		std::cout << "hHess" << std::endl;
		hHess.display();
		std::cout << "Hess" << std::endl;
		Hess.display();
		std::cout << "Mat" << std::endl;
		MatSys.display();
		std::cout << "Vect" << std::endl;
		VectSys.display();*/

		/// STEP 1
		// 1a : avoir PSI(O) - > c'est d�j� fait avant
		//std::cout << "1a" << std::endl;
		
		updateJ();
		//std::cout << "1b" << std::endl;
		updateH();
		//std::cout << "1c" << std::endl;
		computeConstraint();
		//std::cout << "1d" << std::endl;
		setSystem();
		// 1b : solve systeme
		
		correctionEquation();
		
		// 1c : trouver step
		//std::cout << "1f" << std::endl;
		updateStep();
		// 1d : calculer Caff & calculer sigma
		bool sucess = setCenteringParameter();
		if(_saveInterSol && !sucess)
		{
			if (Cgap < CgapBest) {
				_iterBest = _iterGlobal - 1;
				CgapBest = Cgap;
				Xbest = X;
			}
			//break;
		}
		
#ifdef INSTRUMENTATION
		occurencePerBlock.increment(0, 1, 1);
		occurencePerBlock.increment(0, 2, 1);
		occurencePerBlock.increment(0, 3, 1);
		occurencePerBlock.increment(0, 4, 1 );
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

		/// Step 2
		// 2a : calculer mu
		
		updatePerturbedFactor();
		
		// 2b : calculer PSI(*, mu)
		updatePSILlu();
		

#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 5, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		
		/// STEP 3 : solve sytem
		//std::cout << "3a" << std::endl;
		correctionEquation();
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 6, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		
		/// STEP 4 : trouver step
		updateStep();
		
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 7, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		/// STEP 5 : update var
		
		updateVariable();

		
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
		t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
		/// STEP 6
		//std::cout << "6" << std::endl;
		resG = updateRes(_iterGlobal);
			
#ifdef INSTRUMENTATION
		t2 = std::chrono::high_resolution_clock::now();
		timePerBlock.increment(0, 9, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
#endif // INSTRUMENTATION
		//std::cout << iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;

		_iterGlobal++;
	}
	//std::cout << _iterGlobal << " " << resF.get(0, (_iterGlobal - 1) / _stepG) << " " << resF.get(1, (_iterGlobal - 1) / _stepG) << " " << resF.get(2, (_iterGlobal - 1) / _stepG) << std::endl;


#ifdef INSTRUMENTATION	
	occurencePerBlock.increment(0, 5, _iterGlobal);
	occurencePerBlock.increment(0, 6, _iterGlobal);
	occurencePerBlock.increment(0, 7, _iterGlobal);
	occurencePerBlock.increment(0, 8, _iterGlobal);
	occurencePerBlock.increment(0, 9, _iterGlobal / _stepG);

	t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	if (_saveInterSol) {
		if (isnan(Cgap) ||  Cgap<0 || Cgap>CgapBest) {
			//std::cout << "Solution enregistree utilisee  ! " << std::endl;
			X = Xbest;
			//std::cout << _iterBest << " " << resF.get(0, (_iterBest - 1) / _stepG) << " " << resF.get(1, (_iterBest - 1) / _stepG) << " " << resF.get(2, (_iterBest - 1) / _stepG) << std::endl;

		}
	}


	for (int n = 0; n < _nAgent; n++) {
		Pn.set(n + 1, 0, X.get(n, 0));
		Pn.set(n + 2 + _nAgent, 0, X.get(n + _nAgent, 0));
	}

	
	fc = calcFc(&Cost1, &Cost2, &Pn, &tempN2);
	// FB 5
	
	result->setResF(&resF);
	
	

	result->setIter(_iterGlobal);
	

	result->setPn(&Pn);
	
	result->setFc(fc);

#ifdef INSTRUMENTATION
	t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 10, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 10, 1);

	result->setTimeBloc(&timePerBlock, &occurencePerBlock);
#endif // INSTRUMENTATION
	tall = clock() - tall;
	timeOPF = tall;

	result->setTime((float)tall / CLOCKS_PER_SEC);
	
}

void OPFPDIPM::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

	

	Pmin = cas.getPmin();
	Pmax = cas.getPmax();
	Cost2 = cas.getb();
	
	Pn.add(&Pmax, &Pmin); // ou juste projection
	Pn.divide(2);

	for (int n = 1; n < _nAgent; n++) {
		X.set(n, 0, Pn.get(n, 0));
		X.set(n + _nAgent, 0, Pn.get(n + _nAgent, 0));
	}
	


#ifdef INSTRUMENTATION
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 11, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 11, 1);
#endif // INSTRUMENTATION

}

void OPFPDIPM::init(const Simparam& sim, const StudyCase& cas)
{
	// intitilisation des matrixs et variables 
	
	clock_t t = clock();
	//std::cout << "init " << std::endl;

	
	_iterG = sim.getIterG() < 50 ? sim.getIterG() : 50;
	_stepG = sim.getStepG();
	
	_nAgent = cas.getNagent() - 1; // on enl�ve l'agent des pertes
	
	_nBus = cas.getNBus();
	_nLine = cas.getNLine(true); // ne doit pas �tre r�duit ici !!!
	_sizeEq = 2 * (_nBus + 1);
	_sizeInEq = 2 * _nAgent + _nBus + _nLine;
	_sizeVar = 2 * _nAgent + 2 * _nBus;
	_sigma = 0.5;

	//std::cout << _nAgent << " " << _nBus << " " << _nLine << std::endl;
	_nAgentByBus = cas.getNagentByBus();
	_nAgentByBus.increment(0, 0, -1); // don't want the grid agent in this resolution
	
	CoresVoiLin = cas.getCoresVoiLin();
	CoresBusLin = cas.getCoresBusLin();
	CoresLineBus = cas.getCoresLineBus(true);
	nLines = cas.getNLines();
	_CoresBusAgent = cas.getCoresBusAgentLin(); // Cores[n] = b
	BgridLin = cas.getBlin();
	GgridLin = cas.getGlin();
	//std::cout << "BgridLin " << std::endl;
	//BgridLin.display();
	//std::cout << "GgridLin " << std::endl;
	//GgridLin.display();
	//std::cout << "coresBusAgent " << std::endl;
	//_CoresBusAgent.display();
	
	resF = sim.getRes();

	MatrixCPU lowerBound(cas.getLowerBound());
	MatrixCPU upperBound(cas.getUpperBound());
	Pmin = cas.getPmin();
	Pmax = cas.getPmax();
	LowerBound = MatrixCPU(_sizeInEq, 1); //voltage angle, voltage, line...
	UpperBound = MatrixCPU(_sizeInEq, 1);; //voltage angle, voltage, line...
	
	for (int n = 0; n < _nAgent; n++) {
		LowerBound.set(n, 0, Pmin.get(n + 1, 0));
		LowerBound.set(n + _nAgent, 0, Pmin.get(n + 2 + _nAgent, 0));
		UpperBound.set(n , 0, Pmax.get(n + 1, 0));
		UpperBound.set(n + _nAgent, 0, Pmax.get(n + 2 + _nAgent, 0));
	}
	_V0 = upperBound.get(_nBus, 0);
	float eps = 0.01;
	for (int b = 0; b < _nBus; b++) {
		float V = upperBound.get(b + _nBus, 0);
		UpperBound.set(b + 2 * _nAgent, 0, V * V);
		V = lowerBound.get(b + _nBus, 0);
		LowerBound.set(b + 2 * _nAgent, 0, V * V);
		if (b == 0) {
			UpperBound.increment(2 * _nAgent, 0, eps);
			LowerBound.increment(2 * _nAgent, 0, -eps);
		}
	}
	for (int i = 0; i < _nLine; i++) {
		float pij = upperBound.get(2 * _nBus + i, 0);
		UpperBound.set(_nBus + 2 * _nAgent + i, 0, pij);
		pij = lowerBound.get(2 * _nBus + i, 0);
		LowerBound.set(_nBus + 2 * _nAgent + i, 0, pij);
	}
	//std::cout << "Bounds " << std::endl;
	//LowerBound.display();
	//UpperBound.display();
	

	Cost1 = MatrixCPU(cas.geta());
	Cost2 = MatrixCPU(cas.getb());
	Pn = sim.getPn();
	
	if (Pn.max2() == 0) {
		if (_initWithMarket) {
			initWithMarket(sim,cas);
		}
		else if (_initFlatPower) {
			// do nothing
		}
		else {
			Pn.add(&Pmax, &Pmin);
			Pn.divide(2);
		}
		
	}
	
	Pn.project(&Pmin, &Pmax);
	
	Pn.set(0, 0, 0);
	//Pn.display();

	if (_initWithPF) {
		initWithPF(cas);
	}
	else {
		E = MatrixCPU(2 * _nBus, 1);
		for (int i = 0; i < _nBus; i++) {
			E.set(i + _nBus, 0, 1);
		}
	}

	
	dX = MatrixCPU(_sizeVar, 1);
	dY = MatrixCPU(_sizeEq, 1);
	dl = MatrixCPU(_sizeInEq, 1);
	du = MatrixCPU(_sizeInEq, 1);
	dz = MatrixCPU(_sizeInEq, 1);
	dw = MatrixCPU(_sizeInEq, 1);

	Lx = MatrixCPU(_sizeVar, 1);
	Ly = MatrixCPU(_sizeEq, 1);
	Ll = MatrixCPU(_sizeInEq, 1);
	Lu = MatrixCPU(_sizeInEq, 1);
	Lz = MatrixCPU(_sizeInEq, 1);
	Lw = MatrixCPU(_sizeInEq, 1);

	//std::cout << " creation " << std::endl;
	X = MatrixCPU(_sizeVar, 1);
	Xpre = MatrixCPU(_sizeVar, 1);
	
	
	Y = MatrixCPU(_sizeEq, 1);
	l = MatrixCPU(_sizeInEq, 1, valMin);
	u = MatrixCPU(_sizeInEq, 1, valMin);
	z = MatrixCPU(_sizeInEq, 1, valMin);
	w = MatrixCPU(_sizeInEq, 1, -valMin);
	d = MatrixCPU(_sizeInEq, 1);

	dXY = MatrixCPU(_sizeVar + _sizeEq, 1);
	PSI = MatrixCPU(_sizeVar, 1);
	
	tempS1 = MatrixCPU(_sizeInEq, 1);
	tempM1 = MatrixCPU(_sizeVar, 1);
	tempS1bis = MatrixCPU(_sizeInEq, 1);
	tempN2 = MatrixCPU(2 * (_nAgent + 1), 1);

	XHess = MatrixCPU(_sizeVar, _sizeVar);
	hHess = MatrixCPU(_sizeVar, _sizeVar);
	gHess = MatrixCPU(_sizeVar, _sizeVar);
	Hess  = MatrixCPU(_sizeVar, _sizeVar);

	XJac = MatrixCPU(_sizeVar, 1);
	hJac = MatrixCPU(_sizeVar, _sizeEq);
	gJac = MatrixCPU(_sizeVar, _sizeInEq);

	PnBus = MatrixCPU(_nBus * 2, 1);

	MatSys = MatrixCPU(_sizeVar + _sizeEq, _sizeVar + _sizeEq);
	VectSys = MatrixCPU(_sizeVar + _sizeEq, 1);

	

	for (int n = 0; n < _nAgent; n++) {
		X.set(n, 0, Pn.get(n + 1, 0));
		X.set(n + _nAgent, 0, Pn.get(n + 2 + _nAgent, 0));

		XHess.set(n, n, -Cost1.get(n + 1, 0));
		XHess.set(n + _nAgent, n + _nAgent, -Cost1.get(n + 2 + _nAgent, 0));
		l.set(n, 0, min(MAX((X.get(n, 0) - LowerBound.get(n, 0)), valMin), valMax));
		u.set(n, 0, min(MAX((UpperBound.get(n, 0) - X.get(n, 0)), valMin), valMax));
		l.set(n + _nAgent, 0, min(MAX((X.get(n + _nAgent, 0) - LowerBound.get(n + _nAgent, 0)), valMin), valMax));
		u.set(n + _nAgent, 0, min(MAX((UpperBound.get(n + _nAgent, 0) - X.get(n + _nAgent, 0)), valMin), valMax));
	}
	int offset = 2 * _nAgent;
	int indice = 0;
	for (int b = 0; b < _nBus; b++) {
		//float V = _V0*_V0;
		//X.set(b + offset, 0, _V0);// e =1, f = 0;
		float V = E.get(b + _nBus, 0) * E.get(b + _nBus, 0);
		float Vma = UpperBound.get(b + offset, 0);
		float Vmin = LowerBound.get(b + offset, 0);
		V = (Vma - V) * (V > Vma) + (Vmin - V) * (V < Vmin) + V;
		float ei = sqrt(V) * cos(E.get(b, 0));
		float fi = sqrt(V) * sin(E.get(b, 0));
		//std::cout << V << " " << Vmin << " " << Vma << " " << E.get(b, 0) << " " << cos(E.get(b, 0)) << " " << sin(E.get(b, 0)) << " " << ei << " " << fi << std::endl;
		X.set(b + offset, 0, ei);
		X.set(b + _nBus + offset, 0, fi);
		int k = CoresBusLin.get(b, 0);
	    l.set(b + offset, 0, min(MAX((V - Vmin), valMin), valMax));
		u.set(b + offset, 0, min(MAX((Vma - V), valMin), valMax));
		for (int voisin = k + 1; voisin < (k + nLines.get(b, 0)); voisin++) {
			int j = CoresVoiLin.get(voisin, 0);
			if (j > b) {
				float ej = E.get(j + _nBus, 0) * cos(E.get(j, 0));
				float fj = E.get(j + _nBus, 0) * sin(E.get(j, 0));
				float B = BgridLin.get(voisin, 0);
				float G = GgridLin.get(voisin, 0);
				float Pij = (ei * ei + fi * fi - ei * ej - fi * fj) * G + (ei * fj - ej * fi) * B;
				l.set(_nBus + offset + indice , 0, min(MAX((Pij - LowerBound.get(_nBus + offset + indice, 0)), valMin), valMax));
				u.set(_nBus + offset + indice, 0, min(MAX((UpperBound.get(_nBus + offset + indice, 0) - Pij), valMin), valMax));
				indice++;
			}
		}
	}

	//std::cout << " X " << std::endl;
	//X.display();




	
	Xpre.set(&X);
	updateRes(0);
	/*updateJ();
	updateH();
	computeConstraint();
	setSystem();*/

	//std::cout << "-------" << std::endl;
	//GgridLin.display();
	//BgridLin.display();
	//std::cout << "-------" << std::endl;
	//CoresBusLin.display();
	//CoresVoiLin.display();
	//nLines.display();
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	//std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	//std::cout << "---------------------------------------------------------------------------------------" << std::endl;
}

void OPFPDIPM::solveConsensus(float eps, MatrixCPU* PSO)
{
	

	float epsG = eps;


	CgapBest = Cgap;
	Xbest = X;
	float fc = 0;
	float resG = 2 * epsG;

	_iterGlobal = 0;



	while ((_iterGlobal < _iterG) && (resG > epsG)) {
		/// STEP 1
		// 1a : avoir PSI(O) - > c'est d�j� fait avant
		//std::cout << "1a" << std::endl;

		updateJ();

		//std::cout << "1b" << std::endl;

		updateH();

	
		//std::cout << "1c" << std::endl;
		computeConstraint();
		//std::cout << "1d" << std::endl;
		setSystem();
		// 1b : solve systeme
		//std::cout << "1e" << std::endl;
		
		correctionEquation();
		
		// 1c : trouver step
		//std::cout << "1f" << std::endl;
		updateStep();
		// 1d : calculer Caff & calculer sigma
		//std::cout << "1g" << std::endl;
		bool sucess = setCenteringParameter();
		if (_saveInterSol && !sucess)
		{
			if (Cgap < CgapBest) {
				_iterBest = _iterGlobal - 1;
				CgapBest = Cgap;
				Xbest = X;
			}
		}


		/// Step 2
		// 2a : calculer mu
		
		updatePerturbedFactor();
		// 2b : calculer PSI(*, mu)
		updatePSILlu();


		/// STEP 3 : solve sytem
		
		correctionEquation();
		


		/// STEP 4 : trouver step
		//std::cout << "4a" << std::endl;
		updateStep();
		
		/// STEP 5 : update var
		updateVariable();

		/// STEP 6

		resG = updateRes(_iterGlobal);
		//std::cout << iterGlobal << " " << _iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;

		_iterGlobal++;
	}
	std::cout << _iterGlobal << " " << resF.get(0, (_iterGlobal - 1) / _stepG) << " " << resF.get(1, (_iterGlobal - 1) / _stepG) << " " << resF.get(2, (_iterGlobal - 1) / _stepG) << std::endl;

	if (_saveInterSol) {
		if (isnan(Cgap) || Cgap<0 || Cgap>CgapBest) {
			std::cout << "Solution enregistree utilisee  ! " << std::endl;
			X = Xbest;
			std::cout << _iterBest << " " << resF.get(0, (_iterBest - 1) / _stepG) << " " << resF.get(1, (_iterBest - 1) / _stepG) << " " << resF.get(2, (_iterBest - 1) / _stepG) << std::endl;
		}
	}

	for (int n = 0; n < _nAgent; n++) {
		Pn.set(n + 1, 0, X.get(n, 0));
		Pn.set(n + 2 + _nAgent, 0, X.get(n + _nAgent, 0));
	}
	PSO->set(&Pn);
}

void OPFPDIPM::initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO)
{
	// intitilisation des matrixs et variables 

	clock_t t = clock();
	std::cout << "init " << std::endl;
	_rhoSO = rhoSO;

	_iterG = sim.getIterG();
	_stepG = sim.getStepG();

	_nAgent = cas.getNagent() - 1; // on enl�ve l'agent des pertes

	_nBus = cas.getNBus();
	_nLine = cas.getNLine(true); // ne doit pas �tre r�duit ici !!!
	_sizeEq = 2 * (_nBus + 1);
	_sizeInEq = 2 * _nAgent + _nBus + _nLine;
	_sizeVar = 2 * _nAgent + 2 * _nBus;
	_sigma = 0.5;

	std::cout << _nAgent << " " << _nBus << " " << _nLine << std::endl;
	_nAgentByBus = cas.getNagentByBus();
	_nAgentByBus.increment(0, 0, -1); // don't want the grid agent in this resolution

	CoresVoiLin = cas.getCoresVoiLin();
	CoresBusLin = cas.getCoresBusLin();
	CoresLineBus = cas.getCoresLineBus(true);
	nLines = cas.getNLines();
	_CoresBusAgent = cas.getCoresBusAgentLin(); // Cores[n] = b
	BgridLin = cas.getBlin();
	GgridLin = cas.getGlin();
	

	resF = MatrixCPU(3, (_iterG / _stepG) + 1);

	MatrixCPU lowerBound(cas.getLowerBound());
	MatrixCPU upperBound(cas.getUpperBound());
	
	// maibe change to not having constraint on power
	Pmin = MatrixCPU(2 * _nAgent, 1, -1000000); // must not be the real one
	Pmax = MatrixCPU(2 * _nAgent, 1, 1000000); // idem

	LowerBound = MatrixCPU(_sizeInEq, 1); //voltage angle, voltage, line...
	UpperBound = MatrixCPU(_sizeInEq, 1);; //voltage angle, voltage, line...

	for (int n = 0; n < _nAgent; n++) {
		LowerBound.set(n, 0, Pmin.get(n + 1, 0));
		LowerBound.set(n + _nAgent, 0, Pmin.get(n + 2 + _nAgent, 0));
		UpperBound.set(n, 0, Pmax.get(n + 1, 0));
		UpperBound.set(n + _nAgent, 0, Pmax.get(n + 2 + _nAgent, 0));
	}
	_V0 = upperBound.get(_nBus, 0);
	float eps = 0.01;
	for (int b = 0; b < _nBus; b++) {
		float V = upperBound.get(b + _nBus, 0);
		UpperBound.set(b + 2 * _nAgent, 0, V * V);
		V = lowerBound.get(b + _nBus, 0);
		LowerBound.set(b + 2 * _nAgent, 0, V * V);
		if (b == 0) {
			UpperBound.increment(2 * _nAgent, 0, eps);
			LowerBound.increment(2 * _nAgent, 0, -eps);
		}
	}
	for (int i = 0; i < _nLine; i++) {
		float pij = upperBound.get(2 * _nBus + i, 0);
		UpperBound.set(_nBus + 2 * _nAgent + i, 0, pij);
		pij = lowerBound.get(2 * _nBus + i, 0);
		LowerBound.set(_nBus + 2 * _nAgent + i, 0, pij);
	}
	

	Cost1 = MatrixCPU(2 * _nAgent, 1, rhoSO); //must be the consensus rhoPSO, pas important pour l'instant
	Cost2 = MatrixCPU(2 * _nAgent, 1);
	
	Pn = MatrixCPU(2 * _nAgent, 1); // not the real agent

	if (_initWithPF) {
		initWithPF(cas);
	}
	else {
		E = MatrixCPU(2 * _nBus, 1);
		for (int i = 0; i < _nBus; i++) {
			E.set(i + _nBus, 0, 1);
		}
	}


	dX = MatrixCPU(_sizeVar, 1);
	dY = MatrixCPU(_sizeEq, 1);
	dl = MatrixCPU(_sizeInEq, 1);
	du = MatrixCPU(_sizeInEq, 1);
	dz = MatrixCPU(_sizeInEq, 1);
	dw = MatrixCPU(_sizeInEq, 1);

	Lx = MatrixCPU(_sizeVar, 1);
	Ly = MatrixCPU(_sizeEq, 1);
	Ll = MatrixCPU(_sizeInEq, 1);
	Lu = MatrixCPU(_sizeInEq, 1);
	Lz = MatrixCPU(_sizeInEq, 1);
	Lw = MatrixCPU(_sizeInEq, 1);

	//std::cout << " creation " << std::endl;
	X = MatrixCPU(_sizeVar, 1);
	Xpre = MatrixCPU(_sizeVar, 1);


	Y = MatrixCPU(_sizeEq, 1);
	l = MatrixCPU(_sizeInEq, 1, valMin);
	u = MatrixCPU(_sizeInEq, 1, valMin);
	z = MatrixCPU(_sizeInEq, 1, valMin);
	w = MatrixCPU(_sizeInEq, 1, -valMin);
	d = MatrixCPU(_sizeInEq, 1);

	dXY = MatrixCPU(_sizeVar + _sizeEq, 1);
	PSI = MatrixCPU(_sizeVar, 1);

	tempS1 = MatrixCPU(_sizeInEq, 1);
	tempM1 = MatrixCPU(_sizeVar, 1);
	tempS1bis = MatrixCPU(_sizeInEq, 1);
	tempN2 = MatrixCPU(2 * (_nAgent + 1), 1);

	XHess = MatrixCPU(_sizeVar, _sizeVar);
	hHess = MatrixCPU(_sizeVar, _sizeVar);
	gHess = MatrixCPU(_sizeVar, _sizeVar);
	Hess = MatrixCPU(_sizeVar, _sizeVar);

	XJac = MatrixCPU(_sizeVar, 1);
	hJac = MatrixCPU(_sizeVar, _sizeEq);
	gJac = MatrixCPU(_sizeVar, _sizeInEq);

	PnBus = MatrixCPU(_nBus * 2, 1);

	MatSys = MatrixCPU(_sizeVar + _sizeEq, _sizeVar + _sizeEq);
	VectSys = MatrixCPU(_sizeVar + _sizeEq, 1);



	for (int n = 0; n < _nAgent; n++) {
		X.set(n, 0, Pn.get(n + 1, 0));
		X.set(n + _nAgent, 0, Pn.get(n + 2 + _nAgent, 0));

		XHess.set(n, n, -Cost1.get(n + 1, 0));
		XHess.set(n + _nAgent, n + _nAgent, -Cost1.get(n + 2 + _nAgent, 0));
		l.set(n, 0, min(MAX((X.get(n, 0) - LowerBound.get(n, 0)), valMin), valMax));
		u.set(n, 0, min(MAX((UpperBound.get(n, 0) - X.get(n, 0)), valMin), valMax));
		l.set(n + _nAgent, 0, min(MAX((X.get(n + _nAgent, 0) - LowerBound.get(n + _nAgent, 0)), valMin), valMax));
		u.set(n + _nAgent, 0, min(MAX((UpperBound.get(n + _nAgent, 0) - X.get(n + _nAgent, 0)), valMin), valMax));
	}
	int offset = 2 * _nAgent;
	int indice = 0;
	for (int b = 0; b < _nBus; b++) {
		//float V = _V0*_V0;
		//X.set(b + offset, 0, _V0);// e =1, f = 0;
		float V = E.get(b + _nBus, 0) * E.get(b + _nBus, 0);
		float Vma = UpperBound.get(b + offset, 0);
		float Vmin = LowerBound.get(b + offset, 0);
		V = (Vma - V) * (V > Vma) + (Vmin - V) * (V < Vmin) + V;
		float ei = sqrt(V) * cos(E.get(b, 0));
		float fi = sqrt(V) * sin(E.get(b, 0));
		//std::cout << V << " " << Vmin << " " << Vma << " " << E.get(b, 0) << " " << cos(E.get(b, 0)) << " " << sin(E.get(b, 0)) << " " << ei << " " << fi << std::endl;
		X.set(b + offset, 0, ei);
		X.set(b + _nBus + offset, 0, fi);
		int k = CoresBusLin.get(b, 0);
		l.set(b + offset, 0, min(MAX((V - Vmin), valMin), valMax));
		u.set(b + offset, 0, min(MAX((Vma - V), valMin), valMax));
		for (int voisin = k + 1; voisin < (k + nLines.get(b, 0)); voisin++) {
			int j = CoresVoiLin.get(voisin, 0);
			if (j > b) {
				float ej = E.get(j + _nBus, 0) * cos(E.get(j, 0));
				float fj = E.get(j + _nBus, 0) * sin(E.get(j, 0));
				float B = BgridLin.get(voisin, 0);
				float G = GgridLin.get(voisin, 0);
				float Pij = (ei * ei + fi * fi - ei * ej - fi * fj) * G + (ei * fj - ej * fi) * B;
				l.set(_nBus + offset + indice, 0, min(MAX((Pij - LowerBound.get(_nBus + offset + indice, 0)), valMin), valMax));
				u.set(_nBus + offset + indice, 0, min(MAX((UpperBound.get(_nBus + offset + indice, 0) - Pij), valMin), valMax));
				indice++;
			}
		}
	}

	//std::cout << " X " << std::endl;
	//X.display();





	Xpre.set(&X);
	updateRes(0);
	/*updateJ();
	updateH();
	computeConstraint();
	setSystem();*/

	std::cout << "-------" << std::endl;
	//GgridLin.display();
	//BgridLin.display();
	std::cout << "-------" << std::endl;
	//CoresBusLin.display();
	//CoresVoiLin.display();
	//nLines.display();
	//std::cout << "rho " << _rhog << " rhoL " << _rhol << " rho1 " << _rho1 << std::endl;
	std::cout << "fin init temps : " << (float)(clock() - t) / CLOCKS_PER_SEC << std::endl;
	std::cout << "---------------------------------------------------------------------------------------" << std::endl;
}

void OPFPDIPM::updateConsensus(MatrixCPU* Pmarket)
{
	for (int n = 1; n < _nAgent; n++) { // pas l'agent des pertes
		float eta = 0.5 * (Pn.get(n, 0) - Pmarket->get(n, 0));
		etaSO.set(n, 0, etaSO.get(n, 0) + eta);
		eta = 0.5 * (Pn.get(n + _nAgent, 0) - Pmarket->get(n + _nAgent, 0));
		etaSO.set(n + _nAgent, 0, etaSO.get(n + _nAgent, 0) + eta);
	}
	
	
	
	Cost2.add(&Pn, Pmarket);
	Cost2.set(0, 0, 0);
	Cost2.set(_nAgent, 0, 0);
	Cost2.multiply(-0.5);
	Cost2.add(&etaSO);
	Cost2.multiply(_rhoSO);
	

	/*for (int n = 0; n < _nAgent; n++) {
		XHess.set(n, n, -Cost1.get(n + 1, 0));
		XHess.set(n + _nAgent, n + _nAgent, -Cost1.get(n + 2 + _nAgent, 0));
	}*/ // ne change pas
}



float OPFPDIPM::getPLoss()
{
	float Ploss = 0;
	for (int i = 0; i < _nAgent; i++) {
		Ploss -= X.get(i, 0);
	}
	return Ploss;
}

float OPFPDIPM::getQLoss()
{
	float Qloss = 0;
	for (int i = 0; i < _nAgent; i++) {
		Qloss -= X.get(i + _nAgent, 0);
	}
	return Qloss;
}

void OPFPDIPM::updatePerturbedFactor()
{
	_mu = _sigma * Cgap / (2 * _sizeInEq);
}

void OPFPDIPM::correctionEquation()
{
	#ifdef Eigen
		dXY.solveSysEigen(&MatSys, &VectSys);
	#else
		MatrixCPU AD(_sizeVar + _sizeEq, _sizeVar + _sizeEq);
		MatrixCPU PD(_sizeVar + _sizeEq + 1, 1);
		try
		{
			//AD.invertGaussJordan(&MatSys); //inversion matrice Jac
			
			MatSys.LUPFactorization(&AD, &PD);
		}
		catch (const std::exception& e)
		{
			//AD.display();
			std::cout << "error " << e.what() << std::endl;
			exit(-1);
		}

	//dXY.MultiplyMatVec(&AD, &VectSys);
	dXY.solveSys(&AD, &PD, &VectSys);
	#endif
	/**/
	

	for (int i = 0; i < _sizeVar; i++) {
		dX.set(i, 0, dXY.get(i, 0));
	}
	for (int i = 0; i < _sizeEq; i++) {
		dY.set(i, 0, dXY.get(i + _sizeVar, 0));
	}
	tempS1bis.MultiplyMatTransVec(&gJac, &dX);

	// dl
	
	dl.add(&Lz, &tempS1bis);

	// du 
	du.add(&Lw, &tempS1bis);
	du.multiply(-1);

	// dz
	tempS1.multiplyT(&z, &Lz);
	tempS1.add(&Ll);
	tempS1.divideT(&l);

	dz.multiplyT(&z, &tempS1bis);
	dz.divideT(&l);
	dz.add(&tempS1);
	dz.multiply(-1);

	// dw
	tempS1.multiplyT(&w, &Lw);
	tempS1.subtract(&Lu);
	tempS1.divideT(&u);

	dw.multiplyT(&w, &tempS1bis);
	dw.divideT(&u);
	dw.add(&tempS1);



}

void OPFPDIPM::updateStep()
{
	float minP = 1;
	float minD = 1;
	float a, b;
	for (int i = 0; i < _sizeInEq; i++) {
		a = dl.get(i, 0);
		if ( a < 0) {
			b = - l.get(i, 0) / a;
			if ( b < minP) {
				minP = b;
			}
		}
		a = du.get(i, 0);
		if (a < 0) {
			b = -u.get(i, 0) / a;
			if (b < minP) {
				minP = b;
			}
		}
		a = dz.get(i, 0);
		if (a < 0) {
			b = - z.get(i, 0) / a;
			if (b < minD) {
				minD = b;
			}
		}
		a = dw.get(i, 0);
		if (a > 0) {
			b = - w.get(i, 0) / a;// dans le papier il n'y a pas de -, mais il faut bien que cela soit positif...
			if (b < minD) {
				minD = b;
			}
		}
	}
	_stepD = 0.9995 * minD;
	_stepP = 0.9995 * minP;


}

void OPFPDIPM::updateVariable()
{
	dX.multiply(_stepP);
	dl.multiply(_stepP);
	du.multiply(_stepP);

	dY.multiply(_stepD);
	dz.multiply(_stepD);
	dw.multiply(_stepD);

	Xpre.set(&X);

	X.add(&dX);
	//X.set(2 *_nAgent + _nBus, 0, 0); // f0 = 0
	l.add(&dl);
	u.add(&du);
	Y.add(&dY);
	z.add(&dz);
	w.add(&dw);

}

bool OPFPDIPM::setCenteringParameter()
{
	double Caff = 0;
	for (int i = 0; i < _sizeInEq; i++) {
		Caff += (l.get(i, 0) + _stepP * dl.get(i, 0)) * (z.get(i, 0) + _stepD * dz.get(i, 0)) - (u.get(i, 0) + _stepP * du.get(i, 0)) * (w.get(i, 0) + _stepD * dw.get(i, 0));
	}
	
	_sigma = (Caff / Cgap);
	_sigma = _sigma * _sigma * _sigma;

	_sigma = -_sigma * (_sigma < 0) + (1 - _sigma) * (_sigma > 1) + _sigma;
	//std::cout << " Caff " << Caff << " Cgap " << Cgap << " sigma " << _sigma << std::endl;

	if (Cgap < Caff) {
		return false; // fail ?
	}
	else {
		return true; // sucess ???
	}


}

void OPFPDIPM::updateJ()
{
	hJac.set(0.0);
	for (int n = 0; n < _nAgent; n++) {
		XJac.set(n, 0, Cost1.get(n + 1, 0) * X.get(n, 0) + Cost2.get(n + 1, 0));
		XJac.set(n + _nAgent, 0, Cost1.get(n + 2 + _nAgent, 0) * X.get(n + _nAgent, 0) + Cost2.get(n + 2 + _nAgent, 0));
		gJac.set(n, n, 1);
		gJac.set(n + _nAgent, n + _nAgent, 1);
		
		int bus = _CoresBusAgent.get(n + 1, 0);
		//dh1/dp
		hJac.set(n, bus, -1);
		//dh2/dq
		hJac.set(n + _nAgent, _nBus + bus, -1);
	}
	int offset = 2 * _nAgent;
	int indice = 0;
	for (int b = 0; b < _nBus; b++) {
		int k = CoresBusLin.get(b, 0);
		float ei = X.get(offset + b, 0);
		float fi = X.get(offset + _nBus + b, 0);
		
		gJac.set(offset + b, offset + b, 2 * ei);
		gJac.set(offset + _nBus + b, offset + b, 2 *fi);
		
		

		//dh1/dei partiel
		hJac.increment(offset + b, b, 2 * GgridLin.get(k, 0) * ei);
		//dh1/dfi partiel
		hJac.increment(offset + _nBus + b, b, 2 * GgridLin.get(k, 0) * fi);
		//dh2/dei partiel
		hJac.increment(offset + b, b + _nBus, -2 * BgridLin.get(k, 0) * ei);
		//dh2/dfi partiel
		hJac.increment(offset + _nBus + b, b + _nBus, -2 * BgridLin.get(k, 0) * fi);

		for (int voisin = k + 1; voisin < (k + nLines.get(b, 0)); voisin++) {
			int j = CoresVoiLin.get(voisin, 0); // pour les �galit�s -> h c'est tous les bus , pour les in�galit�s -> g on consid�re les lignes
			float ej = X.get(offset + j, 0);
			float fj = X.get(offset + _nBus + j, 0);
			float B = BgridLin.get(voisin, 0);
			float G = GgridLin.get(voisin, 0);
			// ei
			float h = ej * G - B * fj;//h1
			hJac.increment(offset + b, b, h);
			h = -ej * B - fj * G; // h2
			hJac.increment(offset + b, _nBus + b, h);
			h = ei * G - B * fi;
			hJac.increment(offset + j, j, h);
			h = -ei * B - G * fi;
			hJac.increment(offset + j, j + _nBus, h);

			// ej
			h = ei * G + fi * B; // h1
			hJac.set(offset + j, b, h);
			h = -B * ei + G * fi; // h2
			hJac.set(offset + j, _nBus + b, h);

			h = ej * G + fj * B;
			hJac.set(offset + b, j, h);
			h = -B * ej + G * fj;
			hJac.set(offset + b, _nBus + j, h);

			// fi
			h = ej * B + G * fj;
			hJac.increment(offset + _nBus + b, b, h);
			h = ej * G - fj * B;
			hJac.increment(offset + _nBus + b, _nBus + b, h);
			h = ei * B + G * fi;
			hJac.increment(offset + _nBus + j, j, h);
			h = ei * G - B * fi;
			hJac.increment(offset + _nBus + j, j + _nBus, h);

			// fj
			h = -ei * B + fi * G;
			hJac.set(offset + _nBus + j, b, h);
			h = -ei * G - fi * B;
			hJac.set(offset + _nBus + j, _nBus + b, h);

			h = -ej * B + fj * G;
			hJac.set(offset + _nBus + b, j, h);
			h = -ej * G - fj * B;
			hJac.set(offset + _nBus + b, _nBus + j, h);

			if (j > b) {
				// ei
				float g = (2 * ei - ej) * G + fj * B;
				gJac.set(offset + b, offset + _nBus + indice, g);
			
				// ej
				g = - ei* G - fi * B;
				gJac.set(offset + j, offset + _nBus + indice, g);

				// fi
				g = (2 * fi - fj) * G - ej * B;
				gJac.set(offset + _nBus + b, offset + _nBus + indice, g);

				// fj
				g = - fj * G + ei * B;
				gJac.set(offset  + _nBus + j, offset + _nBus + indice, g);

				indice++;
			}
		}
	}


	// bus slack
	hJac.set(offset, 2 * _nBus, 1);
	hJac.set(offset + _nBus, 2 * _nBus + 1, 1);


}

void OPFPDIPM::updateH()
{
	tempS1.set(&z);
	tempS1.divideT(&l);
	d.set(&w);
	d.divideT(&u);
	d.subtract(&tempS1);
	//std::cout << " D " << std::endl;
	//d.display();

	Hess.set(&XHess);

	for (int n = 0; n < 2 * _nAgent; n++) {
		Hess.increment(n, n, d.get(n, 0));
	}



	int offset = 2 * _nAgent;
	int indice = 0;
	for (int b = 0; b < _nBus; b++) {
		int k = CoresBusLin.get(b, 0);
		float y1 = Y.get(b, 0);
		float y2 = Y.get(b + _nBus, 0);
		float d2 = d.get(b + offset, 0);
		float zw2 = z.get(offset + b, 0) + w.get(offset + b, 0);
		double gii = 0; // double car reduction � faire...
		float gdgii  = 0;

		
		// ei
		float hii = 2 * (GgridLin.get(k, 0) * y1 - BgridLin.get(k, 0) * y2);
		hHess.set(offset + b, offset + b, hii);
		
		gdgii = gJac.get(offset + b, offset + b) * d2 * gJac.get(offset + b, offset + b);
		Hess.increment(offset + b, offset + b, gdgii);
		//fi
		hii = 2 * (GgridLin.get(k, 0) * y1 - BgridLin.get(k, 0) * y2);
		hHess.set(offset + b + _nBus, offset + b + _nBus, hii);

		gdgii = gJac.get(offset + b + _nBus, offset + b) * d2 * gJac.get(offset + _nBus + b, offset + b);
		Hess.increment(offset + b, offset + b, gdgii);

		// ei fi 
		gdgii = gJac.get(offset + b, offset + b) * d2 * gJac.get(offset + _nBus + b, offset + b);
		Hess.increment(offset + b, offset + _nBus + b, gdgii);
		Hess.increment(offset + _nBus + b, offset + b, gdgii);

		for (int voisin = k + 1; voisin < (k + nLines.get(b, 0)); voisin++) {
			int j = CoresVoiLin.get(voisin, 0);
			

			float hij = GgridLin.get(voisin, 0) * y1 - BgridLin.get(voisin, 0) * y2;
			hHess.set(offset + b, offset + j, hij);
			hHess.set(offset + j, offset + b, hij);



			if (j > b) { // pour ne pas compter 2 fois chaque ligne 
				float zw3 = z.get(offset + _nBus + indice, 0) + w.get(offset + _nBus + indice, 0);
				float d3 = d.get(offset + _nBus + indice, 0);
				// ei ej
				

				float gij = -GgridLin.get(voisin, 0) * zw3;
				gHess.set(offset + b, offset + j, gij);
				gHess.set(offset + j, offset + b, gij);

				float gdgij = gJac.get(offset + b, offset + _nBus + indice) * d3 * gJac.get(offset + j, offset + _nBus + indice);
				Hess.increment(offset + b, offset + j, gdgij);
				Hess.increment(offset + j, offset + b, gdgij);

				// fi fj 
				hij = GgridLin.get(voisin, 0) * y1 - BgridLin.get(voisin, 0) * y2;
				hHess.set(offset + b + _nBus, offset + j + _nBus, hij);
				hHess.set(offset + j + _nBus, offset + b + _nBus, hij);

				gij = -GgridLin.get(voisin, 0) * zw3;
				gHess.set(offset + b + _nBus, offset + j + _nBus, gij);
				gHess.set(offset + j + _nBus, offset + b + _nBus, gij);

				gdgij = gJac.get(offset + b + _nBus, offset + _nBus + indice) * d3 * gJac.get(offset + j + _nBus, offset + _nBus + indice);
				Hess.increment(offset + b + _nBus, offset + j + _nBus, gdgij);
				Hess.increment(offset + j + _nBus, offset + b + _nBus, gdgij);

				// ei fj 
				hij = - BgridLin.get(voisin, 0) * y1 - GgridLin.get(voisin, 0) * y2;
				hHess.set(offset + b, offset + j + _nBus, hij);
				hHess.set(offset + j + _nBus, offset + b, hij);

				gij = BgridLin.get(voisin, 0) * zw3;
				gHess.set(offset + b, offset + j + _nBus, gij);
				gHess.set(offset + j + _nBus, offset + b, gij);

				gdgij = gJac.get(offset + b, offset + _nBus + indice) * d3 * gJac.get(offset + j + _nBus, offset + _nBus + indice);
				Hess.increment(offset + b, offset + j + _nBus, gdgij);
				Hess.increment(offset + j + _nBus, offset + b, gdgij);

				// fi ej
				hij = BgridLin.get(voisin, 0) * y1 - GgridLin.get(voisin, 0) * y2;
				hHess.set(offset + b + _nBus, offset + j, hij);
				hHess.set(offset + j, offset + b + _nBus, hij);

				gij = -BgridLin.get(voisin, 0) * zw3;
				gHess.set(offset + b + _nBus, offset + j, gij);
				gHess.set(offset + j, offset + b + _nBus, gij);
				
				Hess.increment(offset + b + _nBus, offset + j, gdgij);
				Hess.increment(offset + j, offset + b + _nBus, gdgij);
				/// ei ei or fi fi 
			
				gii += 2 * GgridLin.get(voisin, 0) * zw3;
				gdgii = gJac.get(offset + b, offset + _nBus + indice) * d3 * gJac.get(offset + b, offset + _nBus + indice);
				Hess.increment(offset + b, offset + b, gdgii);
				gdgii = gJac.get(offset + j, offset + _nBus + indice) * d3 * gJac.get(offset + j, offset + _nBus + indice);
				Hess.increment(offset + j, offset + j, gdgii);

				gdgii = gJac.get(offset + _nBus + b, offset + _nBus + indice) * d3 * gJac.get(offset + _nBus + b, offset + _nBus + indice);
				Hess.increment(offset + _nBus + b, offset + _nBus + b, gdgii);
				gdgii = gJac.get(offset + _nBus + j, offset + _nBus + indice) * d3 * gJac.get(offset + _nBus + j, offset + _nBus + indice);
				Hess.increment(offset + _nBus + j, offset + _nBus + j, gdgii);

				gdgii = gJac.get(offset + b, offset + _nBus + indice) * d3 * gJac.get(offset + _nBus + b, offset + _nBus + indice);
				Hess.increment(offset + b, offset + _nBus + b, gdgii);
				Hess.increment(offset + _nBus + b, offset + b, gdgii);
				gdgii = gJac.get(offset + j, offset + _nBus + indice) * d3 * gJac.get(offset + _nBus + j, offset + _nBus + indice);
				Hess.increment(offset + j, offset + _nBus + j, gdgii);
				Hess.increment(offset + _nBus + j, offset + j, gdgii);

				indice++;
			}
			gii += 2 * zw2;
			gHess.set(offset + b, offset + b, gii);
			gHess.set(offset + b + _nBus, offset + b + _nBus, gii);

		}

	}

	Hess.add(&hHess);
	Hess.add(&gHess);


}

void OPFPDIPM::setSystem()
{
	MatSys.setBloc(0, _sizeVar, 0, _sizeVar, &Hess);
	MatSys.setBloc(0, _sizeVar, _sizeVar, _sizeVar + _sizeEq, &hJac);
	for (int i = 0; i < _sizeEq; i++) {
		for (int j = 0; j < _sizeVar; j++) {
			MatSys.set(i + _sizeVar, j, hJac.get(j, i));
		}
	}
	VectSys.setBloc(0, _sizeVar, 0, 1, &PSI);
	VectSys.setBloc(_sizeVar, _sizeVar + _sizeEq, 0, 1, &Ly);
	VectSys.multiply(-1);

}

void OPFPDIPM::computeConstraint()
{
	PnBus.set(0.0);
	Ly.set(0.0);

	for (int n = 0; n < _nAgent; n++) {
		Lz.set(n, 0, X.get(n, 0) - l.get(n, 0) - LowerBound.get(n, 0));
		Lw.set(n, 0, X.get(n, 0) + u.get(n, 0) - UpperBound.get(n, 0));
		Lz.set(n + _nAgent, 0, X.get(n + _nAgent, 0) - l.get(n + _nAgent, 0) - LowerBound.get(n + _nAgent, 0));
		Lw.set(n + _nAgent, 0, X.get(n + _nAgent, 0) + u.get(n + _nAgent, 0) - UpperBound.get(n + _nAgent, 0));

		int bus = _CoresBusAgent.get(n + 1, 0);
		PnBus.increment(bus, 0, X.get(n, 0));
		PnBus.increment(bus + _nBus, 0, X.get(n + _nAgent, 0));

	}
	int offset = 2 * _nAgent;
	int indice = 0;
	for (int b = 0; b < _nBus; b++) {
		int k = CoresBusLin.get(b, 0);
		Ly.increment(b, 0, -PnBus.get(b, 0));
		Ly.increment(b + _nBus, 0, -PnBus.get(b + _nBus, 0));
		float ei = X.get(offset + b, 0);
		float fi = X.get(offset + _nBus + b, 0);
		float V = ei * ei + fi * fi;
		Lz.set(offset + b, 0, V - l.get(offset + b, 0) - LowerBound.get(offset + b, 0));
		Lw.set(offset + b, 0, V + u.get(offset + b, 0) - UpperBound.get(offset + b, 0));

		Ly.increment(b, 0, GgridLin.get(k, 0) * V);
		Ly.increment(b + _nBus, 0, - BgridLin.get(k, 0) * V);
		for (int voisin = k + 1; voisin < (k + nLines.get(b, 0)); voisin++) {
			int j = CoresVoiLin.get(voisin, 0);
			if (j > b) {
				float ej = X.get(offset + j, 0);
				float fj = X.get(offset + _nBus + j, 0);
				float B = BgridLin.get(voisin, 0);
				float G = GgridLin.get(voisin, 0);
				float Pij = (ei * ei + fi * fi - ei * ej - fi * fj) * G + (ei * fj - ej * fi)* B;
				
				float hij = ei * (ej * G - fj * B) + fi * (fj * G + ej * B);
				float hji = ej * (ei * G - fi * B) + fj * (fi * G + ei * B);

				Ly.increment(b, 0, hij);
				Ly.increment(j, 0, hji);

				hij = fi * (ej * G - fj * B) - ei * (fj * G + ej * B);
				hji = fj * (ei * G - fi * B) - ej * (fi * G + ei * B);

				Ly.increment(b + _nBus, 0, hij);
				Ly.increment(j + _nBus, 0, hji);

				//std::cout << " Pij " << Pij << std::endl;

				Lz.set(offset + _nBus + indice, 0, Pij - l.get(offset + _nBus + indice, 0) - LowerBound.get(offset + _nBus + indice, 0));
				Lw.set(offset + _nBus + indice, 0, Pij + u.get(offset + _nBus + indice, 0) - UpperBound.get(offset + _nBus + indice, 0));
				
				indice++;
			}
		}
	}
	Ly.set(2 * _nBus, 0,  X.get(offset, 0) - _V0);
	Ly.set(2 * _nBus + 1, 0, X.get(offset + _nBus, 0));

	Ll.multiplyT(&l, &z);
	Lu.multiplyT(&u, &w);

	PSI.MultiplyMatVec(&hJac, &Y);
	PSI.subtract(&XJac);
	
	tempS1.add(&z, &w);
	tempM1.MultiplyMatVec(&gJac, &tempS1);
	Lx.add(&PSI, &tempM1);
	Lx.multiply(-1);

	tempS1.multiplyT(&w, &Lw);
	tempS1.divideT(&u);

	tempS1bis.multiplyT(&z, &Lz);
	tempS1bis.divideT(&l);

	tempS1.subtract(&tempS1bis);
	

	tempM1.MultiplyMatVec(&gJac, &tempS1);
	PSI.add(&tempM1);

	/*std::cout << "Lx " << std::endl;
	Lx.display();
	std::cout << "Ly " << std::endl;
	Ly.display();
	std::cout << "Lz " << std::endl;
	Lz.display();
	std::cout << "Lw " << std::endl;
	Lw.display();
	std::cout << "Ll " << std::endl;
	Ll.display();
	std::cout << "Lu " << std::endl;
	Lu.display();
	std::cout << "PSI 0 " << std::endl;
	PSI.display();*/


}

void OPFPDIPM::updatePSILlu()
{
	tempS1.multiplyT(&dw, &du);
	tempS1.divideT(&u);

	tempS1bis.multiplyT(&dz, &dl);
	tempS1bis.divideT(&l);

	tempS1bis.subtract(&tempS1);

	for (int i = 0; i < _sizeInEq; i++) {
		float t = tempS1bis.get(i, 0) * _stepP * _stepD;
		float t2 = _mu * (1/u.get(i, 0) - 1/l.get(i, 0));
		tempS1bis.set(i, 0, t - t2);

	}

	tempM1.MultiplyMatVec(&gJac, &tempS1bis);

	for (int i = 0; i < _sizeVar; i++)
	{
		float t = VectSys.get(i, 0);
		VectSys.set(i, 0, t - tempM1.get(i, 0)); // c'est un - car dans le vecteur on a -PSI

	}

	
	Ll.add(-_mu);
	Lu.add(_mu);



}


float OPFPDIPM::updateRes(int indice) 
{
	float resS = 0;
	float resR = 0;
	float resV = 0;
	
	for (int i = 0; i < _sizeInEq; i++) {
		resV += l.get(i, 0) * z.get(i, 0) - u.get(i, 0) * w.get(i, 0);
	}

	Cgap = resV;


	if (indice % _stepG == 0) {

		resS = X.max2(&Xpre);
		resR = Ly.max2();

		resF.set(0, indice / _stepG, resR);
		resF.set(1, indice / _stepG, resS);
		resF.set(2, indice / _stepG, resV);

	}

	return MAX(MAX(resV, resS), resR) * (resV > epsMin);
}




void OPFPDIPM::display() {

	std::cout.precision(3);

	if (_iterGlobal == 0) {
		std::cout << "algorithm not launch" << std::endl;
	}
	else if (_iterGlobal < _iterG) {
		std::cout << "method " << _name << " converged in " << _iterGlobal << " iterations." << std::endl;
		std::cout << "Converged in " << (float)timeOPF / CLOCKS_PER_SEC << " seconds" << std::endl;

	}
	else {
		std::cout << "method " << _name << " not converged in " << _iterGlobal << " iterations." << std::endl;
		std::cout << "time taken " << (float)timeOPF / CLOCKS_PER_SEC << " seconds" << std::endl;
	}
	std::cout << "The power error of this state is (constraint) " << resF.get(0, _iterGlobal / _stepG) << " and convergence " << resF.get(1, _iterGlobal / _stepG) << std::endl;
	std::cout << "===============================================================|" << std::endl;
	std::cout << "      System Summary                                           |" << std::endl;
	std::cout << "===============================================================|" << std::endl;
	std::cout << "Buses            " << _nBus << std::endl;
	std::cout << "Branches         " << _nLine << std::endl;
	std::cout << "Agent            " << _nAgent << std::endl;
	std::cout << "Ploss            " << getPLoss() << std::endl;
	std::cout << "Qloss            " << getQLoss() << std::endl;


	std::cout << std::endl << std::endl;

	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "      Bus Data                                                                                          |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Bus  |      Voltage    |      Power = Generation  + Load       |                  H(X)                 |" << std::endl;
	std::cout << "  #   |      Mag(pu)    |      P (pu)       |        Q (pu )    |       P (pu)      |        Q (pu)     |" << std::endl;
	std::cout << "------|-----------------|-------------------|-------------------|-------------------|-------------------|" << std::endl;


	float seuil = 0.0001;
	int offset = 2 * _nAgent;
	for (int b = 0; b < _nBus; b++) {
		float V = X.get(offset + b, 0) * X.get(offset + b, 0) + X.get(offset + _nBus + b, 0) * X.get(offset + _nBus + b, 0);
		std::cout << std::setw(6) << b << "|" << std::setw(16) << V << " |" << std::setw(19)
			<< PnBus.get(b, 0) << "|" << std::setw(19) << PnBus.get(b + _nBus, 0)
			<< "|" << std::setw(19) << Ly.get(b, 0) << "|" << std::setw(19)
			<< Ly.get(b + _nBus, 0) << "|" << std::endl;

	}
	std::cout << std::endl << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "      Line Data                                                                                         |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Line |    From     |    To      |                            flow                                      |" << std::endl;
	std::cout << "  #   |    Bus      |    Bus     |    P (pu)      |    Q (pu)      |     l (pu)     |     Loss (pu)     |" << std::endl;
	std::cout << "------|-------------|------------|----------------|----------------|----------------|-------------------|" << std::endl;

	int l = 0;
	for (int b = 0; b < _nBus; b++) {
		int k = CoresBusLin.get(b, 0);
		float ei = X.get(offset + b, 0);
		float fi = X.get(offset + _nBus + b, 0);


		for (int voisin = k + 1; voisin < (k + nLines.get(b, 0)); voisin++) {
			int j = CoresVoiLin.get(voisin, 0);
			if (j > b) {

				float ej = X.get(offset + j, 0);
				float fj = X.get(offset + _nBus + j, 0);
				float pij = ei * (ej * GgridLin.get(voisin, 0) - fj * BgridLin.get(voisin, 0)) + fi * (fj * GgridLin.get(voisin, 0) + ej * BgridLin.get(voisin, 0));
				float qij = fi * (ej * GgridLin.get(voisin, 0) - fj * BgridLin.get(voisin, 0)) - ei * (fj * GgridLin.get(voisin, 0) + ej * BgridLin.get(voisin, 0));

				float lij = sqrt((pij * pij + qij * qij) / (ei * ei + fi * fi));
				float lossij = GgridLin.get(voisin, 0) * ((fi - fj) * (fi - fj) + (ei - ej) * (ei - ej));

				std::cout << std::setw(6) << l << "|" << std::setw(12) << b << " |" << std::setw(12)
					<< j << "|" << std::setw(16) << pij
					<< "|" << std::setw(16) << qij << "|" << std::setw(16)
					<< lij << "|" << std::setw(19) << lossij << "|" << std::endl;
				l++;
			}
		}
		
	}
	
	std::cout << std::endl << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "     Constraints                                                                                        |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Bus | Voltage | Voltage | Voltage |        Power Injection          |          Power Injection         |" << std::endl;
	std::cout << "  #  | Mag(pu) | MIN(pu) |  MAX(pu)|  P (pu) | Pmin (pu) | Pmax (pu) |  Q (pu)  | Qmin (pu) | Qmax (pu) |" << std::endl;
	std::cout << "-----|---------|---------|---------|---------|-----------|-----------|----------|-----------|-----------|" << std::endl;
	

	for (int b = 0; b < _nBus; b++) {
		float ei = X.get(offset + b, 0);
		float fi = X.get(offset + _nBus + b, 0);
		float Vi = ei * ei + fi * fi;
		std::cout << std::setw(5) << b << "|" << std::setw(8) << sqrt(Vi) << " |" << std::setw(9)
			<< LowerBound.get(offset + b, 0) << "|" << std::setw(9) << UpperBound.get(offset + b, 0)
			<< "|" << std::setw(9) << PnBus.get(b,0) << "|" << std::setw(11)
			<< "***" << "|" << std::setw(11) << "***" << "|" << std::setw(10) << PnBus.get(b + _nBus, 0)
			<< "|" << std::setw(11) << "***" << "|" << std::setw(11) << "***" << "|" << std::endl;

	}
	std::cout << std::endl << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "      Agent Data                                                                                        |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Agent |  Bus  |  Cost   |  Cost   |        Power Injection          |          Power Injection         |" << std::endl;
	std::cout << "  #    |   #   |  a (pu) |  b (pu) |  P (pu) | Pmin (pu) | Pmax (pu) |  Q (pu)  | Qmin (pu) | Qmax (pu) |" << std::endl;
	std::cout << "-------|-------|---------|---------|---------|-----------|-----------|----------|-----------|-----------|" << std::endl;

	for (int n = 0; n < _nAgent + 1; n++) {
		int b = _CoresBusAgent.get(n, 0);
		std::cout << std::setw(7) << n << "|" << std::setw(7) << b << "|" << std::setw(8) << Cost1.get(n,0) << " |" << std::setw(9)
			<< Cost2.get(n, 0) << "|" << std::setw(9) << Pn.get(n,0) << "|" << std::setw(11)
			<< Pmin.get(n, 0) << "|" << std::setw(11) << Pmax.get(n, 0) << "|" << std::setw(10) << Pn.get(n + _nAgent, 0)
			<< "|" << std::setw(11) << Pmin.get(n + _nAgent + 1, 0) << "|" << std::setw(11) << Pmax.get(n + _nAgent + 1, 0) << "|" << std::endl;
	}


	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "                      END PRINT                                                                         |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;

}

void OPFPDIPM::initWithMarket(const Simparam& sim, const StudyCase& cas)
{
	std::cout << "init Market " << std::endl;
	paramMarketInit.setNAgentLine(_nAgent + 1, 0, true); // il faut remettre l'agent des pertes
	paramMarketInit.setItG(500);
	paramMarketInit.setEpsG(0.001);
	paramMarketInit.setEpsL(0.0001);
	ADMMMarket Market;
	Simparam res(paramMarketInit);
	Market.solveWithMinPower(&res, paramMarketInit, cas);
	Pn = res.getPn();
	//Market.display();

}

void OPFPDIPM::initWithPF(const StudyCase& cas)
{
	std::cout << "init Power Flow " << std::endl;
	if (cas.isRadial()) {
		CPUPFdistPQ PowerFlow;
		
		PowerFlow.init(cas, &Pn);
		PowerFlow.chekcase();
		PowerFlow.solve();
		PowerFlow.calcW(true);
		PowerFlow.calcE();
		//PowerFlow.display2();
		int statu = PowerFlow.getConv();
		if (statu == 1) {
			E = PowerFlow.getE();
		}
		else {
			E = MatrixCPU(2 * _nBus, 1);
			for (int i = 0; i < _nBus; i++) {
				E.set(i + _nBus, 0, 1);
			}
		}
	}
	else {
		CPUPF PowerFlow;
		MatrixCPUD PnD(2 * _nAgent, 1);
		Pn.toMatCPUD(PnD);
		PowerFlow.init(cas, &Pn, &PnD, true);
		PowerFlow.solve();
		//PowerFlow.display2(true);
		int statu = PowerFlow.getConv();
		if (statu == 1) {
			E = PowerFlow.getE();
		}
		else {
			E = MatrixCPU(2 * _nBus, 1);
			for (int i = 0; i < _nBus; i++) {
				E.set(i + _nBus, 0, 1);
			}
		}
		
		
	}
	
	
}

