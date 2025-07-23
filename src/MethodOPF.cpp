#pragma once
#include "../head/MethodOPF.h"


MethodOPF::MethodOPF() : Method()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " OPF Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = "OPF";
	timePerBlock = MatrixCPU(1, 12, 0); // Fb0, Fb11abcd, FB12, Fb2, Fb3, Fb4, Fb5,FB6, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 12, 0); //nb de fois utilis� pendant la simu
}

MethodOPF::MethodOPF(float rho) : Method()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " OPF Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = "OPF";
    _rhol = rho;
	timePerBlock = MatrixCPU(1, 12, 0); // Fb0, Fb11abcd, FB12, Fb2, Fb3, Fb4, Fb5,FB6, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 12, 0); //nb de fois utilis� pendant la simu
}

MethodOPF::~MethodOPF(){
    DELETEA(tempM1);
	DELETEA(tempM);

	DELETEA(X);
	DELETEA(Ypre);
	DELETEA(Y);
	DELETEA(Mu);

	DELETEA(Hinv);
	DELETEA(A);
	DELETEA(Q);

	DELETEA(Childs);
	DELETEA(Chat);
}
	
void MethodOPF::initSimParam(const Simparam& sim){
	
	_rho = sim.getRho();
    _rhoInv = 1 / _rho;
    if (_rhol == 0) {
		_rhol = _rho;
	}

	_iterG = sim.getIterG();
	_stepG = sim.getStepG();

    _iterL = sim.getIterL();
	_stepL = sim.getStepL();
	
	_epsG = sim.getEpsG();
	_epsL = sim.getEpsL();
	

	resF = MatrixCPU(3, (_iterG / _stepG) + 1);
	
	tempN1 = MatrixCPU(_nAgent, 1, 0);
    tempN2 = MatrixCPU(2 * _nAgent, 1);
	tempB2 = MatrixCPU(2 * _nBus, 1);

    Pn = sim.getPn();
}

void MethodOPF::initSize(const StudyCase& cas){
    _nAgent = cas.getNagent();
    _nBus = cas.getNBus();
	_nLine = cas.getNLine(true); // ne doit pas etre reduit ici !!!
    _sizeProbGlob = _nBus;
}

void MethodOPF::initGrid(const StudyCase& cas){

    CoresLineBus = cas.getCoresLineBus(true);

    nChild = MatrixCPU(_nBus, 1);
    Ancestor = MatrixCPU(_nBus, 1, 0); // A_i = bus ant�c�dent de i
	PosChild = MatrixCPU(_nBus, 1, 0); // indice du bus i dans Child[Ai]
	Ancestor.set(0, 0, -1); // the slack bus has no ancestor
	
	ZsRe = cas.getZsRe();
	ZsIm = cas.getZsImag();
	ZsNorm = MatrixCPU(_nLine, 1);
	
	if (!chekcase()) {
		throw std::invalid_argument("not a radial case");
	}

	for (int lold = 0; lold < _nLine; lold++) {
		int l = lold + 1;
		int busTo = l ;
		int busFrom = (int) CoresLineBus.get(lold, 0);
		Ancestor.set(busTo, 0, busFrom);
		nChild.increment(busFrom, 0, 1);
		ZsNorm.set(lold, 0, ZsRe.get(lold, 0) * ZsRe.get(lold, 0) + ZsIm.get(lold, 0) * ZsIm.get(lold, 0));
	}

    MatrixCPU lowerBound(cas.getLowerBound()); //voltage angle, voltage, line...
	MatrixCPU upperBound(cas.getUpperBound()); //voltage angle, voltage, line...
	
    VoltageLimit = MatrixCPU(_nBus, 2); // min, max
	VoltageLimitReal = MatrixCPU(_nBus, 2); // min, max
	if (isCurrentLimited && cas.isCurrentLimit()) {
		FluxLimit = cas.getCurrentLimit();
	}
	else {
		FluxLimit = MatrixCPU(_nBus, 1, 1000); // max
	}

    for (int i = 0; i < _nBus; i++) {
        VoltageLimitReal.set(i, 0, lowerBound.get(_nBus + i, 0));
		VoltageLimitReal.set(i, 1, upperBound.get(_nBus + i, 0));
		VoltageLimit.set(i, 0, lowerBound.get(_nBus + i, 0) * lowerBound.get(_nBus + i, 0) * sqrt((nChild.get(i, 0) + 1) / 2));
		VoltageLimit.set(i, 1, upperBound.get(_nBus + i, 0) * upperBound.get(_nBus + i, 0) * sqrt((nChild.get(i, 0) + 1) / 2));		
    }
}

void MethodOPF::initAgent(const StudyCase& cas){
    _CoresAgentBus = cas.getCoresAgentBusLin();
	_CoresAgentBusBegin = cas.getCoresAgentBusLinBegin();
	_nAgentByBus = cas.getNagentByBus();
    
    PosAgent = MatrixCPU(_nAgent, 1, -1); // indice de l'agent i dans _CoresAgentBus[CoresAgentBegin]
	_CoresBusAgent = cas.getCoresBusAgentLin(); // Cores[n] = b
	
	CoresSoloBusAgent = MatrixCPU(_nBus, 1, -1);
	
	Pmin = cas.getPmin();
	Pmax = cas.getPmax();
	Pbmax = MatrixCPU(2 * _nBus, 1);
	Pbmin = MatrixCPU(2 * _nBus, 1);
	Pb = MatrixCPU(2 * _nBus, 1);

    Cost1 = cas.geta();
	Cost2 = cas.getb();

	
	if (Pn.max2() == 0) {
		Pn.add(&Pmin, &Pmax);
		Pn.divide(2);
	}

    _nAgentByBus.increment(0, 0, -1); // don't want the grid agent in this resolution
	_CoresAgentBusBegin.increment(0, 0, 1); // idem
}


void MethodOPF::allocateTab(){
	// il faut delete pour ne pas faire de fuite memoire
	DELETEA(tempM);
    DELETEA(X);
    DELETEA(Ypre);
    DELETEA(Y);
    DELETEA(Mu);
    DELETEA(Hinv);
    DELETEA(A);
    DELETEA(Q);
    DELETEA(Childs);
    DELETEA(Chat);

    //std::cout << " creation " << std::endl;
    X = new MatrixCPU[_sizeProbGlob];
    Ypre = new MatrixCPU[_sizeProbGlob];
    Y = new MatrixCPU[_sizeProbGlob];
    Mu = new MatrixCPU[_sizeProbGlob];
    
    tempM = new MatrixCPU[_sizeProbGlob];
    Hinv = new MatrixCPU[_sizeProbGlob];
    A = new MatrixCPU[_sizeProbGlob];
    Q = new MatrixCPU[_sizeProbGlob];
    Chat = new MatrixCPU[_sizeProbGlob];
	Childs = new MatrixCPU[_sizeProbGlob];


	for (int i = 0; i< _nBus; i++){
		int sizeOPF = (int) sizeOPFADMM.get(i, 0);
        Childs[i] = MatrixCPU( (int) nChild.get(i, 0), 1);
		X[i] = MatrixCPU(sizeOPF, 1); // {Pi, Qi, vi, li, pi, qi, vAi, (Pci, Qci, lci) for all child Ci}
		Ypre[i] = MatrixCPU(sizeOPF, 1);// Y[i][j] not� dans l'article Yji est ce que i connait sur j
		Y[i] = MatrixCPU(sizeOPF, 1); //Y[i] = {Pi, Qi, vi, li, pi, qi, vAi, (Pci, Qci, lci) for all child Ci}
		Mu[i] = MatrixCPU(sizeOPF, 1);
		A[i] = MatrixCPU(2 + 1*(i>0), sizeOPF);
		Hinv[i] = MatrixCPU(sizeOPF, sizeOPF);
		Q[i] = MatrixCPU(sizeOPF, 1, 0);
		tempM[i] = MatrixCPU(sizeOPF, 1);
	}
	int indice = 0;
	MatrixCPU nChildTemp(_nBus, 1, 0);
	 for (int i = 1; i < _nBus; i++) {
		int Ai = (int) Ancestor.get(i, 0);
		int posChild = (int) nChildTemp.get(Ai, 0);
		Childs[Ai].set(posChild, 0, i);
		PosChild.set(i, 0, posChild);
		nChildTemp.increment(Ai, 0, 1);
    }

}

void MethodOPF::setParam(float rho)
{
	_rho = rho;
}


float MethodOPF::getPLoss(){
    _Ploss = 0;
    switch (losstype)
    {
    case LossType::POWER:
        for (int n = 1; n < _nAgent; n++) {
            int bus = (int) _CoresBusAgent.get(n, 0);
            int In = (int) PosAgent.get(n, 0);

            _Ploss -= X[bus].get(5 + 2 * In, 0);
        }
        break;
    case LossType::CURRENT:
        for (int i = 1; i < _nBus; i++) {
            _Ploss -= X[i].get(3, 0) * ZsRe.get(i - 1, 0);
        }
        break;
    }
    return _Ploss;
}

float MethodOPF::getQLoss(){
    _Qloss = 0;
    switch (losstype)
    {
    case LossType::POWER:
        for (int n = 1; n < _nAgent; n++) {
            int bus = (int) _CoresBusAgent.get(n, 0);
            int In = (int) PosAgent.get(n, 0);

            _Qloss -= X[bus].get(6 + 2 * In, 0);
        }
        break;
    case LossType::CURRENT:
        for (int i = 1; i < _nBus; i++) {
            _Qloss -= X[i].get(3, 0) * ZsIm.get(i - 1, 0);
        }
        break;
    }
    return _Qloss;
}


float MethodOPF::getPLoss2(){
    _Ploss = 0;
    switch (losstype)
    {
    case LossType::POWER:
        for (int i = 1; i < _nBus; i++) {
            _Ploss -= X[i].get(3, 0) * ZsRe.get(i - 1, 0);
        }
        
        break;
    case LossType::CURRENT:
        for (int n = 1; n < _nAgent; n++) {
            int bus = (int) _CoresBusAgent.get(n, 0);
            int In = (int) PosAgent.get(n, 0);

            _Ploss -= X[bus].get(5 + 2 * In, 0);
        }
        break;
    }
    return _Ploss;
}

float MethodOPF::getQLoss2(){
    _Qloss = 0;
    switch (losstype)
    {
    case LossType::POWER:
        for (int i = 1; i < _nBus; i++) {
            _Qloss -= X[i].get(3, 0) * ZsIm.get(i - 1, 0);
        }
        
        break;
    case LossType::CURRENT:
        for (int n = 1; n < _nAgent; n++) {
            int bus = (int) _CoresBusAgent.get(n, 0);
            int In = (int) PosAgent.get(n, 0);

            _Qloss -= X[bus].get(6 + 2 * In, 0);
        }
        break;
    }
    return _Qloss;
}


double MethodOPF::calcFc(MatrixCPUD* cost1, MatrixCPUD* cost2, MatrixCPUD* Pn, MatrixCPUD* tempN1)
{
    tempN1->set(cost1);
    tempN1->multiply(0.5);
    tempN1->multiplyT(Pn);
    tempN1->add(cost2);
    tempN1->multiplyT(Pn);
    return tempN1->sum();
}
float MethodOPF::calcFc()
{
    tempN2.set(&Cost1);
    tempN2.multiply(0.5f);
    tempN2.multiplyT(&Pn);
    tempN2.add(&Cost2);
    tempN2.multiplyT(&Pn);

    return tempN2.sum();
}

int MethodOPF::feasiblePoint() {
	MatrixCPU test(_nBus, 1, -1);
	int counter = 0;
	for (int bus = 0; bus < _nBus; bus++) {

		float Si = X[bus].get(indPi, 0) * X[bus].get(indPi, 0) + X[bus].get(indQi, 0) * X[bus].get(indQi, 0);
		float li = X[bus].get(indli, 0);
		float vi = X[bus].get(indvi, 0);
		float err = Si - li * vi;
		test.set(bus, 0, err);
		if (abs(err) > 0.0001) {
			counter++;
		}
	}
	resF.set(2, (_iterGlobal - 1) / _stepG, test.max2());
	return counter;
}


void MethodOPF::setLagrange(bool lagrange) {
    Lagrange = lagrange;
}
float MethodOPF::getMurhoVar() {
    return _mu;
}
float MethodOPF::getTaurhoVar() {
    return _tau;
}
bool MethodOPF::chekcase(){
    if (_nBus != (_nLine + 1)) {
        std::cout << "wrong number of line " << _nLine << "against " << _nBus << std::endl;
        return false;
    }
    for (int i = 0; i < _nLine; i++) {
        if (CoresLineBus.get(i, 1) != (i + 1)) {
            std::cout << "wrong numerotation of line " << CoresLineBus.get(i, 1) << "against " << (i + 1) << std::endl;
            return false;
        }
        if (CoresLineBus.get(i, 0) > CoresLineBus.get(i, 1)) {
            std::cout << "wrong numerotation of bus " << CoresLineBus.get(i, 0) << "against " << CoresLineBus.get(i, 1) << std::endl;
            return false;
        }
    }
    if (ZsRe.getNLin() == 0  || ZsIm.getNLin() == 0) {
        std::cout << "matrices non defined, ZsRe, Zs Im, Yd" << std::endl;
        ZsRe.display();
        ZsIm.display();
        return false;
    }
    return true;
}

void MethodOPF::updateGlobalProb(){
    
	for (int i = 0; i < _sizeProbGlob; i++) {
		Ypre[i].swap(&Y[i]);
		Y[i].MultiplyMatVec(&Hinv[i], &Q[i]); // solve system by using the inverse
	}
	
	Y[0].set(indvi, 0, 1);
	Y[0].set(indvai, 0, 1);
}

void MethodOPF::updateMu()
{
	for (int i = 0; i < _sizeProbGlob; i++) {
		tempM[i].subtract(&X[i], &Y[i]);
		tempM[i].multiply(_rho);
		Mu[i].add(&tempM[i]);
	}
}

void MethodOPF::updateXWOCurrent()
{
	double x1, x2, x3, x4, c1, c2, c3, c4;
    double lambdaLo, lambdaUp, x3min, x3max, gamma, k2;
	double c1122;
	int nSol = 0;
	int typeSol = 0;
	int BestRoot = 0;
	double bestGamma = -1;
	double p = 0;
	int nRoot = 0;

	for (int i = 1; i < _nBus; i++) {

		bool goodSol = false;
		k2 = sqrt(2.0 / (nChild.get(i, 0) + 1));
		
        c1 = -2 * Chat[i].get(indPi, 0);
        c2 = -2 * Chat[i].get(indQi, 0);
        c3 = -2 * Chat[i].get(indvi, 0) / k2;
        c4 = -2 * Chat[i].get(indli, 0);
        c1122 = c1 * c1 + c2 * c2;
        
        
        x3min = VoltageLimit.get(i, 0);
        x3max = VoltageLimit.get(i, 1);
			
        // case without constraint
        
        x1 = -c1 / 2;
        x2 = -c2 / 2;
        x3 = -c3 / 2;
        x4 = -c4 / 2;
        lambdaUp = 0;
        lambdaLo = 0;
        
        if (x3 < x3min) {
            x3 = x3min;
            lambdaLo = (2 * x3 + c3);
        }
        else if (x3 > x3max) {
            x3 = x3max;
            lambdaUp = -(2 * x3 + c3);
        }
        gamma = k2 * x4 - (x1 * x1 + x2 * x2) / x3; // ce n'est pas vraiment gamma, doit �tre positif
        //std::cout << "x 1 : " << x1 << " " << x2 << " " << x3 * k2 << " " << x4 << " " << (x1 * x1 + x2 * x2) / x3  - k2 * x4 << std::endl;

        if (gamma >= 0) {
            // the solution is good !
            goodSol = true;
        }
        else {
            if (gamma > bestGamma) {
                typeSol = 1;
                bestGamma = gamma;
            }
        }
		
		if (!goodSol) { // cas degenere
			if (c1122 == 0) {
				std::cout << " bus " << i << " : c1= " << c1 << " c2=" << c2 << " c4=" << c4 << " gamma= " << gamma << std::endl;
				x4 = 0;
				goodSol = true;
			}
		}
		// case x3 = x3max lambdaLo = 0
		if (!goodSol) {
			x3 = x3max;

			coefPoly2[0] = 2 * (c4 / (k2 * x3) + 1);
			coefPoly2[1] = 1 / x3;
			coefPoly2[0] = coefPoly2[0] * k2 * k2 / (4 * c1122);
			coefPoly2[1] = coefPoly2[1] * k2 * k2 / (4 * c1122);

			nRoot = resolveRealPolynome3without2term(root2, coefPoly2);
			
			for (int n = 0; n < nRoot; n++) {
				p = root2[n];
								
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
				gamma = (2 * x4 + c4) / k2;
				lambdaUp = -(2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));
				//std::cout << "x2 : " << x1 << " " << x2 << " " << x3 * k2 << " " << x4 << " " << gamma << " " << lambdaUp << std::endl;
				if (gamma >=  0 && lambdaUp >=  0) {
					// the solution is good 
					goodSol = true;
					//nSol = n;
					break;
				}
				if (gamma > bestGamma && lambdaUp > bestGamma) {
					typeSol = 2;
					bestGamma = MYMIN(gamma, lambdaLo);
					BestRoot = n;
				}

			}
			// case x3 = x3min lambdaUp = 0
			if (!goodSol) {
				x3 = x3min;
				 
				coefPoly2[0] = 2 * (c4 / (k2 * x3) + 1);
				coefPoly2[1] = 1 / x3;
				coefPoly2[0] = coefPoly2[0] * k2 * k2 / (4 * c1122);
				coefPoly2[1] = coefPoly2[1] * k2 * k2 / (4 * c1122);

				nRoot = resolveRealPolynome3without2term(root3, coefPoly2);

				for (int n = 0; n < nRoot; n++) {
					p = root3[n];
					//std::cout << "poly " << coefPoly2[0] * p + coefPoly2[1] + p * p * p << std::endl;
					x1 = p * c1 * x3;
					x2 = p * c2 * x3;
					x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
					gamma = (2 * x4 + c4) / k2;
					lambdaLo = (2 * x3 + c3 - gamma * (x1 * x1 + x2 * x2) / (x3 * x3));
					
					//std::cout << "x 3: " << x1 << " " << x2 << " " << x3 * k2 << " " << x4 << " " << gamma << " " << lambdaLo << std::endl;

					if (gamma >= 0 && lambdaLo >= 0) {
						// the solution is good !
						goodSol = true;
						break;
					}
					if (gamma > bestGamma && lambdaLo > bestGamma) {
						typeSol = 3;
						bestGamma = MYMIN(gamma, lambdaLo);
						BestRoot = n;
					}
				}
			}
			// case xmin<x3<xmax lambdaLo = 0 lambdaUp = 0
			if (!goodSol) {
				 
				coefPoly3[0] = c1122 / k2 * (2 * c3 / k2 - c4);
				coefPoly3[1] = (c3 - 2 * c4 / k2);
				coefPoly3[2] = -1;
				coefPoly3[0] = coefPoly3[0] * k2 * k2 / (c1122 * c1122);
				coefPoly3[1] = coefPoly3[1] * k2 * k2 / (c1122 * c1122);
				coefPoly3[2] = coefPoly3[2] * k2 * k2 / (c1122 * c1122);

				nRoot = resvolveRealPolynome4without2term(root4, coefPoly3, Lagrange);

				for (int n = 0; n < nRoot; n++) {
					p = root4[n];
					//std::cout << "poly " <<p * p * p * p + coefPoly3[0] * p*p*p + coefPoly3[1]*p + coefPoly3[2] << std::endl;
					x3 = -(c1122 * p + 2 * c3) / (2 * (c1122 * p * p + 2));
					x1 = p * c1 * x3;
					x2 = p * c2 * x3;
					x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
					gamma = (2 * x4 + c4) / k2;
					//std::cout << "x 4: " << x1 << " " << x2 << " " << x3 * k2 << " " << x4 << " " << gamma << " " << std::endl;

					if (gamma >= 0 && x3 <= x3max && x3 >= x3min) {
						// the solution is good !
						goodSol = true;
						break;
					}if (gamma > bestGamma && (x3max - x3) > bestGamma && (x3 - x3min) > bestGamma) {
						typeSol = 4;
						bestGamma = MYMIN(MYMIN(gamma, (x3max - x3)), (x3 - x3min));
						BestRoot = n;
					}
				}
			}
		}
		if (!goodSol) {
			if (typeSol == 1) {
				// case without constraint
				x1 = -c1 / 2;
				x2 = -c2 / 2;
				x3 = -c3 / 2;
				x3 = (x3max - x3) * (x3 > x3max) + (x3min - x3) * (x3min > x3) + x3;
				x4 = -c4 / 2; // ou  (x1 * x1 + x2 * x2) / (k2 * x3)
			}
			else {
				if (typeSol == 2) {
					x3 = x3max;
					p = root2[BestRoot];
				}
				else if (typeSol == 3) {
					x3 = x3min;
					p = root3[BestRoot];
				}
				else if (typeSol == 4) {
					p = root4[BestRoot];
					x3 = -(c1122 * p + 2 * c3) / (2 * (c1122 * p * p + 2));
					x3 = (x3max - x3) * (x3 > x3max) + (x3min - x3) * (x3min > x3) + x3;
				}
				x1 = p * c1 * x3;
				x2 = p * c2 * x3;
				x4 = (x1 * x1 + x2 * x2) / (x3 * k2);
			}
        }
		
		// X =  {Pi, Qi, vi, li, pi, qi, vAi, (Pci, Qci, lci) for all child Ci}
		X[i].set(indPi, 0, x1);
		X[i].set(indQi, 0, x2);
		X[i].set(indvi, 0, x3* k2);
		X[i].set(indli, 0, x4);
		
		//std::cout << "x F : " << x1 << " " << x2 << " " << x3*k2 << " " << x4 << " " << gamma << std::endl;
		
	}

}


float MethodOPF::updateRes(int indice) 
{
	float resS = 0;
	float resR = 0;
	float resV = 0;
	for (int i = 0; i < _sizeProbGlob; i++) {
		
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
			for (int i = 0; i < _sizeProbGlob; i++) {
				Hinv[i].divide(_tau);
			}
			//std::cout << _iterGlobal << "rho augmente " << _rho << std::endl;
		}
		else if (resS > _mu * resR) {// rho = rho / tau_inc;
			_rho = _rho / _tau;
			for (int i = 0; i < _sizeProbGlob; i++) {
				Hinv[i].multiply(_tau);
			}
			//std::cout << _iterGlobal << "rho diminue " << _rho << std::endl;
		}
	}
	
	return MYMAX(MYMAX(resV, oldrho * resS), resR);
}

float MethodOPF::updateResRhoFixe(int indice)
{
	float resS = 0;
	float resR = 0;
	float resV = 0;
	for (int i = 0; i < _sizeProbGlob; i++) {

		float resTempS = _rho * Y[i].max2(&Ypre[i]);
		float resTempR = Y[i].max2(&X[i]);

		if (resTempS > resS) {
			resS = resTempS;
		}
		if (resTempR > resR) {
			resR = resTempR;
		}

	}
	resF.set(0, indice, resR);
	resF.set(1, indice, resS);
	resF.set(2, indice, resV);

	return MYMAX(MYMAX(resV, resS), resR);
}


void MethodOPF::computePb(){
    Pb.set(0.0);
	for (int i = 0; i < _nBus; i++) {
		int Nb = (int) _nAgentByBus.get(i, 0);
		int begin = (int) _CoresAgentBusBegin.get(i, 0);
		for (int In = 0; In < Nb; In++) {
			int n = (int) _CoresAgentBus.get(In + begin, 0);
			PosAgent.set(n, 0, In);
			Pb.increment(i, 0, Pn.get(n, 0));
			Pb.increment(i + _nBus, 0, Pn.get(n + _nAgent, 0));
		}
	}
}


MatrixCPU MethodOPF::getPb(){
	computePb();
	return Pb;
}
MatrixCPU MethodOPF::getPhi(){
	MatrixCPU Phi(2*_nLine, 1);
	for (int i = 0; i < _nLine; i++)
	{
		Phi.set(i, 0, Y[i + 1].get(indPi,0));
		Phi.set(i + _nLine, 0, Y[i + 1].get(indQi,0));
	}
	return Phi;
}
MatrixCPU MethodOPF::getE(){
	MatrixCPU E(2*_nBus, 1);
	for (int i = 0; i < _nBus; i++)
	{
		E.set(i, 0, Y[i].get(indli,0)); //l_i
		E.set(i + _nBus, 0, Y[i].get(indvi,0)); // v_i
	}
	return E;
}

float MethodOPF::DFSP(int j)
{
	//std::cout << "DFSP " << j << std::endl;
	float p = Pb.get(j, 0);
	for (int i = 0; i < nChild.get(j, 0); i++) {
		int c = (int) Childs[j].get(i, 0);
		p += DFSP(c);
	}
	X[j].set(indPi, 0, p);
	return p;
}
float MethodOPF::DFSQ(int j)
{
	float q = Pb.get(j + _nBus, 0);
	for (int i = 0; i < nChild.get(j, 0); i++) {
		int c = (int) Childs[j].get(i, 0);
		q += DFSQ(c);
	}
	X[j].set(indQi, 0, q);
	return q;
}


void MethodOPF::display(){
    std::cout.precision(3);

    computePb();

	if (_iterGlobal == 0) {
		std::cout << "algorithm not launch" << std::endl;
	}
	else if (_iterGlobal < _iterG) {
		std::cout << "method " << _name << " converged in " << _iterGlobal << " iterations." << std::endl;
		std::cout << "Converged in " << (float) timeOPF / CLOCKS_PER_SEC << " seconds" << std::endl;

	}
	else {
		std::cout << "method " << _name << " not converged in " << _iterGlobal << " iterations." << std::endl;
		std::cout << "time taken " << (float) timeOPF / CLOCKS_PER_SEC << " seconds" << std::endl;
	}
	std::cout << "The power error of this state is (constraint) " << resF.get(0, _iterGlobal / _stepG) << " and convergence " << resF.get(1, _iterGlobal / _stepG) << std::endl;
	std::cout << "===============================================================|" << std::endl;
	std::cout << "      System Summary                                           |" << std::endl;
	std::cout << "===============================================================|" << std::endl;
	std::cout << "Buses            " << _nBus << std::endl;
	std::cout << "Branches         " << _nLine << std::endl;
	std::cout << "Agent            " << _nAgent << std::endl;
	std::cout << "Ploss            " << getPLoss() << " or " << getPLoss2() << std::endl;
	std::cout << "Qloss            " << getQLoss() << " or " << getQLoss2() << std::endl;


	std::cout << std::endl << std::endl;
	
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "      Bus Data                                                                                          |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Bus |    Voltage  |   Power = Generation  + Load    |                Mu voltage and power              |" << std::endl;
	std::cout << "  #  |     Mag(pu) |    P (pu)      |      Q (pu)    |     V (pu)     |      P (pu)    |      Q (pu)    |" << std::endl;
	std::cout << "-----|-------------|----------------|----------------|----------------|----------------|----------------|" << std::endl;

		
	float seuil = 0.0001f;
		
	for (int b = 0; b < _nBus; b++) {
	std::cout << std::setw(5) << b << "|" << std::setw(12) << sqrt(X[b].get(2,0)) << " |" << std::setw(16)
			<< (abs(X[b].get(4, 0)) > seuil) * X[b].get(4, 0) << "|" << std::setw(16) << (abs(X[b].get(5, 0)) > seuil) * X[b].get(5, 0)
			<< "|" << std::setw(16) << Mu[b].get(2, 0) << "|" << std::setw(16)
			<< Mu[b].get(4, 0) << "|" << std::setw(16) << Mu[b].get(5, 0) << "|" << std::endl;

	}
	std::cout << std::endl << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "      Line Data                                                                                         |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Line |    From     |    To      |                           Upstream flow                              |" << std::endl;
	std::cout << "  #   |    Bus      |    Bus     |    P (pu)      |    Q (pu)      |     l (pu)     |     Loss (pu)     |" << std::endl;
	std::cout << "------|-------------|------------|----------------|----------------|----------------|-------------------|" << std::endl;

	for (int l = 0; l < _nLine; l++) {
		int b = l + 1;
		std::cout << std::setw(6) << l << "|" << std::setw(12) << CoresLineBus.get(l, 0) << " |" << std::setw(12)
			<< CoresLineBus.get(l, 1) << "|" << std::setw(16) << X[b].get(0, 0)
			<< "|" << std::setw(16) << X[b].get(1, 0) << "|" << std::setw(16)
			<< X[b].get(3, 0) << "|" << std::setw(19) << X[b].get(3, 0) * ZsRe.get(l, 0) << "|" << std::endl;
	}
	std::cout << std::endl << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "     Constraints                                                                                        |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Bus | Voltage | Voltage | Voltage |        Power Injection          |          Power Injection         |" << std::endl;
	std::cout << "  #  | Mag(pu) | MIN(pu) | MYMAX(pu)|  P (pu) | Pmin (pu) | Pmax (pu) |  Q (pu)  | Qmin (pu) | Qmax (pu) |" << std::endl;
	std::cout << "-----|---------|---------|---------|---------|-----------|-----------|----------|-----------|-----------|" << std::endl;
	

	for (int b = 0; b < _nBus; b++) {
		int nb = (int) _nAgentByBus.get(b, 0);
		std::cout << std::setw(5) << b << "|" << std::setw(8) << sqrt(Y[b].get(2, 0)) << " |" << std::setw(9)
			<< VoltageLimitReal.get(b, 0) << "|" << std::setw(9) << VoltageLimitReal.get(b, 1)
			<< "|" << std::setw(9) << Pb.get(b, 0) << "|" << std::setw(11)
			<< Pbmin.get(b, 0) << "|" << std::setw(11) << Pbmax.get(b,0) << "|" << std::setw(10) << Pb.get(b + _nBus, 0)
			<< "|" << std::setw(11) << Pbmin.get(b + _nBus, 0) << "|" << std::setw(11) << Pbmax.get(b + _nBus, 0) << "|" << std::endl;

	}

    
	std::cout << std::endl << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "      Agent Data                                                                                        |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;
	std::cout << " Agent |  Bus  |  Cost   |  Cost   |        Power Injection          |          Power Injection         |" << std::endl;
	std::cout << "  #    |   #   |  a (pu) |  b (pu) |  P (pu) | Pmin (pu) | Pmax (pu) |  Q (pu)  | Qmin (pu) | Qmax (pu) |" << std::endl;
	std::cout << "-------|-------|---------|---------|---------|-----------|-----------|----------|-----------|-----------|" << std::endl;

	for (int n = 0; n < _nAgent; n++) {
		int b = (int) _CoresBusAgent.get(n, 0);
		std::cout << std::setw(7) << n << "|" << std::setw(7) << b << "|" << std::setw(8) << Cost1.get(n,0) << " |" << std::setw(9)
			<< Cost2.get(n, 0) << "|" << std::setw(9) << Pn.get(n,0) << "|" << std::setw(11)
			<< Pmin.get(n, 0) << "|" << std::setw(11) << Pmax.get(n, 0) << "|" << std::setw(10) << Pn.get(n + _nAgent, 0)
			<< "|" << std::setw(11) << Pmin.get(n + _nAgent, 0) << "|" << std::setw(11) << Pmax.get(n + _nAgent, 0) << "|" << std::endl;
	}


	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "                      END PRINT                                                                         |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;

}
