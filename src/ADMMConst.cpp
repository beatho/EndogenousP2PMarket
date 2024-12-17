#include "../head/ADMMConst.h"



ADMMConst::ADMMConst() : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << " ADMMConst Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
}


ADMMConst::ADMMConst(float rho) : MethodP2P()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "default ADMMConst Constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	_name = NAME;
	_rho = rho;
}

ADMMConst::~ADMMConst()
{
}
void ADMMConst::setParam(float rho)
{
	_rho = rho;
}

void ADMMConst::setTau(float tau)
{
	throw std::domain_error("tau is not define for this method");
}



void ADMMConst::solve(Simparam* result, const Simparam& sim, const StudyCase& cas)
{
#ifdef DEBUG_SOLVE
	cas.display();
	sim.display(1);
#endif // DEBUG_SOLVE
	
	clock_t t =clock();
	// FB 0
	float rho = sim.getRho();
	int iterG = sim.getIterG();
	int iterL = sim.getIterL();
	int stepL = sim.getStepL();
	int stepG = sim.getStepG();
	
	float epsG = sim.getEpsG();
	float epsGC = sim.getEpsGC();
	_ratioEps = epsG / epsGC;
	float epsL = sim.getEpsL();
	int nAgent = sim.getNAgent();
	int nLine = cas.getNLine();
	int nBus = cas.getNBus();
	MatrixCPU BETA(cas.getBeta());
	MatrixCPU connect(cas.getC());
	MatrixCPU a(cas.geta());
	MatrixCPU b(cas.getb());
	MatrixCPU Ub(cas.getUb());
	MatrixCPU Lb(cas.getLb());
	MatrixCPU Pmin(cas.getPmin());
	MatrixCPU Pmax(cas.getPmax());
	MatrixCPU nVoisin(cas.getNvoi());
	MatrixCPU Llimit(cas.getLineLimit());
	MatrixCPU G(cas.getPowerSensi());
	MatrixCPU G2(G);
	G2.multiplyT(&G);

	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);


	MatrixCPU LAMBDA(sim.getLambda());
	MatrixCPU trade(sim.getTrade());
	MatrixCPU Pn(sim.getPn()); // trades sum for each agent
	
	
	MatrixCPU resF(3, (iterG/stepG)+1);
	float fc = 0;


	MatrixCPU MU(nAgent, 1); //  lambda_l/_rho
	MatrixCPU Tlocal(nAgent, nAgent);
	MatrixCPU Tlocal_pre(trade);
	MatrixCPU Tmoy(nAgent,1);
	Tmoy.Moy(&Tlocal_pre, &nVoisin);
	MatrixCPU P(nAgent, 1); // trades mean for each agent
	MatrixCPU Qpart(nLine, nAgent);
	MatrixCPU Qtot(nLine, 1);
	MatrixCPU Kappa1(nLine, 1, 0);
	MatrixCPU Kappa2(nLine, 1, 0);
	MatrixCPU Kappa1_pre(nLine, 1, 0);
	MatrixCPU Kappa2_pre(nLine, 1, 0);


	MatrixCPU alpha(nLine, nAgent);
	

	float rho_p = _rho;
	if (_rho == 0) {
		rho_p = rho;
	}
	float rho1 = sim.getRho1();

	float at1 = rho; // 2*a in the article
	float at2 = rho_p;

	MatrixCPU Ap2(a);
	MatrixCPU Ap1(nVoisin);
	MatrixCPU Ap12(nAgent, 1); // Ap2+Ap1;
	MatrixCPU Bt1(nAgent, nAgent);
	MatrixCPU Bt2(nAgent, nAgent);
	MatrixCPU Bp1(nAgent, 1);
	MatrixCPU Ct(BETA);
	MatrixCPU Cp1(nAgent, 1);
	MatrixCPU Cp2(nAgent, 1);
	MatrixCPU Cp(nAgent, 1);
	MatrixCPU matUb(nAgent, nAgent);
	MatrixCPU matLb(nAgent, nAgent);

	MatrixCPU tempN1(nAgent, 1);
	MatrixCPU temp1N(1, nAgent);
	MatrixCPU tempNN(nAgent, nAgent);
	MatrixCPU tempL1(nLine,1);


	int iterLocal = 0;
	for (int i = 0; i < nAgent; i++) 
	{
		for (int j = 0; j < nAgent;j++) {
			matUb.set(i, j, Ub.get(i, 0));
			matLb.set(i, j, Lb.get(i, 0));
		}
	}
	matUb.multiplyT(&connect);
	matLb.multiplyT(&connect);
	
	
	temp1N.sum(&G2, 1);
	
	
	temp1N.multiply(2*rho1);
	Ap2.addTrans(&temp1N);
	Ap2.multiplyT(&nVoisin);
	Ap2.multiplyT(&nVoisin);
	Ap1.multiply(rho_p);
	Ap12.add(&Ap1, &Ap2);
	Cp1.multiplyT(&b, &nVoisin);

	// FB 1
	Pn.set(&Tmoy);
	Pn.multiplyT(&nVoisin);
	alpha.multiplyTVector(&G, &Pn, 0);
	updateQ(&Qpart, &Qtot, &alpha, nAgent, nLine);
	updateLAMBDA(&LAMBDA, &trade, rho);
	Kappa1_pre.set(&Kappa1);
	Kappa2_pre.set(&Kappa2);
	updateKappa(&Kappa1, &Kappa2, &Llimit, &Qtot);
	updateBt1(&Bt1, &trade, rho, &LAMBDA);
	updateCp2(&Cp2, rho1, &Kappa1, &Kappa2, &G, &tempL1, &Qpart, &nVoisin, nLine, nAgent);
	Cp.add(&Cp1, &Cp2);

	std::cout << "fin init" << std::endl;

	
	

	float resG = 2 * epsG;
	float resL = 2 * epsL;
	int iterGlobal = 0;
	while ((iterGlobal < iterG) && (resG>epsG)) {
		resL = 2 * epsL;
		iterLocal = 0;
		while (iterLocal< iterL && resL>epsL) {
			// FB 2a	
			updateBt2(&Bt2,&Tlocal_pre,&Tmoy,&P,&MU);
			
			updateTl(&Tlocal, at1, at2, &Bt1, &Bt2, &Ct, &matLb, &matUb);
			// FB 2b
			Tmoy.Moy(&Tlocal, &nVoisin); 
			
			// FB 2c
			updateBp1(&Bp1, &MU, &Tmoy);
			updateP(&P, &Ap1, &Ap12, &Bp1, &Cp, &Pmin,&Pmax);
			updateMU(&MU,&Tmoy,&P);

			// FB 3
			resL = calcRes(&Tlocal,&Tlocal_pre,&Tmoy,&P); 
			Tlocal_pre.swap(&Tlocal); 
			iterLocal++;
			
		}
		//std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, iterGlobal / stepG) << " " << resF.get(1, iterGlobal / stepG) << std::endl;
		Tlocal_pre.swap(&Tlocal);
		trade.swap(&Tlocal);
		
		Pn.set(&Tmoy);
		Pn.multiplyT(&nVoisin);
		alpha.multiplyTVector(&G, &Pn, 0);
		updateQ(&Qpart, &Qtot, &alpha, nAgent, nLine);
		updateLAMBDA(&LAMBDA, &trade, rho);
		Kappa1_pre.set(&Kappa1);
		Kappa2_pre.set(&Kappa2);
		updateKappa(&Kappa1, &Kappa2, &Llimit, &Qtot);
		// FB 1b
		updateBt1(&Bt1, &trade, rho, &LAMBDA);
		updateCp2(&Cp2, rho1, &Kappa1, &Kappa2, &G, &tempL1, &Qpart, &nVoisin, nLine, nAgent);
		
		Cp.add(&Cp1, &Cp2);
		

		// FB 4
		resG = updateRes(&resF, &trade, &Tlocal, (iterGlobal/stepG), &Kappa1, &Kappa2, &Kappa1_pre, &Kappa2_pre);
		
		iterGlobal++;
	}
	std::cout << iterGlobal << " " << iterLocal << " " << resL << " " << resF.get(0, (iterGlobal - 1) / stepG) << " " << resF.get(1, (iterGlobal - 1) / stepG) << " " << resF.get(2, (iterGlobal - 1) / stepG) << std::endl;
	// FB 5
	result->setResF(&resF);
	
	result->setLAMBDA(&LAMBDA);
	
	result->setTrade(&trade);
	
	result->setIter(iterGlobal);
	
	Pn.set(&Tmoy);
	Pn.multiplyT(&nVoisin);
	result->setPn(&Pn);
	fc = calcFc(&a, &b, &trade, &Pn, &BETA,&tempN1,&tempNN);
	result->setFc(fc);
	t = clock() - t;
	result->setTime((float)t / CLOCKS_PER_SEC);
	
}

void ADMMConst::updateP0(const StudyCase& cas)
{
	// not used for this method
}

void ADMMConst::init(const Simparam& sim, const StudyCase& cas)
{
	// not used for this method
}


void ADMMConst::updateBt1(MatrixCPU* Bt1, MatrixCPU* trade, float rho, MatrixCPU* LAMBDA)
{
	Bt1->set(trade);
	Bt1->subtractTrans(trade);
	Bt1->multiply(0.5*rho); 
	Bt1->subtract(LAMBDA);
	Bt1->divide(rho);

}

void ADMMConst::updateBt2(MatrixCPU* Bt2, MatrixCPU* Tlocal, MatrixCPU* Tmoy, MatrixCPU* P, MatrixCPU* MU)
{
	Bt2->set(Tlocal);
	Bt2->subtractVector(Tmoy);
	Bt2->addVector(P);
	Bt2->subtractVector(MU);
}

void ADMMConst::updateBp1(MatrixCPU* Bp1, MatrixCPU* MU, MatrixCPU* Tmoy)
{
	Bp1->add(MU, Tmoy);
}

void ADMMConst::updateTl(MatrixCPU* Tlocal, float at1, float at2, MatrixCPU* Bt1, MatrixCPU*Bt2, MatrixCPU* Ct, MatrixCPU* matLb, MatrixCPU* matUb)
{

	float ada = at1 / at2; 
	float apa = at1 + at2;

	Tlocal->set(Bt1);
	Tlocal->multiply(ada);
	Tlocal->add(Bt2);
	Tlocal->multiply(at2);

	Tlocal->subtract(Ct);
	Tlocal->divide(apa); 
	Tlocal->project(matLb, matUb);
}

float ADMMConst::calcRes( MatrixCPU* Tlocal, MatrixCPU* Tlocal_pre, MatrixCPU* Tmoy, MatrixCPU* P)
{
	MatrixCPU temp(*Tlocal);
	temp.subtract(Tlocal_pre);

	MatrixCPU temp2(*Tmoy);
	temp2.subtract(P);
	float d1 = temp.max2();
	float d2 = temp2.max2();
	//std::cout << " ResS " << d1 << " ResR " << d2 << std::endl;

	return d1 * (d1 > d2) + d2 * (d2 >= d1);
}

void ADMMConst::updateP(MatrixCPU* P, MatrixCPU* Ap1, MatrixCPU* Ap12, MatrixCPU* Bp1, MatrixCPU* Cp, MatrixCPU* Pmin, MatrixCPU* Pmax)
{
	P->multiplyT(Ap1, Bp1);
	P->subtract(Cp);
	
	P->divideT(Ap12);
	P->project(Pmin, Pmax);
}

void ADMMConst::updateMU(MatrixCPU* MU, MatrixCPU* Tmoy, MatrixCPU* P)
{
	MU->add(Tmoy);
	MU->subtract(P);
}

void ADMMConst::updateQ(MatrixCPU* Qpart, MatrixCPU* Qtot, MatrixCPU* alpha, int nAgent, int nLine)
{
	for (int l = 0; l < nLine; l++) {
		float qt = 0;
		for (int n = nAgent - 1; n >= 0; n--) {
			qt += alpha->get(l, n);
			if (n > 0) {
				Qpart->set(l, n - 1, qt);
			}
		}
		Qtot->set(l, 0, qt);
	}


}

void ADMMConst::display() {

	std::cout << _name << std::endl;
}
