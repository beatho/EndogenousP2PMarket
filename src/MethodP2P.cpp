#pragma once
#include "../head/MethodP2P.h"
 


MethodP2P::MethodP2P() : Method()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "method constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1abcd, Fb2, Fb3, Fb5, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu
}

MethodP2P::~MethodP2P() 
{
}

void MethodP2P::updateLAMBDA()
{
	tempNN = trade;
	tempNN.addTrans(&trade);
	tempNN.multiply(_rho);
	tempNN.multiply(0.5);
	LAMBDA.add(&LAMBDA, &tempNN);
}
void MethodP2P::updateLambda()
{
	for (int t = 0; t < _nTrade; t++) {
		int k = CoresLinTrans.get(t, 0);
		float lamb = 0.5 * _rhog * (tradeLin.get(t, 0) + tradeLin.get(k, 0));
		LAMBDALin.set(t, 0, LAMBDALin.get(t, 0) + lamb);
	}
}



void MethodP2P::updateKappa()
{
	//
	Kappa1.projectNeg();
	Kappa1.add(&lLimit);
	Kappa1.subtract(&Qtot);
	
	Kappa2.projectNeg();
	Kappa2.add(&lLimit);
	Kappa2.add(&Qtot);
	//
}

void MethodP2P::updateCp2()
{
	tempL1.subtractAbs(&Kappa1, &Kappa2);
	//Cp2->multiplyTrans(G, tempL1, 0);

	float r = 0;
	for (int i = 0; i < _nAgent; ++i)
	{
		r = 0;
		for (int k = 0; k < _nLine; ++k)
		{
			r +=  G.get(k, i) * (tempL1.get(k, 0) + 2 * Qpart.get(k, i));
		}
		Cp2.set(i, 0, r);
	}

	Cp2.multiply(_rho1);
	Cp2.multiplyT(&nVoisin);
}

float MethodP2P::updateResMat(int iter)
{
	tempNN.subtract(&Tlocal, &trade);
	
	float resS = tempNN.max2();
	tempNN.set(&Tlocal);
	tempNN.addTrans(&Tlocal);
	float resR = tempNN.max2();
	
	if (iter > 0 && _tau > 1) {
		if (resR > _mu * resS) {
			_rhog = _tau * _rhog;
			_at1 = _rhog;
			//std::cout << iter << ", rho augmente :" << _rhog << std::endl;
		}
		else if (resS > _mu * resR) {// rho = rho / tau_inc;
			_rhog = _rhog / _tau;
			_at1 = _rhog;
			//std::cout << iter << ", rho diminue :" << _rhog << std::endl;
		}
	}

	resF.set(0, iter, resR);
	resF.set(1, iter, resS);
	

	return resR* (resR > resS) + resS * (resR <= resS);
}
float MethodP2P::updateRes(int iter)
{
	for (int t = 0; t < _nTrade; t++) {
		int k = CoresLinTrans.get(t, 0);
		tempNN.set(t, 0, tradeLin.get(t, 0) + tradeLin.get(k, 0));
	}
	float resR = tempNN.max2();

	
	tempNN.subtract(&Tlocal, &tradeLin);
	float resS = tempNN.max2();

	resF.set(0, iter, resR);
	resF.set(1, iter, resS);

	if (iter > 0 && _tau > 1) {
		if (resR > _mu * resS) {
			_rhog = _tau * _rhog;
			_at1 = _rhog;
			//std::cout << iter << ", rho augmente :" << _rhog << std::endl;
		}
		else if (resS > _mu * resR) {// rho = rho / tau_inc;
			_rhog = _rhog / _tau;
			_at1 = _rhog;
			//std::cout << iter << ", rho diminue :" << _rhog << std::endl;
		}
	}



	return resR * (resR > resS) + resS * (resR <= resS);

}
float MethodP2P::updateResEndo(int iter){
	for (int t = 0; t < _nTrade; t++) {
		int k = CoresLinTrans.get(t, 0);
		tempNN.set(t, 0, tradeLin.get(t, 0) + tradeLin.get(k, 0));
	}
	float resR = tempNN.max2();

	
	tempNN.subtract(&Tlocal, &tradeLin);
	float resS = tempNN.max2();

	// version de l'article
	tempL1.set(&Kappa1);
	tempL2.set(&Kappa2);
	Kappa1_pre.projectNeg();
	Kappa2_pre.projectNeg();
	tempL1.projectNeg();
	tempL2.projectNeg();
	tempL1.subtract(&Kappa1_pre);
	tempL2.subtract(&Kappa2_pre);
	tempL1.multiplyT(&tempL1);
	tempL2.multiplyT(&tempL2);
	tempL1.add(&tempL2);

	float resXf = _ratioEps * sqrt(tempL1.max2());
	resF.set(0, iter, resR);
	resF.set(1, iter, resS);
	resF.set(2, iter, resXf);
	return MYMAX(MYMAX(resXf, resS), resR);
}

float MethodP2P::calcRes()
{
	float d1 = Tlocal.max2(&Tlocal_pre);
	float d2 = P.max2(&Tmoy);

	return d1* (d1 > d2) + d2 * (d2 >= d1);
}

float MethodP2P::calcFc()
{

	float fc = 0;
	tempN1.set(&a);
	tempN1.multiply(0.5);
	tempN1.multiplyT(&Pn);
	tempN1.add(&b);
	tempN1.multiplyT(&Pn);
	
	fc = fc + tempN1.sum();
	tempNN.set(&trade);
	tempNN.multiplyT(&Ct);
	

	fc = fc + tempNN.sum();

	return fc;

}



void MethodP2P::updatePn()
{
	Pn.set(&Tmoy);
	Pn.multiplyT(&nVoisin);
}



void MethodP2P::solveWithMinPower(Simparam* result, const Simparam& sim, const StudyCase& cas)
{
	std::cout << "solveWithMinPower : should not be called" << std::endl;
}



