#include "../head/TestADMMConst.h"

int testADMM()
{
	int n = 1;
	if (!testADMMContruct1()) return n;
	n++;
	if (!testADMMContruct2()) return n;
	n++;
	if (!testADMMContruct3()) return n;
	n++;
	if (!testADMMFunction1()) return n;
	n++;
	if (!testADMMFunction2()) return n;
	n++;
	if (!testADMMFunction3()) return n;
	n++;
	if (!testADMMFunction4()) return n;
	n++;
	if (!testADMMFunction5()) return n;
	n++;
	if (!testADMMFunction6()) return n;
	n++; // 10
	std::cout << "test solve" << std::endl;
	if (!testADMMSolve1()) return n;
	n++;
	if (!testADMMSolve2()) return n;
	n++;
	if (!testADMMSolve3()) return n;
	n++;
	if (!testADMMSolve4()) return n;
	n++;

	return 0;
}

bool testADMMContruct1()
{
	std::cout << "default constructor" << std::endl;
	ADMMConst a;
	return true;
}

bool testADMMContruct2()
{
	float rho = 2;

	std::cout << "param constructor" << std::endl;
	ADMMConst a(rho);
	return true;
}
bool testADMMContruct3()
{
	float rho = 2;

	std::cout << "2 times constructor" << std::endl;
	ADMMConst a;
	a = ADMMConst(rho);
	return true;
}

bool testADMMSolve1()
{
	//solve(Simparam* result, Simparam sim, StudyCase cas);
	StudyCase cas;
	cas.Set2node();

	//cas.display(0);
	//cas.display(1);

	int nAgent = cas.getNagent();
	
	Simparam param(nAgent, 1);
	float epsG = 0.0001f;
	int iterG = 1000;
	/*int iterG = 1;
	int iterL = 5;
	param.setItG(iterG);
	param.setItL(iterL);*/
	Simparam res(param);
	
	
	

	ADMMConst a;

	MatrixCPU Trade(nAgent, nAgent);
	Trade.set(0, 1, -1);
	Trade.set(1, 0, 1);
	a.solve(&res, param, cas);

	MatrixCPU trade = res.getTrade();
	res.display();
	trade.display();
	return trade.isEqual(&Trade,0.001);
}
bool testADMMSolve2()
{
	//solve(Simparam* result, Simparam sim, StudyCase cas);
	StudyCase cas;
	cas.Set4Agent();

	//cas.display(0);
	//cas.display(1);

	int nAgent = cas.getNagent();

	Simparam param(nAgent, 1);
	float epsG = 0.0001f;
	int iterG = 1000;
	param.setEpsG(epsG);
	param.setItG(iterG);
	/*int iterG = 1;
	int iterL = 5;
	param.setItG(iterG);
	param.setItL(iterL);*/
	Simparam res(param);




	ADMMConst a;

	MatrixCPU P(nAgent, 1);
	P.set(0, 0, -11.1764);
	P.set(1, 0, -11.1764);
	P.set(2, 0, 24.033);
	P.set(3, 0, -1.6809);
	a.solve(&res, param, cas);

	MatrixCPU P2 = res.getPn();
	res.display();
	P2.display();
	return P.isEqual(&P2, 0.01);
}

bool testADMMSolve3()
{
	//solve(Simparam* result, Simparam sim, StudyCase cas);
	StudyCase cas;
	cas.Set29node();
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;

	int nAgent = cas.getNagent();

	Simparam param(nAgent, cas.getNLine());
	float epsG = 0.00005f;
	int iterG = 100000;
	param.setEpsG(epsG);
	param.setItG(iterG);
	param.setEpsG(epsG);
	param.setRho(10000);

	Simparam res(param);
	ADMMConst a;
	a.solve(&res, param, cas);

	res.display();


	float fc = -87016;
	float Pn[31] = { -0.927, -4.51, -2.81, -0.795, -0.880, -0.0980, -0.128, -4.07, -3.04, -2.16, -0.593, -3.31, -0.389, -2.71, -2.93, -1.77, -0.490, -2.26, -1.05, -0.138, -2.24, 4.246, 5.19, 3.42, 3.75, 4.41, 2.34, 3.31, 2.57, 3.49, 4.23 };

	MatrixCPU P(31, 1);
	for (int i = 0; i < 31; i++) {
		P.set(i, 0, Pn[i]);
	}
	MatrixCPU P2 = res.getPn();
	P2.display();

	return P2.isEqual(&P, 0.1);
}

bool testADMMSolve4()
{
	//solve(Simparam* result, Simparam sim, StudyCase cas);
	StudyCase cas;
	float lim = 0.8;
	cas.Set2nodeConstraint(lim);
	int nAgent = cas.getNagent();
	Simparam param(nAgent, 1);
	Simparam res(param);


	float value = (1 - lim) * (lim > 1) + lim;

	ADMMConst a;

	MatrixCPU Trade(nAgent, nAgent);
	Trade.set(0, 1, -value);
	Trade.set(1, 0, value);
	a.solve(&res, param, cas);

	MatrixCPU trade = res.getTrade();
	res.display();
	trade.display();
	return trade.isEqual(&Trade, 0.001);
}


bool testADMMFunction1()
{
	/*
void ADMMConst::updateBt1(MatrixCPU* Bt1, MatrixCPU* trade, float rho, MatrixCPU* LAMBDA)
{
	Bt1->set(trade);
	Bt1->subtractTrans(trade);
	Bt1->multiply(0.5*rho); // c'est mieux que de copier une matrix complète
	Bt1->subtract(LAMBDA);
	Bt1->divide(rho);

}*/

	int nAgent = 2;
	float value1 = 2;
	float value2 = 4;
	float value3 = -2;
	float rho = 1.5;
	MatrixCPU Bt1(nAgent, nAgent,value1);
	MatrixCPU Bt11(nAgent, nAgent, value1);
	MatrixCPU trade(nAgent, nAgent,value2);
	MatrixCPU LAMBDA(nAgent, nAgent, value3);
	ADMMConst a;

	a.updateBt1(&Bt1, &trade, rho, &LAMBDA);

	Bt11.set(&trade);
	Bt11.subtractTrans(&trade);
	Bt11.multiply(0.5 * rho); 
	Bt11.subtract(&LAMBDA);
	Bt11.divide(rho);

	return Bt1.isEqual(&Bt11);
}

bool testADMMFunction2()
{
	/*void ADMMConst::updateBt2(MatrixCPU* Bt2, MatrixCPU* Tlocal, MatrixCPU* Tmoy, MatrixCPU* P, MatrixCPU* MU)
{
	Bt2->set(Tlocal);
	Bt2->subtractVector(Tmoy);
	Bt2->addVector(P);
	Bt2->subtractVector(MU);
}*/
	int nAgent = 2;
	float value1 = 0;
	float value2 = 0;
	float value3 = 0;
	float value4 = -4;
	float value5 = 4;
	
	MatrixCPU Bt2(nAgent,nAgent,value1);
	MatrixCPU Bt22(nAgent, nAgent, value1);
	MatrixCPU Tlocal(nAgent, nAgent, value2);
	MatrixCPU Tmoy(nAgent, 1, value3);
	MatrixCPU P(nAgent, 1, value4);
	MatrixCPU MU(nAgent, 1, value5);
	ADMMConst a;
	a.updateBt2(&Bt2, &Tlocal, &Tmoy, &P, &MU);


	Bt22.set(&Tlocal);
	Bt22.subtractVector(&Tmoy);
	Bt22.addVector(&P);
	Bt22.subtractVector(&MU);


	return Bt2.isEqual(&Bt22);
}

bool testADMMFunction3()
{
	/*float ADMMConst::calcRes( MatrixCPU* Tlocal, MatrixCPU* Tlocal_pre, MatrixCPU* Tmoy, MatrixCPU* P)
{
	MatrixCPU temp(*Tlocal);
	temp.subtract(Tlocal_pre);

	MatrixCPU temp2(*Tmoy);
	temp2.subtract(P);
	float d1 = temp.max2();
	float d2 = temp2.max2();

	return d1 * (d1 > d2) + d2 * (d2 >= d1);
}*/
	int nAgent = 4;
	float value1 = 2;
	float value2 = 4;
	float value3 = -2;
	float value4 = 3;
	float d = 0;
	MatrixCPU Tlocal(nAgent, nAgent, value1);
	MatrixCPU Tlocal_pre(nAgent, nAgent, value2);
	MatrixCPU Tmoy(nAgent, 1, value3);
	MatrixCPU P(nAgent, 1, value4);
	ADMMConst a;
	d = a.calcRes(&Tlocal, &Tlocal_pre, &Tmoy, &P);
	
	float d1 = fabs(value1 - value2);
	float d2 = fabs(value3 - value4);
	


	
	return (d==(d1 * (d1 > d2) + d2 * (d2 >= d1)));
}

bool testADMMFunction4()
{
	/*
	void ADMMConst::updateTl(MatrixCPU* Tlocal, float at1, float at2, MatrixCPU* Bt1, MatrixCPU*Bt2, MatrixCPU* Ct, MatrixCPU* matLb, MatrixCPU* matUb)
{

	float ada = at1 / at2; // pourrait être prècalculé
	float apa = at1 + at2;

	Tlocal->set(Bt1);
	Tlocal->multiply(ada);
	Tlocal->add(Bt2);
	Tlocal->multiply(at2);
	Tlocal->subtract(Ct);
	Tlocal->divide(apa); 
	Tlocal->project(matLb, matUb);
}*/

	int nAgent = 4;
	float value1 = 2;
	float value2 = -8;
	float value3 = 1;
	float value4 = 1;
	float value5 = -1;
	float value6 = -30;
	float value7 = 0;
	float value8 = (value1 * value3 + value2 * value4 - value5) / (value3 + value4);
	if (value8 > value7) {
		value8 = value7;
	}
	else if (value8 < value6) {
		value8 = value6;
	}
	MatrixCPU Bt1(nAgent,nAgent,value1);
	float at1 = value3;
	MatrixCPU Bt2(nAgent, nAgent, value2);
	float at2 = value4;
	MatrixCPU Ct(nAgent, nAgent, value5);
	MatrixCPU Lb(nAgent, nAgent, value6);
	MatrixCPU Ub(nAgent, nAgent, value7);
	MatrixCPU Tlocal(nAgent, nAgent);
	MatrixCPU Tlocal2(nAgent, nAgent);
	ADMMConst a;
	a.updateTl(&Tlocal, at1, at2, &Bt1, &Bt2, &Ct, &Lb, &Ub);
	
	float ada = at1 / at2; 
	float apa = at1 + at2;

	Tlocal2.set(&Bt1);
	Tlocal2.multiply(ada);
	Tlocal2.add(&Bt2);
	Tlocal2.multiply(at2);
	Tlocal2.subtract(&Ct);
	Tlocal2.divide(apa);
	Tlocal2.project(&Lb, &Ub);
	
	
	MatrixCPU test(nAgent, nAgent, value8);

	return (Tlocal2.isEqual(&Tlocal) && (Tlocal.isEqual(&test)));
}

bool testADMMFunction5()
{
	/*
	void ADMMConst::updateP(MatrixCPU* P, MatrixCPU* Ap1, MatrixCPU* Ap12, MatrixCPU* Bp1, MatrixCPU* Cp, MatrixCPU* Pmin, MatrixCPU* Pmax)
{
	P->multiplyT(Ap1, Bp1);
	P->subtract(Cp);
	
	P->divideT(Ap12);
	P->project(Pmin, Pmax);
}*/
	int nAgent = 2;
	float value1 = 0.5;
	float value2 = 1;
	float value3 = 1;
	float value5 = 8;
	float value6 = -30;
	float value7 = 0;
	float value8 = (value1 * value2 - value5) / (value3 + value2);
	if (value8 > value7) {
		value8 = value7;
	}
	else if (value8 < value6) {
		value8 = value6;
	}
	MatrixCPU Bp1(nAgent, 1, value1);
	MatrixCPU Ap1(nAgent, 1, value2);
	MatrixCPU Ap2(nAgent, 1, value3);
	MatrixCPU Ap12(nAgent, 1, value3);
	Ap12.add(&Ap1, &Ap2);
	MatrixCPU Cp(nAgent, 1, value5);
	MatrixCPU Lb(nAgent, 1, value6);
	MatrixCPU Ub(nAgent, 1, value7);
	MatrixCPU P(nAgent, 1);
	ADMMConst a;
	a.updateP(&P, &Ap1, &Ap12, &Bp1, &Cp, &Lb, &Ub);
	
	MatrixCPU test(nAgent, 1, value8);
	return test.isEqual(&P);
}

bool testADMMFunction6()
{
	return true;
}
