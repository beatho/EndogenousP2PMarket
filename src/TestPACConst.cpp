#include "../head/TestPACConst.h"

int testPACConst()
{
	int n = 1;
	if (!testPACConstContruct1()) return n;
	n++;
	if (!testPACConstContruct2()) return n;
	n++;
	if (!testPACConstContruct3()) return n;
	n++;
	if (!testPACConstFunction1()) return n;
	n++;
	if (!testPACConstFunction2()) return n;
	n++;
	if (!testPACConstFunction3()) return n;
	n++;
	if (!testPACConstFunction4()) return n;
	n++;
	if (!testPACConstFunction5()) return n;
	n++;
	if (!testPACConstFunction6()) return n;
	n++; // 10
	std::cout << "test solve" << std::endl;
	if (!testPACConstSolve1()) return n;
	n++;
	if (!testPACConstSolve2()) return n;
	n++;
	if (!testPACConstSolve3()) return n;
	n++;
	if (!testPACConstSolve4()) return n;
	n++;

	return 0;
}

bool testPACConstContruct1()
{
	std::cout << "default constructor" << std::endl;
	PACConst a;
	return true;
}

bool testPACConstContruct2()
{
	float rho = 2;

	std::cout << "param constructor" << std::endl;
	PACConst a(rho);
	return true;
}
bool testPACConstContruct3()
{
	float rho = 2;

	std::cout << "2 times constructor" << std::endl;
	PACConst a;
	a = PACConst(rho);
	return true;
}

bool testPACConstSolve1()
{
	//solve(Simparam* result, Simparam sim, StudyCase cas);
	StudyCase cas;
	cas.Set2node();

	cas.display(0);
	//cas.display(1);

	int nAgent = cas.getNagent();
	
	Simparam param(nAgent, 1);
	float epsG = 0.0001f;
	int iterG = 1000;
	param.setEpsG(epsG);
	param.setItG(iterG);
	Simparam res(param);
	
	float lambdaMax = 4.7;
	float lambdaMin = 0.5858;
	/*float rho = 1.5; // 1.5 fonctionne 
	float alpha = 0.9;

	float gammarho = 1 / (rho * lambdaMax * alpha); 
	float gammahatrho = 0.9 * gammarho; 
	float theta = gammahatrho * (gammahatrho <1) + 0.9 * (gammahatrho >= 1);
	float phi = theta;*/

	PACConst a;
	a.setBestRhoGamma(lambdaMax, lambdaMin, cas);
	/*a.setParam(rho);
	a.setGamma(gammarho / rho);
	a.setInitCoef(alpha, phi, theta);
	a.setGammahat(gammahatrho / rho);*/

	MatrixCPU Trade(nAgent, nAgent);
	Trade.set(0, 1, -1);
	Trade.set(1, 0, 1);
	a.solve(&res, param, cas);

	MatrixCPU trade = res.getTrade();
	res.display();
	trade.display();
	return trade.isEqual(&Trade,0.01);
}
bool testPACConstSolve2()
{
	StudyCase cas;
	cas.Set4Agent();

	//cas.display(0);
	//cas.display(1);

	int nAgent = cas.getNagent();

	Simparam param(nAgent, 1);
	float epsG = 0.0001f;
	int iterG = 10000;
	param.setEpsG(epsG);
	param.setItG(iterG);
	Simparam res(param);


	float lambdaMax = 6.4;
	float lambdaMin = 0.1535;
	/*float rho = 0.001; // 0.001 fonctionne en 30 000 itérations


	float gammarho = 1 / (rho * lambdaMax); // 1
	float gammahatrho = 1 * gammarho;*/ //0.8

	PACConst a;
	a.setBestRhoGamma(lambdaMax, lambdaMin, cas);

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

bool testPACConstSolve3()
{
	//solve(Simparam* result, Simparam sim, StudyCase cas);
	StudyCase cas;
	cas.Set29node(); // valeur propre max de 24
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;

	int nAgent = cas.getNagent();

	Simparam param(nAgent, cas.getNLine());
	float epsG = 0.00005f;
	int iterG = 100000;
	param.setEpsG(epsG);
	param.setItG(iterG);

	Simparam res(param);
	float lambdaMax = 24;
	float lambdaMin = 0.12;
	/*float rho = 0.001; // 0.001 fonctionne en 30 000 itérations


	float gammarho = 1 /(rho* lambdaMax); // 1
	float gammahatrho = 1 * gammarho;*/ //0.8
	PACConst a;
	a.setBestRhoGamma(lambdaMax, lambdaMin, cas);
	/*a.setParam(rho);
	a.setGamma(gammarho / rho);
	a.setGammahat(gammahatrho / rho);*/
	a.solve(&res, param, cas);

	res.display();


	float Pn[31] = { -0.927, -4.51, -2.81, -0.795, -0.880, -0.0980, -0.128, -4.07, -3.04, -2.16, -0.593, -3.31, -0.389, -2.71, -2.93, -1.77, -0.490, -2.26, -1.05, -0.138, -2.24, 4.246, 5.19, 3.42, 3.75, 4.41, 2.34, 3.31, 2.57, 3.49, 4.23 };

	MatrixCPU P(31, 1);
	for (int i = 0; i < 31; i++) {
		P.set(i, 0, Pn[i]);
	}
	MatrixCPU P2 = res.getPn();


	return P2.isEqual(&P, 0.1);
}

bool testPACConstSolve4()
{
	//solve(Simparam* result, Simparam sim, StudyCase cas);
	StudyCase cas;
	float lim = 0.8;
	cas.Set2nodeConstraint(lim);
	int nAgent = cas.getNagent();
	float epsG = 0.0001f;
	int iterG = 1000;
	Simparam param(nAgent, 1);
	param.setEpsG(epsG);
	param.setItG(iterG);
	Simparam res(param);

	float lambdaMax = 4.7;
	float lambdaMin = 0.5858;
	/*float rho = 1.5; // 1.5 fonctionne
	float alpha = 0.9;

	float gammarho = 1 / (rho * lambdaMax * alpha);
	float gammahatrho = 0.9 * gammarho;
	float theta = gammahatrho * (gammahatrho <1) + 0.9 * (gammahatrho >= 1);
	float phi = theta;*/

	PACConst a;
	a.setBestRhoGamma(lambdaMax, lambdaMin, cas);
	
	

	float value = (1 - lim) * (lim > 1) + lim;


	MatrixCPU Trade(nAgent, nAgent);
	Trade.set(0, 1, -value);
	Trade.set(1, 0, value);
	a.solve(&res, param, cas);

	MatrixCPU trade = res.getTrade();
	res.display();
	trade.display();
	return trade.isEqual(&Trade, 0.001);
}


bool testPACConstFunction1()
{
	

	int nAgent = 2;
	float value1 = 2;
	float value2 = 4;
	float value3 = -2;
	float rho = 1.5;
	

	return true;
}

bool testPACConstFunction2()
{
	
	int nAgent = 2;
	float value1 = 0;
	float value2 = 0;
	float value3 = 0;
	float value4 = -4;
	float value5 = 4;
	
	


	return true;
}

bool testPACConstFunction3()
{
	/*float PACConstConst::calcRes( MatrixCPU* Tlocal, MatrixCPU* Tlocal_pre, MatrixCPU* Tmoy, MatrixCPU* P)
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
	PACConst a;
	//d = a.calcRes(&Tlocal, &Tlocal_pre, &Tmoy, &P);
	
	float d1 = fabs(value1 - value2);
	float d2 = fabs(value3 - value4);
	


	
	return true; //(d==(d1 * (d1 > d2) + d2 * (d2 >= d1)));
}

bool testPACConstFunction4()
{
	

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
	
	
	MatrixCPU Lb(nAgent, nAgent, value6);
	MatrixCPU Ub(nAgent, nAgent, value7);
	MatrixCPU Tlocal(nAgent, nAgent);
	
	PACConst a;
	
	
	MatrixCPU test(nAgent, nAgent, value8);

	return true;
}

bool testPACConstFunction5()
{
	return true;
}

bool testPACConstFunction6()
{
	return true;
}
