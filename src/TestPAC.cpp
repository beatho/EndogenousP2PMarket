#include "../head/TestPAC.h"

int testPAC()
{
	int n = 1;
	if (!testPACContruct1()) return n;
	n++;
	if (!testPACContruct2()) return n;
	n++;
	if (!testPACContruct3()) return n;
	n++;
	if (!testPACFunction1()) return n;
	n++;
	if (!testPACFunction2()) return n;
	n++;
	if (!testPACFunction3()) return n;
	n++;
	if (!testPACFunction4()) return n;
	n++;
	if (!testPACFunction5()) return n;
	n++;
	if (!testPACFunction6()) return n;
	n++; // 10
	std::cout << "test solve" << std::endl;
	if (!testPACSolve1()) return n;
	n++;
	if (!testPACSolve2()) return n;
	n++;
	if (!testPACSolve3()) return n;
	n++;
	

	return 0;
}

bool testPACContruct1()
{
	std::cout << "default constructor" << std::endl;
	PAC a;
	return true;
}

bool testPACContruct2()
{
	float rho = 2;

	std::cout << "param constructor" << std::endl;
	PAC a(rho);
	return true;
}
bool testPACContruct3()
{
	float rho = 2;

	std::cout << "2 times constructor" << std::endl;
	PAC a;
	a = PAC(rho);
	return true;
}

bool testPACSolve1()
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
	
	float lambdaMax = 4.7f;
	float lambdaMin = 0.5858f;
	/*float rho = 1.5; // 1.5 fonctionne 
	float alpha = 0.9;

	float gammarho = 1 / (rho * lambdaMax * alpha); 
	float gammahatrho = 0.9 * gammarho; 
	float theta = gammahatrho * (gammahatrho <1) + 0.9 * (gammahatrho >= 1);
	float phi = theta;*/

	PAC a;
	std::cout << "test" << std::endl;
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
	return trade.isEqual(&Trade, 0.01f);
}
bool testPACSolve2()
{
	StudyCase cas;
	cas.Set4Agent();

	cas.display();

	int nAgent = cas.getNagent();

	Simparam param(nAgent, 1);
	float epsG = 0.0001f;
	int iterG = 100000;
	param.setEpsG(epsG);
	param.setItG(iterG);
	Simparam res(param);


	float lambdaMax = 6.4f;
	float lambdaMin = 0.1535f;
	/*float rho = 0.001; // 0.001 fonctionne en 30 000 it�rations


	float gammarho = 1 / (rho * lambdaMax); // 1
	float gammahatrho = 1 * gammarho;*/ //0.8

	PAC a;
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
	return P.isEqual(&P2, 0.01f);
}

bool testPACSolve3()
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
	
	float rho = 0.1f; // 0.001 fonctionne en 30 000 it�rations
	PAC a;
	//param.setRho(rho);
	
	float lambdaMax = 23;
	float lambdaMin = 0.1162f;

	float gammarho = 1 /(rho* lambdaMax); // 1
	float gammahatrho = 1 * gammarho; //0.8
	
	a.setBestRhoGamma(lambdaMax, lambdaMin, cas);
	/*
	a.setGamma(gammarho / rho);
	a.setGammahat(gammahatrho / rho);/**/
	a.solve(&res, param, cas);

	res.display();


	float Pn[31] = { -0.927f, -4.51f, -2.81f, -0.795f, -0.880f, -0.0980f, -0.128f, -4.07f, -3.04f, -2.16f, -0.593f, -3.31f, -0.389f, -2.71f, -2.93f, -1.77f, -0.490f, -2.26f, -1.05f, -0.138f, -2.24f, 4.246f, 5.19f, 3.42f, 3.75f, 4.41f, 2.34f, 3.31f, 2.57f, 3.49f, 4.23f };

	MatrixCPU P(31, 1);
	for (int i = 0; i < 31; i++) {
		P.set(i, 0, Pn[i]);
	}
	MatrixCPU P2 = res.getPn();
	P2.display();

	return P2.isEqual(&P, 0.1f);
}

bool testPACSolve4()
{
	//solve(Simparam* result, Simparam sim, StudyCase cas);
	StudyCase cas;
	float lim = 0.8f;
	cas.Set2nodeConstraint(lim);
	int nAgent = cas.getNagent();
	Simparam param(nAgent, 1);
	Simparam res(param);


	float value = (1 - lim) * (lim > 1) + lim;

	PAC a;

	MatrixCPU Trade(nAgent, nAgent);
	Trade.set(0, 1, -value);
	Trade.set(1, 0, value);
	a.solve(&res, param, cas);

	MatrixCPU trade = res.getTrade();
	res.display();
	trade.display();
	return trade.isEqual(&Trade, 0.001f);
}


bool testPACFunction1()
{
	

	int nAgent = 2;
	float value1 = 2;
	float value2 = 4;
	float value3 = -2;
	float rho = 1.5f;
	

	return true;
}

bool testPACFunction2()
{
	
	int nAgent = 2;
	float value1 = 0;
	float value2 = 0;
	float value3 = 0;
	float value4 = -4;
	float value5 = 4;
	
	


	return true;
}

bool testPACFunction3()
{
	/*float PACConst::calcRes( MatrixCPU* Tlocal, MatrixCPU* Tlocal_pre, MatrixCPU* Tmoy, MatrixCPU* P)
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
	PAC a;
	//d = a.calcRes(&Tlocal, &Tlocal_pre, &Tmoy, &P);
	
	float d1 = fabs(value1 - value2);
	float d2 = fabs(value3 - value4);
	


	
	return true; //(d==(d1 * (d1 > d2) + d2 * (d2 >= d1)));
}

bool testPACFunction4()
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
	
	PAC a;
	
	
	MatrixCPU test(nAgent, nAgent, value8);

	return true;
}

bool testPACFunction5()
{
	return true;
}

bool testPACFunction6()
{
	return true;
}
