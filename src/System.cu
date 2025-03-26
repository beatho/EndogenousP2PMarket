#include "../head/System.h"
#include "../head/System.cuh"


// To DO MatrixCPU genererP0(path,dateMonth)
// To DO updateCas(P0)

System::System() {
#ifdef DEBUG_CONSTRUCTOR
	std::cout << " system constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR

	
	_simparam = Simparam(_case.getNagent(), _case.getNLine());
	_methode = nullptr;
	_result = new Simparam(_simparam);

}


System::System(float rho, int iterMaxGlobal, int iterMaxLocal, float epsGlobal, float epsLocal, std::string nameMethode, int nAgent, float P, float dP, float a, float da, float b, float db)
{
	
	int _nAgent = nAgent;
	if (nAgent == 0) {
		_case = StudyCase();
	}
	else {
		_case = StudyCase(nAgent,  P,  dP,  a,  da,  b,  db);
	}
	setMethod(nameMethode);
	

	_simparam = Simparam(rho, iterMaxGlobal, iterMaxLocal, epsGlobal, epsLocal, _nAgent);
	_result = new Simparam(_simparam);
	
}
System::~System()
{
	DELETEB(_methode);
	DELETEB(_methodePF)
	DELETEB(_methodePFGPU);
	DELETEB(_result);
}

Simparam System::solve()
{
	_methode->solve(_result, _simparam, _case);
 	return *_result;
}

Simparam System::solvePF()
{
	if(usePFGPU){
		MatrixGPU PQ = MatrixGPU(_case.getPobj(), 1);
		MatrixGPUD PQD = MatrixGPUD(_case.getPobjD(), 1);

		_methodePFGPU->init(_case, &PQ, &PQD, useDoublePF);
		_methodePFGPU->solve();
		_methodePFGPU->display(true);
	} else{
		MatrixCPU PQ = _case.getPobj();
		MatrixCPUD PQD = _case.getPobjD();

		_methodePF->init(_case, &PQ, &PQD, useDoublePF);
		_methodePF->solve();
		_methodePF->display(true);
	}


    return *_result;
}

ResultInterface* System::solve(ResultInterface* res, ParamInterface* param, StudyCaseInterface* caseInter, bool AC){
	
	
	if(AC){
		_case.SetACStudyCaseFromInterface(caseInter);
	}else{
		_case.SetDCStudyCaseFromInterface(caseInter);
	}
	_simparam.setFromInterface(param, AC);
	*_result = _simparam;
	
	_methode->solve(_result, _simparam, _case);
 	
	_result->convertToResultInterface(res);

	return res;
}
ResultInterface* System::solvePF(ResultInterface* res, ParamInterface* param, StudyCaseInterface* caseInter){
	
	MatrixCPU sizes(param->getSize());
	_case.SetACStudyCaseFromInterface(caseInter);
	if(usePFGPU){
		MatrixGPU PQ = MatrixGPU(_case.getPobj(), 1);
		MatrixGPUD PQD = MatrixGPUD(_case.getPobjD(), 1);

		_methodePFGPU->init(_case, &PQ, &PQD, useDoublePF);
		_methodePFGPU->solve();
		_methodePFGPU->display(true);
		
		res->setvarPhysic(_methodePFGPU->getW(), MatrixCPU(2*sizes.get(0, nLineP_ind), 1), _methodePFGPU->getE());
		res->setResult(_methodePFGPU->getIter(), 1, _methodePFGPU->getTime(), _methodePFGPU->getRes(), MatrixCPU(0,0));

	} else{
		MatrixCPU PQ = _case.getPobj();
		MatrixCPUD PQD = _case.getPobjD();

		_methodePF->init(_case, &PQ, &PQD, useDoublePF);
		_methodePF->solve();
		_methodePF->display(true);

		res->setvarPhysic(_methodePF->getW(), MatrixCPU(2*sizes.get(0, 2), 1), _methodePF->getE());
		res->setResult(_methodePF->getIter(), 1, _methodePF->getTime(), _methodePF->getRes(), MatrixCPU(0,0));
	}


    return res;
}


void System::solveIntervalle(std::string path, MatrixCPU* interval, int nCons, int nGen)
{
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	if (interval->getNCol() != 2 || interval->getNLin() != 4) {
		throw std::invalid_argument("interval must be 4*2, year, month, day, hour");
	}
	
	int Nsimu = getNbSimu(interval);

	std::cout << "Simulation count " <<Nsimu << std::endl;

	
	bool fin = false;
	bool bissextile = false;
	int year = interval->get(0, 0);
	if ((year%4 == 0 && year%100 != 0) || (year % 400 == 0)) {
		bissextile = true;
	}
	int month = interval->get(1, 0);
	int day = interval->get(2, 0);
	int hour = interval->get(3, 0);
	int dayl = dayMonth[month-1];
	if (month == 2 && bissextile) {
		dayl = 29;
	}

		
	int Nhour = 24 * dayl;
	std::string date = generateMonth(year, month);
	

	MatrixCPU P0Global(nCons, Nhour);
	
	MatrixCPU P0(nCons, 1);
	int indiceP0 = (day-1) * 24 + hour;
	
	StudyCase cas;
	_temps = MatrixCPU(1,Nsimu);
	_iter = MatrixCPU(1, Nsimu);
	_conv= MatrixCPU(1, Nsimu,-1);
	_fc = MatrixCPU(1, Nsimu);
	_ResR = MatrixCPU(1, Nsimu);
	_ResS = MatrixCPU(1, Nsimu);
	_ResX = MatrixCPU(1, Nsimu);
	int indice = 0;
	int stepG = _simparam.getStepG();
	std::string nameP0;
	
	nameP0 = path + "/load/Month/";
	generateP0(&P0Global, nameP0, date);
	P0Global.getBloc(&P0, 0, nCons, indiceP0, indiceP0 + 1);
	//P0.display();
	indiceP0 = indiceP0 + 1;
	cas.setReduce(true);
	
	cas.SetEuropeP0(path, &P0, true);
	if (_simparam._lineLimitMin) {
		cas.setLineLimitMin(_simparam._lineLimitMin); // on met une valeur minimale pour verifier si cela aide la convergence
	}
	else if (_simparam.offsetConstraint) {
		cas.setLineLimitRelaxation(_simparam.offsetConstraint);
	}
	
	//cas.display();
	 
	setStudyCase(cas);
	 
	//resetMethod();
	//display(1);
	float epsG = _simparam.getEpsG();
	clock_t t = 0;
	std::cout << "-";
	/*MatrixGPU Sensi(cas.getPowerSensi(), 1);
	MatrixGPU g(cas.getLineLimit());
	g.saveCSV("PowerLineEuropeHugeOffsetJunePowerTech.csv", mode, 1);
	g.transferGPU();*/
	while (!fin) {
		
		if (indice == (Nsimu - 1)) {
			fin = true;
		}
		
		/// simulation
		
		t = clock();
		solve();
		t = clock() - t;
		std::cout << "-";
		_temps.set(0, indice, (float)t / CLOCKS_PER_SEC);

		/*if (cas.getNLine() != 0) {
			MatrixGPU Pn(_result->getPn(), 1);
			g.multiply(&Sensi, &Pn);
			g.transferCPU();
			g.saveCSV("PowerLineEuropeOffsetJune.csv", mode, 1);
			g.transferGPU();
		}*/

		_iter.set(0, indice, _result->getIter());
		_fc.set(0, indice, _result->getFc());
		int iter = _result->getIter();
		MatrixCPU Res(_result->getRes());
		//Res.saveCSV("ResRhoRhoRho.csv", std::fstream::in | std::fstream::out | std::fstream::app);
		float resR = Res.get(0, (iter - 1) / stepG);
		float resS = Res.get(1, (iter - 1) / stepG);
		float resX = Res.get(2, (iter - 1) / stepG);
		_ResR.set(0, indice, resR);
		_ResS.set(0, indice, resS);
		_ResX.set(0, indice, resX);
		if (resR <= epsG && resS <= epsG && resX <= epsG) {
			_conv.set(0, indice, 1);
		}
		else {
			_conv.set(0, indice, 0);
		}
		if (!fin) {
			t = clock();
			/*_simparam.setTrade(&(res.getTrade()));
			_simparam.setLAMBDA(&(res.getLambda()));
			_simparam.setMU(&(res.getMU()));
			_simparam.setPn(&res.getPn());*/
			hour = hour + 1;
			if (hour == 24) {
				std::cout << day << "\n";
				day = day + 1;
				hour = 0;
				if (day > dayl) {
					day = 1;
					month = month + 1;
					if (month > 12) {
						year = year + 1;
						month = 1;
						if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)) {
							bissextile = true;
						}
					}
					dayl = dayMonth[month - 1];
					if (month == 2 && bissextile) {
						dayl = 29;
					}
					Nhour = 24 * dayl;
					date = generateMonth(year, month);
					P0Global.setSize(nCons, Nhour);
					generateP0(&P0Global, nameP0, date);
					indiceP0 = 0;
				}
			}
			P0Global.getBloc(&P0, 0, nCons, indiceP0, indiceP0 + 1);
			
			cas.UpdateP0(&P0);
			setStudyCase(cas);
			if (_simparam._warmstart) {
				UpdateP0();
			}
			
			indiceP0 = indiceP0 + 1;
			indice = indice + 1;
		}
		
	}
	std::cout << std::endl;
	
	std::cout << "times sum : " << _temps.sum() << std::endl;
	//std::cout << "simulation count : " << Nsimu << std::endl;
	/*std::cout << "iter :" << std::endl;
	_iter.display();
	//std::cout << "times :" << std::endl;
	//_temps.display();
	//std::cout << "conv :" << std::endl;
	//_conv.display();
	//std::cout << "fc : " << std::endl;
	//_fc.display();
	std::cout << "Res_R : " << std::endl;
	_ResR.display();
	std::cout << "Res_S : " << std::endl;
	_ResS.display();
	std::cout << "Res_X : " << std::endl;
	_ResX.display();*/
	
}


void System::solveIntervalle(std::string path, std::string name, MatrixCPU* interval)
{

	if (interval->getNCol() != 2 || interval->getNLin() != 4) {
		throw std::invalid_argument("interval must be 4*2, year, month, day, hour");
	}
	int Nsimu = getNbSimu(interval);

	std::cout << "Simulation count " << Nsimu << std::endl;


	bool fin = false;
	bool bissextile = false;
	int year = interval->get(0, 0);
	if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)) {
		bissextile = true;
	}
	int month = interval->get(1, 0);
	int day = interval->get(2, 0);
	int hour = interval->get(3, 0);
	int dayl = dayMonth[month - 1];
	if (month == 2 && bissextile) {
		dayl = 29;
	}


	int Nhour = 24 * dayl;
	std::string date = generateMonth(year, month);


	
	int indiceP0 = (day - 1) * 24 + hour;

	StudyCase cas;
	
	_temps = MatrixCPU(1, Nsimu);
	_iter = MatrixCPU(1, Nsimu);
	_conv = MatrixCPU(1, Nsimu, -1);
	_fc = MatrixCPU(1, Nsimu);
	_ResR = MatrixCPU(1, Nsimu);
	_ResS = MatrixCPU(1, Nsimu);
	_ResX = MatrixCPU(1, Nsimu);
	int indice = 0;
	int stepG = _simparam.getStepG();
	std::string nameP0;

	nameP0 = path + "/load/Month/";

	int nCons = getNFileline(nameP0 + date + ".txt"); // long...
		
	MatrixCPU P0Global(nCons, Nhour);
	MatrixCPU P0(nCons, 1);
	generateP0(&P0Global, nameP0, date);
	P0Global.getBloc(&P0, 0, nCons, indiceP0, indiceP0 + 1);
	//P0.display();
	indiceP0 = indiceP0 + 1;
	cas.SetStudyCase(path, name, &P0, false);
	cas.setReduce(true);
	//cas.display();

	//cas.display(2);
	setStudyCase(cas);
	try
	{
		_methode->setBestParam(cas);
	}
	catch (const std::exception&)
	{
		// juste une m�thode o� ce n'est pas implement�
	}
		
	
	//_case.display();
	//resetMethod(); grrrrrrrrrrrrrrrrrrrrrrrrr probl�me avec les matrices qui reste sur GPU
	//display(1);
	float epsG = _simparam.getEpsG();
	clock_t t = 0; 
	std::cout << "-";
	while (!fin) {

		if (indice == (Nsimu - 1)) {
			fin = true;
		}

		/// simulation

		t = clock();
		Simparam res;
		try
		{
			res = solve();
		}
		catch (const std::exception& e)
		{
			std::cout << e.what() << std::endl;
			return;
		}
		
		t = clock() - t;
		std::cout << "-";
		_temps.set(0, indice, (float)t / CLOCKS_PER_SEC);
		/*std::cout << "Echange entre les agents " << std::endl;
		//displayTradesAgent();
		std::cout << "flux dans les lignes " << std::endl;
		if (cas.getNLine() != 0) {
			MatrixCPU Sensi(cas.getPowerSensi());
			MatrixCPU Pn(res.getPn());
			MatrixCPU g(cas.getNLine(), 1, 0);
			g.multiply(&Sensi, &Pn);
			//cas.displayLineCores(&g);
		}
		std::cout << "---------------------- - " <<std::endl;*/
		//MatrixCPU lambda(res.getLambda());
		
		//lambda.display();

		_iter.set(0, indice, res.getIter());
		_fc.set(0, indice, res.getFc());
		int iter = res.getIter();
		MatrixCPU Res(res.getRes());
		//Res.saveCSV("ResRhoRhoRho.csv", std::fstream::in | std::fstream::out | std::fstream::app);
		float resR = Res.get(0, (iter - 1) / stepG);
		float resS = Res.get(1, (iter - 1) / stepG);
		float resX = Res.get(2, (iter - 1) / stepG);
		_ResR.set(0, indice, resR);
		_ResS.set(0, indice, resS);
		_ResX.set(0, indice, resX);
		if (resR <= epsG && resS <= epsG) {
			_conv.set(0, indice, 1);
		}
		else {
			_conv.set(0, indice, 0);
		}
		if (!fin) {
			t = clock();
			hour = hour + 1;
			if (hour == 24) {
				std::cout << "\n";
				day = day + 1;
				hour = 0;
				if (day > dayl) {
					day = 1;
					month = month + 1;
					if (month > 12) {
						year = year + 1;
						month = 1;
						if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)) {
							bissextile = true;
						}
					}
					dayl = dayMonth[month - 1];
					if (month == 2 && bissextile) {
						dayl = 29;
					}
					Nhour = 24 * dayl;
					date = generateMonth(year, month);
					P0Global.setSize(nCons, Nhour);
					generateP0(&P0Global, nameP0, date);
					indiceP0 = 0;
				}
			}
			P0Global.getBloc(&P0, 0, nCons, indiceP0, indiceP0 + 1);
			_case.UpdateP0(&P0);
			//setStudyCase(cas);
			if (_simparam._warmstart) {
				UpdateP0();
			}
			indiceP0 = indiceP0 + 1;
			indice = indice + 1;
		}

	}
	std::cout << std::endl;

	std::cout << "times sum : " << _temps.sum() << std::endl;
	//std::cout << "simulation count : " << Nsimu << std::endl;
	std::cout << "iter :" << std::endl;
	_iter.display();
	//std::cout << "times :" << std::endl;
	//_temps.display();
	//std::cout << "conv :" << std::endl;
	//_conv.display();
	//std::cout << "fc : " << std::endl;
	//_fc.display();
	/*std::cout << "Res_R : " << std::endl;
	_ResR.display();
	std::cout << "Res_S : " << std::endl;
	_ResS.display();
	std::cout << "Res_X : " << std::endl;
	_ResX.display();*/

}


void System::solveIntervalle(std::string path, int begin, int end, int chosenAgentGen)
{
	if (begin > end || begin < 0 || end > 60 * 24) {
		throw std::invalid_argument("begin and end must be withinn a day with one minute step");
	}
	int Nsimu = end - begin + 1;
	std::cout << "Simulation count " << Nsimu << std::endl;

	StudyCase cas;
	
	cas.SetEuropeTestFeeder(path, chosenAgentGen, begin);
	_simparam.setNAgentLine(cas.getNagent(), cas.getNLine(), cas.isAC());
	_result->setNAgentLine(cas.getNagent(), cas.getNLine(), cas.isAC());

	if (_methode->_name == "ADMMConst1" || _methode->_name == "ADMMGPUConst4") {
		std::cout << "ajout du DC" << std::endl;
		cas.genDCGridFromAC(); // pour utiliser les m�thodes DC
		cas.setReduce(true);
	}
	_temps = MatrixCPU(1, Nsimu);
	_iter = MatrixCPU(1, Nsimu);
	_conv = MatrixCPU(1, Nsimu, -1);
	_fc = MatrixCPU(1, Nsimu);
	_ResR = MatrixCPU(1, Nsimu);
	_ResS = MatrixCPU(1, Nsimu);
	_ResX = MatrixCPU(1, Nsimu);

	float epsG = _simparam.getEpsG();
	float epsGC = _simparam.getEpsGC();
	int stepG = _simparam.getStepG();
	clock_t t = 0;// clock();
	std::cout << "-";
	for (int iter = 0; iter < Nsimu; iter++) {
		t = clock();
		_methode->solve(_result,_simparam, cas);
		t = clock() - t;
		std::cout << "-";
		_temps.set(0, iter, (float) t / CLOCKS_PER_SEC);
		int iterFinal = _result->getIter();
		_iter.set(0, iter, iterFinal);
		_fc.set(0, iter, _result->getFc());
		MatrixCPU Res(_result->getRes());
		//Res.saveCSV("ResRhoRhoRho.csv", std::fstream::in | std::fstream::out | std::fstream::app);
		float resR = Res.get(0, (iterFinal - 1) / stepG);
		float resS = Res.get(1, (iterFinal - 1) / stepG);
		float resX = Res.get(2, (iterFinal - 1) / stepG);
		_ResR.set(0, iter, resR);
		_ResS.set(0, iter, resS);
		_ResX.set(0, iter, resX);
		if (resR <= epsG && resS <= epsG && resX <= epsGC) {
			_conv.set(0, iter, 1);
		}
		else {
			_conv.set(0, iter, 0);
		}
		if (iter < Nsimu - 1) {
			if ((iter + 1) % 60 == 0) {
				std::cout << std::endl;
			}
			cas.nextStepPobj();
			_methode->updateP0(cas);
		}
	}
	std::cout << std::endl;

	std::cout << "times sum : " << _temps.sum() << std::endl;
	//std::cout << "simulation count : " << Nsimu << std::endl;
	std::cout << "iter :" << std::endl;
	_iter.display();
	//std::cout << "times :" << std::endl;
	//_temps.display();
	//std::cout << "conv :" << std::endl;
	//_conv.display();
	//std::cout << "fc : " << std::endl;
	//_fc.display();
	/*std::cout << "Res_R : " << std::endl;
	_ResR.display();
	std::cout << "Res_S : " << std::endl;
	_ResS.display();
	std::cout << "Res_X : " << std::endl;
	_ResX.display();*/
}


void System::UpdateP0()
{
	_methode->updateP0(_case);
}

void System::resetMethod()
{
	
	setMethod(_methode->_name);
}

void System::resetParam()
{
	
	_simparam.setNAgentLine(0, 0, false);
}

void System::removeLink(int i, int j)
{
	_case.removeLink(i, j);
}

void System::addLink(int i, int j)
{
	_case.addLink(i, j);
}

Agent System::removeAgent(int agent)
{
	return _case.removeAgent(agent);
}

void System::restoreAgent(Agent& agent, bool all)
{
	_case.restoreAgent(agent, all);
}

void System::setBestRho(float rhoMax, bool rhoVar, float rhoTest)
{
	int nAgent = _case.getNagent();
	float rhoMin = 0.01;
	if (rhoMax == 0) {
		rhoMax = 0.2 * nAgent;
	}
	float epsRho = 0.01;
	float dRhoMin = 0.005;

	


	float rho_a = rhoMin;
	float rho_b = rhoMax;
	int multiplieur = 1;
	while ((rho_b - rho_a) > epsRho) {
		float dRho = multiplieur * (rho_a + rho_b) / 100;
		dRho = dRho > dRhoMin ? dRho : dRhoMin;
		float rho_x = 0.5 * (rho_a + rho_b);

		float rho_c = rho_x - (dRho / 2);
		float rho_d = rho_x + (dRho / 2);
		setRho(rho_c);
		*_result = solve();
		int iter_c = _result->getIter();
		int iterL_c = _result->getIterLTot();

		setRho(rho_d);
		*_result = solve();
		int iter_d = _result->getIter();
		int iterL_d = _result->getIterLTot();

		if (iter_d < iter_c) {
			rho_a = rho_d;
			multiplieur = 1;
		}
		else if (iter_d > iter_c) {
			rho_b = rho_c;
			multiplieur = 1;
		}
		else {
			std::cout << iter_d << " " << iter_c << " " << iterL_d << " " << iterL_c << std::endl;
			multiplieur++;
		}

	}

	float rho_x = 0.5 * (rho_a + rho_b);
	setRho(rho_x);
	*_result = solve();
	int iter_x = _result->getIter();
	int iterL_x = _result->getIterLTot();
	float time_x = _result->getTime();

	float rho_the = nAgent * 0.05;
	setRho(rho_the);
	*_result = solve();
	int iter_the = _result->getIter();
	int iterL_the = _result->getIterLTot();
	float time_the = _result->getTime();
	int iter_test = iter_the + iter_x; 
	int iterL_test = 0;
	float time_test = 0;
	if (rhoTest != 0 && rhoTest != rho_the && rhoTest != rho_x) {
		setRho(rhoTest);
		*_result = solve();
		iter_test = _result->getIter();
		iterL_test = _result->getIterLTot();
		time_test = _result->getTime();
	}

	float rhoBest = 0;
	if (iter_x < iter_the) {
		if (iter_x < iter_test) {
			rhoBest = rho_x;
		}
		else if (iter_x > iter_test) {
			rhoBest = rhoTest;
		}
		else {
			if (iterL_x <= iterL_test) {
				rhoBest = rho_x;
			}
			else {
				rhoBest = rhoTest;
			}
		}
	} else if (iter_x > iter_the) {
		if (iter_the < iter_test) {
			rhoBest = rho_the;
		}
		else if (iter_the > iter_test) {
			rhoBest = rhoTest;
		}
		else {
			if (iterL_the <= iterL_test) {
				rhoBest = rho_the;
			}
			else {
				rhoBest = rhoTest;
			}
		}
	}
	else {
		if (iter_x > iter_test) {
			rhoBest = rhoTest;	
		}
		else {
			if (iterL_the <= iterL_x) {
				rhoBest = rho_the;
			}
			else {
				rhoBest = rho_x;
			}
		}
	}
	setRho(rhoBest);
	std::cout << "Best rho find is " << rhoBest << std::endl;
	
}

void System::setStudyCase(const StudyCase& cas)
{
	_case = cas;
	if (cas.getNagent() != _simparam.getNAgent() || cas.getNLine() != _simparam.getNLine())
	{
		std::cout << "wrong number of agent or branch, simparam and result update " << std::endl;
		_simparam.setNAgentLine(cas.getNagent(), cas.getNLine(), cas.isAC());
		_result->setNAgentLine(cas.getNagent(), cas.getNLine(), cas.isAC());
		//std::cout << "end update" << std::endl;
	}

}

void System::setStudyCase(std::string fileName)
{
	_case.SetACFromFile(fileName);
	if (_case.getNagent() != _simparam.getNAgent() || _case.getNLine() != _simparam.getNLine())
	{
		std::cout << "wrong number of agent or branch, simparam and result update " << std::endl;
		_simparam.setNAgentLine(_case.getNagent(), _case.getNLine(), _case.isAC());
		_result->setNAgentLine(_case.getNagent(), _case.getNLine(), _case.isAC());
		//std::cout << "end update" << std::endl;
	}

}



void System::setSimparam(const Simparam& param)
{
	_simparam = param;
	DELETEB(_result);
	_result = new Simparam(_simparam);
	if (_case.getNagent() != param.getNAgent())
	{
		std::cout << "wrong number of agent, simparam and result update" << std::endl;
		std::cout << "if it is not wanted change the study case before doing that" << std::endl;
		_simparam.setNagent(_case.getNagent());
		_result->setNagent(_case.getNagent());
	}

}
void System::setMethod(std::string nameMethode) {
	
	DELETEB(_methode);
	useOPF = false;
	if (!nameMethode.compare(sADMMMarket)) {
		_methode = new ADMMMarket;
	}
	else if ((!nameMethode.compare(sADMMMarketMP))) {
		_methode = new ADMMMarketOpenMP;
	}else if ((!nameMethode.compare(sADMMMarketGPU))) {
		_methode = new ADMMMarketGPU;
	}
	else if (!nameMethode.compare(sADMMConst)) {
		_methode = new ADMMConst();
	}
	else if ((!nameMethode.compare(sADMMConst1))) {
		_methode = new ADMMConst1;
	}
	#ifdef OSQP
	else if ((!nameMethode.compare(sOSQPConst))) {
		_methode = new OSQPConst;
	}
	else if ((!nameMethode.compare(sADMMGPUConstCons))) {
		_methode = new ADMMGPUConstCons;
	}
	else if ((!nameMethode.compare(sDCOPFOSQP))) {
		useOPF = true;
		_methode = new DCOPFOSQP;
	}
	#endif
	else if ((!nameMethode.compare(sADMMGPUConst1))) {
		_methode = new ADMMGPUConst1;
	}
	else if ((!nameMethode.compare(sADMMGPUConst1T))) {
		_methode = new ADMMGPUConst1T;
	}
	else if ((!nameMethode.compare(sADMMGPUConst2))) {
		_methode = new ADMMGPUConst2;
	}
	else if ((!nameMethode.compare(sADMMGPUConst3))) {
		_methode = new ADMMGPUConst3;
	}
	else if ((!nameMethode.compare(sADMMGPUConst4))) {
		_methode = new ADMMGPUConst4;
	}
	else if ((!nameMethode.compare(sADMMGPUConst5))) {
		_methode = new ADMMGPUConst5;
	}
	else if ((!nameMethode.compare(sADMMGPUConstCons2))) {
		_methode = new ADMMGPUConstCons2;
	}
	else if ((!nameMethode.compare(sADMMGPUConstCons3))) {
		_methode = new ADMMGPUConstCons3;
	}
	else if ((!nameMethode.compare(sADMMACConst1))) {
		_methode = new ADMMACConst1;
	}
	else if ((!nameMethode.compare(sPAC))) {
		_methode = new PAC;
	}
	else if ((!nameMethode.compare(sPACConst))) {
		_methode = new PACConst;
	}
	else if ((!nameMethode.compare(sOPFADMM))) {
		_methode = new OPFADMM;
		useOPF = true;
	}
	else if ((!nameMethode.compare(sOPFADMMGPU))) {
		_methode = new OPFADMMGPU;
		useOPF = true;
	}else if ((!nameMethode.compare(sOPFADMM2))) {
		_methode = new OPFADMM2;
		useOPF = true;
	}
	else if ((!nameMethode.compare(sOPFADMMGPU2))) {
		_methode = new OPFADMMGPU2;
		useOPF = true;
	}
	else if ((!nameMethode.compare(sEndoMarketCons))){
		_methode = new MarEndoCons;
	}
	else if ((!nameMethode.compare(sEndoMarketConsGPU))){
		_methode = new MarEndoConsGPU;
	}
	else if ((!nameMethode.compare(sEndoMarketDirect))){
		_methode = new MarketEndoDirect;
	}
	else if ((!nameMethode.compare(sEndoMarketDirectGPU))){
		_methode = new MarketEndoDirectGPU;
	}	
	else {
		std::cout << "unknonwn method " << nameMethode << " !" << std::endl;
	}
}

void System::setMethodPF(std::string nameMethode, bool isDouble)
{
	useDoublePF = isDouble;
	if (!nameMethode.compare(sNR)) {
		_methodePF = new CPUPF;
	}
	else if ((!nameMethode.compare(sNRGPU))) {
		_methodePFGPU = new GPUPF;
		usePFGPU = true;
	}else if ((!nameMethode.compare(sGS))) {
		_methodePF = new CPUPFGS;
	}
	else if ((!nameMethode.compare(sGSGPU))) {
		_methodePFGPU = new GPUPFGS;
		usePFGPU = true;
	}
	else if ((!nameMethode.compare(sDistPQ))) {
		_methodePF = new CPUPFdistPQ;
	}
	else if ((!nameMethode.compare(sDistPQGPU))) {
		_methodePFGPU = new GPUPFdistPQ;
		usePFGPU = true;
	}else {
		std::cout << "unknonwn method " << nameMethode << " !" << std::endl;
	}
}
void System::setRho(float rho) {
	_simparam.setRho(rho);
	_result->setRho(rho);
}
void System::setRho1(float rho1)
{
	_simparam.setRho1(rho1);
	_result->setRho1(rho1);
}
void System::setRhoL(float rho)
{
	if (!useOPF) {
		((MethodP2P*) _methode)->setParam(rho);
	}
}

void System::setIter(int iterG, int iterL) {
	_simparam.setItG(iterG);
	_simparam.setItL(iterL);
	_result->setItG(iterG);
	_result->setItL(iterL);
}

void System::setItIntern(int iter)
{
	_simparam.setItIntern(iter);
}

void System::setStep(int stepG, int stepL)
{
	_simparam.setStep(stepG, stepL);
	_result->setStep(stepG, stepL);
}

void System::setStep(int stepG, int stepL, int stepIntern)
{
	_simparam.setStep(stepG, stepL, stepIntern);
}

void System::setEpsG(float epsG)
{
	_simparam.setEpsG(epsG);
}

void System::setEpsGC(float epsgC)
{
	_simparam.setEpsGC(epsgC);
}

void System::setEpsIntern(float eps)
{
	_simparam.setEpsIntern(eps);
}

void System::setEpsL(float epsL)
{
	_simparam.setEpsL(epsL);
}

void System::setTrade(MatrixCPU* trade)
{
	_simparam.setTrade(trade);
	MatrixCPU Pn(trade->getNLin(), 1);
	Pn.sum(trade);
	_simparam.setPn(&Pn);
}
void System::setLineLimitMin(float lineMin) {
	_simparam._lineLimitMin = lineMin;
}

void System::setConstraintRelaxation(float factor)
{
	float epsGC = _simparam.getEpsGC();
	_simparam.offsetConstraint = epsGC * factor;
	//_case.setLineLimitRelaxation(epsGC);
}

void System::setWarmStart(bool warmstart)
{
	_simparam._warmstart = warmstart;
}

MatrixCPU System::getRes() const {
	return _result->getRes();
}
MatrixCPU System::getTrade() const {
	return _result->getTrade();
}

MatrixCPU System::getTemps() const
{
	return _temps;
}

MatrixCPU System::getIter() const
{
	return _iter;
}

MatrixCPU System::getConv() const
{
	return _conv;
}

MatrixCPU System::getFc() const
{
	return _fc;
}

MatrixCPU System::getResR() const
{
	return _ResR;
}

MatrixCPU System::getResS() const
{
	return _ResS;
}



MatrixCPU System::getResX() const
{
	return _ResX;
}

int System::getNTrade() const
{
	return (_case.getNvoi()).sum();
}


void System::display(int type) // type=0 result, type =1 simparam & methode, type=2 case
{
	if (type==1) {
		std::cout << "Simparam : " << std::endl;
		_simparam.display(1);
		
	}
	else if (type == 2)
	{
		std::cout << "Case : " << std::endl;
		_case.display();
		_case.display(1);
		_case.display(2);
	}
	else {
		std::cout << "Method : ";
		_methode->display();
		_result->display();
	}
}

void System::displayTradesAgent()
{
	int N = _case.getNagent();
	int nCons = _case.getNCons();
	int nGen = N - nCons;
	MatrixCPU Trade(_result->getTrade());
	MatrixCPU Pn(_result->getPn());
	for (int n = 0; n < nCons; n++) {
		
		std::cout << "Agent consomateur n " << n << " Puissance echangee " << Pn.get(n,0) *_case.getSbase() << std::endl;
		for (int i = nCons; i < N; i++) {
			if (Trade.get(i, n) > 0.0001) { // positif car on regarde ce que le g�n�rateur vend
				std::cout << "         achete " << Trade.get(i, n) * _case.getSbase() << " MWh au generateur " << i << std::endl;
			}
		}
	}
}
void System::displayTime(std::string fileName) const
{
	_result->displayTime(fileName);
}

void System::setMethod(Method* method)
{
	//DELETEB(_methode); hum cela ne lui appartient pas...
	_methode = method;
}

MatrixCPU System::getPn() const
{
	return _result->getPn();
}

int System::getNFileline(std::string nameFile)
{
	int number_of_lines = 0;
	std::string line;
	std::ifstream myfile(nameFile);

	while (std::getline(myfile, line))
		++number_of_lines;
	return number_of_lines;
}

int System::getNbSimu(MatrixCPU* interval) const
{
	int year1, year2, month1, month2, day1, day2, hour1, hour2;
	year1 = interval->get(0, 0);
	year2 = interval->get(0, 1);
	month1 = interval->get(1, 0);
	month2 = interval->get(1, 1);
	day1 = interval->get(2, 0);
	day2 = interval->get(2, 1);
	hour1 = interval->get(3, 0);
	hour2 = interval->get(3, 1);
	// verifier que c'est possible
	if (year1 > year2) {
		throw std::invalid_argument("date1 must be before date2 (year)");
	}
	else if (year1 == year2) {
		if (month1 > month2) {
			throw std::invalid_argument("date1 must be before date2 (month)");
		}
		else if (month1 == month2) {
			if (day1 > day2) {
				throw std::invalid_argument("date1 must be before date2 (day)");
			}
			else if (day1 == day2) {
				if (hour1 > hour2) {
					throw std::invalid_argument("date1 must be before date2 (hour)");
				}
			}
		}
	}

	int m[12] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 };
	int dayMonth[12] = { 31, 28, 31, 30, 31, 30 , 31, 31, 30, 31, 30, 31 };

	int yearref = 2012;
	int dy = (year1 - yearref);
	int N1 = (dy * 365 + m[month1 - 1] + 1 + day1 - 1) * 24 + hour1;
	N1 = N1 + dy / 4 - dy / 100 + dy / 400;
	if ((dy % 4 == 0 && dy % 100 != 0) || (dy % 400 == 0)) {
		if (month1 < 3) {
			N1 = N1 - 24;
		}
	}

	dy = year2 - yearref;
	int N2 = (dy * 365 + m[month2 - 1] + 1 + day2 - 1) * 24 + hour2;
	N2 = N2 + dy / 4 - dy / 100 + dy / 400;

	if ((dy % 4 == 0 && dy % 100 != 0) || (dy % 400 == 0)) {
		if (month2 < 3) {
			N2 = N2 - 24;
		}
	}
	return N2 - N1 + 1;
}

int System::getNagent() const
{
	return _case.getNagent();
}



std::string System::generateDate(int year, int month, int day, int hour)
{
	std::string smonth;
	std::string sday;
	std::string shour;
	if (month < 10) {
		smonth = "0" + std::to_string(month);
	}
	else {
		smonth = std::to_string(month);
	}
	if (day < 10) {
		sday = "0" + std::to_string(day);
	}
	else {
		sday = std::to_string(day);
	}
	if (hour < 10) {
		shour = "0" + std::to_string(hour);
	}
	else {
		shour = std::to_string(hour);
	}



	std::string d = std::to_string(year) + "-" + smonth + "-" + sday + " " + shour + "-00-00";

	return d;
}
std::string System::generateMonth(int year, int month)
{
	std::string smonth;
	std::string sday;
	std::string shour;
	if (month < 10) {
		smonth = "0" + std::to_string(month);
	}
	else {
		smonth = std::to_string(month);
	}
	
	std::string d = std::to_string(year) + "-" + smonth;

	return d;
}

void System::generateP0(MatrixCPU* P0, std::string path, std::string month) {
	
	P0->setFromFile(path + month + ".txt", 1);
}
