
#include "../head/main.cuh"
#include "../head/main.h"

// pour l'agent des pertes de Q, lui permettre de consommer et de vendre peut poser probl�me 
// car il peut faire les 2 alors que l'on aimerait qu'il soit inactif, en vrai il fait juste intermediaire mais bon...


// tester sans r�actif et r�sistance nulle
// trouver o� j'ai cass� ? 



// DC - endogene
//R�solution directe centralis� et d�centralis�
//Consensus avec DC - OPF
//Comparaison avec powerTech

// Resoudre probleme:
// init avec march� d�j� resolu ou d�j� bien avanc�
// faire que y impose valeur march� mais pas inverse
// 
// 


// Cas Italy : le bus 136 est reli� par la suisse, il n'y a donc pas de ligne !!!! 


// mingw32 - make.exe

// Simulation

int main2() {
#ifdef DEBUG_TEST
	std::cout << "test Agent err =" << testAgent() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	std::cout << "test StudyCase err =" << testStudyCase() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	std::cout << "test Simparam err =" << testSimparam() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	std::cout << "test ADMMConst err =" << testADMM() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	std::cout << "test ADMMGPU err =" << testADMMGPUConst1() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	std::cout << "test ADMMGPU5 err =" << testADMMGPU5() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	std::cout << "test ADMMGPU6 err =" << testADMMGPU6() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	std::cout << "test ADMMGPUConst2 err =" << testADMMGPUConst2() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	std::cout << "test ADMMGPUConst3 err =" << testADMMGPUConst3() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	std::cout << "test systeme err =" << testSysteme() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	std::cout << "test systeme err =" << result << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	std::cout << "test Matrix err =" << testMatrix() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	std::cout << "test MatrixGPU err =" << testMGPU() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	std::cout << "test OSQP err =" << testOSQP() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	testADMMGPUtemp();
	std::cout << "-------------------------------------------------------- " << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
#endif // DEBUG_TEST
	srand(time(nullptr));
	std::cout.precision(6);
	//std::cout << "test StudyCase err =" << testStudyCase() << std::endl;
	//std::cout << "test Utilities err =" << testUtilities() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	//std::cout << "test Matrix err =" << testMatrix() << std::endl;
	//std::cout << "test MatrixGPU err =" << testMGPU() << std::endl;
	//std::cout << "test ADMMConst err =" << testADMM() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	try
	{
		SimuStatMarketEndoArticle();
		//compareCPUGPU();
		//std::cout << "test PAC err =" << testPAC() << std::endl;
	}
	catch (const std::exception& e)
	{
		std::cout << e.what() << std::endl;
	}
	
	std::cout << "-------------------------------------------------------- " << std::endl;
	//std::cout << "test PACConst err =" << testPACConst() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	//std::cout << "test ADMM1Const err =" << testADMM1() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	//std::cout << "test ADMMGPU err =" << testADMMGPUConst1() << std::endl;
	std::cout << "-------------------------------------------------------- " << std::endl;
	//std::cout << testCalculResX(1) << std::endl;
	//std::cout << "test systeme err =" << testSysteme() << std::endl;
	try
	{
		//testMarket();
		//testCPUPF();
		//testOPF();
		//testMarketEndo();

		//SimuTemporalWOConstraint("Italy");
		//SimuTemporalTestFeederEndo();
		//SimuTemporalTestFeederEndoAll();
		//SimuTemporal("Europe");
		//SimuStatMarketEndoACAgent();


		//SimuStatMarketEndo();
		//SimuStatMarketEndoAC();
		//SimuCompare();
		
		//SimuTemporalLlimit("Europe");
		//comparaisonArticle();

		//SimuStatOPFCompare();
		//SimuCompareISGT();
		//
		//SimuCompareParra()
		
		//SimuStatMarketEndoACAgent();
		//SimuSensiStudyCase();
		//SimuTemporalWOConstraint("Europe");
		//SimuStatMarketEndoGrid();
		//std::cout << testCalculVGS(2) / BILLION << " s" << std::endl;
		//SimuTemporalTestFeeder();
		//std::cout << testCalculPnShared(0, 512, 1) / BILLION << " s" << std::endl;
		//std::cout << testCalculChat(0, 64, 9) / BILLION << " s" << std::endl;
		
	}
	catch (const std::exception& e)
	{
		std::cout << e.what() << std::endl;
	}
	

	//SimuStatPFTransport();
	//testADMMGPUtemp();
	//SimuStatPFSGE();
	
	//testCPUPF();
	//testADMMACConst();
	//getErrorSensiPwerLine();
	//SimuStatPowerTech();

	/*StudyCase cas;
	cas.SetACFromFile("case39");
	cas.display();*/
	//res.display();

	//comparaisonArticle();
	//testExemple();
	//SimulationTempStepOpti();
	//SimuTemporalLlimit("Europe");
	//SimuTemporalConvergence("Europe");
	//SimuTemporal();
	//SimuTemporal("France");
	//std::cout << testCalculLAMBDABt1(1) / 1000000000 << " s" << std::endl;
	//std::cout << testCalculLAMBDABt1(1) / 1000000000 << " s" << std::endl;
	//std::cout << testCalculQpart(7) / 1000000000 << " s" << std::endl;
	//std::cout << testCalculQpart(2) / 1000000000 << " s" << std::endl;
	// testCalculChat
	//SimuTemporalRho();
	
	  
	
	return 0;
}




void SimulationTempStepOpti()
{
	// date 0->4 1 janvier 2012
	// date 0->4 1 juin 2013
	// date 10-15 20 octobre 2014

	std::string path = "data/";
	std::string fileName = "SimutemporalStep5a.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	
	std::string methode = "ADMMGPUConst5a";
	
	float epsG = 0.001f;
	float epsGC = 0.1f;
	float epsL = 0.00001f;



	MatrixCPU interval(4, 2);
	interval.set(0, 0, 2013);
	interval.set(0, 1, 2013);
	interval.set(1, 0, 6);
	interval.set(1, 1, 6);
	interval.set(2, 0, 10);
	interval.set(2, 1, 10);
	interval.set(3, 0, 0);
	interval.set(3, 1, 4);

	int iterG = 80000;
	int iterL = 1500;

	//int nCons = 25;
	//int nGen = 25;
	int nGen = 969;
	int nCons = 1494;
	int nAgent = nGen + nCons;



	float rho = 12; //rhoAgent * nAgent; //0.056 * nAgent;
	float rho1 = 0.0003; //0.0004;
	float rhoL = 12;

	
	const int nStep = 9;
	int Steps[nStep] = { 1, 2, 5, 10, 25, 50, 100, 200, 500 };
	MatrixCPU StepLocals(1, nStep);
	MatrixCPU StepGlobal(1, nStep);
	for (int i = 0; i < nStep; i++) {
		StepLocals.set(0, i, Steps[i]);
		StepGlobal.set(0, i, Steps[i]);
	}

	System sys;
	sys.setIter(iterG, iterL);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	sys.setEpsGC(epsGC);
	sys.setMethod(methode);
	sys.setRho(rho);
	sys.setRhoL(rhoL);
	sys.setRho1(rho1);


	MatrixCPU Param(1, 18);
	Param.set(0, 0, nAgent);
	Param.set(0, 1, epsGC);
	Param.set(0, 2, epsG);
	Param.set(0, 3, epsL);
	Param.set(0, 4, iterG);
	Param.set(0, 5, iterL);
	Param.set(0, 6, 0);
	Param.set(0, 7, rho);
	Param.set(0, 8, rhoL);
	Param.set(0, 9, rho1);
	Param.set(0, 10, interval.get(0, 0));
	Param.set(0, 11, interval.get(1, 0));
	Param.set(0, 12, interval.get(2, 0));
	Param.set(0, 13, interval.get(3, 0));
	Param.set(0, 14, interval.get(0, 1));
	Param.set(0, 15, interval.get(1, 1));
	Param.set(0, 16, interval.get(2, 1));
	Param.set(0, 17, interval.get(3, 1));
	Param.saveCSV(fileName, mode);
	
	StepLocals.saveCSV(fileName, mode);
	StepGlobal.saveCSV(fileName, mode);


	for (int local = 0; local < 1; local++) {
		for (int global = 0; global < nStep; global++) {
			sys.setStep(StepGlobal.get(0, global), StepLocals.get(0, local));
			
		
			clock_t t = clock();
			sys.solveIntervalle(path, &interval, nCons, nGen);
			t = clock() - t;

			std::cout << "temps simulation : " << t / CLOCKS_PER_SEC << std::endl;


			MatrixCPU temps(sys.getTemps());
			MatrixCPU iter(sys.getIter());
			MatrixCPU conv(sys.getConv());
			MatrixCPU fc(sys.getFc());
			//MatrixCPU ResR(sys.getResR());
			//MatrixCPU ResS(sys.getResS());

			temps.display();
			iter.display();

			temps.saveCSV(fileName, mode);
			iter.saveCSV(fileName, mode);
			fc.saveCSV(fileName, mode);
			//ResR.saveCSV(fileName, mode);
			//ResS.saveCSV(fileName, mode);
			conv.saveCSV(fileName, mode);

			std::cout << "-------------------------------------------------------- " << std::endl;
			sys.resetParam();
		}
	}



}

void testExemple()
{
	System sys;
	std::string method = "ADMMGPUConstCons5";
	//std::string method = "ADMMConst";
	sys.setMethod(method);
	std::string path = "data/";
	
	float epsG = 0.0001f;
	float epsL = 0.00005f;
	int iterG = 20000;
	int iterL = 3000;
	int stepG = 100;
	int stepL = 5;
	sys.setStep(stepG, stepL);
	sys.setIter(iterG, iterL);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	StudyCase cas;
	cas.Set4nodeBis(path);
	cas.setReduce(false);
	cas.display(0);
	float rho = 1.5;
	float rho1 = 0.01;
	float rhol = 2;
	sys.setStudyCase(cas);
	sys.setRho(rho);
	sys.setRho1(rho1);
	sys.setRhoL(rhol);
	Simparam res = sys.solve();
	std::cout << "---------------------------------------------------" << std::endl;
	res.display(0);
	
	MatrixCPU Sensi(cas.getPowerSensi());
	MatrixCPU Pn(res.getPn());
	MatrixCPU g(cas.getNLine(), 1, 0);
	std::cout << "----------------" << std::endl;
	std::cout << "Puissance des agents " << std::endl;
	Pn.display();
	std::cout << "----------------" << std::endl;
	//std::cout << "Echange entre les agents " << std::endl;
	//sys.displayTradesAgent();
	std::cout << "----------------" << std::endl;
	std::cout << "flux dans les lignes " << std::endl;
	if (cas.getNLine() != 0) {
		g.multiply(&Sensi, &Pn);
		cas.displayLineCores(&g);
	}
	///MatrixCPU residual(res.getRes());
}

void testADMMGPUtemp()
{
	System sys;
	float epsG = 0.00001f;
	float epsL = 0.000001f;
	int iterG = 5000;
	int iterL = 1500;
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	sys.setStep(10, 10);
	sys.setIter(iterG, iterL);
	sys.setRho1(0.8); //10
	const int nMethode = 5; // 12
	//std::string methodes[nMethode] = { "ADMMConst", "ADMMConst","OSQPConst", "ADMMGPUConst1", "ADMMGPUConst1T", "ADMMGPUConst2", "ADMMGPUConst3", "ADMMGPUConst4", "ADMMGPUConstCons", "PAC", "PACConst", "DCOPFOSQP"}; // ADMMGPUConstCons2 ADMMGPUConstCons3
	std::string methodes[nMethode] = { "ADMMGPUConst1", "ADMMGPUConst1T", "ADMMGPUConst2", "ADMMGPUConst3", "ADMMGPUConst4" };
	std::cout << "--------------------CAS SANS CONTRAINTE----------------------------- " << std::endl;
	float fc = -87016;
	float Pn[31] = { -0.927, -4.51, -2.81, -0.795, -0.880, -0.0980, -0.128, -4.07, -3.04, -2.16, -0.593, -3.31, -0.389, -2.71, -2.93, -1.77, -0.490, -2.26, -1.05, -0.138, -2.24, 4.246, 5.19, 3.42, 3.75, 4.41, 2.34, 3.31, 2.57, 3.49, 4.23 };

	MatrixCPU P(31, 1);
	for (int i = 0; i < 31; i++) {
		P.set(i, 0, Pn[i]);
	}

	MatrixCPU result(nMethode, 1, -1);
	MatrixCPU temps(nMethode, 1, -1);
	StudyCase cas1;
	cas1.Set29node();
	cas1.display();
	sys.setStudyCase(cas1);

	for (int i = 0; i < nMethode-1; i++) {
		std::string methode = methodes[i];
		std::cout << methode << std::endl;
		if (!methode.compare("ADMMGPUConst4")) {
			sys.setEpsL(epsL / 10);
		}
		else {
			sys.setEpsL(epsL);
		}
		sys.setMethod(methode);
		sys.setRho(10000);
		Simparam res = sys.solve();
		
		res.display();
		
		MatrixCPU P2 = res.getPn();
		P2.display();
			
		int test = P2.isEqual(&P, 0.1) + 2 * (fabs(fc - res.getFc())<= 1);
		result.set(i, 0, test);
		temps.set(i, 0, res.getTime());
		std::cout << "-------------------------------------------------------- " << std::endl;
		
	}

	
	std::cout << "--------------------CAS AVEC CONTRAINTE----------------------------- " << std::endl;
	MatrixCPU result2(nMethode, 1, -1);
	MatrixCPU temps2(nMethode, 1, -1);
	StudyCase cas;
	float lim = 0.8;
	cas.Set2nodeConstraint(lim);
	sys.setStudyCase(cas);
	
	float value = (1 - lim) * (lim > 1) + lim;
	float fc1 = -0.960524;
	MatrixCPU P1(2, 1);
	P1.set(0, 0, -value);
	P1.set(1, 0, value);
	for (int i = 0 ; i < nMethode-1; i++) {
		std::string methode = methodes[i];
		if (!methode.compare("ADMMGPUConst4")) {
			sys.setEpsL(epsL / 10);
		}
		else {
			sys.setEpsL(epsL);
		}
		std::cout << methode << std::endl;
		sys.setMethod(methode);
		sys.setRho(1.5);
		Simparam res = sys.solve();
		res.display();
		MatrixCPU P2 = res.getPn();
		P2.display();
		int test = P2.isEqual(&P1, 0.01) + 2 * (fabs(fc1 - res.getFc()) <= 0.1);
		result2.set(i, 0, test);
		temps2.set(i, 0, res.getTime());
		std::cout << "-------------------------------------------------------- " << std::endl;
	}
	std::cout << "--------------------CAS SANS GAMMA POUR COMPARER AVEC OPF----------------------------- " << std::endl;

	MatrixCPU result3(nMethode, 1, -1);
	MatrixCPU temps3(nMethode, 1, -1);
	cas.Set39Bus();
	cas.setLineLimit(26, 200);
	cas.setReduce(true);
	cas.genBetaUniforme(0);
	sys.setStudyCase(cas);
	float rho = 1000;
	float epsGC = 0.01f;
	epsG = 0.01f;
	epsL = 0.0001f;
	//float rho1 = 0.00028;
	float rho1 = 0.5;
	float rhoL = 1.5;
	sys.setRho1(rho1);
	sys.setRho(rho);
	sys.setEpsGC(epsGC);
	float fcRef = 0;
	float errM = 0;
	MatrixCPU PnRef;
	MatrixCPU err;
	for (int i = 0; i < nMethode-1; i++) {
		std::string methode = methodes[i];
		std::cout << methode << std::endl;
		if (!methode.compare("ADMMGPUConst4")) {
			sys.setEpsL(epsL / 10);
		}
		else {
			sys.setEpsL(epsL);
		}
		sys.setMethod(methode);
		//sys.display(2);
		Simparam res = sys.solve();
		res.display();
		MatrixCPU P2 = res.getPn();

		if (i == 0) { // ce sera la ref
			fcRef = res.getFc();
			PnRef = res.getPn();
			err = MatrixCPU(PnRef.getNLin(), 1, 1);
			std::cout << "fcRef " << fcRef << std::endl;
			result3.set(i, 0, 3);
			temps3.set(i, 0, res.getTime());
		}
		else {
			fc = res.getFc();
			MatrixCPU P2 = res.getPn();
			err.RelativeEror(&PnRef, &P2);
			errM = err.max2();
			std::cout << "fc " << fc << " err " << errM << std::endl;
			//err.display();
			int test = (errM < 0.01) + 2 * (fabs((fc - fcRef) / fcRef) <= 0.01);
			result3.set(i, 0, test);
			temps3.set(i, 0, res.getTime());
		}


		std::cout << "-------------------------------------------------------- " << std::endl;
		sys.resetParam();
	}


	std::cout << "--------------------GRAND CAS SANS GAMMA POUR COMPARER AVEC OPF----------------------------- " << std::endl;
	std::string name = "France";
	std::string path = "data/" + name + "/";
	const int nMethodeBig = 3;
	std::string methodesBig[nMethodeBig] = { "ADMMConst", "ADMMGPUConst4", "DCOPFOSQP"};
	MatrixCPU result4(nMethodeBig, 1, -1);
	MatrixCPU temps4(nMethodeBig, 1, -1);


	MatrixCPU interval(4, 2);
	interval.set(0, 0, 2013);
	interval.set(0, 1, 2013);
	interval.set(1, 0, 6);
	interval.set(1, 1, 6);
	interval.set(2, 0, 1);
	interval.set(2, 1, 1);
	interval.set(3, 0, 0);
	interval.set(3, 1, 0);

	rho = 125;
	rhoL = rho;
	rho1 = 0.0000001;
	epsG = 0.001f;
	epsGC = 1.0f;
	epsL = 0.0001f;
	sys.setEpsGC(epsGC);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	for (int i = 0; i < 0; i++) {
		std::string methode = methodesBig[i];
		std::cout << methode << std::endl;
		sys.setMethod(methode);
		sys.setRho(rho);
		sys.setRho1(rho1);
		sys.setRhoL(rhoL);
		std::cout << "----------------Simu---------------------------- " << std::endl;
		// Debut simu

		clock_t t = clock();

		sys.solveIntervalle(path, name, &interval);
		t = clock() - t;

	
		MatrixCPU fcMat(sys.getFc());
		if (i == 0) { // ce sera la ref
			fcRef = fcMat.get(0, 0);
			PnRef = sys.getPn();
			err = MatrixCPU(PnRef.getNLin(), 1, 1);
			std::cout << "fcRef " << fcRef << std::endl;
			result4.set(i, 0, 3);
			temps4.set(i, 0, (float)t / CLOCKS_PER_SEC);
		}
		else {
			fc = fcMat.get(0, 0);
			MatrixCPU P2 = sys.getPn();
			err.RelativeEror(&PnRef, &P2);
			errM = err.max2();
			std::cout << "fc " << fc << " err " << errM <<std::endl;
			//err.display();
			int test = ( errM < 0.01) + 2 * (fabs((fc - fcRef) / fcRef) <= 0.01);
			result4.set(i, 0, test);
			temps4.set(i, 0, (float) t / CLOCKS_PER_SEC);
		}


		std::cout << "-------------------------------------------------------- " << std::endl;
		sys.resetParam();
	}

	std::cout << "-------------------- RESULTATS ----------------------------- " << std::endl;

	result.display();
	temps.display();
	result2.display();
	temps2.display();
	result3.display();
	temps3.display();
	result4.display();
	temps4.display();

}

std::string generateDate(int year, int month, int day, int hour)
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
		
		

	std::string d = std::to_string(year) + "-" + smonth + "-" + sday + " " + shour +"-00-00";
	
	return d;
}

float pow10(int n)
{
	float v = 1;
	for (int i = 0; i < n; i++) {
		v = v * 10;
	}
	return v;
}

int getNFileline(std::string nameFile)
{
	int number_of_lines = 0;
	std::string line;
	std::ifstream myfile(nameFile);

	while (std::getline(myfile, line))
		++number_of_lines;
	return number_of_lines;
}

void getErrorSensiPwerLine()
{
	int _nGen = 969;
	int _nCons = 1494;
	int _nAgent = _nGen + _nCons;
	StudyCase cas;
	std::string path = "data/";
	std::string fileName3 = path + "SensiPowerReduceEurope.txt";
	std::string fileName4 = path + "lineLimitReduceEurope.txt";
	int _nLineConstraint = getNFileline(fileName4);
	MatrixCPU _SensiPowerReduce = MatrixCPU(_nLineConstraint, _nAgent);
	MatrixCPU _lineLimitsReduce = MatrixCPU(_nLineConstraint, 1);
	_SensiPowerReduce.setFromFile(fileName3);
	_lineLimitsReduce.setFromFile(fileName4);

	MatrixCPU errorSensi(_nLineConstraint, 1);

	for (int l = 0; l < _nLineConstraint; l++) {
		float sum = 0;
		for (int n = 0; n < _nAgent; n++) {
			sum += abs(_SensiPowerReduce.get(l, n));
		}
		errorSensi.set(l, 0, sum);
	}
	std::cout << "resolution of the PowerLine " << errorSensi.max2() << std::endl;
	std::cout << "Puissance max ligne " << _lineLimitsReduce.max2() << std::endl;

}

/* Simu marche*/

void comparaisonArticle()
{
	System sys;
	std::string fileName = "CompareArticle2.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;

	std::string method = "ADMMConst";
	//std::string method = "ADMMGPUConst1T";
	sys.setMethod(method);

	float epsG = 0.00001f;
	float epsGC = 0.001f;
	float epsL = 0.000005f;
	int   iterG = 5000;
	int   iterL = 10000;
	float stepG = 1;
	float stepL = 1;
	float rho = 10;
	//float rho1 = 0.00028;
	float rho1 = 5;
	float rhoL = rho;
	float distCost = 10;
	/*MatrixCPU MatDistCost = MatrixCPU(4, 4, distCost);
	for (int i = 0; i < 4; i++) {
		MatDistCost.set(i, 3, distCost / 2); // moins de co�t pour les �changes avec la zone 4
		MatDistCost.set(3, i, distCost / 2); // moins de co�t pour les �changes avec la zone 4
	}*/
	//MatDistCost.set(3, 3, distCost / 2); // moins de co�t pour les �changes dans la zone 4 uniquement

	sys.setIter(iterG, iterL);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	sys.setEpsGC(epsGC);
	sys.setStep(stepG, stepL);
	StudyCase cas;
	cas.Set39Bus("data/", false);




	// cas.genBetaUniforme(10);

	//cas.setDistance(1, "data/Distance39.txt");
	//cas.genBetaDistance(distCost);

	//cas.genBetaDistanceByZone(&MatDistCost);
	MatrixCPU beta(cas.getBeta());
	//beta.display();
	//std::cout << "set Line limit " << std::endl;
	cas.setLineLimit(26, 200);
	cas.setReduce(true);
	MatrixCPU Sensi(cas.getPowerSensi());
	Sensi.display();
	//cas.display(1);
	sys.setStudyCase(cas);
	sys.setRho(rho);
	sys.setRho1(rho1);
	sys.setRhoL(rhoL);


	MatrixCPU Param(1, 13);
	Param.set(0, 0, cas.getNagent());
	Param.set(0, 1, epsGC);
	Param.set(0, 2, epsG);
	Param.set(0, 3, epsL);
	Param.set(0, 4, iterG);
	Param.set(0, 5, iterL);
	Param.set(0, 6, stepG);
	Param.set(0, 7, stepL);
	Param.set(0, 8, rho);
	Param.set(0, 9, rhoL);
	Param.set(0, 10, rho1);
	Param.set(0, 11, distCost);
	Param.set(0, 12, 2);
	std::cout << " Save Y:1, N:0" << std::endl;
	int inputUser = 0;
	std::cin >> inputUser;
	if (inputUser) {
		Param.saveCSV(fileName, mode);
	}

	MatrixCPU Result(1, 7);



	std::cout << "Solve :" << std::endl;
	Simparam res = sys.solve();
	res.display();
	int iter = res.getIter();
	float fc = res.getFc();

	MatrixCPU Res(res.getRes());
	//Res.saveCSV("ResRhoRhoRho.csv", std::fstream::in | std::fstream::out | std::fstream::app);
	float resR = Res.get(0, (iter - 1) / stepG);
	float resS = Res.get(1, (iter - 1) / stepG);
	float resX = Res.get(2, (iter - 1) / stepG);
	float temps = res.getTime();
	std::cout << "----------------" << std::endl;
	Result.set(0, 0, iter);
	Result.set(0, 1, temps);
	Result.set(0, 2, fc);
	Result.set(0, 3, resR);
	Result.set(0, 4, resS);
	Result.set(0, 5, resX);


	std::cout << "----------------" << std::endl;

	std::cout << "Resultat de simulation fc " << fc << std::endl;
	std::cout << "temps de simulation " << temps << std::endl;
	std::cout << "Nombre d'it�ration " << iter << std::endl;
	std::cout << " Residus R,S,X : " << resR << ", " << resS << ", " << resX << std::endl;
	std::cout << "----------------" << std::endl;
	std::cout << "Echange entre les agents " << std::endl;
	sys.displayTradesAgent();
	std::cout << "----------------" << std::endl;
	std::cout << "flux dans les lignes " << std::endl;
	if (cas.getNLine() != 0) {
		MatrixCPU Sensi(cas.getPowerSensi());

		MatrixCPU Pn(res.getPn());
		if (inputUser) {
			Pn.saveCSV(fileName, mode, 1);
		}
		MatrixCPU g(cas.getNLine(), 1, 0);
		g.multiply(&Sensi, &Pn);
		cas.displayLineCores(&g);
		Result.set(0, 6, g.get(0, 0));
	}/**/
	if (inputUser) {
		Result.saveCSV(fileName, mode);
	}

}

void SimuStatADMMGPU()
{
	//std::string fileName = "ADMMGPU5_article_Release.csv";
	//std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	std::string methode = "ADMMGPUConst3";
	float rhoAgent = 0.05;
	float P = 1000;
	float dP = 400;
	float a = 0.07;
	float da = 0.02;
	float b = 50;
	float db = 20;
	int nNAgent = 1;
	int nAgentMax = 3000;
	int offset = 0;
	int nSimu = 50;
	int iterGlobal = 1500;
	int iterLocal = 700;
	int stepG = 5;
	int stepL = 5;
	float epsG = 0.001f;
	float epsL = 0.0001f;

	MatrixCPU Param(1, 10);
	Param.set(0, 0, nAgentMax);
	Param.set(0, 1, nNAgent);
	Param.set(0, 2, rhoAgent);
	Param.set(0, 3, epsG);
	Param.set(0, 4, epsL);
	Param.set(0, 5, iterGlobal);
	Param.set(0, 6, iterLocal);
	Param.set(0, 7, stepG);
	Param.set(0, 8, stepL);
	Param.set(0, 9, nSimu);
	//Param.saveCSV(fileName, mode);

	MatrixCPU Agents(1, nNAgent);
	MatrixCPU temps(nSimu, nNAgent, -1);
	MatrixCPU iters(nSimu, nNAgent, -1);
	MatrixCPU fcs(nSimu, nNAgent, -1);
	MatrixCPU ResF(2, iterGlobal / stepG);
	MatrixCPU ResR(nSimu, nNAgent, -1);
	MatrixCPU ResS(nSimu, nNAgent, -1);
	System sys;
	sys.setIter(iterGlobal, iterLocal);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	sys.setStep(stepG, stepL);
	sys.setMethod(methode);

	for (int agent = 0; agent < nNAgent; agent++) {
		std::cout << "--------- --------- --------- --------- ----------" << std::endl;

		int agents = (agent + 1) * (nAgentMax - offset) / nNAgent + offset;
		std::cout << agents << std::endl;
		std::cout << "--------- --------- --------- --------- ----------" << std::endl;
		Agents.set(0, agent, agents);

		for (int j = 0; j < nSimu; j++) {
			StudyCase cas(agents, P, dP, a, da, b, db);
			sys.setStudyCase(cas);
			std::cout << "-";
			float rho = rhoAgent * agents;
			sys.setRho(rho);
			clock_t t = clock();
			Simparam res = sys.solve();
			float temp = clock() - t;
			int iter = res.getIter();
			temps.set(j, agent, temp / CLOCKS_PER_SEC);
			iters.set(j, agent, iter);
			fcs.set(j, agent, res.getFc());
			ResF = res.getRes();
			ResR.set(j, agent, ResF.get(0, (iter - 1) / stepG));
			ResS.set(j, agent, ResF.get(1, (iter - 1) / stepG));

		}
		std::cout << std::endl;
	}
	Agents.display();
	temps.display();
	iters.display();
	float temptotal = temps.sum();
	std::cout << "temps total " << temptotal << " temps moyen " << temptotal / nSimu << std::endl;
	/*Agents.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	iters.saveCSV(fileName, mode);
	fcs.saveCSV(fileName, mode);
	ResR.saveCSV(fileName, mode);
	ResS.saveCSV(fileName, mode);*/
}

void SimuTemporalLlimit(std::string name) {
	int type = 0;
	std::string path;
	int nGen;
	int nCons;
	int nAgent;
	//int nLine;
	if (!name.compare("Europe")) {
		type = 1;
		path = "data/";
		nGen = 969;
		nCons = 1494;
		nAgent = nGen + nCons;
		//nLine = 2156;
	}
	else {
		path = "data/" + name + "/";
	}
	std::string fileName = "SimutemporalLimitGPU" + name + ".csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nMethode = 1;
	std::string methodes[nMethode] = { "ADMMGPUConst4" };
	//std::string method = "ADMMGPUConst1";
	//std::string method = "ADMMConst";

	//std::string methodes[nMethode] = { "ADMMGPUConst1","ADMMGPUConst2","ADMMGPUConst3", "ADMMGPUConst3" };
	float epsG = 0.001f;
	float epsGC = 0.1f;
	float epsL = 0.0001f;
	int warmstart = 0;

	MatrixCPU interval(4, 2);
	interval.set(0, 0, 2013);
	interval.set(0, 1, 2013);
	interval.set(1, 0, 6);
	interval.set(1, 1, 6);
	interval.set(2, 0, 3);
	interval.set(2, 1, 3);
	interval.set(3, 0, 0);
	interval.set(3, 1, 23);

	int iterG = 20000;
	int iterL = 5000;//500;


	const int nLimit = 5; // 
	MatrixCPU LilimitTab(1, nLimit);
	
	for (int i = 0; i < nLimit; i++) {
		LilimitTab.set(0, i, 300 * i);
	}

	float rho = 12; //rhoAgent * nAgent; //0.056 * nAgent;
	float rho1 = 0.00037; //0.0004;
	float rhoL = 12;

	float stepG = 10;
	float stepL = 1;


	System sys;
	StudyCase cas;
	sys.setIter(iterG, iterL);
	sys.setStep(stepG, stepL);
	sys.setEpsG(epsG);
	sys.setEpsGC(epsGC);
	sys.setEpsL(epsL);
	sys.setWarmStart(warmstart == 1);


	MatrixCPU Param(1, 21);
	Param.set(0, 0, nAgent);
	Param.set(0, 1, epsGC);
	Param.set(0, 2, epsG);
	Param.set(0, 3, epsL);
	Param.set(0, 4, iterG);
	Param.set(0, 5, iterL);
	Param.set(0, 6, stepG);
	Param.set(0, 7, stepL);
	Param.set(0, 8, rho);
	Param.set(0, 9, rho1);
	Param.set(0, 10, rhoL);
	Param.set(0, 11, interval.get(0, 0));
	Param.set(0, 12, interval.get(1, 0));
	Param.set(0, 13, interval.get(2, 0));
	Param.set(0, 14, interval.get(3, 0));
	Param.set(0, 15, interval.get(0, 1));
	Param.set(0, 16, interval.get(1, 1));
	Param.set(0, 17, interval.get(2, 1));
	Param.set(0, 18, interval.get(3, 1));
	Param.set(0, 19, nLimit);
	Param.set(0, 20, warmstart);
	
	
	std::cout << " Save Y:1, N:0" << std::endl;
	int inputUser;
	std::cin >> inputUser;
	if (inputUser) {
		
	}


	for (int i = 0; i < nMethode; i++) {
		std::string methode = methodes[i];
		sys.setMethod(methode);
		sys.setRho(rho);
		sys.setRho1(rho1);
		sys.setRhoL(rhoL);
		for (int j = 0; j < nLimit; j++) {
			float limit = LilimitTab.get(0, j);
			
			sys.setLineLimitMin(limit);

			std::cout << "----------------Simu---------------------------- " << std::endl;
			// Debut simu
			clock_t t = clock();
			if (type) {
				sys.solveIntervalle(path, &interval, nCons, nGen);
			}
			else {
				sys.solveIntervalle(path, name, &interval);
			}

			t = clock() - t;

			std::cout << "temps simulation : " << t / CLOCKS_PER_SEC << std::endl;

			MatrixCPU iter(sys.getIter());
			MatrixCPU ResR(sys.getResR());
			MatrixCPU ResS(sys.getResS());
			MatrixCPU ResX(sys.getResX());
			MatrixCPU tempTab(sys.getTemps());

			if (inputUser) {
				if(j==0 && i == 0){
					Param.set(0, 0, sys.getNagent());
					Param.saveCSV(fileName, mode);
					LilimitTab.saveCSV(fileName, mode);
				}
				
				iter.saveCSV(fileName, mode);
				ResR.saveCSV(fileName, mode);
				ResS.saveCSV(fileName, mode);
				ResX.saveCSV(fileName, mode);
				tempTab.saveCSV(fileName, mode);
				iter.display();
				ResR.display();
				ResS.display();
				ResX.display();
				tempTab.display();
			}
			else {
				iter.display();
				ResR.display();
				ResS.display();
				ResX.display();
				tempTab.display();
			}
			std::cout << "-------------------------------------------------------- " << std::endl;
			sys.resetParam();
			sys.resetMethod();
		}
		
	}
}

void SimuTemporalConvergence(std::string name) {
	int type = 0;
	std::string path;
	int nGen;
	int nCons;
	//int nAgent;
	//int nLine;
	if (!name.compare("Europe")) {
		type = 1;
		path = "data/";
		nGen = 969;
		nCons = 1494;
		//nAgent = nGen + nCons;
		//nLine = 2156;
	}
	else {
		path = "data/" + name + "/";
	}
	std::string fileName = "SGESimutemporalConvergenceRho" + name + ".csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nMethode = 1;
	std::string methodes[nMethode] = { "ADMMGPUConst4" };
	//std::string method = "ADMMGPUConst1";
	//std::string method = "ADMMConst";

	//std::string methodes[nMethode] = { "ADMMGPUConst1","ADMMGPUConst2","ADMMGPUConst3", "ADMMGPUConst3" };
	float epsG = 0.001f;
	float epsGC = 1.0f;
	float epsL = 0.00001f;
	int warmstart = 1;
	float offset = 0;

	MatrixCPU interval(4, 2);
	interval.set(0, 0, 2013);
	interval.set(0, 1, 2013);
	interval.set(1, 0, 6);
	interval.set(1, 1, 6);
	interval.set(2, 0, 3);
	interval.set(2, 1, 3);
	interval.set(3, 0, 10);
	interval.set(3, 1, 19);

	int iterG = 100000;
	int iterL = 10000;//500;


	float rho = 10; //rhoAgent * nAgent; //0.056 * nAgent;
	float rho1 = 0.0005; //0.00037;
	float rhoL = 10;

	float stepG = 10;
	float stepL = 1;


	System sys;
	StudyCase cas;
	sys.setIter(iterG, iterL);
	sys.setStep(stepG, stepL);
	sys.setEpsG(epsG);
	sys.setEpsGC(epsGC);
	sys.setEpsL(epsL);
	sys.setWarmStart(warmstart == 1);
	sys.setConstraintRelaxation(offset);

	MatrixCPU Param(1, 21);
	Param.set(0, 0, 0);
	Param.set(0, 1, epsGC);
	Param.set(0, 2, epsG);
	Param.set(0, 3, epsL);
	Param.set(0, 4, iterG);
	Param.set(0, 5, iterL);
	Param.set(0, 6, stepG);
	Param.set(0, 7, stepL);
	Param.set(0, 8, rho);
	Param.set(0, 9, rho1);
	Param.set(0, 10, rhoL);
	Param.set(0, 11, interval.get(0, 0));
	Param.set(0, 12, interval.get(1, 0));
	Param.set(0, 13, interval.get(2, 0));
	Param.set(0, 14, interval.get(3, 0));
	Param.set(0, 15, interval.get(0, 1));
	Param.set(0, 16, interval.get(1, 1));
	Param.set(0, 17, interval.get(2, 1));
	Param.set(0, 18, interval.get(3, 1));
	Param.set(0, 19, offset);
	Param.set(0, 20, warmstart);


	std::cout << " Save Y:1, N:0" << std::endl;
	int inputUser;
	std::cin >> inputUser;
	if (inputUser) {
		Param.saveCSV(fileName, mode);
	}


	for (int i = 0; i < nMethode; i++) {
		std::string methode = methodes[i];
		sys.setMethod(methode);
		sys.setRho(rho);
		sys.setRho1(rho1);
		sys.setRhoL(rhoL);
		
			

		std::cout << "----------------Simu---------------------------- " << std::endl;
		// Debut simu
		clock_t t = clock();
		if (type) {
			sys.solveIntervalle(path, &interval, nCons, nGen);
		}
		else {
			sys.solveIntervalle(path, name, &interval);
		}

		t = clock() - t;

		std::cout << "temps simulation : " << t / CLOCKS_PER_SEC << std::endl;

		MatrixCPU iter(sys.getIter());
		MatrixCPU ResR(sys.getResR());
		MatrixCPU ResS(sys.getResS());
		MatrixCPU ResX(sys.getResX());
		MatrixCPU tempTab(sys.getTemps());
		MatrixCPU fcTab(sys.getFc());

		if (inputUser) {
			tempTab.saveCSV(fileName, mode);
			iter.saveCSV(fileName, mode);
			ResR.saveCSV(fileName, mode);
			ResS.saveCSV(fileName, mode);
			ResX.saveCSV(fileName, mode);
			fcTab.saveCSV(fileName, mode);
			iter.display();
			ResR.display();
			ResS.display();
			ResX.display();
			tempTab.display();
		}
		else {
			iter.display();
			ResR.display();
			ResS.display();
			ResX.display();
			tempTab.display();
		}
		std::cout << "-------------------------------------------------------- " << std::endl;
		sys.resetParam();
		sys.resetMethod();
		

	}
}

void SimuTemporal(std::string name)
{
	int type = 0;
	std::string path;
	int nGen;
	int nCons;
	int nAgent = 0;
	//int nLine;
	if (!name.compare("Europe")) {
		type = 1;
		nGen = 969;
		nCons = 1494;
		nAgent = nGen + nCons;
		//nLine = 2156;
	}
	
	path = "data/" + name + "/";
	

	std::string fileName = "SimutemporalDCEndoGPURho" + name + "Stock.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	//const int nMethode = 4;
	//std::string methodes[nMethode] = { "ADMMConst", "OSQPConst","ADMMGPUConst1T","ADMMGPUConst4" };
	const int nMethode = 1;
	std::string methodes[nMethode] = {"ADMMGPUConst4" };
	//std::string methodes[nMethode] = { "ADMMGPUConst1","ADMMGPUConst2","ADMMGPUConst3", "ADMMGPUConst4" };
	//std::string methodes[nMethode] = { "ADMMConst", "ADMMGPUConst1T","ADMMGPUConst3", "ADMMGPUConst4"};
	float epsG = 0.01f;
	float epsGC = 1.0f;
	float epsL = 0.00001f;
	float offset = 1;


	MatrixCPU interval(4, 2);
	interval.set(0, 0, 2012);
	interval.set(0, 1, 2012);
	interval.set(1, 0, 11);
	interval.set(1, 1, 11);
	interval.set(2, 0, 1);
	interval.set(2, 1, 1);
	interval.set(3, 0, 0);
	interval.set(3, 1, 23);

	int iterG = 30000;
	int iterL = 2000;



	
	float rho = 125;
	float rho1 = 0.001;
	//float rho1 = 0.0000001;
	//float rho1 = 125;
	//float rhoL = 125;

	float stepG = 5;
	float stepL = 1;

	System sys;

	sys.setIter(iterG, iterL);
	sys.setStep(stepG, stepL);
	sys.setEpsGC(epsGC);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	if (offset)
	{
		sys.setConstraintRelaxation(offset);
	}


	MatrixCPU Param(1, 18);
	Param.set(0, 0, nAgent);
	Param.set(0, 1, epsGC);
	Param.set(0, 2, epsG);
	Param.set(0, 3, epsL);
	Param.set(0, 4, iterG);
	Param.set(0, 5, iterL);
	Param.set(0, 6, stepG);
	Param.set(0, 7, stepL);
	Param.set(0, 8, rho);
	Param.set(0, 9, rho1);
	Param.set(0, 10, interval.get(0, 0));
	Param.set(0, 11, interval.get(1, 0));
	Param.set(0, 12, interval.get(2, 0));
	Param.set(0, 13, interval.get(3, 0));
	Param.set(0, 14, interval.get(0, 1));
	Param.set(0, 15, interval.get(1, 1));
	Param.set(0, 16, interval.get(2, 1));
	Param.set(0, 17, interval.get(3, 1));
	Param.saveCSV(fileName, mode);


	for (int i = 0; i < nMethode; i++) {
		std::string methode = methodes[i];
		sys.setMethod(methode);
		sys.setRho(rho);
		sys.setRho1(rho1);
		std::cout << "----------------Simu---------------------------- " << std::endl;
		// Debut simu

		clock_t t = clock();

		if (type) {
			sys.solveIntervalle(path, &interval, nCons, nGen);
		}
		else {
			sys.solveIntervalle(path, name, &interval);
		}
		t = clock() - t;

		std::cout << "calculation time : " << (float)t / CLOCKS_PER_SEC << std::endl;

		MatrixCPU temps(sys.getTemps());
		MatrixCPU iter(sys.getIter());
		MatrixCPU conv(sys.getConv());
		MatrixCPU fc(sys.getFc());
		MatrixCPU ResR(sys.getResR());
		MatrixCPU ResS(sys.getResS());
		MatrixCPU ResX(sys.getResX());

		temps.display();
		ResR.display();
		ResS.display();
		ResX.display();
		iter.display();
		temps.saveCSV(fileName, mode);
		iter.saveCSV(fileName, mode);
		fc.saveCSV(fileName, mode);
		ResR.saveCSV(fileName, mode);
		ResS.saveCSV(fileName, mode);
		ResX.saveCSV(fileName, mode);
		conv.saveCSV(fileName, mode);

		std::cout << "-------------------------------------------------------- " << std::endl;
		/*sys.displayTime(fileName);*/
		sys.resetParam();
	}

}


void SimuTemporalWOConstraint(std::string name)
{
	int type = 0;
	std::string path;
	int nGen;
	int nCons;
	int nAgent = 0;
	//int nLine;
	if (!name.compare("Europe")) {
		type = 1;
		nGen = 969;
		nCons = 1494;
		nAgent = nGen + nCons;
		//nLine = 2156;
	}
	path = "data/" + name + "/";
	
	std::string fileName;
#ifdef INSTRUMENTATION
	fileName = "SimutemporalFB" + name + ".csv";
#else
	fileName = "SimuTemporalMarketrho" + name + "Stock.csv";
#endif // INSTRUMENTATION

	
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	//const int nMethode = 8;
	//std::string methodesName[nMethode] = { "ADMM", "ADMMMP","ADMMGPU", "PAC","PACMP","PACGPU", "OSQP","OSQPCen"};
	const int nMethode = 1;
	std::string methodesName[nMethode] = {"ADMMMP"};

	Method* methodes[nMethode];
	//methodes[0] = new ADMMMarket;
	//methodes[0] = new ADMMMarketOpenMP;
	//methodes[0] = new ADMMMarketGPU;
	//methodes[0] = new PAC;
	//methodes[0] = new PACOpenMP;
	//methodes[0] = new PACGPU;
	//methodes[0] = new OSQP;
	//methodes[0] = new OSQPCentralized2;




	float epsG = 0.01f;
	float epsL = 0.005f;

	MatrixCPU interval(4, 2);
	interval.set(0, 0, 2013);
	interval.set(0, 1, 2013);
	interval.set(1, 0, 3);
	interval.set(1, 1, 3);
	interval.set(2, 0, 1);
	interval.set(2, 1, 1);
	interval.set(3, 0, 0);
	interval.set(3, 1, 3);

	int iterG = 10000;
	int iterL = 1000;


	float stepG = 2;
	float stepL = 2;
	float rho = 10;
	float rhoL = 10;


	MatrixCPU Param(1, 18);
	Param.set(0, 0, nAgent);
	Param.set(0, 1, 0);
	Param.set(0, 2, epsG);
	Param.set(0, 3, epsL);
	Param.set(0, 4, iterG);
	Param.set(0, 5, iterL);
	Param.set(0, 6, stepG);
	Param.set(0, 7, stepL);
	Param.set(0, 8, rho);
	Param.set(0, 9, rhoL);
	Param.set(0, 10, interval.get(0, 0));
	Param.set(0, 11, interval.get(1, 0));
	Param.set(0, 12, interval.get(2, 0));
	Param.set(0, 13, interval.get(3, 0));
	Param.set(0, 14, interval.get(0, 1));
	Param.set(0, 15, interval.get(1, 1));
	Param.set(0, 16, interval.get(2, 1));
	Param.set(0, 17, interval.get(3, 1));
	

	System sys;

	sys.setIter(iterG, iterL);
	sys.setStep(stepG, stepL);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);



	for (int i = 0; i < nMethode; i++) {
		std::cout << "----------------Simu---------------------------- " << std::endl;
		sys.setMethod(methodes[i]);

		sys.setRho(rho);
		//sys.setRhoL(rhoL);
		// Debut simu

		clock_t t = clock();

		if (type) {
			sys.solveIntervalle(path, &interval, nCons, nGen);
		}
		else {
			sys.solveIntervalle(path, name, &interval);
		}
		t = clock() - t;

		std::cout << "calculation time : " << (float)t / CLOCKS_PER_SEC << std::endl;

		MatrixCPU temps(sys.getTemps());
		MatrixCPU iter(sys.getIter());
		MatrixCPU conv(sys.getConv());
		MatrixCPU fc(sys.getFc());
		MatrixCPU ResR(sys.getResR());
		MatrixCPU ResS(sys.getResS());

		std::cout << "-------------------------------------------------------- " << std::endl;
		if (i == 0) {
			int nAgent = sys.getNagent();
			Param.set(0, 0, nAgent);
			//Param.saveCSV(fileName, mode);
		}
#ifdef INSTRUMENTATION
		sys.displayTime(fileName);
#else
		temps.saveCSV(fileName, mode);
		iter.saveCSV(fileName, mode);
		conv.saveCSV(fileName, mode);
		fc.saveCSV(fileName, mode);
		ResR.saveCSV(fileName, mode);
		ResS.saveCSV(fileName, mode);
		iter.display();
#endif // INSTRUMENTATION

		
		sys.resetParam();
	}
	
	sys.setMethod(nullptr); /**/
	for (int i = 0; i < nMethode; i++) {
		DELETEB(methodes[i]);
	}
}


void SimuTemporalRho(std::string name)
{
	int type = 0;
	std::string path;
	int nGen;
	int nCons;
	//int nAgent;
	//int nLine;
	if (!name.compare("Europe")) {
		type = 1;
		path = "data/";
		nGen = 969;
		nCons = 1494;
		//nAgent = nGen + nCons;
		//nLine = 2156;
	}
	else {
		path = "data/" + name + "/";
	}
	std::string fileName = "SimutemporalRho" + name + ".csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nMethode = 1;
	std::string methodes[nMethode] = { "ADMMGPUConst1T" };
	//std::string method = "ADMMGPUConst1";
	//std::string method = "ADMMConst";

	//std::string methodes[nMethode] = { "ADMMGPUConst1","ADMMGPUConst2","ADMMGPUConst3", "ADMMGPUConst3" };
	float epsG = 0.01f;
	float epsGC = 0.1f;
	float epsL = 0.005f;

	MatrixCPU interval(4, 2);
	interval.set(0, 0, 2013);
	interval.set(0, 1, 2013);
	interval.set(1, 0, 6);
	interval.set(1, 1, 6);
	interval.set(2, 0, 1);
	interval.set(2, 1, 1);
	interval.set(3, 0, 0);
	interval.set(3, 1, 0);

	int iterG = 20000;
	int iterL = 2000;//500;


	const int nRho = 7; // rhoAgent = 1 semble une bonne valeur
	MatrixCPU rhoTab(2, nRho, 7);
	rhoTab.set(1, 0, 0.0002);

	for (int i = 1; i < nRho; i++) {
		rhoTab.set(0, i, rhoTab.get(0, i - 1) + 1);
		rhoTab.set(1, i, rhoTab.get(1, i - 1) + 0.0001);
	}

	float stepG = 10;
	float stepL = 5;

	System sys;
	StudyCase cas;
	sys.setIter(iterG, iterL);
	sys.setStep(stepG, stepL);
	sys.setEpsG(epsG);
	sys.setEpsGC(epsGC);
	sys.setEpsL(epsL);


	MatrixCPU Param(1, 16);
	Param.set(0, 0, 1);
	Param.set(0, 1, epsG);
	Param.set(0, 2, epsGC);
	Param.set(0, 3, epsL);
	Param.set(0, 4, iterG);
	Param.set(0, 5, iterL);
	Param.set(0, 6, stepG);
	Param.set(0, 7, stepL);
	Param.set(0, 8, interval.get(0, 0));
	Param.set(0, 9, interval.get(1, 0));
	Param.set(0, 10, interval.get(2, 0));
	Param.set(0, 11, interval.get(3, 0));
	Param.set(0, 12, interval.get(0, 1));
	Param.set(0, 13, interval.get(1, 1));
	Param.set(0, 14, interval.get(2, 1));
	Param.set(0, 15, interval.get(3, 1));
	std::cout << " Save Y:1, N:0" << std::endl;
	int inputUser;
	std::cin >> inputUser;
	if (inputUser) {
		Param.saveCSV(fileName, mode);
		rhoTab.saveCSV(fileName, mode);
	}


	for (int i = 0; i < nMethode; i++) {
		std::string methode = methodes[i];
		sys.setMethod(methode);
		for (int k = 0; k < nRho; k++) {
			float rho = rhoTab.get(0, k);
			sys.setRho(rho);
			for (int j = 0; j < nRho; j++) {
				float rho1 = rhoTab.get(1, j);
				float rhoL = 1 * MYMAX(rho, rho1);
				sys.setRho1(rho1);
				sys.setRhoL(rhoL);

				std::cout << "----------------Simu---------------------------- " << std::endl;
				// Debut simu
				clock_t t = clock();
				if (type) {
					sys.solveIntervalle(path, &interval, nCons, nGen);
				}
				else {
					sys.solveIntervalle(path, name, &interval);
				}

				t = clock() - t;

				std::cout << "temps simulation : " << t / CLOCKS_PER_SEC << std::endl;

				MatrixCPU iter(sys.getIter());
				MatrixCPU ResR(sys.getResR());
				MatrixCPU ResS(sys.getResS());
				MatrixCPU ResX(sys.getResX());

				if (inputUser) {
					iter.saveCSV(fileName, mode);
					ResR.saveCSV(fileName, mode);
					ResS.saveCSV(fileName, mode);
					ResX.saveCSV(fileName, mode);
				}


				std::cout << "-------------------------------------------------------- " << std::endl;
				sys.resetParam();
				sys.resetMethod();
			}
		}
	}

}

void SimuTemporal()
{
	std::string path = "data/";
	//std::string fileName2 = "SimutemporalFBTimeOffset.csv";
	std::string fileName = "SimutemporalCompareMethodeJDD.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nMethode = 4;
	std::string methodes[nMethode] = { "ADMMConst", "OSQPConst","ADMMGPUConst1T","ADMMGPUConst4" };
	//const int nMethode = 1;
	//std::string methodes[nMethode] = { "ADMMGPUConst1T"};
	//std::string methodes[nMethode] = { "ADMMGPUConst1","ADMMGPUConst1T","ADMMGPUConst2", "ADMMGPUConst3", "ADMMGPUConst4", "ADMMGPUConst5", "ADMMGPUConst5a" };
	//std::string methodes[nMethode] = { "ADMMConst", "ADMMGPUConst1T" };
	float epsG = 0.01f;
	float epsGC = 5.0f;
	float epsL = 0.001f;
	float lineLimMin = 0; // MW
	int offset = 2;
	int warmstart = 1;

	MatrixCPU interval(4, 2);
	interval.set(0, 0, 2013);
	interval.set(0, 1, 2013);
	interval.set(1, 0, 6);
	interval.set(1, 1, 6);
	interval.set(2, 0, 1);
	interval.set(2, 1, 1);
	interval.set(3, 0, 0);
	interval.set(3, 1, 23);

	int iterG = 200000;
	int iterL = 1000;//500;


	int nGen = 969;
	int nCons = 1494;
	int nAgent = nGen + nCons;



	float rho = 12; //rhoAgent * nAgent; //0.056 * nAgent;
	float rho1 = 0.00037; //0.0004;
	float rhoL = 12;

	float stepG = 10;
	float stepL = 1;

	System sys;
	StudyCase cas;
	sys.setIter(iterG, iterL);
	sys.setStep(stepG, stepL);
	sys.setEpsGC(epsGC);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	sys.setLineLimitMin(lineLimMin);
	sys.setWarmStart(warmstart == 1);
	if (offset)
	{
		sys.setConstraintRelaxation(offset);
	}


	MatrixCPU Param(1, 22);
	Param.set(0, 0, nAgent);
	Param.set(0, 1, epsGC);
	Param.set(0, 2, epsG);
	Param.set(0, 3, epsL);
	Param.set(0, 4, iterG);
	Param.set(0, 5, iterL);
	Param.set(0, 6, stepG);
	Param.set(0, 7, stepL);
	Param.set(0, 8, rho);
	Param.set(0, 9, rho1);
	Param.set(0, 10, rhoL);
	Param.set(0, 11, interval.get(0, 0));
	Param.set(0, 12, interval.get(1, 0));
	Param.set(0, 13, interval.get(2, 0));
	Param.set(0, 14, interval.get(3, 0));
	Param.set(0, 15, interval.get(0, 1));
	Param.set(0, 16, interval.get(1, 1));
	Param.set(0, 17, interval.get(2, 1));
	Param.set(0, 18, interval.get(3, 1));
	Param.set(0, 19, lineLimMin);
	Param.set(0, 20, warmstart);
	Param.set(0, 21, offset);
	Param.saveCSV(fileName, mode);
	//Param.saveCSV(fileName2, mode);


	for (int i = 0; i < nMethode; i++) {
		std::string methode = methodes[i];
		sys.setMethod(methode);
		if (!methode.compare("ADMMGPUConst4")) {
			sys.setEpsL(epsL / 10);
		}

		sys.setRho(rho);
		sys.setRho1(rho1);
		sys.setRhoL(rhoL);
		std::cout << "----------------Simu---------------------------- " << std::endl;
		// Debut simu

		clock_t t = clock();

		sys.solveIntervalle(path, &interval, nCons, nGen);
		t = clock() - t;

		std::cout << "calculation time : " << (float)t / CLOCKS_PER_SEC << std::endl;


		MatrixCPU temps(sys.getTemps());
		MatrixCPU iter(sys.getIter());
		MatrixCPU conv(sys.getConv());
		MatrixCPU fc(sys.getFc());
		MatrixCPU ResR(sys.getResR());
		MatrixCPU ResS(sys.getResS());
		MatrixCPU ResX(sys.getResX());

		temps.display();
		ResR.display();
		ResS.display();
		ResX.display();
		iter.display();
		temps.saveCSV(fileName, mode);
		iter.saveCSV(fileName, mode);
		fc.saveCSV(fileName, mode);
		ResR.saveCSV(fileName, mode);
		ResS.saveCSV(fileName, mode);
		ResX.saveCSV(fileName, mode);
		conv.saveCSV(fileName, mode);/**/

		std::cout << "-------------------------------------------------------- " << std::endl;
		//sys.displayTime(fileName2);
		sys.resetParam();
	}

}

void SimuCompare()
{
	std::string fileName = "ComparaisonMarket2_All_400.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	
	const int nMethode = 8;
	std::vector<int> indices = { 0, 1, 2, 3, 4, 5, 6, 7 };
	std::string methodesName[nMethode] = { "OSQPCentalized","PAC", "PACOpenMP", "PACGPU", "OSQP", "ADMMMarket","ADMMMarketOpenMP", "ADMMMarketGPU" };
	Method* methodes[nMethode];
	//methodes[0] = new OSQPCentralized2;
	methodes[1] = new PAC;
	methodes[2] = new PACOpenMP;
	methodes[3] = new PACGPU;
	methodes[4] = new PACGPU;
	methodes[5] = new ADMMMarket;
	methodes[6] = new ADMMMarketOpenMP;
	methodes[7] = new ADMMMarketGPU;


	if (nMethode != indices.size()) {
		throw std::domain_error("not good number of methods");
	}
	// for random case
	float P = 100;
	float dP = 10;
	float a = 1;
	float da = 0.1;
	float b = 2;
	float db = 1;
	float gamma = 2;
	float dGamma = 1;
	float propCons = 0.3;
	float propPro = 0;	
	float propGenNFle = 0.1f;

	float Q = 10;
	float dQ = 2;
	/*float P = 100;
	float dP = 20;
	float Q = 10;
	float dQ = 2;
	float a = 0.07; // pour P et Q
	float da = 0.02;
	float b = 10;
	float db = 4;
	float gamma = 8;
	float dGamma = 2;
	float propCons = 0.375f;
	float propPro = 0;
	float propGenNFle = 0.125f;*/

	bool AC = false;
	

	int nNAgent = 1;
	int nAgentMax = 400;
	int offset = 0;
	float rho = 1;
	int nSimu = 50;
	int iterGlobal = 100;
	int iterLocal = 1000;
	int stepG = 10;
	int stepL = 10;
	float epsG = 0.1f;
	float epsL = 0.01f;

	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;
	

	MatrixCPU Param(1, 25);
	Param.set(0, 0, nAgentMax);
	Param.set(0, 1, nNAgent);
	Param.set(0, 2, rho);
	Param.set(0, 3, epsG);
	Param.set(0, 4, epsL);
	Param.set(0, 5, iterGlobal);
	Param.set(0, 6, iterLocal);
	Param.set(0, 7, stepG);
	Param.set(0, 8, stepL);
	Param.set(0, 9, nSimu);
	Param.set(0, 10, nMethode);
	Param.set(0, 11, P);
	Param.set(0, 12, dP);
	Param.set(0, 13, Q);
	Param.set(0, 14, dQ);
	Param.set(0, 15, a);
	Param.set(0, 16, da);
	Param.set(0, 17, b);
	Param.set(0, 18, db);
	Param.set(0, 19, gamma);
	Param.set(0, 20, dGamma);
	Param.set(0, 21, propCons);
	Param.set(0, 22, propPro);
	Param.set(0, 23, propGenNFle);
	Param.set(0, 24, AC);

	Param.saveCSV(fileName, mode);

	MatrixCPU Agents(1, nNAgent);
	MatrixCPU temps(nMethode * nSimu, nNAgent, nanf(""));
	MatrixCPU iters(nMethode * nSimu, nNAgent, nanf(""));
	MatrixCPU fcs(nMethode * nSimu, nNAgent, nanf(""));
	MatrixCPU ResF(2, iterGlobal / stepG);
	MatrixCPU ResR(nMethode * nSimu, nNAgent, nanf(""));
	MatrixCPU ResS(nMethode * nSimu, nNAgent, nanf(""));
	System sys;
	StudyCase cas;
	sys.setIter(iterGlobal, iterLocal);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	sys.setStep(stepG, stepL);


	for (int agent = 0; agent < nNAgent; agent++) {
		std::cout << "--------- --------- --------- --------- ----------" << std::endl;

		int agents = (agent + 1) * (nAgentMax - offset) / nNAgent + offset;
		std::cout << agents << std::endl;

		Agents.set(0, agent, agents);

		for (int j = 0; j < nSimu; j++) {
			if (AC) {
				cas = StudyCase(agents, P, dP, P, dP, a, da, a, da, b, db, gamma, dGamma, propCons, propGenNFle, propPro);
			}
			else
			{
				cas = StudyCase(agents, P, dP, a, da, b, db, gamma, dGamma, propCons, propGenNFle, propPro);
			}
			sys.setStudyCase(cas);
			std::cout << "-";
			//float rho = rhoAgent * agents;
			sys.setRho(rho);
			std::random_shuffle(indices.begin(), indices.end());
			for (int i = 0; i < nMethode; i++) {
				if (indices[i] > 0) {

					if (indices[i] > 0 && indices[i] < 4) {
						methodes[indices[i]]->setBestParam(cas);
					}

					sys.setMethod(methodes[indices[i]]);


					//clock_t t = clock();
					t1 = std::chrono::high_resolution_clock::now();
					Simparam res = sys.solve();
					t2 = std::chrono::high_resolution_clock::now();
					//clock_t temp = clock() - t;
					int iter = res.getIter();
					temps.set(indices[i] * nSimu + j, agent, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / BILLION);
					iters.set(indices[i] * nSimu + j, agent, iter);
					fcs.set(indices[i] * nSimu + j, agent, res.getFc());
					ResF = res.getRes();
					ResR.set(indices[i] * nSimu + j, agent, ResF.get(0, (iter - 1) / stepG));
					ResS.set(indices[i] * nSimu + j, agent, ResF.get(1, (iter - 1) / stepG));
				}
				//sys.resetParam();
			}
		}
		std::cout << std::endl;
	}
	Agents.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	iters.saveCSV(fileName, mode);
	fcs.saveCSV(fileName, mode);
	ResR.saveCSV(fileName, mode);
	ResS.saveCSV(fileName, mode);

	Agents.display();
	temps.display();
	iters.display();
	fcs.display();
	ResR.display();
	ResS.display();/**/

	for (int i = 0; i < nMethode; i++) {
		DELETEB(methodes[i]);
	}
	sys.setMethod(nullptr);
}

void SimuCompareAll()
{
	std::string fileName = "ComparaisonMarket_Release_50_200.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nMethode = 7;
	std::vector<int> indices = { 0, 1, 2, 3, 4, 5, 6 };
	std::string methodesName[nMethode] = { "OSQP", "PAC", "PACOpenMP", "PACGPU", "ADMMMarket","ADMMMarketOpenMP", "ADMMMarketGPU" };
	Method* methodes[nMethode];
	//methodes[0] = new OSQP;
	methodes[1] = new PAC;
	methodes[2] = new PACOpenMP;
	methodes[3] = new PACGPU;
	methodes[4] = new ADMMMarket;
	methodes[5] = new ADMMMarketOpenMP;
	methodes[6] = new ADMMMarketGPU;


	if (nMethode != indices.size()) {
		throw std::domain_error("not good number of methods");
	}

	float P = 100;
	float dP = 20;
	float Q = 10;
	float dQ = 2;
	float a = 0.07; // pour P et Q
	float da = 0.02;
	float b = 10;
	float db = 4;
	float gamma = 8;
	float dGamma = 2;
	float propCons = 0.5f;
	float propPro = 0;
	float propGenNFle = 0.125f;

	int nNAgent = 4;
	int nAgentMax = 200;
	int offset = 0;
	float rhoAgent = 0.25;
	int nSimu = 50;
	int iterGlobal = 4000;
	int iterLocal = 1000;
	int stepG = 10;
	int stepL = 10;
	float epsG = 0.01f;
	float epsL = 0.001f;

	MatrixCPU Param(1, 24);
	Param.set(0, 0, nAgentMax);
	Param.set(0, 1, nNAgent);
	Param.set(0, 2, rhoAgent);
	Param.set(0, 3, epsG);
	Param.set(0, 4, epsL);
	Param.set(0, 5, iterGlobal);
	Param.set(0, 6, iterLocal);
	Param.set(0, 7, stepG);
	Param.set(0, 8, stepL);
	Param.set(0, 9, nSimu);
	Param.set(0, 10, nMethode);
	Param.set(0, 11, P);
	Param.set(0, 12, dP);
	Param.set(0, 13, Q);
	Param.set(0, 14, dQ);
	Param.set(0, 15, a);
	Param.set(0, 16, da);
	Param.set(0, 17, b);
	Param.set(0, 18, db);
	Param.set(0, 19, gamma);
	Param.set(0, 20, dGamma);
	Param.set(0, 21, propCons);
	Param.set(0, 22, propPro);
	Param.set(0, 23, propGenNFle);

	Param.saveCSV(fileName, mode);

	MatrixCPU Agents(1, nNAgent);
	MatrixCPU temps(nMethode * nSimu, nNAgent, -1);
	MatrixCPU iters(nMethode * nSimu, nNAgent, -1);
	MatrixCPU fcs(nMethode * nSimu, nNAgent, -1);
	MatrixCPU ResF(2, iterGlobal / stepG);
	MatrixCPU ResR(nMethode * nSimu, nNAgent, -1);
	MatrixCPU ResS(nMethode * nSimu, nNAgent, -1);
	System sys;
	sys.setIter(iterGlobal, iterLocal);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	sys.setStep(stepG, stepL);


	for (int agent = 0; agent < nNAgent; agent++) {
		std::cout << "--------- --------- --------- --------- ----------" << std::endl;

		int agents = (agent + 1) * (nAgentMax - offset) / nNAgent + offset;
		std::cout << agents << std::endl;

		Agents.set(0, agent, agents);

		for (int j = 0; j < nSimu; j++) {
			StudyCase cas(agents, P, dP, Q, dQ, a, da, a, da, b, db, gamma, dGamma, propCons, propGenNFle, propPro);
			sys.setStudyCase(cas);
			std::cout << "-";
			float rho = rhoAgent * agents;
			sys.setRho(rho);
			std::random_shuffle(indices.begin(), indices.end());
			for (int i = 0; i < nMethode; i++) {
				if (indices[i] > 0 && indices[i] < 4) {
					methodes[indices[i]]->setBestParam(cas);
				}
				sys.setMethod(methodes[indices[i]]);
				clock_t t = clock();
				Simparam res = sys.solve();
				float temp = clock() - t;
				int iter = res.getIter();
				temps.set(indices[i] * nSimu + j, agent, temp / CLOCKS_PER_SEC);
				iters.set(indices[i] * nSimu + j, agent, iter);
				fcs.set(indices[i] * nSimu + j, agent, res.getFc());
				ResF = res.getRes();
				ResR.set(indices[i] * nSimu + j, agent, ResF.get(0, (iter - 1) / stepG));
				ResS.set(indices[i] * nSimu + j, agent, ResF.get(1, (iter - 1) / stepG));
				sys.resetParam();
			}
		}
		std::cout << std::endl;
	}
	Agents.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	iters.saveCSV(fileName, mode);
	fcs.saveCSV(fileName, mode);
	ResR.saveCSV(fileName, mode);
	ResS.saveCSV(fileName, mode);

	for (int i = 0; i < nMethode; i++) {
		DELETEB(methodes[i]);
	}
	sys.setMethod(nullptr);
}

void SimuCompareParra()
{
	std::string fileName = "ComparaisonMarketACPara_300b.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nMethode = 4;
	std::vector<int> indices = { 0, 1, 2, 3};
	std::string methodesName[nMethode] = { "PACOpenMP", "PACGPU", "ADMMMarketOpenMP", "ADMMMarketGPU" };
	Method* methodes[nMethode];
	 
	methodes[0] = new PACOpenMP;
	methodes[1] = new PACGPU;
	methodes[2] = new ADMMMarketOpenMP;
	methodes[3] = new ADMMMarketGPU;


	if (nMethode != indices.size()) {
		throw std::domain_error("not good number of methods");
	}

	float P = 100;
	float dP = 20;
	float Q = 10;
	float dQ = 2;
	float a = 0.07; // pour P et Q
	float da = 0.02;
	float b = 10;
	float db = 4;
	float gamma = 8;
	float dGamma = 2;
	float propCons = 0.375f;
	float propPro = 0;
	float propGenNFle = 0.125f;
	bool AC = true;


	int nNAgent = 1;
	int nAgentMax = 300;
	int offset = 0;
	float rhoAgent = 0.025;
	int nSimu = 50;
	int iterGlobal = 50;
	int iterLocal = 1000;
	int stepG = 10;
	int stepL = 10;
	float epsG = 0.01f;
	float epsL = 0.001f;

	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;
	

	MatrixCPU Param(1, 25);
	Param.set(0, 0, nAgentMax);
	Param.set(0, 1, nNAgent);
	Param.set(0, 2, rhoAgent);
	Param.set(0, 3, epsG);
	Param.set(0, 4, epsL);
	Param.set(0, 5, iterGlobal);
	Param.set(0, 6, iterLocal);
	Param.set(0, 7, stepG);
	Param.set(0, 8, stepL);
	Param.set(0, 9, nSimu);
	Param.set(0, 10, nMethode);
	Param.set(0, 11, P);
	Param.set(0, 12, dP);
	Param.set(0, 13, Q);
	Param.set(0, 14, dQ);
	Param.set(0, 15, a);
	Param.set(0, 16, da);
	Param.set(0, 17, b);
	Param.set(0, 18, db);
	Param.set(0, 19, gamma);
	Param.set(0, 20, dGamma);
	Param.set(0, 21, propCons);
	Param.set(0, 22, propPro);
	Param.set(0, 23, propGenNFle);
	Param.set(0, 24, AC);

	Param.saveCSV(fileName, mode);

	MatrixCPU Agents(1, nNAgent);
	MatrixCPU temps(nMethode * nSimu, nNAgent, nanf(""));
	MatrixCPU iters(nMethode * nSimu, nNAgent, nanf(""));
	MatrixCPU fcs(nMethode * nSimu, nNAgent, nanf(""));
	MatrixCPU ResF(2, iterGlobal / stepG);
	MatrixCPU ResR(nMethode * nSimu, nNAgent, nanf(""));
	MatrixCPU ResS(nMethode * nSimu, nNAgent, nanf(""));
	System sys;
	StudyCase cas;
	sys.setIter(iterGlobal, iterLocal);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	sys.setStep(stepG, stepL);



	for (int agent = 0; agent < nNAgent; agent++) {
		std::cout << "--------- --------- --------- --------- ----------" << std::endl;

		int agents = (agent + 1) * (nAgentMax - offset) / nNAgent + offset;
		std::cout << agents << std::endl;

		Agents.set(0, agent, agents);

		for (int j = 0; j < nSimu; j++) {
			if (AC) {
				cas = StudyCase(agents, P, dP, P, dP, a, da, a, da, b, db, gamma, dGamma, propCons, propGenNFle, propPro);
			}
			else
			{
				cas = StudyCase(agents, P, dP, a, da, b, db, propCons, propPro);
			}
			sys.setStudyCase(cas);
			std::cout << "-";
			float rho = rhoAgent * agents;
			sys.setRho(rho);
			std::random_shuffle(indices.begin(), indices.end());
			for (int i = 0; i < nMethode; i++) {
				//if (indices[i] > 2) {
					if (indices[i] < 2) {
						methodes[indices[i]]->setBestParam(cas);
					}
					sys.setMethod(methodes[indices[i]]);
					t1 = std::chrono::high_resolution_clock::now();
					Simparam res = sys.solve();
					t2 = std::chrono::high_resolution_clock::now();
					//clock_t temp = clock() - t;
					int iter = res.getIter();

					temps.set(indices[i] * nSimu + j, agent, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / BILLION);
					iters.set(indices[i] * nSimu + j, agent, iter);
					fcs.set(indices[i] * nSimu + j, agent, res.getFc());
					ResF = res.getRes();
					ResR.set(indices[i] * nSimu + j, agent, ResF.get(0, (iter - 1) / stepG));
					ResS.set(indices[i] * nSimu + j, agent, ResF.get(1, (iter - 1) / stepG));
					//
				//}
			}
		}
		std::cout << std::endl;
	}
	Agents.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	iters.saveCSV(fileName, mode);
	fcs.saveCSV(fileName, mode);
	ResR.saveCSV(fileName, mode);
	ResS.saveCSV(fileName, mode);

	for (int i = 0; i < nMethode; i++) {
		DELETEB(methodes[i]);
	}
	sys.setMethod(nullptr);
}

void SimuStatPowerTech()
{
	std::string fileName = "PowerTechSizeEvolutionReduceOffsetCPU2.csv";
	std::string fileName2 = "PowerTechSizeEvolutionReduceOffsetCPUFB2.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	std::string methode = "ADMMConst";
	std::string path = "data/";
	float rhoAgent = 10;
	float rhoLine = 0.0005;

	float Pconso = 60;
	float dPconso = 50;
	float Propcons = 0.60;
	float bProd = 20;
	float dbProd = 18;
	float Pprod = 300;
	float dPprod = 250;
	float gamma = 4;
	float dgamma = 2;
	float limit = 1000;
	float dlimit = 300;

	int nNAgent = 1;
	int nAgentMax = 600;
	int nNLine = 10;
	int nLineMax = 900;
	int offsetAgent = 0;
	int offsetLine = 0;
	int nSimu = 20;

	int offset = 2;
	int iterGlobal = 50000;
	int iterLocal = 5000;
	int stepG = 10;
	int stepL = 1;
	float epsG = 0.01f;
	float epsL = 0.001f;
	float epsGC = 1.0f;
	MatrixCPU Param(1, 18);
	Param.set(0, 0, nAgentMax);
	Param.set(0, 1, nNAgent);
	Param.set(0, 2, rhoAgent);
	Param.set(0, 3, nLineMax);
	Param.set(0, 4, nNLine);
	Param.set(0, 5, rhoLine);
	Param.set(0, 6, epsG);
	Param.set(0, 7, epsGC);
	Param.set(0, 8, epsL);
	Param.set(0, 9, iterGlobal);
	Param.set(0, 10, iterLocal);
	Param.set(0, 11, stepG);
	Param.set(0, 12, stepL);
	Param.set(0, 13, nSimu);
	Param.set(0, 14, offset);
	Param.set(0, 15, Pconso);
	Param.set(0, 16, Pprod);
	Param.set(0, 17, limit);
	Param.saveCSV(fileName, mode);
	//Param.saveCSV(fileName2, mode);

	MatrixCPU Agents(1, nNAgent);
	MatrixCPU Lines(1, nNLine);

	MatrixCPU temps(nNLine * nSimu, nNAgent, -1);
	MatrixCPU iters(nNLine * nSimu, nNAgent, -1);
	MatrixCPU fcs(nNLine * nSimu, nNAgent, -1);
	MatrixCPU ResF(3, iterGlobal / stepG);
	MatrixCPU ResR(nNLine * nSimu, nNAgent, -1);
	MatrixCPU ResS(nNLine * nSimu, nNAgent, -1);
	MatrixCPU ResX(nNLine * nSimu, nNAgent, -1);

	System sys;
	sys.setIter(iterGlobal, iterLocal);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	sys.setStep(stepG, stepL);
	sys.setMethod(methode);


	/*
	*	void genGridFromFile(std::string path, bool alreadyDefine=true);
		void genAgents(int nAgent, float propCons, float Pconso, float dPconso,  float bProd, float dbProd, float Pprod, float dPpord, float gamma, float dGamma); // gamma = -1 distance ?
		void genLinkGridAgent();
		void computeSensiPower();
		void genLineLimit(int nLine, float limit, float dlimit );
	*/
	for (int line = 0; line < nNLine; line++) {
		std::cout << "--------- --------- --------- --------- ----------" << std::endl;

		int lines = line * (nLineMax - offsetLine) / MYMAX((nNLine - 1), 1) + offsetLine;
		std::cout << lines << std::endl;
		//std::cout << "--------- --------- --------- --------- ----------" << std::endl;
		Lines.set(0, line, lines);
		for (int agent = 0; agent < nNAgent; agent++) {
			//std::cout << "--------- --------- --------- --------- ----------" << std::endl;

			int agents = (agent + 1) * (nAgentMax - offsetAgent) / nNAgent + offsetAgent;
			std::cout << agents << std::endl;
			std::cout << "-|-|-|-|-|-|-|-|-  -|-|-|-|-|-|-|-|-|" << std::endl;
			Agents.set(0, agent, agents);

			for (int j = 0; j < nSimu; j++) {
				StudyCase cas;
				std::cout << "-";
				cas.genGridFromFile(path);
				cas.genAgents(agents, Propcons, Pconso, dPconso, bProd, dbProd, Pprod, dPprod, gamma, dgamma); // gamma = -1 distance ?
				cas.genLinkGridAgent();
				cas.genLineLimit(lines, limit, dlimit);
				cas.setReduce(true);
				cas.setLineLimitRelaxation(epsGC);
				sys.setStudyCase(cas);
				//cas.display();

				//float rho = rhoAgent * agents;
				//float rho1 = rhoLine * lines;
				sys.setRho(rhoAgent);
				sys.setRho1(rhoLine);
				std::cout << "|";
				clock_t t = clock();
				Simparam res = sys.solve();
				float temp = clock() - t;

				int iter = res.getIter();
				temps.set(line * nSimu + j, agent, temp / CLOCKS_PER_SEC);
				iters.set(line * nSimu + j, agent, iter);
				fcs.set(line * nSimu + j, agent, res.getFc());
				ResF = res.getRes();
				ResR.set(line * nSimu + j, agent, ResF.get(0, (iter - 1) / stepG));
				ResS.set(line * nSimu + j, agent, ResF.get(1, (iter - 1) / stepG));
				ResX.set(line * nSimu + j, agent, ResF.get(2, (iter - 1) / stepG));
				sys.displayTime(fileName2);
				sys.resetMethod();
			}
			std::cout << std::endl;
		}
	}
	Agents.display();
	temps.display();
	iters.display();
	float temptotal = temps.sum();
	std::cout << "temps total " << temptotal << " temps moyen " << temptotal / (nSimu * nNLine * nNAgent) << std::endl;

	Lines.saveCSV(fileName, mode);
	Agents.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	iters.saveCSV(fileName, mode);
	fcs.saveCSV(fileName, mode);
	ResR.saveCSV(fileName, mode);
	ResS.saveCSV(fileName, mode);
	ResX.saveCSV(fileName, mode);/**/
}


void SimuSensiStudyCase()
{
	std::string fileName = "ComparaisonMarkeSensiAgent_100_4.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;

	const int nMethode = 1;
	std::string methodesName[nMethode] = { "ADMMMarketGPU" };
	Method* methodes = new ADMMMarketGPU;
	
	int nAgent = 300;

	float P0min = 0;
	float P0max = 100;
	float aMin = 0.1;
	float aMax = 10;
	float gammaMin = 0;
	float gammaMax = 100;
	float propConsoMin = 0.1;
	float propConsoMax = 0.9;
	float borneMin = 10;
	float borneMax = 10;


	float rho = 100;
	float rhoL = rho;
	int nSimu = 1000;

	int iterGlobal = 10000;
	int iterLocal = 10000;
	int stepG = 10;
	int stepL = 10;
	float epsG = 0.0001f;
	float epsL = 0.00001f;

	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;
	

	MatrixCPU Param(1, 20);
	Param.set(0, 0, nAgent);
	Param.set(0, 1, rho);
	Param.set(0, 2, rhoL);
	Param.set(0, 3, epsG);
	Param.set(0, 4, epsL);
	Param.set(0, 5, iterGlobal);
	Param.set(0, 6, iterLocal);
	Param.set(0, 7, stepG);
	Param.set(0, 8, stepL);
	Param.set(0, 9, nSimu);
	Param.set(0, 10, P0min);
	Param.set(0, 11, P0max);
	Param.set(0, 12, aMin);
	Param.set(0, 13, aMax);
	Param.set(0, 14, propConsoMin);
	Param.set(0, 15, propConsoMax);
	Param.set(0, 16, gammaMin);
	Param.set(0, 17, gammaMax);
	Param.set(0, 18, borneMin);
	Param.set(0, 19, borneMax);


	Param.saveCSV(fileName, mode);

	MatrixCPU Rho(1, 1, rho);
	MatrixCPU temps(nMethode * nSimu, 1, nanf(""));
	MatrixCPU iters(nMethode * nSimu, 1, nanf(""));
	
	System sys;
	StudyCase cas;
	sys.setIter(iterGlobal, iterLocal);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	sys.setStep(stepG, stepL);
	sys.setMethod(methodes);

	
	for (int j = 0; j < nSimu; j++) {
		if (j % 50 == 0) {
			std::cout << std::endl;
		}
		cas.genAgentsFullRandom(nAgent, aMin, aMax, P0min, P0max, gammaMin, gammaMax, propConsoMin, propConsoMax, borneMin, borneMax);
		sys.setStudyCase(cas);
		std::cout << "-";
		
		sys.setRho(rho);
		
		t1 = std::chrono::high_resolution_clock::now();
		Simparam res = sys.solve();
		t2 = std::chrono::high_resolution_clock::now();
		
		int iter = res.getIter();
		temps.set(j, 0, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / BILLION);
		iters.set(j, 0, iter);
	
		
		//sys.resetParam();
		
	}
	
	Rho.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	iters.saveCSV(fileName, mode);
	

	
	DELETEA(methodes);
	
	sys.setMethod(nullptr);
}


/* Simu PF */

void SimuStatPFSGE() {
	std::string fileName = "SGEPosterPowerFlowStat500.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;

	float Pconso = 0.005;
	float dPconso = 0.001;
	float Propcons = 0.50;
	float PropGen = 0.5;
	float bProd = 5;
	float dbProd = 2;
	float Pprod = 0.005;
	float dPprod = 0.001;
	float gamma = 4;
	float dgamma = 2;
	float dQ = 0.005;
	float length = 0.001;
	float dlength = 0.0005;

	int nNBus = 1;
	int nBusMax = 500;
	const int offsetBus = 0;
	int nSimu = 50;

	float factorAgent[] = { 0.2, 0.5, 1, 2 };

	float factorBuses[] = { 2, 0.3, 0.3 , 1 };
	int nCasAgent = 4;
	int nCasBuses = 4;
	int nCas = nCasAgent * nCasBuses;
	int nMethod = 5;


	MatrixCPU Param(1, 13);
	Param.set(0, 0, nBusMax);
	Param.set(0, 1, nNBus);
	Param.set(0, 2, Pconso);
	Param.set(0, 3, dPconso);
	Param.set(0, 4, Pprod);
	Param.set(0, 5, dPprod);
	Param.set(0, 6, dQ);
	Param.set(0, 7, length);
	Param.set(0, 8, dlength);
	Param.set(0, 9, nSimu);
	Param.set(0, 10, nCasAgent);
	Param.set(0, 11, nCasBuses);
	Param.set(0, 12, nMethod);

	Param.saveCSV(fileName, mode);

	MatrixCPU Buses(1, nNBus);

	MatrixCPU temps(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU temps2(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU iters(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU PResult(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU QResult(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU Residuals(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU ConvResult(nCas * nMethod, nNBus * nSimu, -3);
	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;

	for (int bus = 0; bus < nNBus; bus++) {
		std::cout << "--------- --------- --------- --------- ----------" << std::endl;
		int buses = offsetBus > 0 ? bus * (nBusMax - offsetBus) / MYMAX(nNBus - 1, 1) + offsetBus : (bus + 1) * nBusMax / nNBus;
		Buses.set(0, bus, buses);
		for (int casBus = 0; casBus < nCasBuses; casBus++) {
			int nDeep = 0;
			int nBranch = 0;
			switch (casBus)
			{
			case 0:
				nDeep = sqrt(factorBuses[0] * buses);
				nBranch = sqrt(factorBuses[0] * buses);
				break;
			case 1:
				nDeep = sqrt(factorBuses[1] * buses);
				nBranch = buses;
				break;
			case 2:
				nDeep = buses;
				nBranch = sqrt(factorBuses[2] * buses);
				break;
			case 3:
				nDeep = factorBuses[3] * buses;
				nBranch = factorBuses[3] * buses;
				break;
			default:
				std::cout << "error casBus =" << casBus << " not defined" << std::endl;
				return;
			}
			for (int casAgent = 0; casAgent < nCasAgent; casAgent++) {
				int agents = factorAgent[casAgent] * buses;
				for (int simu = 0; simu < nSimu; simu++) {
					CPUPF PF;
					GPUPF PFG;
					CPUPFGS PFGS;
					GPUPFGS PFGSG;
					CPUPFdist PFDist;

					int i = (casBus * nCasAgent + casAgent) * nMethod;
					int j = bus * nSimu + simu;

					StudyCase cas;
					std::cout << "-";
					cas.genGridBT(buses, nBranch, nDeep, length, dlength);
					cas.genAgentsAC(agents, Propcons, PropGen, Pconso, dPconso, bProd, dbProd, dQ, Pprod, dPprod, gamma, dgamma);
					cas.genLinkGridAgent();


					MatrixCPU PQ = cas.getPobj();
					MatrixGPU PQG = MatrixGPU(cas.getPobj(), 1);

					//std::cout << "*****************************************************************************" << std::endl;
					t1 = std::chrono::high_resolution_clock::now();
					PFDist.init(cas, &PQ);
					if (!PFDist.chekcase()) {
						std::cout << "pas reseau de distribution" << std::endl;
						return;
					}



					PFDist.solve();
					PFDist.calcW(true);
					t2 = std::chrono::high_resolution_clock::now();

					PResult.set(i, j, PFDist.getP0());
					QResult.set(i, j, PFDist.getQ0());
					temps.set(i, j, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
					temps2.set(i, j, PFDist.getTime());
					iters.set(i, j, PFDist.getIter());
					Residuals.set(i, j, PFDist.getRes());
					ConvResult.set(i, j, PFDist.getConv());
					i++;
					int conv = PFDist.getConv();
					//std::cout << "*****************************************************************************" << std::endl;


					t1 = std::chrono::high_resolution_clock::now();
					PFG.init(cas, &PQG);
					PFG.solve();
					t2 = std::chrono::high_resolution_clock::now();

					PResult.set(i, j, PFG.getP0());
					QResult.set(i, j, PFG.getQ0());
					temps.set(i, j, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
					temps2.set(i, j, PFG.getTime());
					iters.set(i, j, PFG.getIter());
					Residuals.set(i, j, PFG.getRes());
					ConvResult.set(i, j, PFG.getConv());
					i++;
					//std::cout << "*****************************************************************************" << std::endl;

					t1 = std::chrono::high_resolution_clock::now();


					PF.init(cas, &PQ);
					PF.solve();
					t2 = std::chrono::high_resolution_clock::now();

					PResult.set(i, j, PF.getP0());
					QResult.set(i, j, PF.getQ0());
					temps.set(i, j, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
					temps2.set(i, j, PF.getTime());
					iters.set(i, j, PF.getIter());
					Residuals.set(i, j, PF.getRes());
					ConvResult.set(i, j, PF.getConv());
					i++;

					//std::cout << "*****************************************************************************" << std::endl;
					if (conv != -1) {
						t1 = std::chrono::high_resolution_clock::now();
						PFGS.init(cas, &PQ);

						PFGS.solve();
						t2 = std::chrono::high_resolution_clock::now();

						PResult.set(i, j, PFGS.getP0());
						QResult.set(i, j, PFGS.getQ0());
						temps.set(i, j, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
						temps2.set(i, j, PFGS.getTime());
						iters.set(i, j, PFGS.getIter());
						Residuals.set(i, j, PFGS.getRes());
						ConvResult.set(i, j, PFGS.getConv());
						i++;
						//std::cout << "*****************************************************************************" << std::endl;

						t1 = std::chrono::high_resolution_clock::now();
						PFGSG.init(cas, &PQG);
						PFGSG.solve();
						t2 = std::chrono::high_resolution_clock::now();

						PResult.set(i, j, PFGSG.getP0());
						QResult.set(i, j, PFGSG.getQ0());
						temps.set(i, j, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
						temps2.set(i, j, PFGSG.getTime());
						iters.set(i, j, PFGSG.getIter());
						Residuals.set(i, j, PFGSG.getRes());
						ConvResult.set(i, j, PFGSG.getConv());
						i++;
					}
				}
			}
			std::cout << std::endl;
		}
		//std::cout << "-|-|-|-|-|-|-|-|-  -|-|-|-|-|-|-|-|-|" << std::endl;

	}
	/*Agents.display();
	temps.display();
	iters.display();*/
	ConvResult.display();
	float temptotal = temps2.sum();
	std::cout << "temps total " << temptotal << " temps moyen " << temptotal / (nSimu * nNBus * nCas * nMethod) << std::endl;

	Buses.saveCSV(fileName, mode);

	PResult.saveCSV(fileName, mode);
	QResult.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	temps2.saveCSV(fileName, mode);
	iters.saveCSV(fileName, mode);
	Residuals.saveCSV(fileName, mode);
	ConvResult.saveCSV(fileName, mode);

}

void SimuStatPFTransport() {
	std::string fileName = "PowerFlowStatTransport.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;

	float Pconso = 0.50;
	float dPconso = 0.01;
	float Propcons = 0.50;
	float PropGen = 0.5;
	float bProd = 5;
	float dbProd = 2;
	float Pprod = 0.50;
	float dPprod = 0.01;
	float gamma = 4;
	float dgamma = 2;
	float dQ = 0.1;
	float length = 0.005;
	float dlength = 0.0005;

	int nNBus = 1;
	int nBusMax = 50;
	const int offsetBus = 0;
	int nSimu = 3;

	float factorAgent[] = { 0.2, 0.5, 1, 2 };

	float factorLine[] = { 2, 1, 0.5 , 0.5 };
	float factorDLine[] = { 1, 0, 0.25 , 0,5 };

	int nCasAgent = 4;
	int nCasBuses = 4;
	int nCas = nCasAgent * nCasBuses;
	int nMethod = 4;


	MatrixCPU Param(1, 13);
	Param.set(0, 0, nBusMax);
	Param.set(0, 1, nNBus);
	Param.set(0, 2, Pconso);
	Param.set(0, 3, dPconso);
	Param.set(0, 4, Pprod);
	Param.set(0, 5, dPprod);
	Param.set(0, 6, dQ);
	Param.set(0, 7, length);
	Param.set(0, 8, dlength);
	Param.set(0, 9, nSimu);
	Param.set(0, 10, nCasAgent);
	Param.set(0, 11, nCasBuses);
	Param.set(0, 12, nMethod);

	Param.saveCSV(fileName, mode);

	MatrixCPU Buses(1, nNBus);

	MatrixCPU temps(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU temps2(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU iters(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU PResult(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU QResult(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU Residuals(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU ConvResult(nCas * nMethod, nNBus * nSimu, -3);
	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;

	for (int bus = 0; bus < nNBus; bus++) {
		std::cout << "--------- --------- --------- --------- ----------" << std::endl;
		int buses = offsetBus > 0 ? bus * (nBusMax - offsetBus) / MYMAX(nNBus - 1, 1) + offsetBus : (bus + 1) * nBusMax / nNBus;
		Buses.set(0, bus, buses);
		for (int casBus = 0; casBus < nCasBuses; casBus++) {
			int nLines = 0;
			int dnLines = 0;
			switch (casBus)
			{
			case 0:
				nLines = factorLine[0] * buses;
				dnLines = factorDLine[0];
				break;
			case 1:
				nLines = factorLine[1] * (buses * (buses - 1)) / 2;
				dnLines = factorDLine[1];
				break;
			case 2:
				nLines = factorLine[2] * (buses * sqrt(buses)) / 2;
				dnLines = factorDLine[2] * sqrt(buses);
				break;
			case 3:
				nLines = factorLine[3] * (buses * sqrt(buses)) / 2;
				dnLines = factorDLine[3] * sqrt(buses);
				break;
			default:
				std::cout << "error casBus =" << casBus << " not defined" << std::endl;
				return;
			}
			for (int casAgent = 0; casAgent < nCasAgent; casAgent++) {
				int agents = factorAgent[casAgent] * buses;
				for (int simu = 0; simu < nSimu; simu++) {
					CPUPF PF;
					GPUPF PFG;
					CPUPFGS PFGS;
					GPUPFGS PFGSG;


					int i = (casBus * nCasAgent + casAgent) * nMethod;
					int j = bus * nSimu + simu;

					StudyCase cas;
					std::cout << "-";
					cas.genGridHTB(buses, nLines, dnLines, length, dlength);
					cas.genAgentsAC(agents, Propcons, PropGen, Pconso, dPconso, bProd, dbProd, dQ, Pprod, dPprod, gamma, dgamma);
					cas.genLinkGridAgent();


					MatrixCPU PQ = cas.getPobj();
					MatrixGPU PQG = MatrixGPU(cas.getPobj(), 1);

					//std::cout << "*****************************************************************************" << std::endl;

					t1 = std::chrono::high_resolution_clock::now();


					PF.init(cas, &PQ);
					PF.solve();
					t2 = std::chrono::high_resolution_clock::now();

					PResult.set(i, j, PF.getP0());
					QResult.set(i, j, PF.getQ0());
					temps.set(i, j, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
					temps2.set(i, j, PF.getTime());
					iters.set(i, j, PF.getIter());
					Residuals.set(i, j, PF.getRes());
					ConvResult.set(i, j, PF.getConv());
					i++;
					int conv = PF.getConv();

					//std::cout << "*****************************************************************************" << std::endl;


					t1 = std::chrono::high_resolution_clock::now();
					PFG.init(cas, &PQG);
					PFG.solve();
					t2 = std::chrono::high_resolution_clock::now();

					PResult.set(i, j, PFG.getP0());
					QResult.set(i, j, PFG.getQ0());
					temps.set(i, j, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
					temps2.set(i, j, PFG.getTime());
					iters.set(i, j, PFG.getIter());
					Residuals.set(i, j, PFG.getRes());
					ConvResult.set(i, j, PFG.getConv());
					i++;


					//std::cout << "*****************************************************************************" << std::endl;
					if (conv != -1) {
						t1 = std::chrono::high_resolution_clock::now();
						PFGS.init(cas, &PQ);

						PFGS.solve();
						t2 = std::chrono::high_resolution_clock::now();

						PResult.set(i, j, PFGS.getP0());
						QResult.set(i, j, PFGS.getQ0());
						temps.set(i, j, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
						temps2.set(i, j, PFGS.getTime());
						iters.set(i, j, PFGS.getIter());
						Residuals.set(i, j, PFGS.getRes());
						ConvResult.set(i, j, PFGS.getConv());
						i++;
						//std::cout << "*****************************************************************************" << std::endl;

						t1 = std::chrono::high_resolution_clock::now();
						PFGSG.init(cas, &PQG);
						PFGSG.solve();
						t2 = std::chrono::high_resolution_clock::now();

						PResult.set(i, j, PFGSG.getP0());
						QResult.set(i, j, PFGSG.getQ0());
						temps.set(i, j, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
						temps2.set(i, j, PFGSG.getTime());
						iters.set(i, j, PFGSG.getIter());
						Residuals.set(i, j, PFGSG.getRes());
						ConvResult.set(i, j, PFGSG.getConv());
						i++;
					}
				}
			}
			std::cout << std::endl;
		}
		//std::cout << "-|-|-|-|-|-|-|-|-  -|-|-|-|-|-|-|-|-|" << std::endl;

	}
	/*Agents.display();
	temps.display();
	iters.display();*/
	ConvResult.display();
	float temptotal = temps2.sum();
	std::cout << "temps total " << temptotal << " temps moyen " << temptotal / (nSimu * nNBus * nCas * nMethod) << std::endl;

	Buses.saveCSV(fileName, mode);

	PResult.saveCSV(fileName, mode);
	QResult.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	temps2.saveCSV(fileName, mode);
	iters.saveCSV(fileName, mode);
	Residuals.saveCSV(fileName, mode);
	ConvResult.saveCSV(fileName, mode);

}

void SimuStatPFCompare() {
	std::string fileName = "ComparaisonPFCPU_Release_1500.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nMethode = 4;
	std::vector<int> indices = { 0, 1, 2, 3};
	std::string methodesName[nMethode] = { "Curr", "NR", "GS", "BackPQ" };

	CPUPF* methodes[nMethode];
	methodes[0] = new CPUPFdist;
	methodes[1] = new CPUPF;
	methodes[2] = new CPUPFGS;
	methodes[3] = new CPUPFdistPQ;


	if (nMethode != indices.size()) {
		throw std::domain_error("not good number of methods");
	}


	float Pconso = 0.05;
	float dPconso = 0.02;
	float Propcons = 0.5;
	float PropGen = 0.5;
	float bProd = 5;
	float dbProd = 2;
	float Pprod = 0.05;
	float dPprod = 0.02;
	float gamma = 4;
	float dgamma = 2;
	float dQ = 0.05;
	float length = 0.001;
	float dlength = 0.0005;

	int nNBus = 1;
	int nBusMax = 1500;
	const int offsetBus = 0;
	int nSimu = 50;
	int million  = 1000000;

	float factorAgent = 0.75;

	
	int nCasAgent = 1;
	int nCasBuses = 1;
	int nCas = nCasAgent * nCasBuses;



	MatrixCPU Param(1, 13);
	Param.set(0, 0, nBusMax);
	Param.set(0, 1, nNBus);
	Param.set(0, 2, Pconso);
	Param.set(0, 3, dPconso);
	Param.set(0, 4, Pprod);
	Param.set(0, 5, dPprod);
	Param.set(0, 6, dQ);
	Param.set(0, 7, length);
	Param.set(0, 8, dlength);
	Param.set(0, 9, nSimu);
	Param.set(0, 10, nCasAgent);
	Param.set(0, 11, nCasBuses);
	Param.set(0, 12, nMethode);

	Param.saveCSV(fileName, mode);

	MatrixCPU Buses(1, nNBus);
	MatrixCPU temps(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU temps2(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU iters(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU PResult(nMethode * nSimu, nNBus, nanf(""));

	MatrixCPU QResult(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU Residuals(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU ConvResult(nMethode * nSimu, nNBus, nanf(""));
    StudyCase cas;

	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;

	for (int bus = 0; bus < nNBus; bus++) {
		std::cout << "--------- --------- --------- --------- ----------" << std::endl;
		int buses = offsetBus > 0 ? bus * (nBusMax - offsetBus) / MYMAX(nNBus - 1, 1) + offsetBus : (bus + 1) * nBusMax / nNBus;
		Buses.set(0, bus, buses);
		
		int nDeep = 2 * buses;
		int nBranch = buses;
		int j = bus;	
		
		int agents = factorAgent * buses;
		for (int simu = 0; simu < nSimu; simu++) {
			
			std::cout << "-";
			cas.genGridBT(buses, nBranch, nDeep, length, dlength);
			cas.genAgentsAC(agents, Propcons, PropGen, Pconso, dPconso, bProd, dbProd, dQ, Pprod, dPprod, gamma, dgamma);
			cas.genLinkGridAgent();
			
			MatrixCPU PQ = cas.getPobj();
			//PQ.display();
			std::random_shuffle(indices.begin(), indices.end());
			for (int i = 0; i < nMethode; i++) {
				int k = indices[i] * nSimu + simu;
				//std::cout <<" methode  " << k <<std::endl;
				t1 = std::chrono::high_resolution_clock::now();
				methodes[indices[i]]->init(cas, &PQ);
				methodes[indices[i]]->solve();
				t2 = std::chrono::high_resolution_clock::now();

				PResult.set(k, j, methodes[indices[i]]->getP0());
				QResult.set(k, j, methodes[indices[i]]->getQ0());
				temps.set(k, j, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()/million);
				temps2.set(k, j, methodes[indices[i]]->getTime());
				iters.set(k, j, methodes[indices[i]]->getIter());
				Residuals.set(k, j, methodes[indices[i]]->getRes());
				ConvResult.set(k, j, methodes[indices[i]]->getConv());
			}
		}
		
		std::cout << std::endl;
		//std::cout << "-|-|-|-|-|-|-|-|-  -|-|-|-|-|-|-|-|-|" << std::endl;

	}
	
	temps.display();
	/*iters.display();
	PResult.display();
	*/
	ConvResult.display();
	
	float temptotal = temps2.sum();
	std::cout << "temps total " << temptotal << " temps moyen " << temptotal / (nSimu * nNBus * nCas * nMethode) << std::endl;

	Buses.saveCSV(fileName, mode);

	PResult.saveCSV(fileName, mode);
	QResult.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	temps2.saveCSV(fileName, mode);
	iters.saveCSV(fileName, mode);
	Residuals.saveCSV(fileName, mode);
	ConvResult.saveCSV(fileName, mode);

}
 
/* Simu OPF*/

void SimuTemporalTestFeeder() {
	std::string path = "data/ACGrid/";
	std::string fileName2 = "SimutemporalFBTestFeeder.csv";
	std::string fileName = "SimutemporalTestFeederAll.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nMethode = 4;
	std::string methodesName[nMethode] = { "OPFADMM", "OPFADMMGPU", "OPFADMM2", "OPFADMMGPU2"};
	Method* methodes[nMethode];
	methodes[0] = new OPFADMM;
	methodes[1] = new OPFADMMGPU;
	methodes[2] = new OPFADMM2;
	methodes[3] = new OPFADMMGPU2;

	bool saveTime = true;


	float epsG = 0.0005f;
	float epsL = 0.00001f;

	int iterG = 30000;
	int iterL = 10000;//500;

	int chosenAgenGen = 0;
	int nAgent = 57;
	if (chosenAgenGen != 0) {
		nAgent = 57; // may change
	}

	int begin = 0;
	int end = 60 * 1 - 1;

	float rho = 2;


	float stepG = 100;
	float stepL = 1;

	System sys;
	sys.setIter(iterG, iterL);
	sys.setStep(stepG, stepL);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);

	MatrixCPU Param(1, 11);
	Param.set(0, 0, chosenAgenGen);
	Param.set(0, 1, nAgent);
	Param.set(0, 2, epsG);
	Param.set(0, 3, epsL);
	Param.set(0, 4, iterG);
	Param.set(0, 5, iterL);
	Param.set(0, 6, stepG);
	Param.set(0, 7, stepL);
	Param.set(0, 8, rho);
	Param.set(0, 9, begin);
	Param.set(0, 10, end);

	if (saveTime) {
		Param.saveCSV(fileName, mode);
	}
#ifdef INSTRUMENTATION
	Param.saveCSV(fileName2, mode);
#endif // INSTRUMENTAION



	for (int i = 0; i < 1; i++) {
		
		sys.setMethod(methodes[i]);
		sys.setRho(rho);
		std::cout << "----------------Simu---------------------------- " << std::endl;
		// Debut simu
		std::cout << "methode " << methodesName[i] << std::endl;
		clock_t t = clock();

		sys.solveIntervalle(path, begin, end, chosenAgenGen);
		t = clock() - t;

		std::cout << "calculation time : " << (float)t / CLOCKS_PER_SEC << std::endl;


		MatrixCPU temps(sys.getTemps());
		MatrixCPU iter(sys.getIter());
		MatrixCPU conv(sys.getConv());
		MatrixCPU fc(sys.getFc());
		MatrixCPU ResR(sys.getResR());
		MatrixCPU ResS(sys.getResS());
		MatrixCPU ResX(sys.getResX());

		temps.display();
		ResR.display();
		ResS.display();
		ResX.display();
		iter.display();

		if (saveTime) {
			temps.saveCSV(fileName, mode);
			iter.saveCSV(fileName, mode);
			fc.saveCSV(fileName, mode);
			ResR.saveCSV(fileName, mode);
			ResS.saveCSV(fileName, mode);
			ResX.saveCSV(fileName, mode);
			conv.saveCSV(fileName, mode);/**/
		}


		std::cout << "-------------------------------------------------------- " << std::endl;
#ifdef INSTRUMENTATION
		sys.displayTime(fileName2);
#endif // INSTRUMENTATION


		sys.resetParam();
	}
	for (int i = 0; i < nMethode; i++) {
		DELETEB(methodes[i]);
	}
	sys.setMethod(nullptr);

}


void SimuStatOPF() {
	std::string fileName = "OptimalPowerFlowStat.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nMethode = 5;
	std::vector<int> indices = { 0, 1, 2, 3, 4 };
	std::string methodesName[nMethode] = { "PDIPM", "OPFADMM", "OPFADMM2","OPFADMMGPU", "OPFADMMGPU2" };
	Method* methodes[nMethode];
	//methodes[0] = new OPFPDIPM;
	methodes[1] = new OPFADMM;
	methodes[2] = new OPFADMM2;
	methodes[3] = new OPFADMMGPU;
	methodes[4] = new OPFADMMGPU2;

	float Pconso = 0.005;
	float dPconso = 0.001;
	float Propcons = 0.4;
	float PropGen = 0.25;
	float bProd = 5;
	float dbProd = 2;
	float Pprod = 0.005;
	float dPprod = 0.001;
	float gamma = 0;
	float dgamma = 0;
	float dQ = 0.005;
	float length = 0.001;
	float dlength = 0.0005;

	int nNBus = 1;
	int nBusMax = 50;
	const int offsetBus = 0;
	int nSimu = 10;
	float rhoAgent = 0.25;
	
	int iterGlobal = 4000;
	int iterLocal = 1000;
	int stepG = 10;
	int stepL = 10;
	float epsG = 0.01f;
	float epsL = 0.001f;


	float factorAgent[] = { 0.2, 0.5, 1, 2, 5};

	float factorBuses[] = { 2, 0.3, 0.3 , 1 };
	int nCasAgent = 5;
	int nCasBuses = 4;
	int nCas = nCasAgent * nCasBuses;
	int nMethod = 5;


	MatrixCPU Param(1, 13);
	Param.set(0, 0, nBusMax);
	Param.set(0, 1, nNBus);
	Param.set(0, 2, Pconso);
	Param.set(0, 3, dPconso);
	Param.set(0, 4, Pprod);
	Param.set(0, 5, dPprod);
	Param.set(0, 6, dQ);
	Param.set(0, 7, length);
	Param.set(0, 8, dlength);
	Param.set(0, 9, nSimu);
	Param.set(0, 10, nCasAgent);
	Param.set(0, 11, nCasBuses);
	Param.set(0, 12, nMethod);

	Param.saveCSV(fileName, mode);

	MatrixCPU Buses(1, nNBus);

	MatrixCPU temps(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU temps2(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU iters(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU PResult(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU QResult(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU Residuals(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU Residuals2(nCas * nMethod, nNBus * nSimu, -3);
	MatrixCPU ConvResult(nCas * nMethod, nNBus * nSimu, -3);
	System sys;
	sys.setIter(iterGlobal, iterLocal);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	sys.setStep(stepG, stepL);

	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;

	for (int bus = 0; bus < nNBus; bus++) {
		std::cout << "--------- --------- --------- --------- ----------" << std::endl;
		int buses = offsetBus > 0 ? bus * (nBusMax - offsetBus) / MYMAX(nNBus - 1, 1) + offsetBus : (bus + 1) * nBusMax / nNBus;
		Buses.set(0, bus, buses);
		for (int casBus = 0; casBus < nCasBuses; casBus++) {
			int nDeep = 0;
			int nBranch = 0;
			switch (casBus)
			{
			case 0:
				nDeep = sqrt(factorBuses[0] * buses);
				nBranch = sqrt(factorBuses[0] * buses);
				break;
			case 1:
				nDeep = sqrt(factorBuses[1] * buses);
				nBranch = buses;
				break;
			case 2:
				nDeep = buses;
				nBranch = sqrt(factorBuses[2] * buses);
				break;
			case 3:
				nDeep = factorBuses[3] * buses;
				nBranch = factorBuses[3] * buses;
				break;
			default:
				std::cout << "error casBus =" << casBus << " not defined" << std::endl;
				return;
			}
			for (int casAgent = 0; casAgent < nCasAgent; casAgent++) {
				int agents = factorAgent[casAgent] * buses;
				for (int simu = 0; simu < nSimu; simu++) {
					int i = (casBus * nCasAgent + casAgent) * nMethod;
					int j = bus * nSimu + simu;
					StudyCase cas;
					
					cas.genGridBT(buses, nBranch, nDeep, length, dlength);
					std::cout << "*";
					cas.genAgentsAC(agents, Propcons, PropGen, Pconso, dPconso, bProd, dbProd, dQ, Pprod, dPprod, gamma, dgamma);
					std::cout << "|";
					cas.genLinkGridAgent();
					std::cout << "/";
					sys.setStudyCase(cas);
					std::cout << "-";
					float rho = rhoAgent * agents;
					sys.setRho(rho);
					std::random_shuffle(indices.begin(), indices.end());
					for (int m = 0; m < nMethode; m++) {
						sys.setMethod(methodes[indices[m]]);

						t1 = std::chrono::high_resolution_clock::now();
						Simparam res = sys.solve();
						MatrixCPU Pn = res.getPn();
						MatrixCPU Residuals = res.getRes();
						int iter = res.getIter();
						t2 = std::chrono::high_resolution_clock::now();

						PResult.set(i + indices[m], j, Pn.get(1,0));
						QResult.set(i + indices[m], j, Pn.get(1 + agents, 0));
						temps.set(i + indices[m], j, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
						temps2.set(i + indices[m], j, res.getTime());
						iters.set(i + indices[m], j, iter);
						Residuals.set(i + indices[m], j, Residuals.get(0, (iter - 1) / stepG));
						Residuals2.set(i + indices[m], j, Residuals.get(1, (iter - 1) / stepG));
					}
				}
			}
			std::cout << std::endl;
		}
		//std::cout << "-|-|-|-|-|-|-|-|-  -|-|-|-|-|-|-|-|-|" << std::endl;

	}
	for (int i = 0; i < nMethode; i++) {
		DELETEB(methodes[i]);
	}
	sys.setMethod(nullptr);
	/*Agents.display();
	temps.display();
	iters.display();*/
	ConvResult.display();
	float temptotal = temps2.sum();
	std::cout << "temps total " << temptotal << " temps moyen " << temptotal / (nSimu * nNBus * nCas * nMethod) << std::endl;

	Buses.saveCSV(fileName, mode);

	PResult.saveCSV(fileName, mode);
	QResult.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	temps2.saveCSV(fileName, mode);
	iters.saveCSV(fileName, mode);
	Residuals.saveCSV(fileName, mode);
	ConvResult.saveCSV(fileName, mode);

}

void SimuStatOPFCompare() {
	std::string fileName = "ComparaisonOPFCPU_Release_300_v2_bonus.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	/*const int nMethode = 5;
	std::vector<int> indices = { 0, 1, 2, 3, 4 };
	std::string methodesName[nMethode] = { "PDIPM", "OPFADMM", "OPFADMM2","OPFADMMGPU", "OPFADMMGPU2" };*/
	const int nMethode = 3;
	std::vector<int> indices = { 0, 1, 2 };
	std::string methodesName[nMethode] = { "PDIPM", "OPFADMM", "OPFADMM2" };
	MethodOPF* methodes[nMethode];
	//methodes[0] = new OPFPDIPM;
	methodes[1] = new OPFADMM;
	methodes[2] = new OPFADMM2;
	//methodes[3] = new OPFADMMGPU;
	//methodes[4] = new OPFADMMGPU2;


	if (nMethode != indices.size()) {
		throw std::domain_error("not good number of methods");
	}

	// Market
	float Pconso = 0.05;
	float dPconso = 0.02;
	float Propcons = 0.5;
	float PropNFleGen = 0.25;
	float PropGen = 1 - PropNFleGen - Propcons;
	float bProd = 5;
	float dbProd = 2;
	float Pprod = 0.05;
	float dPprod = 0.02;
	float gamma = 4; // inutile car pas de trade
	float dgamma = 2; // idem
	float dQ = 0.05;
	float length = 0.001;
	float dlength = 0.0005;
	float factorAgent = 0.75;


	// cases
	int nNBus = 1;
	int nBusMax = 300;
	const int offsetBus = 0;
	int nSimu = 2;
	int million = 1000000;
	int nCasAgent = 1;
	int nCasBuses = 1;
	//int nCas = nCasAgent * nCasBuses;
	

	// simulation
	int iterGlobal = 50000;
	int iterLocal = 1000;
	int stepG = 10;
	int stepL = 10;
	float epsG = 0.01f;
	float epsL = 0.001f;
	float rho = 2;



	MatrixCPU Param(1, 22);
	Param.set(0, 0, nBusMax);
	Param.set(0, 1, nNBus);
	Param.set(0, 2, length);
	Param.set(0, 3, dlength);
	Param.set(0, 4, Pconso);
	Param.set(0, 5, dPconso);
	Param.set(0, 6, Pprod);
	Param.set(0, 7, dPprod);
	Param.set(0, 8, dQ);
	Param.set(0, 9, Propcons);
	Param.set(0, 10, PropGen);
	Param.set(0, 11, nSimu);
	Param.set(0, 12, nCasAgent);
	Param.set(0, 13, nCasBuses);
	Param.set(0, 14, nMethode);
	Param.set(0, 15, rho);
	Param.set(0, 16, epsG);
	Param.set(0, 17, epsL);
	Param.set(0, 18, iterGlobal);
	Param.set(0, 19, iterLocal);
	Param.set(0, 20, stepG);
	Param.set(0, 21, stepL);



	Param.saveCSV(fileName, mode);


	MatrixCPU Buses(1, nNBus);
	MatrixCPU temps(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU iters(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU PResult(nMethode * nSimu, nNBus, nanf(""));

	MatrixCPU QResult(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU ResR(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU ResS(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU ResV(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU Fc(nMethode * nSimu, nNBus, nanf(""));
	System sys;
	StudyCase cas;
	Simparam res;
	sys.setIter(iterGlobal, iterLocal);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	sys.setStep(stepG, stepL);
	sys.setRho(rho);

	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;

	for (int bus = 0; bus < nNBus; bus++) {
		std::cout << "--------- --------- --------- --------- ----------" << std::endl;
		int buses = offsetBus > 0 ? bus * (nBusMax - offsetBus) / MYMAX(nNBus - 1, 1) + offsetBus : (bus + 1) * nBusMax / nNBus;
		Buses.set(0, bus, buses);

		int nDeep = 2 * buses;
		int nBranch = buses;
		int j = bus;

		int agents = factorAgent * buses;
		for (int simu = 0; simu < nSimu; simu++) {

			std::cout << "-";
			cas.genGridBT(buses, nBranch, nDeep, length, dlength);
			cas.genAgentsAC(agents, Propcons, PropGen, Pconso, dPconso, bProd, dbProd, dQ, Pprod, dPprod, gamma, dgamma);
			cas.genLinkGridAgent();

			sys.setStudyCase(cas);
			std::random_shuffle(indices.begin(), indices.end());
			for (int i = 0; i < nMethode; i++) {
				int k = indices[i] * nSimu + simu;
				//std::cout << " Simu  " << simu << " methode  " << i << std::endl;
				//if (indices[i] > 0) {
					sys.setMethod(methodes[indices[i]]);
					t1 = std::chrono::high_resolution_clock::now();
					res = sys.solve();
					t2 = std::chrono::high_resolution_clock::now();
					MatrixCPU Pn = res.getPn();
			
					MatrixCPU ResF = res.getRes();
					int iter = res.getIter();
				
					PResult.set(k, j, Pn.get(1, 0));
					QResult.set(k, j, Pn.get(agents + 2, 0));
					temps.set(k, j, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
					iters.set(k, j, iter );
					ResR.set(k, j, ResF.get(0, (iter - 1) / stepG));
					ResS.set(k, j, ResF.get(1, (iter - 1) / stepG));
					ResV.set(k, j, ResF.get(2, (iter - 1) / stepG));
					Fc.set(k, j, res.getFc());
					sys.resetParam();
				//}
			}
		}

		std::cout << std::endl;
		//std::cout << "-|-|-|-|-|-|-|-|-  -|-|-|-|-|-|-|-|-|" << std::endl;

	}

		
	Buses.saveCSV(fileName, mode);

	PResult.saveCSV(fileName, mode);
	QResult.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	iters.saveCSV(fileName, mode);
	ResR.saveCSV(fileName, mode);
	ResS.saveCSV(fileName, mode);
	ResV.saveCSV(fileName, mode);
	Fc.saveCSV(fileName, mode);


	for (int i = 0; i < nMethode; i++) {
		DELETEB(methodes[i]);
	}
	sys.setMethod(nullptr);
}

void SimuCompareISGT() {
	std::string fileName = "ComparaisonOPFISGT200.csv";
	std::string fileName2 = "Residuals400bis";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	/*const int nMethode = 5;
	std::vector<int> indices = { 0, 1, 2, 3, 4 };
	std::string methodesName[nMethode] = { "PDIPM", "OPFADMM", "OPFADMM2","OPFADMMGPU", "OPFADMMGPU2" };*/
	const int nMethode = 4;
	std::vector<int> indices = { 0, 1, 2, 3};
	std::string methodesName[nMethode] = { "OPFADMM", "OPFADMM2","OPFADMMGPU", "OPFADMMGPU2" };
	Method* methodes[nMethode];
	
	methodes[0] = new OPFADMM;
	methodes[1] = new OPFADMM2;
	methodes[2] = new OPFADMMGPU;
	methodes[3] = new OPFADMMGPU2;


	if (nMethode != indices.size()) {
		throw std::domain_error("not good number of methods");
	}

	// Market
	float Pconso = 0.05;
	float dPconso = 0.02;
	float Propcons = 0.5;
	float PropNFleGen = 0.25;
	float PropGen = 1 - PropNFleGen - Propcons;
	float bProd = 5;
	float dbProd = 2;
	float Pprod = 0.05;
	float dPprod = 0.02;
	float gamma = 4; // inutile car pas de trade
	float dgamma = 2; // idem
	float dQ = 0.05;
	float length = 0.001;
	float dlength = 0.0005;
	float factorAgent = 0.75;


	// cases
	int nNBus = 1;
	int nBusMax = 400;
	const int offsetBus = 0;
	int nSimu = 50;
	int million = 1000000;
	int nCasAgent = 1;
	int nCasBuses = 1;
	//int nCas = nCasAgent * nCasBuses;


	// simulation
	int iterGlobal = 3000;
	int iterLocal = 2000;
	int stepG = 1;
	int stepL = 1;
	float epsG = 0.01f;
	float epsL = 0.00001f;
	float rho = 2;



	MatrixCPU Param(1, 22);
	Param.set(0, 0, nBusMax);
	Param.set(0, 1, nNBus);
	Param.set(0, 2, length);
	Param.set(0, 3, dlength);
	Param.set(0, 4, Pconso);
	Param.set(0, 5, dPconso);
	Param.set(0, 6, Pprod);
	Param.set(0, 7, dPprod);
	Param.set(0, 8, dQ);
	Param.set(0, 9, Propcons);
	Param.set(0, 10, PropGen);
	Param.set(0, 11, nSimu);
	Param.set(0, 12, nCasAgent);
	Param.set(0, 13, nCasBuses);
	Param.set(0, 14, nMethode);
	Param.set(0, 15, rho);
	Param.set(0, 16, epsG);
	Param.set(0, 17, epsL);
	Param.set(0, 18, iterGlobal);
	Param.set(0, 19, iterLocal);
	Param.set(0, 20, stepG);
	Param.set(0, 21, stepL);



	//	Param.saveCSV(fileName, mode);


	MatrixCPU Buses(1, nNBus);
	MatrixCPU temps(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU iters(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU PResult(nMethode * nSimu, nNBus, nanf(""));

	MatrixCPU QResult(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU ResR(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU ResS(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU ResV(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU CountRelax(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU Fc(nMethode * nSimu, nNBus, nanf(""));
	System sys;
	StudyCase cas;
	Simparam res;
	sys.setIter(iterGlobal, iterLocal);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	sys.setStep(stepG, stepL);
	sys.setRho(rho);

	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;

	for (int bus = 0; bus < nNBus; bus++) {
		std::cout << "--------- --------- --------- --------- ----------" << std::endl;
		int buses = offsetBus > 0 ? bus * (nBusMax - offsetBus) / MYMAX(nNBus - 1, 1) + offsetBus : (bus + 1) * nBusMax / nNBus;
		Buses.set(0, bus, buses);

		int nDeep = 2 * buses;
		int nBranch = buses;
		int j = bus;

		int agents = factorAgent * buses;
		std::cout << " nBus " << buses << " nAgent " << agents << std::endl;
		for (int simu = 0; simu < nSimu; simu++) {

			std::cout << "-";
			cas.genGridBT(buses, nBranch, nDeep, length, dlength);
			cas.genAgentsAC(agents, Propcons, PropGen, Pconso, dPconso, bProd, dbProd, dQ, Pprod, dPprod, gamma, dgamma);
			cas.genLinkGridAgent();

		
			sys.setStudyCase(cas);
			std::random_shuffle(indices.begin(), indices.end());
			for (int i = 0; i < nMethode; i++) {
				int k = indices[i] * nSimu + simu;
				//std::cout << " Simu  " << simu << " methode  " << indices[i] << std::endl;
				
				sys.setMethod(methodes[indices[i]]);
				t1 = std::chrono::high_resolution_clock::now();
				res = sys.solve();
				t2 = std::chrono::high_resolution_clock::now();
				//int count = methodes[indices[i]]->feasiblePoint();
				MatrixCPU Pn = res.getPn();
				MatrixCPU ResF = res.getRes();
				int iter = res.getIter();
				ResF.saveCSV(fileName2 + methodesName[indices[i]] + ".csv");
				std::cout << indices[i] << " " << iter << ", " << (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million ;
				//PResult.set(k, j, Pn.get(1, 0));
				//QResult.set(k, j, Pn.get(agents + 2, 0));
				//temps.set(k, j, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
				iters.set(k, j, iter);
				//ResR.set(k, j, ResF.get(0, (iter - 1) / stepG));
				//ResS.set(k, j, ResF.get(1, (iter - 1) / stepG));
				//ResV.set(k, j, ResF.get(2, (iter - 1) / stepG));
				//CountRelax.set(k, j, count);
				//Fc.set(k, j, res.getFc());
				//sys.resetParam();
				
			}
			std::cout << std::endl;
		}

		std::cout << std::endl;
		//std::cout << "-|-|-|-|-|-|-|-|-  -|-|-|-|-|-|-|-|-|" << std::endl;

	}


	Buses.saveCSV(fileName, mode);

	PResult.saveCSV(fileName, mode);
	QResult.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	iters.saveCSV(fileName, mode);
	ResR.saveCSV(fileName, mode);
	ResS.saveCSV(fileName, mode);
	//ResV.saveCSV(fileName, mode);
	Fc.saveCSV(fileName, mode);/**/
	//CountRelax.saveCSV(fileName, mode);
	iters.display();


	for (int i = 0; i < nMethode; i++) {
		DELETEB(methodes[i]);
	}
	sys.setMethod(nullptr);
}



/* Simu Marche Endo*/

void SimuStatMarketEndo() {
	std::string fileName = "ComparaisonMarketEndo_All_150.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	/*const int nMethode = 5;
	std::vector<int> indices = { 0, 1, 2, 3, 4 };
	std::string methodesName[nMethode] = { "PDIPM", "OPFADMM", "OPFADMM2","OPFADMMGPU", "OPFADMMGPU2" };*/
	//const int nMethode = 4;
	//std::vector<int> indices = { 0, 1, 2, 3 };
	//std::string methodesName[nMethode] = { "EndoDirect", "EndoConsensus","DC-EndoPF", "AC-EndoPF" };
	
	const int nMethode = 8;
	std::vector<int> indices = { 0, 1, 2, 3, 4, 5, 6, 7 };
	std::string methodesName[nMethode] = { "EndoDirect", "EndoConsensus","DC-EndoPF", "AC-EndoPF", "EndoDirectGPU", "EndoConsensusGPU","DC-EndoPFGPU", "AC-EndoPFGPU" };
	
	Method* methodes[nMethode];

	bool save = true;

	methodes[0] = new MarketEndoDirect;
	methodes[1] = new MarEndoCons;
	methodes[2] = new ADMMConst;
	methodes[3] = new EndoPF;
	methodes[4] = new MarketEndoDirectGPU;
	methodes[5] = new MarEndoConsGPU;
	methodes[6] = new ADMMGPUConst4;
	methodes[7] = new EndoPFGPU;

	if (nMethode != indices.size()) {
		throw std::domain_error("not good number of methods");
	}

	// Market
	float Pconso = 0.05;
	float dPconso = 0.01;
	float Propcons = 0.5;
	float PropNFleGen = 0.25;
	float PropGen = 1 - PropNFleGen - Propcons;
	float bProd = 1;
	float dbProd = 0.1;
	float Pprod = 0.01;
	float dPprod = 0.005;
	float gamma = 1; // 
	float dgamma = 0.2; //  
	float dQ = 0.01;
	float length = 0.001;
	float dlength = 0.0005;
	float factorAgent = 1;
	
	
	// cases
	int nNBus = 1;
	int nBusMax = 150;
	const int offsetBus = 0;
	int nSimu = 50;
	int million = 1000000;
	int nCasAgent = 1;
	int nCasBuses = 1;
	//int nCas = nCasAgent * nCasBuses;


	// simulation
	int iterGlobal = 10000;
	int iterLocal = 2000;
	int stepG = 1;
	int stepL = 1;
	float epsG = 0.01f;
	float epsL = 0.001f;
	float rho = 10;
	float rho1 = rho/2;
	float epsGC = 0.1;

	


	MatrixCPU Param(1, 22);
	Param.set(0, 0, nBusMax);
	Param.set(0, 1, nNBus);
	Param.set(0, 2, length);
	Param.set(0, 3, dlength);
	Param.set(0, 4, Pconso);
	Param.set(0, 5, dPconso);
	Param.set(0, 6, Pprod);
	Param.set(0, 7, dPprod);
	Param.set(0, 8, dQ);
	Param.set(0, 9, Propcons);
	Param.set(0, 10, PropGen);
	Param.set(0, 11, nSimu);
	Param.set(0, 12, nCasAgent);
	Param.set(0, 13, nCasBuses);
	Param.set(0, 14, nMethode);
	Param.set(0, 15, rho);
	Param.set(0, 16, epsG);
	Param.set(0, 17, epsL);
	Param.set(0, 18, iterGlobal);
	Param.set(0, 19, iterLocal);
	Param.set(0, 20, stepG);
	Param.set(0, 21, stepL);


	if (save) {
		Param.saveCSV(fileName, mode);
	}
	


	MatrixCPU Buses(1, nNBus);
	MatrixCPU temps(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU iters(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU PResult(nMethode * nSimu, nNBus, nanf(""));

	MatrixCPU QResult(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU ResR(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU ResS(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU ResV(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU Fc(nMethode * nSimu, nNBus, nanf(""));
	System sys;
	
	/*sys.setIter(iterGlobal, iterLocal);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	sys.setStep(stepG, stepL);
	sys.setRho(rho);*/

	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;

	for (int bus = 0; bus < nNBus; bus++) {
		std::cout << "--------- --------- --------- --------- ----------" << std::endl;
		int buses = offsetBus > 0 ? bus * (nBusMax - offsetBus) / MYMAX(nNBus - 1, 1) + offsetBus : (bus + 1) * nBusMax / nNBus;
		Buses.set(0, bus, buses);

		int nDeep = 2 * buses;
		int nBranch = buses;
		int j = bus;

		int agents = factorAgent * buses;
		//std::cout << " cas avec " << buses << " bus et " << agents << "agents et une profondeur de " << nDeep << " et largeur de " << nBranch << std::endl;

		for (int simu = 0; simu < nSimu; simu++) {

			std::cout << "-";
		
			
			
			StudyCase cas;
			cas.genGridBT(buses, nBranch, nDeep, length, dlength);
			cas.genAgentsAC(agents, Propcons, PropGen, Pconso, dPconso, bProd, dbProd, dQ, Pprod, dPprod, gamma, dgamma);
			
			
			cas.genLinkGridAgent();
			cas.genDCGridFromAC(); // pour utiliser les m�thodes DC
			cas.setReduce(true);
			
			Simparam param(cas.getNagent(), cas.getNLine(true), true);
			param.setEpsL(epsL);
			param.setEpsG(epsG);//0.0001f
			param.setItG(iterGlobal); //500000
			param.setItL(iterLocal);
			param.setStep(stepG, stepL);
			param.setRho(rho);
			param.setRho1(rho1);
			param.setEpsGC(epsGC);

			Simparam res(param);
			//sys.setStudyCase(cas);
			std::random_shuffle(indices.begin(), indices.end());
			for (int i = 0; i < nMethode; i++) {
				int k = indices[i] * nSimu + simu;
				//std::cout << " Simu  " << simu << " methode  " << i << " : " << methodesName[indices[i]] << std::endl;
				//if (indices[i] < 3) {
				
				//sys.setMethod(methodes[indices[i]]);
				t1 = std::chrono::high_resolution_clock::now();
				methodes[indices[i]]->solve(&res, param, cas);
				t2 = std::chrono::high_resolution_clock::now();
				
				MatrixCPU Pn = res.getPn();
				
				MatrixCPU ResF = res.getRes();
				int iter = res.getIter();
				//std::cout << "enregistrement de " << Pn.get(1, 0) << " en pose " << k << " methode " << indices[i] << std::endl;
				PResult.set(k, j, Pn.get(1, 0));
				QResult.set(k, j, Pn.get(agents + 2, 0));
				temps.set(k, j, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
				iters.set(k, j, iter);
				ResR.set(k, j, ResF.get(0, (iter - 1) / stepG));
				ResS.set(k, j, ResF.get(1, (iter - 1) / stepG));
				ResV.set(k, j, ResF.get(2, (iter - 1) / stepG));

				Fc.set(k, j, res.getFc());
				sys.resetParam();
				//}
			}
		}

		std::cout << std::endl;
		//std::cout << "-|-|-|-|-|-|-|-|-  -|-|-|-|-|-|-|-|-|" << std::endl;

	}

	if (save) {
		Buses.saveCSV(fileName, mode);

		PResult.saveCSV(fileName, mode);
		QResult.saveCSV(fileName, mode);
		temps.saveCSV(fileName, mode);
		iters.saveCSV(fileName, mode);
		ResR.saveCSV(fileName, mode);
		ResS.saveCSV(fileName, mode);
		ResV.saveCSV(fileName, mode);
		Fc.saveCSV(fileName, mode);
	}
	
	
	iters.display();
	PResult.display();

	for (int i = 0; i < nMethode; i++) {
		DELETEB(methodes[i]);
	}
	sys.setMethod(nullptr);
}

void SimuStatMarketEndoAC() {
	std::string fileName = "ComparaisonMarketEndoACAll_800.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;

	const int nMethode = 6;
	std::vector<int> indices = { 0, 1, 2, 3, 4, 5 };
	std::string methodesName[nMethode] = { "EndoDirect", "EndoConsensus", "AC EndoPF", "EndoDirectGPU", "EndoConsensusGPU", "AC EndoPFGPU"};
	Method* methodes[nMethode];
	methodes[0] = new MarketEndoDirect;
	methodes[1] = new MarEndoCons;
	methodes[2] = new EndoPF;
	methodes[3] = new MarketEndoDirectGPU;
	methodes[4] = new MarEndoConsGPU;
	methodes[5] = new EndoPFGPU;
	bool save = true;

	

	if (nMethode != indices.size()) {
		throw std::domain_error("not good number of methods");
	}

	// Market
	float Pconso = 0.05;
	float dPconso = 0.01;
	float Propcons = 0.5;
	float PropNFleGen = 0.25;
	float PropGen = 1 - PropNFleGen - Propcons;
	float bProd = 1;
	float dbProd = 0.1;
	float Pprod = 0.01;
	float dPprod = 0.005;
	float gamma = 1; // 
	float dgamma = 0.2; //  
	float dQ = 0.01;
	float length = 0.001;
	float dlength = 0.0005;
	float factorAgent = 1;


	// cases
	int nNBus = 1;
	int nBusMax = 800;
	const int offsetBus = 0;
	int nSimu = 50;
	//int million = 1000000;
	int nCasAgent = 1;
	int nCasBuses = 1;
	//int nCas = nCasAgent * nCasBuses;


	// simulation
	int iterGlobal = 5000;
	int iterLocal = 5000;
	int stepG = 1;
	int stepL = 1;
	float epsG = 0.1f;
	float epsL = 0.001f;
	float rho = 10;
	float rho1 = rho / 2;
	float epsGC = 0.05;




	MatrixCPU Param(1, 23);
	Param.set(0, 0, nBusMax);
	Param.set(0, 1, nNBus);
	Param.set(0, 2, length);
	Param.set(0, 3, dlength);
	Param.set(0, 4, Pconso);
	Param.set(0, 5, dPconso);
	Param.set(0, 6, Pprod);
	Param.set(0, 7, dPprod);
	Param.set(0, 8, dQ);
	Param.set(0, 9, Propcons);
	Param.set(0, 10, PropGen);
	Param.set(0, 11, nSimu);
	Param.set(0, 12, nCasAgent);
	Param.set(0, 13, nCasBuses);
	Param.set(0, 14, nMethode);
	Param.set(0, 15, rho);
	Param.set(0, 16, epsG);
	Param.set(0, 17, epsGC);
	Param.set(0, 18, epsL);
	Param.set(0, 19, iterGlobal);
	Param.set(0, 20, iterLocal);
	Param.set(0, 21, stepG);
	Param.set(0, 22, stepL);


	if (save) {
		Param.saveCSV(fileName, mode);
	}



	MatrixCPU Buses(1, nNBus);
	MatrixCPU temps(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU iters(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU PResult(nMethode * nSimu, nNBus, nanf(""));

	MatrixCPU QResult(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU ResR(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU ResS(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU ResV(nMethode * nSimu, nNBus, nanf(""));
	MatrixCPU Fc(nMethode * nSimu, nNBus, nanf(""));
	//System sys;

	/*sys.setIter(iterGlobal, iterLocal);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	sys.setStep(stepG, stepL);
	sys.setRho(rho);*/

	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;

	for (int bus = 0; bus < nNBus; bus++) {
		std::cout << "--------- --------- --------- --------- ----------" << std::endl;
		int buses = offsetBus > 0 ? bus * (nBusMax - offsetBus) / MYMAX(nNBus - 1, 1) + offsetBus : (bus + 1) * nBusMax / nNBus;
		Buses.set(0, bus, buses);

		int nDeep = 2 * buses;
		int nBranch = buses;
		int j = bus;

		int agents = factorAgent * buses;
		//std::cout << " cas avec " << buses << " bus et " << agents << "agents et une profondeur de " << nDeep << " et largeur de " << nBranch << std::endl;

		for (int simu = 0; simu < nSimu; simu++) {
		

			std::cout << "-";



			StudyCase cas;
			cas.genGridBT(buses, nBranch, nDeep, length, dlength);
			cas.genAgentsAC(agents, Propcons, PropGen, Pconso, dPconso, bProd, dbProd, dQ, Pprod, dPprod, gamma, dgamma);


			cas.genLinkGridAgent();
			cas.genDCGridFromAC(); // pour utiliser les m�thodes DC
			cas.setReduce(true);

			Simparam param(cas.getNagent(), cas.getNLine(true), true);
			param.setEpsL(epsL);
			param.setEpsG(epsG);//
			param.setItG(iterGlobal); //1000
			param.setItL(iterLocal);
			param.setStep(stepG, stepL);
			param.setRho(rho);
			param.setRho1(rho1);
			param.setEpsGC(epsGC);

			Simparam res(param);
			//sys.setStudyCase(cas);
			std::random_shuffle(indices.begin(), indices.end());
			for (int i = 0; i < nMethode; i++) {
				int k = indices[i] * nSimu + simu;
				//std::cout << " Simu  " << simu << " methode  " << indices[i] << " : " << methodesName[indices[i]] << std::endl;
				//if (indices[i] < 3) {

				//sys.setMethod(methodes[indices[i]]);
				t1 = std::chrono::high_resolution_clock::now();
				methodes[indices[i]]->solve(&res, param, cas);
				t2 = std::chrono::high_resolution_clock::now();

				MatrixCPU Pn = res.getPn();

				MatrixCPU ResF = res.getRes();
				int iter = res.getIter();
				//std::cout << "enregistrement de " << Pn.get(1, 0) << " en pose " << k << " methode " << indices[i] << std::endl;
				PResult.set(k, j, Pn.get(1, 0));
				QResult.set(k, j, Pn.get(agents + 2, 0));
				temps.set(k, j, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / BILLION);
				iters.set(k, j, iter);
				ResR.set(k, j, ResF.get(0, (iter - 1) / stepG));
				ResS.set(k, j, ResF.get(1, (iter - 1) / stepG));
				ResV.set(k, j, ResF.get(2, (iter - 1) / stepG));

				Fc.set(k, j, res.getFc());
				
				//CHECK_LAST_CUDA_ERROR();
				//sys.resetParam();
				//}
			}
		}
		std::cout << std::endl;
		//std::cout << "-|-|-|-|-|-|-|-|-  -|-|-|-|-|-|-|-|-|" << std::endl;

	}


	if (save) {
		Buses.saveCSV(fileName, mode);

		PResult.saveCSV(fileName, mode);
		QResult.saveCSV(fileName, mode);
		temps.saveCSV(fileName, mode);
		iters.saveCSV(fileName, mode);
		ResR.saveCSV(fileName, mode);
		ResS.saveCSV(fileName, mode);
		ResV.saveCSV(fileName, mode);
		Fc.saveCSV(fileName, mode);
	}


	temps.display();
	iters.display();
	for (int i = 0; i < nMethode; i++) {
		DELETEB(methodes[i]);
	}

	
	//sys.setMethod(nullptr);
}


void SimuStatMarketEndoACAgent() {
	std::string fileName = "ComparaisonMarketEndoAllBus100_Agent600.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;

	const int nMethode = 6;
	std::vector<int> indices = { 0, 1, 2, 3, 4, 5 };
	std::string methodesName[nMethode] = { "EndoDirect", "EndoConsensus", "AC EndoPF", "EndoDirectGPU", "EndoConsensusGPU", "AC EndoPFGPU" };
	Method* methodes[nMethode];
	methodes[0] = new MarketEndoDirect;
	methodes[1] = new MarEndoCons;
	methodes[2] = new EndoPF;
	methodes[3] = new MarketEndoDirectGPU;
	methodes[4] = new MarEndoConsGPU;
	methodes[5] = new EndoPFGPU;
	bool save = true;



	if (nMethode != indices.size()) {
		throw std::domain_error("not good number of methods");
	}
	
	//float factorAgent = 1;
	int nBus = 100;


	// cases
	int nNAgent = 1;
	int nAgentMax = 600;
	float factorPower = ((float) nBus / nAgentMax);


	// Market
	float Pconso = 0.05 * factorPower;
	float dPconso = 0.01 * factorPower;
	float Propcons = 0.5;
	float PropNFleGen = 0.25;
	float PropGen = 1 - PropNFleGen - Propcons;
	float bProd = 1;
	float dbProd = 0.1;
	float Pprod = 0.01 * factorPower;
	float dPprod = 0.005 * factorPower;
	float gamma = 1; // 
	float dgamma = 0.2; //  
	float dQ = 0.01 * factorPower;
	float length = 0.001;
	float dlength = 0.0005;



	const int offsetAgent = 0;
	int nSimu = 50;
	int nCasAgent = 1;
	int nCasBuses = 1;
	//int nCas = nCasAgent * nCasBuses;


	// simulation
	int iterGlobal = 1000;
	int iterLocal = 3000;
	int stepG = 10;
	int stepL = 1;
	float epsG = 0.1f;
	float epsL = 0.01f;
	float rho = 50;
	float rho1 = rho / 2;
	float epsGC = 1;




	MatrixCPU Param(1, 22);
	Param.set(0, 0, nAgentMax);
	Param.set(0, 1, nNAgent);
	Param.set(0, 2, length);
	Param.set(0, 3, dlength);
	Param.set(0, 4, Pconso);
	Param.set(0, 5, dPconso);
	Param.set(0, 6, Pprod);
	Param.set(0, 7, dPprod);
	Param.set(0, 8, dQ);
	Param.set(0, 9, Propcons);
	Param.set(0, 10, PropGen);
	Param.set(0, 11, nSimu);
	Param.set(0, 12, nCasAgent);
	Param.set(0, 13, nCasBuses);
	Param.set(0, 14, nMethode);
	Param.set(0, 15, rho);
	Param.set(0, 16, epsG);
	Param.set(0, 17, epsL);
	Param.set(0, 18, iterGlobal);
	Param.set(0, 19, iterLocal);
	Param.set(0, 20, stepG);
	Param.set(0, 21, stepL);


	if (save) {
		Param.saveCSV(fileName, mode);
	}



	MatrixCPU Agents(1, nNAgent);
	MatrixCPU temps(nMethode * nSimu, nNAgent, nanf(""));
	MatrixCPU iters(nMethode * nSimu, nNAgent, nanf(""));
	MatrixCPU PResult(nMethode * nSimu, nNAgent, nanf(""));

	MatrixCPU QResult(nMethode * nSimu, nNAgent, nanf(""));
	MatrixCPU ResR(nMethode * nSimu, nNAgent, nanf(""));
	MatrixCPU ResS(nMethode * nSimu, nNAgent, nanf(""));
	MatrixCPU ResV(nMethode * nSimu, nNAgent, nanf(""));
	MatrixCPU Fc(nMethode * nSimu, nNAgent, nanf(""));
	//System sys;

	/*sys.setIter(iterGlobal, iterLocal);
	sys.setEpsG(epsG);
	sys.setEpsL(epsL);
	sys.setStep(stepG, stepL);
	sys.setRho(rho);*/

	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;

	for (int agent = 0; agent < nNAgent; agent++) {
		std::cout << "--------- --------- --------- --------- ----------" << std::endl;
		int agents = offsetAgent > 0 ? agent * (nAgentMax - offsetAgent) / MYMAX(nNAgent - 1, 1) + offsetAgent : (agent + 1) * nAgentMax / nNAgent;
		Agents.set(0, agent, agents);

		int nDeep = 2 * nBus;
		int nBranch = nBus;
		int j = agent;

		
		//std::cout << " cas avec " << buses << " bus et " << agents << "agents et une profondeur de " << nDeep << " et largeur de " << nBranch << std::endl;

		for (int simu = 0; simu < nSimu; simu++) {


			std::cout << "-";



			StudyCase cas;
			cas.genGridBT(nBus, nBranch, nDeep, length, dlength);
			cas.genAgentsAC(agents, Propcons, PropGen, Pconso, dPconso, bProd, dbProd, dQ, Pprod, dPprod, gamma, dgamma);


			cas.genLinkGridAgent();
			
			Simparam param(cas.getNagent(), cas.getNLine(true), true);
			param.setEpsL(epsL);
			param.setEpsG(epsG);//0.0001f
			param.setItG(iterGlobal); //500000
			param.setItL(iterLocal);
			param.setStep(stepG, stepL);
			param.setRho(rho);
			param.setRho1(rho1);
			param.setEpsGC(epsGC);

			Simparam res(param);
			//sys.setStudyCase(cas);
			std::random_shuffle(indices.begin(), indices.end());
			for (int i = 0; i < nMethode; i++) {
				int k = indices[i] * nSimu + simu;
				//std::cout << " Simu  " << simu << " methode  " << i << " : " << methodesName[indices[i]] << std::endl;
				//if (indices[i] < 3) {

				//sys.setMethod(methodes[indices[i]]);
				t1 = std::chrono::high_resolution_clock::now();
				methodes[indices[i]]->solve(&res, param, cas);
				t2 = std::chrono::high_resolution_clock::now();

				MatrixCPU Pn = res.getPn();

				MatrixCPU ResF = res.getRes();
				int iter = res.getIter();
				//std::cout << "enregistrement de " << Pn.get(1, 0) << " en pose " << k << " methode " << indices[i] << std::endl;
				PResult.set(k, j, Pn.get(1, 0));
				QResult.set(k, j, Pn.get(agents + 2, 0));
				temps.set(k, j, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / BILLION);
				iters.set(k, j, iter);
				ResR.set(k, j, ResF.get(0, (iter - 1) / stepG));
				ResS.set(k, j, ResF.get(1, (iter - 1) / stepG));
				ResV.set(k, j, ResF.get(2, (iter - 1) / stepG));

				Fc.set(k, j, res.getFc());

				//CHECK_LAST_CUDA_ERROR();
				//sys.resetParam();
				//}
			}
		}
		std::cout << std::endl;
		//std::cout << "-|-|-|-|-|-|-|-|-  -|-|-|-|-|-|-|-|-|" << std::endl;

	}


	if (save) {
		Agents.saveCSV(fileName, mode);

		PResult.saveCSV(fileName, mode);
		QResult.saveCSV(fileName, mode);
		temps.saveCSV(fileName, mode);
		iters.saveCSV(fileName, mode);
		ResR.saveCSV(fileName, mode);
		ResS.saveCSV(fileName, mode);
		ResV.saveCSV(fileName, mode);
		Fc.saveCSV(fileName, mode);
	}


	temps.display();
	iters.display();
	for (int i = 0; i < nMethode; i++) {
		DELETEB(methodes[i]);
	}


	//sys.setMethod(nullptr);
}

void SimuTemporalTestFeederEndo() {
	std::string path = "data/ACGrid/";
	std::string fileName2 = "SimutemporalFBTestFeederShort.csv";
	std::string fileName = "SimutemporalTestFeederEndo.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nMethode = 4;
	std::string methodesName[nMethode] = { "EndoDirect", "EndoConsensus", "EndoDirect GPU", "EndoConsensus GPU"};
	Method* methodes[nMethode];
	methodes[0] = new MarketEndoDirect;
	methodes[1] = new MarEndoCons;
	methodes[2] = new MarketEndoDirectGPU;
	methodes[3] = new MarEndoConsGPU;
	bool saveTime = true; //true  false


	float epsG = 0.05f;
	float epsGC = 0.005f;
	float epsL = 0.0005f;

	int iterG = 20000;
	int iterL = 5000;//500;

	int chosenAgenGen = 0;
	int nAgent = 57;
	if (chosenAgenGen != 0) {
		nAgent = 57; // may change
	}

	int begin = 0;
	int end = 30 - 1; //60 * 2 - 1; 30 - 1

	float rho = 10;


	float stepG = 1;
	float stepL = 1;

	

	System sys;
	sys.setIter(iterG, iterL);
	sys.setStep(stepG, stepL);
	sys.setEpsG(epsG);
	sys.setEpsGC(epsG);
	sys.setEpsL(epsL);

	MatrixCPU Param(1, 12);
	Param.set(0, 0, nAgent);
	Param.set(0, 1, epsGC);
	Param.set(0, 2, epsG);
	Param.set(0, 3, epsL);
	Param.set(0, 4, iterG);
	Param.set(0, 5, iterL);
	Param.set(0, 6, stepG);
	Param.set(0, 7, stepL);
	Param.set(0, 8, rho);
	Param.set(0, 9, chosenAgenGen);
	Param.set(0, 10, begin);
	Param.set(0, 11, end);
	

	if (saveTime) {
		//Param.saveCSV(fileName, mode);
	}
#ifdef INSTRUMENTATION
	//Param.saveCSV(fileName2, mode);
#endif // INSTRUMENTAION



	for (int i = 2; i < 3; i++) {//nMethode

		sys.setMethod(methodes[i]);
		sys.setRho(rho);
		std::cout << "----------------Simu---------------------------- " << std::endl;
		// Debut simu
		std::cout << "methode " << methodesName[i] << std::endl;
		clock_t t = clock();

		sys.solveIntervalle(path, begin, end, chosenAgenGen);
		t = clock() - t;

		std::cout << "calculation time : " << (float)t / CLOCKS_PER_SEC << std::endl;


		MatrixCPU temps(sys.getTemps());
		MatrixCPU iter(sys.getIter());
		MatrixCPU conv(sys.getConv());
		MatrixCPU fc(sys.getFc());
		MatrixCPU ResR(sys.getResR());
		MatrixCPU ResS(sys.getResS());
		MatrixCPU ResX(sys.getResX());

		temps.display();
		ResR.display();
		ResS.display();
		ResX.display();
		iter.display();

		if (saveTime) {
			temps.saveCSV(fileName, mode);
			iter.saveCSV(fileName, mode);
			fc.saveCSV(fileName, mode);
			ResR.saveCSV(fileName, mode);
			ResS.saveCSV(fileName, mode);
			ResX.saveCSV(fileName, mode);
			conv.saveCSV(fileName, mode);/**/
		}


		std::cout << "-------------------------------------------------------- " << std::endl;
#ifdef INSTRUMENTATION
		sys.displayTime(fileName2);
#endif // INSTRUMENTATION


		sys.resetParam();
	}
	for (int i = 0; i < nMethode; i++) {
		DELETEB(methodes[i]);
	}
	sys.setMethod(nullptr);

}

void SimuTemporalTestFeederEndoAll()  {
	std::string path = "data/ACGrid/";
	std::string fileName2 = "SimutemporalFBTestFeederShort.csv";
	std::string fileName = "SimutemporalTestFeederEndo.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nMethode = 8;
	std::string methodesName[nMethode] = { "EndoDirect", "EndoConsensus", "Ac EndoPF", "Dc EndoPF", "EndoDirect GPU", "EndoConsensus GPU", "AC EndoPF GPU", "DC EndoPF GPU" };
	Method* methodes[nMethode];
	methodes[0] = new MarketEndoDirect;
	methodes[1] = new MarEndoCons;
	methodes[2] = new EndoPF;
	methodes[3] = new ADMMConst;
	methodes[4] = new MarketEndoDirectGPU;
	methodes[5] = new MarEndoConsGPU;
	methodes[6] = new EndoPFGPU;
	methodes[7] = new ADMMGPUConst4;
	bool saveTime = false; //true  false


	float epsG = 0.05f;
	float epsGC = 0.005f;
	float epsL = 0.0005f;

	int iterG = 1000;
	int iterL = 5000;//500;

	int chosenAgenGen = 0;
	int nAgent = 57;
	if (chosenAgenGen != 0) {
		nAgent = 57; // may change
	}

	int begin = 0;
	int end = 5 - 1; //60 * 2 - 1; 30 - 1

	float rho = 10;


	float stepG = 1;
	float stepL = 1;



	System sys;
	sys.setIter(iterG, iterL);
	sys.setStep(stepG, stepL);
	sys.setEpsG(epsG);
	sys.setEpsGC(epsG);
	sys.setEpsL(epsL);

	MatrixCPU Param(1, 12);
	Param.set(0, 0, nAgent);
	Param.set(0, 1, epsGC);
	Param.set(0, 2, epsG);
	Param.set(0, 3, epsL);
	Param.set(0, 4, iterG);
	Param.set(0, 5, iterL);
	Param.set(0, 6, stepG);
	Param.set(0, 7, stepL);
	Param.set(0, 8, rho);
	Param.set(0, 9, chosenAgenGen);
	Param.set(0, 10, begin);
	Param.set(0, 11, end);


	if (saveTime) {
		//Param.saveCSV(fileName, mode);
	}
#ifdef INSTRUMENTATION
	Param.saveCSV(fileName2, mode);
#endif // INSTRUMENTAION



	for (int i = 6; i < 7; i++) {//nMethode

		sys.setMethod(methodes[i]);
		sys.setRho(rho);
		std::cout << "----------------Simu---------------------------- " << std::endl;
		// Debut simu
		std::cout << "methode " << methodesName[i] << std::endl;
		clock_t t = clock();

		sys.solveIntervalle(path, begin, end, chosenAgenGen);
		t = clock() - t;

		std::cout << "calculation time : " << (float)t / CLOCKS_PER_SEC << std::endl;


		MatrixCPU temps(sys.getTemps());
		MatrixCPU iter(sys.getIter());
		MatrixCPU conv(sys.getConv());
		MatrixCPU fc(sys.getFc());
		MatrixCPU ResR(sys.getResR());
		MatrixCPU ResS(sys.getResS());
		MatrixCPU ResX(sys.getResX());

		temps.display();
		ResR.display();
		ResS.display();
		ResX.display();
		iter.display();

		if (saveTime) {
			temps.saveCSV(fileName, mode);
			iter.saveCSV(fileName, mode);
			fc.saveCSV(fileName, mode);
			ResR.saveCSV(fileName, mode);
			ResS.saveCSV(fileName, mode);
			ResX.saveCSV(fileName, mode);
			conv.saveCSV(fileName, mode);/**/
		}


		std::cout << "-------------------------------------------------------- " << std::endl;
#ifdef INSTRUMENTATION
		sys.displayTime(fileName2);
#endif // INSTRUMENTATION


		sys.resetParam();
	}
	for (int i = 0; i < nMethode; i++) {
		DELETEB(methodes[i]);
	}
	sys.setMethod(nullptr);

}


void SimuStatMarketEndoGrid() {
	std::string fileName = "ComparaisonMarketEndoGrid_100.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nMethod = 4;
	std::vector<int> indices = { 0, 1, 2, 3 };
	std::string methodesName[nMethod] = { "EndoDirect", "EndoConsensus", "EndoDirectGPU", "EndoConsensusGPU" };
	Method* methodes[nMethod];
	methodes[0] = new MarketEndoDirect;
	methodes[1] = new MarEndoCons;
	methodes[2] = new MarketEndoDirectGPU;
	methodes[3] = new MarEndoConsGPU;


	// Market
	float Pconso = 0.05;
	float dPconso = 0.01;
	float Propcons = 0.5;
	float PropNFleGen = 0.25;
	float PropGen = 1 - PropNFleGen - Propcons;
	float bProd = 1;
	float dbProd = 0.1;
	float Pprod = 0.01;
	float dPprod = 0.005;
	float gamma = 1; // 
	float dgamma = 0.2; //  
	float dQ = 0.01;
	float length = 0.001;
	float dlength = 0.0005;


	// cases
	int nNBus = 1;
	int nBusMax = 100;
	const int offsetBus = 0;
	int nSimu = 50;
	//int million = 1000000;
	int nCasAgent = 4;
	int nCasBuses = 4;
	int nCas = nCasAgent * nCasBuses;

	// simulation
	int iterGlobal = 1000;
	int iterLocal = 1000;
	int stepG = 1;
	int stepL = 1;
	float epsG = 0.1f;
	float epsL = 0.01f;
	float rho = 50;
	float rho1 = rho / 2;
	float epsGC = 1;
	
	

	float factorAgent[] = { 0.5, 1, 2, 5 };
	//float factorBuses[] = { 2, 0.3, 0.3 , 1 }; Normal, Line, Balance, OneStep





	MatrixCPU Param(1, 22);
	Param.set(0, 0, nBusMax);
	Param.set(0, 1, nNBus);
	Param.set(0, 2, length);
	Param.set(0, 3, dlength);
	Param.set(0, 4, Pconso);
	Param.set(0, 5, dPconso);
	Param.set(0, 6, Pprod);
	Param.set(0, 7, dPprod);
	Param.set(0, 8, dQ);
	Param.set(0, 9, Propcons);
	Param.set(0, 10, PropGen);
	Param.set(0, 11, nSimu);
	Param.set(0, 12, nCasAgent);
	Param.set(0, 13, nCasBuses);
	Param.set(0, 14, nMethod);
	Param.set(0, 15, rho);
	Param.set(0, 16, epsG);
	Param.set(0, 17, epsL);
	Param.set(0, 18, iterGlobal);				 
	Param.set(0, 19, iterLocal);
	Param.set(0, 20, stepG);
	Param.set(0, 21, stepL);


	Param.saveCSV(fileName, mode);

	MatrixCPU Buses(1, nNBus);

	MatrixCPU temps(nCas * nMethod, nNBus * nSimu, nanf(""));
	MatrixCPU iters(nCas * nMethod, nNBus * nSimu, nanf(""));
	MatrixCPU PResult(nCas * nMethod, nNBus * nSimu, nanf(""));
	MatrixCPU QResult(nCas * nMethod, nNBus * nSimu, nanf(""));
	MatrixCPU ResR(nCas * nMethod, nNBus * nSimu, nanf(""));
	MatrixCPU ResS(nCas * nMethod, nNBus * nSimu, nanf(""));
	MatrixCPU ResV(nCas * nMethod, nNBus * nSimu, nanf(""));
	MatrixCPU Fc(nCas * nMethod, nNBus * nSimu, nanf(""));

	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;


	

	for (int bus = 0; bus < nNBus; bus++) {
		std::cout << "--------- --------- --------- --------- ----------" << std::endl;
		int buses = offsetBus > 0 ? bus * (nBusMax - offsetBus) / MYMAX(nNBus - 1, 1) + offsetBus : (bus + 1) * nBusMax / nNBus;
		Buses.set(0, bus, buses);
		for (int casBus = 0; casBus < nCasBuses; casBus++) {
			int nBranch = buses;
			int nDeep = buses;
			
			for (int casAgent = 0; casAgent < nCasAgent; casAgent++) {
				int agents = factorAgent[casAgent] * buses;
				float factorPower = ((float) buses / agents);
				for (int simu = 0; simu < nSimu; simu++) {

					int i = (casBus * nCasAgent + casAgent) * nMethod;
					int j = bus * nSimu + simu;

					StudyCase cas;
					std::cout << "-";
					cas.genGridBTSpecial(buses, nBranch, nDeep, length, dlength, (RadialType) casBus);
					cas.genAgentsAC(agents, Propcons, PropGen, Pconso * factorPower, dPconso * factorPower, bProd, dbProd, dQ * factorPower, Pprod * factorPower, dPprod * factorPower, gamma, dgamma);
					cas.genLinkGridAgent();

					
				
					
					Simparam param(cas.getNagent(), cas.getNLine(true), true);
					param.setEpsL(epsL);
					param.setEpsG(epsG);//0.0001f
					param.setItG(iterGlobal); //500000
					param.setItL(iterLocal);
					param.setStep(stepG, stepL);
					param.setRho(rho);
					param.setRho1(rho1);
					param.setEpsGC(epsGC);

					Simparam res(param);
					//sys.setStudyCase(cas);
					std::random_shuffle(indices.begin(), indices.end());


					for (int method = 0; method < nMethod; method++) {
						int k = i + indices[method];
				
						t1 = std::chrono::high_resolution_clock::now();
						methodes[indices[method]]->solve(&res, param, cas);
						t2 = std::chrono::high_resolution_clock::now();

						MatrixCPU Pn = res.getPn();

						MatrixCPU ResF = res.getRes();
						int iter = res.getIter();
						//std::cout << "enregistrement de " << Pn.get(1, 0) << " en pose " << k << " " << j << " methode " << indices[method] << std::endl;
						PResult.set(k, j, Pn.get(1, 0));
						QResult.set(k, j, Pn.get(agents + 2, 0));
						temps.set(k, j, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / BILLION);
						iters.set(k, j, iter);
						ResR.set(k, j, ResF.get(0, (iter - 1) / stepG));
						ResS.set(k, j, ResF.get(1, (iter - 1) / stepG));
						ResV.set(k, j, ResF.get(2, (iter - 1) / stepG));

						Fc.set(k, j, res.getFc());

						//CHECK_LAST_CUDA_ERROR();
						//sys.resetParam();
						//}
					}

				}
			}
			std::cout << std::endl;
		}
		//std::cout << "-|-|-|-|-|-|-|-|-  -|-|-|-|-|-|-|-|-|" << std::endl;

	}
	/*Agents.display();
	temps.display();
	iters.display();*/
	
	float temptotal = temps.sum();
	std::cout << "temps total " << temptotal << " temps moyen " << temptotal / (nSimu * nNBus * nCas * nMethod) << std::endl;

	Buses.saveCSV(fileName, mode);

	PResult.saveCSV(fileName, mode);
	QResult.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
	iters.saveCSV(fileName, mode);
	ResR.saveCSV(fileName, mode);
	ResS.saveCSV(fileName, mode);
	ResV.saveCSV(fileName, mode);
	Fc.saveCSV(fileName, mode);

}


void SimuStatMarketEndoArticle(){
	std::string fileName = "ComparaisonMarketEndoAgentArticle.csv";
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nMethod = 4;
	std::vector<int> indices = { 0, 1, 2, 3 };
	std::string methodesName[nMethod] = { "EndoDirect", "EndoConsensus", "EndoDirectGPU", "EndoConsensusGPU" };
	Method* methodes[nMethod];
	methodes[0] = new MarketEndoDirect;
	methodes[1] = new MarEndoCons;
	methodes[2] = new MarketEndoDirectGPU;
	methodes[3] = new MarEndoConsGPU;


	// Market
	float Pconso = 0.5f;
	float dPconso = 0.1f;
	float Propcons = 0.5f;
	float PropNFleGen = 0.25f;
	float PropGen = 1 - PropNFleGen - Propcons;
	float bProd = 1;
	float dbProd = 0.1f;
	float Pprod = 2;
	float dPprod = 0.1f;
	float gamma = 2; // 
	float dgamma = 1; //  
	float dQ = 0.005;
	


	// cases
	int nSimu = 50;
	const int nNAgent = 6;
	int tabNagent[nNAgent] = {10, 20, 50, 100, 120, 150};
	//int million = 1000000;


	// simulation
	int stepG = 1;
	int stepL = 1;
	int stepIntern = 1;

	if(stepG < stepIntern){
		stepG = stepIntern;
	}

	int iterL = 5000;
	int iterG  = 10000;
	int iterIntern = 5000;

	float epsL = 0.0001f;
	float epsG = 0.001f;
	float epsGC = 0.0005f;
	float epsIntern = 0.001f;

	float rhoInit = 1; // 1 pour cas 2 noeuds, 5 pour cas9, cas 10, 10 cas 69

	MatrixCPU Param(1, 22);
	/*Param.set(0, 0, nBusMax);
	Param.set(0, 1, nNBus);
	Param.set(0, 2, length);
	Param.set(0, 3, dlength);*/
	Param.set(0, 4, Pconso);
	Param.set(0, 5, dPconso);
	Param.set(0, 6, Pprod);
	Param.set(0, 7, dPprod);
	Param.set(0, 8, dQ);
	Param.set(0, 9, Propcons);
	Param.set(0, 10, PropGen);
	Param.set(0, 11, nSimu);
	//Param.set(0, 12, nCasAgent);
	//Param.set(0, 13, nCasBuses);
	Param.set(0, 14, nMethod);
	Param.set(0, 15, rhoInit);
	Param.set(0, 16, epsG);
	Param.set(0, 17, epsL);
	Param.set(0, 18, iterG);				 
	Param.set(0, 19, iterL);
	Param.set(0, 20, stepG);
	Param.set(0, 21, stepL);


	Param.saveCSV(fileName, mode);

	MatrixCPU Agents(1, nNAgent);
	MatrixCPU temps(nMethod * nSimu, nNAgent, -1);
	MatrixCPU iters(nMethod * nSimu, nNAgent, -1);
	
	StudyCase cas;
	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;
	cas.SetACFromFile("case85");
	for (int agent = 0; agent < 2; agent++) {
		std::cout << " N agent : " << tabNagent[agent] << std::endl;
		std::cout << "--------- --------- --------- --------- ----------" << std::endl;
		int agents = tabNagent[agent];
		//Agents.set(0, agent, agents);
		for (int simu = 0; simu < nSimu; simu++) {

			int i = simu * nMethod;
			int j = agent;

			
			std::cout << "-";
			
			cas.genAgentsAC(agents, Propcons, PropNFleGen, Pconso, dPconso, bProd, dbProd, dQ, Pprod, dPprod, gamma, dgamma);
			cas.genLinkGridAgent();
			
			Simparam param(cas.getNagent(), cas.getNLine(true), true);
			param.setEpsL(epsL);
			param.setEpsG(epsG);//0.005f FB
			param.setEpsGC(epsGC); //0.001f FB
			param.setEpsIntern(epsIntern);
			param.setItG(iterG); //20000 500000
			param.setItL(iterL);
			param.setItIntern(iterIntern);
			param.setStep(stepG, stepL, stepIntern);
			
			//for (int i = 0; i < 5; i++) {
			param.setRho(rhoInit);
			Simparam res(param);
			//sys.setStudyCase(cas);
			std::random_shuffle(indices.begin(), indices.end());


			for (int method = 0; method < nMethod; method++) {
				int k = i + indices[method];
		
				t1 = std::chrono::high_resolution_clock::now();
				methodes[indices[method]]->solve(&res, param, cas);
				t2 = std::chrono::high_resolution_clock::now();

				int iter = res.getIter();
						
				temps.set(k, j, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / BILLION);
				iters.set(k, j, iter);
						
			}
		}
		std::cout << std::endl;
		//std::cout << "-|-|-|-|-|-|-|-|-  -|-|-|-|-|-|-|-|-|" << std::endl;

	}
	/*Agents.display();
	temps.display();
	iters.display();*/
	
	float temptotal = temps.sum();
	
	Agents.saveCSV(fileName, mode);

	
	temps.saveCSV(fileName, mode);
	iters.saveCSV(fileName, mode);
	
}




/* Test fonctionnel */

void testCPUPF()
{
	// donc fonctionne pour 3,9, 14, 30, 57 : m�me resultat
	// fonctionne pour 85 et c'est un reseau de distribution ! (mais pas de prod dedans ...)
	// Fonctionne pour 69
	// ne fonctionne pas pour 118
	// 

	// ne fonctione pas en float pour 141, le noeud 85 pose probleme, a des imp�dance 1000fois plus grande que les autres (probleme de pr�cision num�rique)
	// en gros, m�me en partant des m�me tensions, on ne trouve pas la m�me puissance r�active � un noeud (85) 
	// Poutant je calcule en double ! Peut �tre s�parer les calculs associ� � G et B (comme les G ont les m�mes ordre de grandeur entre eux)
	// Il faut stocker G et B, et E en double pour que cela marche !!! (oulala le passage sur GPU va �tre pas tr�s efficace).
	// 
	// 
	// diverge compeltement pour 300 ou 145  !!!!

	// Le cas 39 converge vers la bonne solution si proche de la solution, une autre solution si on est plus loin (il faut passer Seidel en double !)
	// et matlab est d'accord si on passe tous les noeuds en PQ !
	// idem pour le cas ieee30, qui converge vers une autre solution (rem : il y a diff�rents niveau de tension avec transfo dans ce cas !)

	bool methodeToCompute[] = { true, true, true, true, true, true, true, true, true, true, true };
	int nMethode = 11; //bool setDouble = true;
	CPUPF PF;
	CPUPFGS PFGS;

	GPUPF PFG;
	GPUPFGS PFGSG;

	CPUPF PF2;
	CPUPFGS PFGS2;

	GPUPF PFG2;
	GPUPFGS PFGSG2;
	
	CPUPFdist PF3;
	CPUPFdistPQ PF3PQ;
	GPUPFdistPQ PF3PQG;
	std::cout << "--------------------------------------------------------" << std::endl;
	bool setE = true; // not true
	bool setSol = false;
	bool setSolve = true;
	bool displayAll = true;
	bool setError = false;
	bool save = false;
	
		

	// 
	// 136ma n'a pas les bonnes impedances ...
	StudyCase cas2;
	int choseCase = 0;
	std::string fileName = "TimeByBlockPF";
	std::string chosenCase = "";
	float Power = 0;
	switch (choseCase)
	{
	case 0:
		chosenCase = "case10ba";// 9, 30 57   118 300 | radial: case10ba case4_dist  case69 case85
		cas2.SetACFromFile(chosenCase); //case_ACTIVSg2000	
		cas2.display();
		break;
	case 1:
		chosenCase = "RandRadial";
		cas2.genGridBT(100, 80, 200, 0.001, 0.0005);
		cas2.genAgentsAC(50, 0.5, 0.5, 1, 0.5, 1, 0.1, 0.5, 1, 0.5, 1, 0.2);
		cas2.genLinkGridAgent();
		break;
	case 2:
		chosenCase = "RandHTB";

		cas2.genGridHTB(10, 20, 1, 0.01, 0.0005);
		cas2.genAgentsAC(16, 0.5, 0.5, Power, 0, 1, 0, 0.5, Power, 0, 1, 0.2);
		cas2.genLinkGridAgent();
		break;
	case 3:
		chosenCase = "EuropeTestFeeder";
		cas2.SetEuropeTestFeeder();
		cas2.display();
		break;
	case 4:
		chosenCase = "2node";
		cas2.SetAC2node();
		cas2.display();
		break;
	case 5:
		chosenCase = "3node";
		cas2.SetAC3Bus();
		cas2.display();
		break;
	default:
		std::cout << "unknown choice, case 0" << std::endl;
		break;
	}



	fileName += chosenCase + ".csv";

	int million = 1000000;
	
	MatrixCPU results(4, nMethode, nanf(""));
	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;

	
#ifdef INSTRUMENTATION
	save = false;
#endif // INSTRUMENTATION

	double epsi = 0;
	double err = 0.1;
	//cas2.display();	
	//std::cout << " *** setPQ *** " << std::endl;
	MatrixCPU PQ2 = cas2.getPobj();
	MatrixCPUD PQ2D = cas2.getPobjD();


	if (choseCase == 4) {
		PQ2.set(1, 0, -1);
		PQ2.set(2, 0, 1);
		PQ2D.set(1, 0, -1);
		PQ2D.set(2, 0, 1);
		PQ2.set(4, 0, -0.9);
		PQ2.set(5, 0, 0.9);
		PQ2D.set(4, 0, -0.9);
		PQ2D.set(5, 0, 0.9);
	}

	//PQ2.display();

	MatrixGPU PQ2G = MatrixGPU(PQ2, 1);
	MatrixGPUD PQ2GD = MatrixGPUD(PQ2D, 1);
	MatrixCPUD sol = cas2.getSolPF();
	MatrixCPUD voltInit = cas2.getVoltageInitD();

	int nBus = cas2.getNBus();
	MatrixCPU E(2 * nBus, 1);
	MatrixCPUD ED(2 * nBus, 1);
	int method = 0;
	MatrixCPUD W(2 * nBus, 1);
	
	if (setE) {
		//std::cout << " *** setE *** " << std::endl;
		for (int i = 0; i < nBus; i++) {
			if (setSol) {
				E.set(i + nBus, 0, sol.get(i, 0) + epsi);
				E.set(i, 0, sol.get(i, 1) + epsi);
				ED.set(i + nBus, 0, sol.get(i, 0) + epsi);
				ED.set(i, 0, sol.get(i, 1) + epsi);
				W.set(i, 0, sol.get(i, 2));
				W.set(i + nBus, 0, sol.get(i, 3));
			}
			else {
				E.set(i + nBus, 0, voltInit.get(i + nBus, 0) + epsi);
				E.set(i, 0, voltInit.get(i, 0) + epsi);
				ED.set(i + nBus, 0, voltInit.get(i + nBus, 0) + epsi);
				ED.set(i, 0, voltInit.get(i, 0) + epsi);
			}
			if (setError) {
				epsi = err;
			}
		}
		//E.display();
	}

	if (nBus > 150) {
		displayAll = false;
	}
	MatrixGPU EG = MatrixGPU(E, 1);
	MatrixGPUD EGD = MatrixGPUD(ED, 1);



	std::cout << "******************************* Newton Raphson CPU ************************************" << std::endl;
	if (methodeToCompute[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		PF.init(cas2, &PQ2);
		if (setE) {
			PF.setE(&E);
		}
		if (setSolve) {
			PF.solve();
		}
		else {
			PF.calcW(true);
		}
		if (save) {
			PF.saveTimeBlock(fileName);
		}
		t2 = std::chrono::high_resolution_clock::now();
		PF.display(displayAll);
		results.set(0, method, PF.getP0());
		results.set(1, method, PF.getQ0());
		results.set(2, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(3, method, PF.getIter());
		
		
	}
	method++;
	if (methodeToCompute[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		PF2.init(cas2, &PQ2, &PQ2D, true);
		if (setE) {
			PF2.setE(&ED);
		}
		if (setSolve) {
			PF2.solve();
		}
		else {
			PF2.calcW(true);
		}
		if (save) {
			PF2.saveTimeBlock(fileName);
		}

		t2 = std::chrono::high_resolution_clock::now();
		PF2.display(displayAll);
		results.set(0, method, PF2.getP0());
		results.set(1, method, PF2.getQ0());
		results.set(2, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(3, method, PF2.getIter());
		/**/
	}
	method++;

	std::cout << "*****************************Newton Raphson GPU*******************************" << std::endl;
	if (methodeToCompute[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		PFG.init(cas2, &PQ2G);
		if (setE) {
			PFG.setE(&EG);
		}
		if (setSolve) {
			PFG.solve();
		}
		else {
			PFG.calcW(true);
		}
		if (save) {
			PFG.saveTimeBlock(fileName);
		}
		t2 = std::chrono::high_resolution_clock::now();
		PFG.display(displayAll);
		results.set(0, method, PFG.getP0());
		results.set(1, method, PFG.getQ0());
		results.set(2, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(3, method, PFG.getIter());
	}
	method++;

	if (methodeToCompute[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		PFG2.init(cas2, &PQ2G, &PQ2GD, true);
		if (setE) {
			PFG2.setE(&EGD);
		}
		if (setSolve) {
			PFG2.solve();
		}
		else {
			PFG2.calcW(true);
		}
		if (save) {
			PFG2.saveTimeBlock(fileName);
		}

		t2 = std::chrono::high_resolution_clock::now();
		PFG2.display(displayAll);
		results.set(0, method, PFG2.getP0());
		results.set(1, method, PFG2.getQ0());
		results.set(2, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(3, method, PFG2.getIter());
		
	}
	method++;/**/


	std::cout << "*******************************Gauss Seidel CPU*********************************" << std::endl;
	if (methodeToCompute[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		PFGS.init(cas2, &PQ2);
		if (setE) {
			PFGS.setE(&E);
		}
		if (setSolve) {
			PFGS.solve();
		}
		else {
			PFGS.calcW(true);
		}
		if (save) {
			PFGS.saveTimeBlock(fileName);
		}

		PFGS.calcE();
		t2 = std::chrono::high_resolution_clock::now();
		PFGS.display(displayAll);
		results.set(0, method, PFGS.getP0());
		results.set(1, method, PFGS.getQ0());
		results.set(2, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(3, method, PFGS.getIter());
	}
	method++;
	if (methodeToCompute[method]) {
		
		t1 = std::chrono::high_resolution_clock::now();
		PFGS2.init(cas2, &PQ2, &PQ2D, true);
		if (setE) {
			PFGS2.setE(&ED);
		}
		if (setSolve) {
			PFGS2.solve();
		}
		else {
			PFGS2.calcW(true);
		}
		if (save) {
			PFGS2.saveTimeBlock(fileName);
		}
		PFGS2.calcE();
		t2 = std::chrono::high_resolution_clock::now();
		PFGS2.display(displayAll);
		results.set(0, method, PFGS2.getP0());
		results.set(1, method, PFGS2.getQ0());
		results.set(2, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(3, method, PFGS2.getIter());
			
		
	}
	method++;/**/
	std::cout << "*******************************Gauss Seidel GPU**************************************" << std::endl;
	
	if (methodeToCompute[method]) {

		t1 = std::chrono::high_resolution_clock::now();
		PFGSG.init(cas2, &PQ2G);
		if (setE) {
			PFGSG.setE(&EG);
		}
		if (setSolve) {
			PFGSG.solve();
		}
		else {
			PFGSG.calcW(true);
		}
		if (save) {
			PFGSG.saveTimeBlock(fileName);
		}

		PFGSG.calcE();
		t2 = std::chrono::high_resolution_clock::now();
		PFGSG.display(displayAll);
		results.set(0, method, PFGSG.getP0());
		results.set(1, method, PFGSG.getQ0());
		results.set(2, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(3, method, PFGSG.getIter());
	}
	method++;/**/
	if (methodeToCompute[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		PFGSG2.init(cas2, &PQ2G, &PQ2GD, true);
		if (setE) {
			PFGSG2.setE(&EGD);
		}
		if (setSolve) {
			PFGSG2.solve();
		}
		else {
			PFGSG2.calcW(true);
		}
		if (save) {
			PFGSG2.saveTimeBlock(fileName);
		}
		PFGSG2.calcE();
		t2 = std::chrono::high_resolution_clock::now();
		PFGSG2.display(displayAll);
		results.set(0, method, PFGSG2.getP0());
		results.set(1, method, PFGSG2.getQ0());
		results.set(2, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(3, method, PFGSG2.getIter());
	}
	method++;
	

	std::cout << "*****************************Current Flow ********************************" << std::endl;
	
	if (methodeToCompute[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		PF3.init(cas2, &PQ2);
		if (PF3.chekcase()) {
			if (setE) {
				PF3.setE(&E);
			}
			if (setSolve) {
				PF3.solve();
			}
			if (save) {
				PF3.saveTimeBlock(fileName);
			}
			PF3.calcE();
			PF3.calcW(true);
			t2 = std::chrono::high_resolution_clock::now();
			PF3.display(displayAll);
			results.set(0, method, PF3.getP0());
			results.set(1, method, PF3.getQ0());
			results.set(2, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
			results.set(3, method, PF3.getIter());
		}
		else {
			std::cout << "pas reseau de distribution" << std::endl;
		}/**/
	}
	method++;

	std::cout << "*****************************Back-For PQ ********************************" << std::endl;


	if (methodeToCompute[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		PF3PQ.init(cas2, &PQ2);
		if (PF3PQ.chekcase()) {
			if (setE) {
				PF3PQ.setE(&E);
			}
			if (setSolve) {
				PF3PQ.solve();
			}
			if (save) {
				PF3PQ.saveTimeBlock(fileName);
			}
			PF3PQ.calcE();
			PF3PQ.calcW(true);
			t2 = std::chrono::high_resolution_clock::now();
			PF3PQ.display(displayAll);
			results.set(0, method, PF3PQ.getP0());
			results.set(1, method, PF3PQ.getQ0());
			results.set(2, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
			results.set(3, method, PF3PQ.getIter());
		}
		else {
			std::cout << "pas reseau de distribution" << std::endl;
		}
	}
	method++;

	std::cout << "*****************************Back-For PQ GPU ********************************" << std::endl;


	if (methodeToCompute[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		PF3PQG.init(cas2, &PQ2G);
		if (PF3PQG.chekcase()) {
			if (setE) {
				PF3PQG.setE(&EG);
			}
			if (setSolve) {
				PF3PQG.solve();
			}
			if (save) {
				PF3PQG.saveTimeBlock(fileName);
			}
			PF3PQG.calcE();
			PF3PQG.calcW(true);
			t2 = std::chrono::high_resolution_clock::now();
			PF3PQG.display(displayAll);
			results.set(0, method, PF3PQG.getP0());
			results.set(1, method, PF3PQG.getQ0());
			results.set(2, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
			results.set(3, method, PF3PQG.getIter());
		}
		else {
			std::cout << "pas reseau de distribution" << std::endl;
		}
	}
	method++;
	std::cout << "NRC    NRDC    NRG    NRGD    GSC    GSCD    GSG    GSGD    CuC    PQ    PQGPU" << std::endl;

	results.display();


}

void testOPF(int choseCase, std::string chosenCase)
{
	StudyCase cas;
	int million = 1000000;
	//int choseCase = 2;// rxRX
	std::string fileName = "OPFISGTResidualscase85.csv"; //"TimeByBlockPF";
	//std::string chosenCase = "";
	
	int nMethode = 6;
	bool methodeToCompute[] = { true, true, true, true, true, true };
	int nMethodeReal = 6;

	OPFADMM opfADMM;
	OPFADMM2 opfADMM2;
	OPFADMMCons opfADMMCons;
	OPFADMMGPU opfADMMGPU;
	OPFADMMGPU2 opfADMMGPU2;
	OPFADMMConsGPU opfADMMConsGPU;

	bool saveResiduals = false; //true false
	bool doubleSolve   = false;
	MatrixCPU results(5, nMethode, -1);
	int method = 0;
	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;

	float rhoInit = 1; // 1 pou cas 2 ou 10 50 cas 69
	

	switch (choseCase)
	{
	case 0:
		cas.SetAC2node();
		cas.display();
		break;
	case 1:
		//chosenCase = "case85"; //case10ba case4_dist case85 case69 ?, pas radial : case30 case300
		cas.SetACFromFile(chosenCase);
		cas.display();
		break;
	case 2:
		cas.SetEuropeTestFeeder();
		cas.display();
		break;
	case 4:
		cas.SetAC3Bus();
		cas.display();
		break;
	case 5:
		//chosenCase = "RandRadial";
		cas.genGridBT(30, 80, 200, 0.001, 0.0005);
		cas.genAgentsAC(30, 0.5, 0.5, 1, 0.5, 1, 0.1, 0.5, 1, 0.5, 1, 0.2);
		cas.genLinkGridAgent();
		cas.display();
		break;
	default:
		throw std::invalid_argument("unknown choseCase");
		break;
	}

	Simparam param(cas.getNagent(), cas.getNLine(true), true);
	float epsL = 0.001f;
	float epsG = 0.01f;
	int stepG = 1;
	int stepL = 5;
	int iterG = 10000;
	int iterL = 10000;
	param.setEpsL(epsL);
	param.setEpsG(epsG);
	param.setItG(iterG); //500000  10000
	param.setItL(iterL);
	param.setStep(stepG, stepL);
	param.setRho(rhoInit);
	 

	Simparam res(param);

	int nAgent = cas.getNagent();
	MatrixCPU Pn;
	bool radial = cas.isRadial();

	MatrixCPU Param(1, 13);
	Param.set(0, 0, choseCase);
	Param.set(0, 1, nAgent);
	Param.set(0, 2, cas.getNBus());
	Param.set(0, 3, nMethodeReal);
	Param.set(0, 4, rhoInit);
	Param.set(0, 5, epsG);
	Param.set(0, 6, epsL);
	Param.set(0, 7, iterG);
	Param.set(0, 8, iterL);
	Param.set(0, 9, stepG);
	Param.set(0, 10, stepL);
	Param.set(0, 11, opfADMM.getMurhoVar());
	Param.set(0, 12, opfADMM.getTaurhoVar());

	std::cout << "*********************************** OPFADMM ************************************" << std::endl;

	
	/*param.display(1);*/
	if (radial && methodeToCompute[method])
	{
		t1 = std::chrono::high_resolution_clock::now();
		opfADMM.solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		if(doubleSolve) opfADMM.solve(&res, param, cas);
		Pn = res.getPn();
		if (saveResiduals) {
			Param.saveCSV(fileName);
			MatrixCPU Residuals(res.getRes());
			Residuals.saveCSV(fileName);
		}
		

		Pn.display();
		// 
		//opfADMM.display();

		results.set(0, method, Pn.get(0, 0));
		results.set(1, method, Pn.get(nAgent, 0));
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());
		res.display();
	}


	
	
	method++;
	std::cout << "********************************** OPFADMM 2 ***********************************" << std::endl;

	
	/*param.display(1);*/
	if (radial && methodeToCompute[method])
	{
		t1 = std::chrono::high_resolution_clock::now();
		opfADMM2.solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		if(doubleSolve) opfADMM2.solve(&res, param, cas);
		Pn = res.getPn();
		Pn.display();
		//opfADMM2.display();
		if (saveResiduals) {
			MatrixCPU Residuals(res.getRes());
			Residuals.saveCSV(fileName);
		}

		results.set(0, method, Pn.get(0, 0));
		results.set(1, method, Pn.get(nAgent, 0));
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());
		res.display();
	}

	
	method++;

	std::cout << "*****************************OPFADMM  Cons****************************************" << std::endl;


	if (radial && methodeToCompute[method])
	{
		t1 = std::chrono::high_resolution_clock::now();
		opfADMMCons.solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		if(doubleSolve) opfADMMCons.solve(&res, param, cas);
		Pn = res.getPn();
		Pn.display();
		//opfADMMCons.display();
		if (saveResiduals) {
			MatrixCPU Residuals(res.getRes());
			Residuals.saveCSV(fileName);
		}

		results.set(0, method, Pn.get(0, 0));
		results.set(1, method, Pn.get(nAgent , 0));
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

		res.display();
	}

	

	method++;

	std::cout << "*********************************OPFADMM GPU**********************************" << std::endl;

	//param.setEpsL(epsL / 20);
	if (radial && methodeToCompute[method])
	{
		t1 = std::chrono::high_resolution_clock::now();
		opfADMMGPU.solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		if(doubleSolve) opfADMMGPU.solve(&res, param, cas);
		Pn = res.getPn();
		Pn.display();
		//opfADMMGPU.display();
		if (saveResiduals) {
			MatrixCPU Residuals(res.getRes());
			Residuals.saveCSV(fileName);
		}

		results.set(0, method, Pn.get(0, 0));
		results.set(1, method, Pn.get(nAgent, 0));
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

		res.display();
		//opfADMMGPU.feasiblePoint();
	}


	CHECK_LAST_CUDA_ERROR();
	method++;


	std::cout << "*****************************OPFADMM 2 GPU****************************************" << std::endl;
	
	if (radial && methodeToCompute[method])
	{
		t1 = std::chrono::high_resolution_clock::now();
		opfADMMGPU2.solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		if(doubleSolve) opfADMMGPU2.solve(&res, param, cas);
		Pn = res.getPn();
		//Pn.display();
		if (saveResiduals) {
			MatrixCPU Residuals(res.getRes());
			Residuals.saveCSV(fileName);
		}
		//opfADMMGPU.display();

		results.set(0, method, Pn.get(0, 0));
		results.set(1, method, Pn.get(nAgent, 0));
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

		res.display();
		opfADMMGPU2.feasiblePoint();
	}

	CHECK_LAST_CUDA_ERROR();

	method++;
	std::cout << "*****************************OPFADMM Cons GPU****************************************" << std::endl;

	res = param;
	if (radial && methodeToCompute[method])
	{
		t1 = std::chrono::high_resolution_clock::now();
		opfADMMConsGPU.solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		if(doubleSolve) opfADMMConsGPU.solve(&res, param, cas);
		Pn = res.getPn();
		//Pn.display();
		//opfADMMConsGPU.display();
		if (saveResiduals) {
			MatrixCPU Residuals(res.getRes());
			Residuals.saveCSV(fileName);
		}

		results.set(0, method, Pn.get(0, 0));
		results.set(1, method, Pn.get(nAgent, 0));
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

		res.display();
	}


	CHECK_LAST_CUDA_ERROR();
	method++;

	

	std::cout << "*****************************************************************************" << std::endl;
	std::cout << "OPFADMM -" << " OPFADMM 2 " << " OPFADMM Cons -" << " OPFADMMGPU -" << " OPFADMM 2 -" << " OPFADMM Cons GPU" << std::endl;
	results.display();


}

void testADMMACConst()
{
	StudyCase cas;
	cas.SetAC3Bus();

	System sys;
	sys.setStudyCase(cas);
	sys.setMethod("ADMMACConst1");

	float epsG = 0.0001;
	float epsGC = 0.001;
	float epsL = 0.000001;
	sys.setEpsG(epsG);
	sys.setEpsGC(epsGC);
	sys.setEpsL(epsL);

	float rhoG = 60;
	float rho1 = 60;
	sys.setRho(rhoG);
	sys.setRho1(rho1);

	int iterL = 20000;
	int iterG = 20000;
	sys.setIter(iterG, iterL);

	Simparam res = sys.solve();
}

void testMarket(int choseCase, std::string chosenCase, bool AC, int sizeN)
{
	StudyCase cas;
	int million = 1000000;
	std::string fileName = "MarketTimeByBlock"; //MarketTimeByBlock
		
	int nMethode = 8;
	bool methodeToSimule[10] = {false, true, true, true, false, false, false, false, false, false };

	int method = 0;
	bool doubleSolve = (sizeN==1);
	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;

	// for random case
	float P = 100;
	float dP = 10;
	float a = 1;
	float da = 0.1;

	float b = 2;
	float db = 1;
	int agents = sizeN;
	float beta = 2;
	float dBeta = 1;
	float propCons = 0.4;
	float propPro = 0;

	//bool AC = true;// true false

	// for european case
	int nCons = 1494;
	std::string path = "data/";
	int Nhour = 24 * 31; // janvier
	MatrixCPU P0Global(nCons, Nhour);
	MatrixCPU P0(nCons, 1);
	std::string nameP0 = path + "Europe/load/Month/2012-01.txt";
	//P0Global.setFromFile(nameP0, 1);
	//P0Global.getBloc(&P0, 0, nCons, 0, 1); // 1ere heure de l'ann�e
	//int nStep = 1; // WIP
	//int nSimu = 1; // WIP

	float rho = 1;
	//float factor = 0.9;


	
	switch (choseCase)
	{
	case 0:
		if (AC) {
			cas.SetAC2node();
			//chosenCase = "caseAC2node";
		}
		else {
			cas.Set2node();
			nMethode += 2;
			cas.setReduce(true);
			//chosenCase = "case2node";
		}
		//cas.display();
		break;
	case 1:
		
		if (AC) {
			//chosenCase = "case10ba"; //case10ba case4_dist case85 case118
			cas.SetACFromFile(chosenCase);
		}
		else
		{
			cas.Set29node();
			nMethode += 2;
			cas.setReduce(true);
			//chosenCase = "case29node";
		}
		//cas.display();
		break;
	case 2:
		if (AC) {
			cas.SetEuropeTestFeeder();
			//chosenCase = "caseTestFeeder";
		}
		else {
			cas.SetEuropeP0WithoutConstraint(path, &P0);
			nMethode += 2;
			cas.setReduce(true);
			//chosenCase = "caseEurope";
		}
		//cas.display();
		break;
	case 3:
		if (AC) {
			cas.SetAC3Bus();
			//chosenCase = "caseAC3node";
		}
		else
		{
			cas.Set3Bus();
			nMethode += 2;
			cas.setReduce(true);
			//chosenCase = "case3node";
		}
		//cas.display();
		break;
	case 4:
		if (AC) {
			//cas = StudyCase(agents, P, dP, P, dP, a, da, a, da, b, db, beta, dBeta, propCons, 0, propPro);
			//chosenCase = "randomAC" + std::to_string(agents);
			cas.genGridBT(2, 2, 2, 0.01, 0.005);
			cas.genAgentsAC(sizeN, 0.5, 0.25, 0.5, 0.1, 1, 0.1, 0.005, 2, 0.1, 2, 1);
			cas.genLinkGridAgent();
			}
		else
		{
			cas = StudyCase(agents, P, dP, a, da, b, db, beta, dBeta, propCons, 0, propPro);
			nMethode += 2;
			cas.setReduce(true);
			//chosenCase = "random" + std::to_string(agents);
		}
		//cas.display();
		break;
	default:
		throw std::invalid_argument("unknown choseCase");
		break;
	}

	cas.display();
	fileName += chosenCase + ".csv";
	MatrixCPU results(5, nMethode, -1);
	//OSQPCentralized2 osqpCen;
	ADMMMarket admmMarket;
	ADMMMarketOpenMP admmMarketOpenMP;
	ADMMMarketGPU admmMarketGPU;
	//OSQP osqp;
	PAC pac;
	PACGPU pacGPU;
	PACOpenMP pacOpenMP;
	ADMMConst admmMarketEndo;
	ADMMGPUConst4 admmMarketEndoGPU;

	Simparam param(cas.getNagent(), cas.getNLine(true), cas.getNLine(), AC);
	float epsG = 0.00001f;
	float epsL = 0.000001f; //float epsL = 0.0005f;
	//float epsG = 0.001f;
	int iterGlobal = 100000; // 100000
	int iterLocal = 2000;
	int stepG = 1;
	int stepL = 5;

	param.setEpsL(epsL);
	param.setEpsG(epsG);
	param.setItG(iterGlobal); //500000
	param.setItL(iterLocal);
	param.setStep(stepG, stepL);

	MatrixCPU Param(1, 10);
	Param.set(0, 0, cas.getNagent());
	Param.set(0, 1, rho);
	Param.set(0, 2, epsG);
	Param.set(0, 3, epsL);
	Param.set(0, 4, iterGlobal);
	Param.set(0, 5, iterLocal);
	Param.set(0, 6, stepG);
	Param.set(0, 7, stepL);
	Param.set(0, 8, nMethode);
	Param.set(0, 9, AC);

#ifdef INSTRUMENTATION
	Param.saveCSV(fileName);
#endif // INSTRUMENTATION

	param.setRho(rho);
	
	Simparam res(param);
	int nAgent = cas.getNagent();
	MatrixCPU Pn;
	/*param.display(1);*/
	#ifdef OSQP
	std::cout << "**************************     OSQP Cen      ****************************************" << std::endl;

	if (!AC && methodeToSimule[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		osqpCen.solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		Pn = res.getPn();
		Pn.display();


		results.set(0, method, Pn.get(1, 0));
		if (AC) {
			results.set(1, method, Pn.get(nAgent + 1, 0));
		}
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION
		res.display();
	}
#endif

	method++;


	std::cout << "**************************     ADMMMarket      ****************************************" << std::endl;

	if (methodeToSimule[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		admmMarket.solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		if(doubleSolve) admmMarket.solve(&res, param, cas);
		Pn = res.getPn();
		Pn.display();


		results.set(0, method, Pn.get(1, 0));
		if (AC) {
			results.set(1, method, Pn.get(nAgent + 1, 0));
		}
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION
		res.display();
	}
	

	method++;
	std::cout << "***************************** admmMarketOpenMP     ****************************************" << std::endl;

	if (methodeToSimule[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		admmMarketOpenMP.solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		if(doubleSolve) admmMarketOpenMP.solve(&res, param, cas);
		Pn = res.getPn();
		Pn.display();

		results.set(0, method, Pn.get(1, 0));
		if (AC) {
			results.set(1, method, Pn.get(nAgent + 1, 0));
		}
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION

		res.display();
	}
	method++;
	if (method > (nMethode + 1)) {
		throw std::invalid_argument("nMethod is too small");
	}
	std::cout << "********************************** admmMarketGPU **************************************" << std::endl;

	if (methodeToSimule[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		admmMarketGPU.solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		if(doubleSolve) admmMarketGPU.solve(&res, param, cas);
		Pn = res.getPn();
		Pn.display();

		results.set(0, method, Pn.get(1, 0));
		if (AC) {
			results.set(1, method, Pn.get(nAgent + 1, 0));
		}
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION
		res.display();
		
	}
	method++;
	if (method > (nMethode + 1)) {
		throw std::invalid_argument("nMethod is too small");
	}
	CHECK_LAST_CUDA_ERROR();
	std::cout << "********************************** OSQP  *****************************************" << std::endl;
#ifdef OSQP
	if (methodeToSimule[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		osqp.solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		Pn = res.getPn();
		Pn.display();

		results.set(0, method, Pn.get(1, 0));
		if (AC) {
			results.set(1, method, Pn.get(nAgent + 1, 0));
		}
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION

		res.display();
	}
#endif
	method++;
	if (method > (nMethode)) {
		throw std::invalid_argument("nMethod is too small");
	}

	std::cout << "******************************       PAC     ***************************************" << std::endl;

	if (methodeToSimule[method]) {
		pac.setBestRhoGammaHeuristic(cas);
		t1 = std::chrono::high_resolution_clock::now();
		pac.solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();

		if(doubleSolve) pac.solve(&res, param, cas);

		Pn = res.getPn();
		Pn.display();

		results.set(0, method, Pn.get(1, 0));
		if (AC) {
			results.set(1, method, Pn.get(nAgent + 1, 0));
		}
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION

		res.display();
	}
	method++;
	if (method > (nMethode)) {
		throw std::invalid_argument("nMethod is too small");
	}
	std::cout << "******************************       PAC OpenMP    ***************************************" << std::endl;

	if (methodeToSimule[method]) {
		pacOpenMP.setBestRhoGammaHeuristic(cas);
		t1 = std::chrono::high_resolution_clock::now();
		pacOpenMP.solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();

		if(doubleSolve) pacOpenMP.solve(&res, param, cas);

		Pn = res.getPn();
		Pn.display();

		results.set(0, method, Pn.get(1, 0));
		if (AC) {
			results.set(1, method, Pn.get(nAgent + 1, 0));
		}
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION

		res.display();
	}
	method++;
	if (method > (nMethode)) {
		throw std::invalid_argument("nMethod is too small");
	}

	std::cout << "******************************       PAC GPU     ***************************************" << std::endl;

	if (methodeToSimule[method]) {
		pacGPU.setBestRhoGammaHeuristic(cas);
		t1 = std::chrono::high_resolution_clock::now();
		pacGPU.solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();

		if(doubleSolve) pacGPU.solve(&res, param, cas);

		Pn = res.getPn();
		Pn.display();

		results.set(0, method, Pn.get(1, 0));
		if (AC) {
			results.set(1, method, Pn.get(nAgent + 1, 0));
		}
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION

		res.display();
	}
	method++;
	if (method > (nMethode)) {
		throw std::invalid_argument("nMethod is too small");
	}
	CHECK_LAST_CUDA_ERROR();
	std::cout << "********************************  admmMarketEndo  ***************************************" << std::endl;

	if (!AC && methodeToSimule[method]) {

		t1 = std::chrono::high_resolution_clock::now();
		admmMarketEndo.solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();

		if(doubleSolve) admmMarketEndo.solve(&res, param, cas);

		Pn = res.getPn();
		Pn.display();
		results.set(0, method, Pn.get(1, 0));
		if (AC) {
			results.set(1, method, Pn.get(nAgent + 1, 0));
		}
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION

		res.display();
		
		
		if (method + 1 > (nMethode)) {
			throw std::invalid_argument("nMethod is too small");
		}
	}
	method++;
	
	std::cout << "**********************************   admmMarketEndoGPU    *******************************************" << std::endl;

	if (!AC && methodeToSimule[method]) {

		t1 = std::chrono::high_resolution_clock::now();
		admmMarketEndoGPU.solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();

		if(doubleSolve) admmMarketEndoGPU.solve(&res, param, cas);

		Pn = res.getPn();
		Pn.display();

		results.set(0, method, Pn.get(1, 0));
		if (AC) {
			results.set(1, method, Pn.get(nAgent + 1, 0));
		}
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION

		res.display();
		
		if (method + 1 > (nMethode)) {
			throw std::invalid_argument("nMethod is too small");
		}

	}
	method++;
	CHECK_LAST_CUDA_ERROR();
	std::cout << "*****************************************************************************" << std::endl;
	std::cout << "OSQPCen,  ADMMMarket , ADMMMarketOpenMP, ADMMMarketGPU, OSQP, PAC, PAC OpenMp, PACGPU ";
	if (!AC) {
		std::cout << " admmMarketEndo , admmMarketEndoGPU ";
	}
	std::cout << std::endl;


	results.display();
}

void testMarketEndo(int choseCase, std::string chosenCase, int sizeN, int sizeB)
{
	MarEndoCons* marketEndoCPU = new MarEndoCons;
	MarEndoConsGPU* marketEndoGPU = new MarEndoConsGPU;

	MarketEndoDirect* marketEndoDirCPU = new MarketEndoDirect;
	MarketEndoDirectGPU* marketEndoDirGPU = new MarketEndoDirectGPU;

	EndoPF* endoPF = new EndoPF;
	EndoPFGPU* endoPFGPU = new EndoPFGPU;

	ADMMConst* endoDCPF = new ADMMConst;
	ADMMGPUConst4* endoDCPFGPU = new ADMMGPUConst4;

	StudyCase cas;
	int nMethode = 8;
	bool methodeToSimule[8] = { true, true, false, false, true, true, false, false }; ///false

	int million = 1000000;
	bool doubleSolve =  (sizeN==1); // true  false
	std::string fileName = "TimeByBlockMarketEndo";
	//std::string chosenCase = "";
	int offsetAgent = 0; // which agent we kept the value to compare results
	
	std::cout << "choseCase " << choseCase << "chosenCase " << chosenCase << "doubleSolve " << doubleSolve << " N " << sizeN << " B " << sizeB <<std::endl;

	int stepG = 1;
	int stepL = 1;
	int stepIntern = 1;

	if(stepG < stepIntern){
		stepG = stepIntern;
	}

	int iterL = 7000;
	int iterG  = 10000;
	int iterIntern = 5000;

	float epsL = 0.0001f;
	float epsG = 0.001f;
	float epsGC = 0.001f;
	float epsIntern = 0.0001f;

	float rhoInit = 1; // 1 pour cas 2 noeuds, 5 pour cas9, cas 10, 10 cas 69

	MatrixCPU results(5, nMethode, nanf(""));
	int method = 0;
	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;

	switch (choseCase)
	{
	case 0:
		//chosenCase = "case10ba";// 9, 30 57  118 case30  // radial case10ba case4_dist case85 case69
		cas.SetACFromFile(chosenCase); //case_ACTIVSg2000	
		cas.display();
		break;
	case 1:
		//chosenCase = "RandRadial";
		cas.genGridBT(sizeB, 0.8 * sizeB, 2*sizeB, 0.01, 0.005);
		cas.genAgentsAC(sizeN, 0.5, 0.25, 0.5, 0.1, 1, 0.1, 0.005, 2, 0.1, 0, 0);
		cas.genLinkGridAgent();
		cas.display();
		
		break;
	case 2:
		//chosenCase = "RandHTB";
		cas.SetACFromFile(chosenCase);
		cas.genAgentsAC(sizeN, 0.5, 0.25, 0.5, 0.1, 1, 0.1, 0.005, 2, 0.1, 0, 0);
		cas.genLinkGridAgent();
		cas.display();
		break;
	case 3:
		//chosenCase = "EuropeTestFeeder";
		cas.SetEuropeTestFeeder();
		cas.display();
		break;
	case 4:
		//chosenCase = "2node";
		cas.SetAC2node();
		cas.display();
		break;
	case 5:
		//chosenCase = "3node";
		cas.SetAC3Bus();
		cas.display();
		cas.display(1);
		
		break;
	case 6:
		//chosenCase = "case30";// 9, 30 57 69 118  // radial case10ba case4_dist case85
		cas.SetACFromFileSimplify(chosenCase); //case_ACTIVSg2000	
		//cas.display();
		break;
	default:
		//std::cout << "unknown choice, case 0" << std::endl;
		return;
	}
	if (cas.isAC() && (methodeToSimule[3] || methodeToSimule[7])) {
		std::cout << "ajout du DC" << std::endl;
		cas.genDCGridFromAC(); // pour utiliser les m�thodes DC
		cas.setReduce(true);
		//cas.display(1);
		//cas.display(2);
	}
	fileName += chosenCase + ".csv";
	Simparam param(cas.getNagent(), cas.getNLine(true), true);
	
	param.setEpsL(epsL);
	param.setEpsG(epsG);//0.005f FB
	param.setEpsGC(epsGC); //0.001f FB
	param.setEpsIntern(epsIntern);
	param.setItG(iterG); //20000 500000
	param.setItL(iterL);
	param.setItIntern(iterIntern);
	param.setStep(stepG, stepL, stepIntern);
	
	//for (int i = 0; i < 5; i++) {
	param.setRho(rhoInit);
	Simparam res(param);
	int nAgent = cas.getNagent();
	MatrixCPU Pn;
	bool radial = cas.isRadial();


	std::cout << "********************* Market Endogen CPU Direct ************************************" << std::endl;


	param.setRho(rhoInit);
	if (radial && methodeToSimule[method])
	{
		t1 = std::chrono::high_resolution_clock::now();
		marketEndoDirCPU->solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		if(doubleSolve) marketEndoDirCPU->solve(&res, param, cas);
		
		MatrixCPU Res = res.getRes();
		//Res.saveCSV("ResCPU.csv");
		Pn = res.getPn();
		Pn.display();
		//opfADMM.display();

		results.set(0, method, Pn.get(offsetAgent, 0));
		results.set(1, method, Pn.get(nAgent + offsetAgent, 0));
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION
		res.display();
	}

	method++;

	std::cout << "********************* Market Endogen CPU  ************************************" << std::endl;


	param.setRho(rhoInit);
	if (radial && methodeToSimule[method])
	{
		t1 = std::chrono::high_resolution_clock::now();
		marketEndoCPU->solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		if (doubleSolve) marketEndoCPU->solve(&res, param, cas);

		
		Pn = res.getPn();
		Pn.display();
		//opfADMM.display();

		results.set(0, method, Pn.get(offsetAgent, 0));
		results.set(1, method, Pn.get(nAgent + offsetAgent, 0));
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());
#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION
		res.display();
	}

	method++;


	std::cout << "********************* Endo with PF CPU  ************************************" << std::endl;

	param.setRho(rhoInit*20);
	if (methodeToSimule[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		endoPF->solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		if (doubleSolve) endoPF->solve(&res, param, cas);

		Pn = res.getPn();
		Pn.display();
		//opfADMM.display();

		results.set(0, method, Pn.get(offsetAgent, 0));
		results.set(1, method, Pn.get(nAgent + offsetAgent, 0));
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION

		res.display();
	}
	

	method++;


	std::cout << "********************* Endo with DC-PF CPU  ************************************" << std::endl;
	param.setRho(rhoInit);
	if (methodeToSimule[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		endoDCPF->solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		if (doubleSolve) endoDCPF->solve(&res, param, cas);

		Pn = res.getPn();
		Pn.display();
		//opfADMM.display();

		results.set(0, method, Pn.get(offsetAgent, 0));
		//results.set(1, method, Pn.get(nAgent + 1, 0));
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION

		res.display();
	}


	method++;
	std::cout << "********************* Market Endogen GPU Direct ************************************" << std::endl;

	CHECK_LAST_CUDA_ERROR();
	param.setRho(rhoInit);
	if (radial && methodeToSimule[method])
	{
		t1 = std::chrono::high_resolution_clock::now();
		marketEndoDirGPU->solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		if (doubleSolve) marketEndoDirGPU->solve(&res, param, cas);
		MatrixCPU Res = res.getRes();
		//Res.saveCSV("ResGPU.csv");
		Pn = res.getPn();
		Pn.display();
		//opfADMM.display();

		results.set(0, method, Pn.get(offsetAgent, 0));
		results.set(1, method, Pn.get(nAgent + offsetAgent, 0));
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION

		res.display();
	}
	CHECK_LAST_CUDA_ERROR();

	method++;
	std::cout << "********************* Market Endogen GPU  ************************************" << std::endl;


	param.setRho(rhoInit);
	if (radial && methodeToSimule[method])
	{
		t1 = std::chrono::high_resolution_clock::now();
		marketEndoGPU->solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		if (doubleSolve) marketEndoGPU->solve(&res, param, cas);

		Pn = res.getPn();
		Pn.display();
		//marketEndoGPU.display();

		results.set(0, method, Pn.get(offsetAgent, 0));
		results.set(1, method, Pn.get(nAgent + offsetAgent, 0));
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION

		res.display();
	}

	method++;
	CHECK_LAST_CUDA_ERROR();
	
	std::cout << "********************* Endo with PF GPU  ************************************" << std::endl;

	param.setRho(rhoInit * 20);
	if (methodeToSimule[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		endoPFGPU->solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		if (doubleSolve) endoPFGPU->solve(&res, param, cas);

		Pn = res.getPn();
		Pn.display();
		//opfADMM.display();

		results.set(0, method, Pn.get(offsetAgent, 0));
		results.set(1, method, Pn.get(nAgent + offsetAgent, 0));
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION
		
		res.display();
	}


	method++;
	CHECK_LAST_CUDA_ERROR();
	
	
	std::cout << "********************* Endo with DC-PF GPU  ************************************" << std::endl;

	param.setRho(rhoInit);
	if (methodeToSimule[method]) {
		t1 = std::chrono::high_resolution_clock::now();
		endoDCPFGPU->solve(&res, param, cas);
		t2 = std::chrono::high_resolution_clock::now();
		if (doubleSolve) endoDCPFGPU->solve(&res, param, cas);

		Pn = res.getPn();
		Pn.display();
		//opfADMM.display();

		results.set(0, method, Pn.get(offsetAgent, 0));
		//results.set(1, method, Pn.get(nAgent + 1, 0));
		results.set(2, method, res.getFc());
		results.set(3, method, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / million);
		results.set(4, method, res.getIter());

#ifdef INSTRUMENTATION
		res.displayTime(fileName);
#endif // INSTRUMENTATION

		res.display();
	}


	method++;
	CHECK_LAST_CUDA_ERROR();
	
	std::cout << "*****************************************************************************" << std::endl;
	std::cout << "EndoDirect - EndoCons - ACEndo - DCEndo - EndoDirectGPU - EndoConsGPU - ACEndoGPU - DCEndoGPU" << std::endl;
	
	
	
	results.display();

	DELETEB(endoDCPF);
	DELETEB(endoDCPFGPU);
	DELETEB(endoPF);
	DELETEB(endoPFGPU);
	DELETEB(marketEndoDirGPU);
	DELETEB(marketEndoDirCPU);
	DELETEB(marketEndoCPU);
	DELETEB(marketEndoGPU);


}






/*


*/