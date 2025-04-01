#include "../head/TestSysteme.h"





int testSysteme()
{
	int n = 1;
	if (!testSysConstruct1()) return n;
	n++;
	if (!testSysConstruct2()) return n;
	n++;
	if (!testSysSolve()) return n;


	return 0;
}

bool testSysConstruct1()
{
	std::cout << "default constructor" << std::endl;
	System sys;

	return true;
}


bool testSysConstruct2()
{
	int nAgent = 4;
	float P = 10;
	float dP = 2;
	float a = 0.06f;
	float da = 0.01f;
	float b = 50;
	float db = 10;
	float rho = 1;
	int iterMaxGlobal = 500;
	int iterMaxLocal = 300;
	float epsGlobal = 0.001f;
	float epsLocal = 0.0001f;
	std::string name = "ADMMConst";
	std::cout << "param constructor" << std::endl;
	System sys(rho, iterMaxGlobal, iterMaxLocal, epsGlobal, epsLocal, name, nAgent, P, dP, a, da, b, db);
	//sys.display();
	return true;
}

bool testSysSolve()
{
	std::cout << "simulation test" << std::endl;
	StudyCase cas;
	cas.Set2node();
	
	

	System sys;
	sys.setStudyCase(cas);
	sys.display(1);
	sys.display();
	Simparam res = sys.solve();
	res.display();


	return true;
}
