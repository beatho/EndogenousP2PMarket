 #include "../head/TestSimparam.h"

int testSimparam()
{
	int n = 1;
	if (!testSimConstruc1()) return n;
	n++;
	if (!testSimConstruc2()) return n;
	n++;
	if (!testSimConstruc3()) return n;
	n++;
	if (!testSimConstruc4()) return n;
	n++;
	if (!testIterLtot()) return n;
	n++;

	return 0;
}

/*Simparam() {
		std::cout << "constructeur simparam" << std::endl;
	};
	Simparam(int nAgent);
	Simparam(float rho, int iterMaxGlobal, int iterMaxLocal, float epsGlobal, float epsLocal, int _nAgent);
	*/
bool testSimConstruc1()
{
	std::cout << "default constructor" << std::endl;
	Simparam sim;
	sim.display(1);
	return true;
}
bool testSimConstruc2()
{
	int nAgent = 5;
	std::cout << "agent count constructor" << std::endl;
	Simparam sim(nAgent);
	sim.display(1);
	return true;
}

bool testSimConstruc3()
{
	int nAgent = 5;
	float rho = 1;
	int iterMaxGlobal = 500;
	int iterMaxLocal = 300;
	float epsGlobal = 0.001f;
	float epsLocal = 0.0001f;
	
	std::cout << "param constructor" << std::endl;
	Simparam sim(rho, iterMaxGlobal, iterMaxLocal, epsGlobal, epsLocal, nAgent);
	sim.display(1);
	return true;
}
bool testSimConstruc4()
{
	int nAgent = 5;
	std::cout << "2 times constructor" << std::endl;
	Simparam sim;
	sim = Simparam(nAgent);
	sim.display(1);
	return true;
}

bool testIterLtot()
{
	Simparam res;
	int value = 300;
	res.setItLTot(value);

	int value2 = res.getIterLTot();


	return (value==value2);
}
