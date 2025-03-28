#include "../head/TestStudyCase.h"

int testStudyCase()
{
	int n = 1;
	if (!testSCConstru()) return n;
	n++;
	if (!testSCConnect()) return n;
	n++;
	if (!testSC39Bus()) return n;
	n++;
	if (!testSCUpdateP0()) return n;
	n++;
	if (!testRemoveLink()) return n;
	n++;
	if (!testAddLink()) return n;
	n++;
	if (!testGenRandom()) return n;
	n++;
	
	
	return 0;
}



bool testSCConstru()
{
	/*  StudyCase();
	StudyCase(int nAgent, float P, float dP, float a, float da, float b, float db);
	*/
	int nAgent = 4;
	float P = 10;
	float dP = 2;
	float a = 0.06f;
	float da = 0.f;
	float b = 50;
	float db = 10;

	std::cout << "default constructor" << std::endl;
	StudyCase s1;
	std::cout << " constructor " << std::endl;
	StudyCase s2( nAgent,  P,  dP,  a,  da,  b, db);
	
	std::cout << "29 nodes " << std::endl;
	StudyCase s3;
	s3.Set29node();
	
	std::cout << "2-times constructor" << std::endl;
	StudyCase s4;
	s4 = StudyCase(nAgent, P, dP, a, da, b, db);
	
	std::cout << "2 nodes " << std::endl;
	StudyCase s5;
	s5.Set2node();
	

	return true;
}

bool testSCConnect()
{
	int nPro = 1;
	int nGen = 2;
	int nCons = 3;
	int nAgent = nPro + nGen + nCons;
	MatrixCPU connec(nAgent, nAgent);
	for (int i = 0; i < nCons;i++) {
		for (int j = nCons; j < nAgent;j++) {
			connec.set(i, j, 1);
		}
	}
	for (int i = nCons; i < nCons + nGen; i++) {
		for (int j = 0; j < nCons;j++) {
			connec.set(i, j, 1);
		}
		for (int j = nCons + nGen; j < nAgent;j++) {
			connec.set(i, j, 1);
		}
	}
	for (int i = nCons + nGen; i < nAgent;i++) {
		for (int j = 0; j < nCons; j++) {
			connec.set(i, j, 1);
		}
		for (int j = nCons; j < nCons + nGen;j++) {
			connec.set(i, j, 1);
		}
	}
	StudyCase s2(nAgent,1,0,1,0,1,0,(float) nCons/nAgent, (float) nPro/nAgent);
	MatrixCPU connec2 = s2.getC();
	

	return connec2.isEqual(&connec);

}

bool testSC39Bus()
{
	StudyCase cas;
	std::string path = "test/";
	cas.Set39Bus(path);
	cas.display();
	std::cout << "-------------------------------------------------------- " << std::endl;
	return true;
}


bool testSCUpdateP0()
{
	StudyCase cas;
	std::string path = "test/";
	std::string date =  "2012-03";
	std::string date2 = "2012-01";
	std::string nameP0 = path + date + ".txt";
	std::string nameP02 = path + date2 + ".txt";
	
	MatrixCPU P0(1494, 1);
	P0.setFromFile(nameP0, 1);
	cas.SetEuropeP0(path, &P0, true);
		
	//cas.saveCSV("StudyCaseEurope.csv", false);
	P0.setFromFile(nameP02,1);
	cas.UpdateP0(&P0);
	//cas.display();
	return true;
}

bool testRemoveLink() {
	int nAgent = 4; // 2 cons, 2 prod
	float P = 10;
	float dP = 2;
	float a = 0.06f;
	float da = 0.01f;
	float b = 50;
	float db = 10;

	StudyCase s(nAgent, P, dP, a, da, b, db);
	MatrixCPU Connect(s.getC());
	MatrixCPU Voisin(s.getNvoi());
	Connect.set(0, 2, 0);
	Connect.set(2, 0, 0);
	Voisin.set(0, 0, 1);
	Voisin.set(2, 0, 1);


	try
	{
		s.removeLink(0, 1);
	}
	catch (const std::exception&)
	{
		try
		{
			s.removeLink(0, 4);
		}
		catch (const std::exception&)
		{
			s.removeLink(0, 2);
			MatrixCPU Connect2(s.getC());
			MatrixCPU Voisin2(s.getNvoi());
			s.display();
			Agent a = s.getAgent(0);
			Agent a2 = s.getAgent(2);
			
			a.display();
			a2.display();
			std::cout << "-------------------------------------------------------- " << std::endl;
			return (Connect.isEqual(&Connect2) && Voisin.isEqual(&Voisin2));

		}
		return false;
	}
	return false;

}

bool testAddLink() {
	int nAgent = 4; // 2 cons, 2 prod
	float P = 10;
	float dP = 2;
	float a = 0.06f;
	float da = 0.01f;
	float b = 50;
	float db = 10;

	StudyCase s(nAgent, P, dP, a, da, b, db);
	MatrixCPU Connect(s.getC());
	MatrixCPU Voisin(s.getNvoi());
	s.removeLink(0, 2);

	try
	{
		s.addLink(0, 1);
	}
	catch (const std::exception&)
	{
		try
		{
			s.addLink(0, 4);
		}
		catch (const std::exception&)
		{
			s.addLink(0, 2);
			MatrixCPU Connect2(s.getC());
			MatrixCPU Voisin2(s.getNvoi());
			std::cout << "-------------------------------------------------------- " << std::endl;
			return (Connect.isEqual(&Connect2) && Voisin.isEqual(&Voisin2));

		}
		return false;
	}
	return false;

}

bool testGenRandom()
{
	int nAgent = 12;
	float propCons = 0.25f;
	float Pconso = 10;
	float dPconso = 1;
	float bProd = 2;
	float dbProd = 0.1f;
	float Pprod = 20;
	float dPprod = 5;
	float Gamma = 15;
	float dGamma = 3;
	std::string path = "data/";
	int nLine = 30;
	float limit = 100;
	float dlimit = 10;

	StudyCase cas;
	cas.genAgents(nAgent, propCons, Pconso, dPconso, bProd, dbProd, Pprod, dPprod, Gamma, dGamma);
	cas.display();
	cas.genGridFromFile(path,true);
	cas.setReduce(true);
	std::cout << "Gen Link :" << std::endl;
	cas.genLinkGridAgent();
	//cas.display();
	std::cout << "Original constraints :" << std::endl;
	MatrixCPU Llimit(cas.getLineLimit());
	Llimit.display();
	std::cout << "Modified constraints :" << std::endl;
	cas.genLineLimit(nLine, limit, dlimit);
	
	cas.display(1);
	
	std::cout << "-------------------------------------------------------- " << std::endl;
	return true;
}
