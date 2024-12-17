#include "../head/TestAgent.h"

int testAgent()
{
	int n = 1;
	if (!testAConstruc()) return n;
	n++;
	if (!testASet()) return n;
	n++;
	if (!testAupdateP0()) return n;
	return 0;
}

bool testAConstruc()
{
	int id = 1;
	float pLim1 = -10;
	float pLim2 = -1;
	float cost1 = 0.06;
	float cost2 = 0.5;
	int nVoisin = 0;
	MatrixCPU connec(1, 1);
	int nAgent = 1; 
	int type = 1;

	std::cout << "default constructor " << std::endl;
	Agent a1;
	a1.display();
	std::cout << "constructor " << std::endl;
	Agent a2(id, pLim1, pLim2, cost1, cost2, nVoisin, &connec, nAgent, type);
	a2.display();
	std::cout << "copy constructor" << std::endl;
	Agent a3(a2);
	a3.display();
	std::cout << "2-times constructor" << std::endl;
	Agent a4;
	a4 = Agent(id, pLim1, pLim2, cost1, cost2, nVoisin, &connec, nAgent, type);
	a4.display();

	return true;
}
bool testASet()
{
	int id = 1;
	float pLim1 = -10;
	float pLim2 = -1;
	float cost1 = 0.06;
	float cost2 = 0.5;
	int nVoisin = 0;
	MatrixCPU connec(0, 0);
	int nAgent = 1;
	int type = 1;

	Agent a1;
	std::cout << "setter " << std::endl;
	a1.setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &connec, nAgent, type);
	a1.display();
	return true;
}

bool testAupdateP0()
{
	int id = 1;
	float pLim1 = -10;
	float pLim2 = -1;
	float cost1 = 0.06;
	float cost2 = 0.5;
	int nVoisin = 0;
	MatrixCPU connec(0, 0);
	int nAgent = 1;
	int type = 1;
	float P0 = 10;
	float dP = 0.2;

	Agent a1;
	a1.setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &connec, nAgent, type);
	a1.display();
	a1.updateP0(P0, dP);
	a1.display();

	float Pmin = -(1 + dP) * P0;
	float Pmax = -(1 - dP) * P0;
	float a = 1;
	float b = P0 * a;

	if (Pmin != a1.getPmin()) return false;
	if (Pmax != a1.getPmax()) return false;
	if (a != a1.getA()) return false;
	if (b != a1.getB()) return false;
	if (0 != a1.getUb()) return false;
	if (Pmin != a1.getLb()) return false;

	return true;
}
