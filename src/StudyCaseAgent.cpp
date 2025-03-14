#include "../head/StudyCaseAgent.h"




void StudyCaseAgent::setMatFromFile(const std::string& path, const std::string& date, MatrixCPU* Pgen, MatrixCPU* P0, MatrixCPU* costGen)
{
	std::string namePgen = "/PowerMaxCountry.csv";
	std::string nameP0 = "/load/Country_" + date + ".txt";
	std::string namecostGen = "/CoefPoly.csv";
	Pgen->setFromFile(path + namePgen);
	P0->setFromFile(path + nameP0,1);
	costGen->setFromFile(path + namecostGen);

}

void StudyCaseAgent::setGenFromFile(const std::string& path, MatrixCPU* Pgen, MatrixCPU* costGen, MatrixCPU* BusGen)
{
	MatrixCPU Generator(_nGen, 3);
	Generator.setFromFile(path, 1);
	
	for (int i = 0; i < _nGen; i++) {
		
		costGen->set(i, 0, Generator.get(i, 0));
		Pgen->set(i, 0, Generator.get(i, 1));
		BusGen->set(i, 0, Generator.get(i, 2));
	}
}

void StudyCaseAgent::initCaseFromPobj()
{
	int offset = 0;
	if (_AC) {
		offset = 1;
	}

	float dP = 0.1; // 10% of flexibility
	int nVoisin;
	float pLim1, pLim2, cost1, cost2, P0;
	float qLim1, qLim2, cost1Q, cost2Q, Q0, S0;
	
	int step = temporalStep;
	if (step > _PobjTemp.getNLin()) {
		std::cout << _PobjTemp.getNLin() << " " << step << std::endl;
		throw std::invalid_argument("temporel step too hight Pobj not defined");
	}

	for (int id = offset; id < _nCons; id++)
	{
		// P
		
		S0 = _factor.get(id - offset,0) * _PobjTemp.get(step, id - offset);
		//std::cout << S0 << " " << _factor.get(id, 0) << " " << _Pobj.get(id, step) << std::endl;
		if (_AC) {
			P0 = S0 * _PF.get(id - offset, 0);
		}
		else {
			P0 = S0;
		}
		_Pobj.set(id, 0, P0);
		pLim1 = -(1 + dP) * P0;
		pLim2 = -(1 - dP) * P0;
		cost1 = 1;
		cost2 = P0 * cost1;

		nVoisin = _nGen + _nPro;
		(_agents[id]).setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);
		_Ub.set(id, 0, 0);
		_Lb.set(id, 0, pLim1);

		_Pmin.set(id, 0, pLim1);
		_Pmax.set(id, 0, pLim2);
		_a.set(id, 0, cost1);
		_b.set(id, 0, cost2);
		_nVoisin.set(id, 0, nVoisin);

		if (_AC) {
			// Q
			Q0 = S0 * sqrt(1 - _PF.get(id - offset, 0) * _PF.get(id - offset, 0));
			_Pobj.set(id + _nAgent, 0, Q0);
			float randomFloat = rand1();
			//int signe = -1; // 2 * (randomFloat > 0.8) - 1; // 4 chance sur 5 d'�tre inductif (m�me signe que P)
			int signe = 2 * (id % 5) - 1; // un sur 5 ?
			Q0 = signe * Q0; 
			qLim1 = Q0 * (1 - dP + 2 * dP * (Q0 < 0));
			qLim2 = Q0 * (1 + dP - 2 * dP * (Q0 < 0));

			if (qLim1 == 0) {
				qLim1 = -1;
			} if (qLim2 == 0) {
				qLim2 = 1;
			}
			cost1Q = 1;
			cost2Q = Q0 * cost1Q;
			nVoisin = _nAgent - 1;
			_Ub.set(id + _nAgent, 0, qLim2 * (qLim2 > 0));
			_Lb.set(id + _nAgent, 0, qLim1 * (qLim1 < 0));

			_Pmin.set(id + _nAgent, 0, qLim1);
			_Pmax.set(id + _nAgent, 0, qLim2);
			_a.set(id + _nAgent, 0, cost1Q);
			_b.set(id + _nAgent, 0, cost2Q);
			_nVoisin.set(id + _nAgent, 0, nVoisin);
		}
	}
}

float StudyCaseAgent::rand1() const
{
	float a = (float)(rand()) / ((float)(RAND_MAX));
	return a;
}

int StudyCaseAgent::randab(int a, int b) const
{
	return a + (rand() % (b - a));
}

float StudyCaseAgent::randabfloat(float min, float max) const
{
	return min + rand1() * (max - min);
}

void StudyCaseAgent::initMat()
{
	DELETEA(_agents);
	_agents = new Agent[_nAgent];
	_a = MatrixCPU(_nAgent, 1);
	_b = MatrixCPU(_nAgent, 1);
	_Ub = MatrixCPU(_nAgent, 1);
	_Lb = MatrixCPU(_nAgent, 1);
	_Pobj = MatrixCPU(_nAgent, 1);
	_PobjD = MatrixCPUD(_nAgent, 1);
	_Pmin = MatrixCPU(_nAgent, 1);
	_Pmax = MatrixCPU(_nAgent, 1);
	_nVoisin = MatrixCPU(_nAgent, 1);

	_connect = MatrixCPU(_nAgent, _nAgent);
	_BETA = MatrixCPU(_nAgent, _nAgent);

	genConnec(&_connect);
}

void StudyCaseAgent::initMatAC()
{
	DELETEA(_agents);
	_AC = true;
	_agents = new Agent[_nAgent];
	_a = MatrixCPU(2 * _nAgent, 1, 0);
	_b = MatrixCPU(2 * _nAgent, 1);
	_Ub = MatrixCPU(2 * _nAgent, 1);
	_Lb = MatrixCPU(2 * _nAgent, 1);
	_Pobj = MatrixCPU(2 * _nAgent, 1);
	_PobjD = MatrixCPUD(2 * _nAgent, 1);
	_Pmin = MatrixCPU(2 * _nAgent, 1);
	_Pmax = MatrixCPU(2 * _nAgent, 1);
	_nVoisin = MatrixCPU(2 * _nAgent, 1);

//  Pour puissance active seulement
	_BETA = MatrixCPU(_nAgent, _nAgent);
	_connect = MatrixCPU(_nAgent, _nAgent);
	genConnec(&_connect); 
}

int StudyCaseAgent::getNFileline(std::string nameFile)
{
	int number_of_lines = 0;
	std::string line;
	std::ifstream myfile(nameFile);

	while (std::getline(myfile, line))
		++number_of_lines;
	return number_of_lines;
}



StudyCaseAgent::StudyCaseAgent()
{
	 _nAgent = 0;
	 _nPro = 0;
	 _nGen = 0;
	 _nCons = 0;
	 _nGenNFle = 0;
	 _timeInit = 0;
	 
}

StudyCaseAgent::StudyCaseAgent(int nAgent, float P, float dP, float a, float da, float b, float db, float propCons, float propPro)
{
	clock_t t = clock();
	srand(time(nullptr));
	if (propCons < 0 || propPro < 0 || (1 - propCons - propPro) < 0) {
		throw std::invalid_argument("propCons and propPro are proportion <1 and >0");
	}	
	_nAgent = nAgent;
	_nCons = nAgent * propCons;
	_nPro = nAgent * propPro;
	_nGen = nAgent - _nCons - _nPro;

	
	float pLim1,pLim2,cost1,cost2;
	initMat();
	genBetaUniforme(0);

	int nVoisin;
	
	for (int id = 0; id < nAgent; id++)
	{
		if (id < _nCons) { // consumer
			float P0 = -P + dP * 2 * (rand1() - 0.5);
			pLim1 = 1.1 * P0;
			pLim2 = 0.9 * P0;
			cost1 = a + da * 2 * (rand1() - 0.5);
			cost2 = -P0 * a;
			nVoisin = _nGen + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, nAgent, 1);
			_Ub.set(id, 0, 0);
			_Lb.set(id, 0, pLim1);
		}
		else if (id < (_nCons + _nGen)) { // generator
			pLim1 = 0;
			pLim2 = P + dP * 2 * (rand1() - 0.5);
			cost1 = a + da * 2 * (rand1() - 0.5);
			cost2 = b + db * 2 * (rand1() - 0.5);
			nVoisin = _nCons + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, nAgent, 2);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, 0);
		}
		else { // prosumer
			
			pLim1 = -P + dP * 2 * (rand1() - 0.5);
			pLim2 = P + dP * 2 * (rand1() - 0.5);
			cost1 = a + da * 2 * (rand1() - 0.5);
			cost2 = b + db * 2 * (rand1() - 0.5);
			nVoisin = _nGen + _nCons;
			
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, nAgent, 3);
			
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, pLim1);
		}
		_a.set(id, 0, cost1);
		_b.set(id, 0, cost2);

		_Pmin.set(id, 0, pLim1);
		_Pmax.set(id, 0, pLim2);
		_nVoisin.set(id, 0, nVoisin);
	}
	//
	
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;


}

StudyCaseAgent::StudyCaseAgent(int nAgent, float P, float dP, float Q, float dQ, float a, float da, float aQ, float daQ, float b, float db, float Gamma, float dGamma, float propCons, float propGenNFle, float propPro)
{
	srand(time(nullptr));
	if (propCons < 0 || propPro < 0 || propGenNFle < 0 || (1 - propCons - propPro - propGenNFle) < 0) {
		throw std::invalid_argument("propCons or propPro or propGenNFle are proportion <1 and >0");
	} if (dP > P || db > b || dQ > Q || dGamma > Gamma) {
		throw std::invalid_argument("variation must be smaller than average");
	}
	if (P < 0) {
		P = -P;
	}
	_AC = true;
	_nAgent = nAgent + 1;
	_nCons = nAgent * propCons + 1;
	_nPro = nAgent * propPro;
	_nGenNFle = nAgent * propGenNFle;
	_nGen = _nAgent - _nCons - _nPro;
	_nGenFle = _nGen - _nGenNFle;
	


	float pLim1, pLim2, cost1, cost2, costQ1, costQ2, qLim1, qLim2, P0, Q0, a0, aQ0;
	
	initMatAC();
	
	genBetaRandom(Gamma, dGamma);

	int nVoisin;
	
	(_agents[0]).setAgent(0, 0, 0, 0, 0, _nGen + _nPro, &_connect, _nAgent, 1);
	
	_Lb.set(0, 0, -10000 * P); // pour ne pas avoir besoin de le modifier
	_Ub.set(_nAgent, 0, 10000 * P); // idem
	_Lb.set(_nAgent, 0, -10000 * P);
	_nVoisin.set(0, 0, _nGen + _nPro);
	_nVoisin.set(_nAgent, 0, _nAgent - 1);

	for (int id = 1; id < _nAgent; id++)
	{
		if (id < _nCons) { // consumer
			float P0 = P + dP * 2 * (rand1() - 0.5);
			Q0 = dQ * 2 * (rand1() - 0.5);
			a0 = a + da * 2 * (rand1() - 0.5);
			aQ0 = a + da * 2 * (rand1() - 0.5);
			pLim1 = -1.1 * P0 / _Sbase;
			pLim2 = -0.9 * P0 / _Sbase;
			qLim1 = Q0 / _Sbase * (1.05 - 0.1 * (Q0 > 0));
			qLim2 = Q0 / _Sbase * (0.95 + 0.1 * (Q0 > 0));
			cost1 = a0 * (_Sbase * _Sbase);
			cost2 = P0 * a0 * _Sbase;
			costQ1 = aQ0 * (_Sbase * _Sbase);
			costQ2 = -Q0 * _Sbase * aQ0;
			nVoisin = _nGen + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);
			_Pobj.set(id, 0, -P0 / _Sbase);
			_PobjD.set(id, 0, -P0 / _Sbase);
			_Ub.set(id, 0, 0);
			_Lb.set(id, 0, pLim1);
		}
		else if (id < _nCons + _nGenNFle) {
			P0 = P + dP * 2 * (rand1() - 0.5);
			Q0 = dQ * 2 * (rand1() - 0.5);
			a0 = a + da * 2 * (rand1() - 0.5);
			aQ0 = a + da * 2 * (rand1() - 0.5);
			pLim1 = 0.9 * P0 / _Sbase;
			pLim2 = 1.1 * P0 / _Sbase;
			qLim1 = -dQ;
			qLim2 = dQ;
			cost1 = a0 * (_Sbase * _Sbase);
			cost2 = - P0 * _Sbase * a0;
			costQ1 = aQ0 * (_Sbase * _Sbase);
			costQ2 = -Q0 * _Sbase * aQ0;
			nVoisin = _nCons + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);
			_Pobj.set(id, 0, P0 / _Sbase);
			_PobjD.set(id, 0, P0 / _Sbase);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, 0);
		}
		else if (id < _nCons + _nGen)// generator
		{ 
			P0 = P + dP * 2 * (rand1() - 0.5);
			Q0 = dQ * 2 * (rand1() - 0.5);
			a0 = a + da * 2 * (rand1() - 0.5);
			aQ0 = a + da * 2 * (rand1() - 0.5);
			pLim1 = 0;
			pLim2 = P0 / _Sbase;
			qLim1 = Q0 / _Sbase * (1.05 - 0.1 * (Q0 > 0));
			qLim2 = Q0 / _Sbase * (0.95 + 0.1 * (Q0 > 0));
			
			cost1 = a0 * (_Sbase * _Sbase);
			cost2 = b * _Sbase + db * 2 * (rand1() - 0.5) * _Sbase;
			costQ1 = aQ0 * (_Sbase * _Sbase);
			costQ2 = -Q0 * _Sbase * aQ0;
			nVoisin = _nCons + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 2);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, 0);
			_Pobj.set(id, 0, P0 / _Sbase);
			_PobjD.set(id, 0, P0 / _Sbase);
		}
		else { // prosumer
			P0 = dP * 2 * (rand1() - 0.5); // positive or negative
			Q0 = dQ * 2 * (rand1() - 0.5); // positive or negative
			a0 = a + da * 2 * (rand1() - 0.5);
			aQ0 = a + da * 2 * (rand1() - 0.5);

			pLim1 = -dP;
			pLim2 = dP;
			qLim1 = -dQ;
			qLim2 = dQ;
			cost1 = a0 * (_Sbase * _Sbase);
			cost2  = -P0 * _Sbase * a0;
			costQ1 = aQ0 * (_Sbase * _Sbase);
			costQ2 = -Q0 * _Sbase * aQ0;
			nVoisin = _nGen + _nCons + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);
			_Pobj.set(id, 0, P0 / _Sbase);
			_PobjD.set(id, 0, P0 / _Sbase);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, pLim1);
		}
		
		_PobjD.set(id + _nAgent, 0, Q0 / _Sbase);
		_Pobj.set(id + _nAgent, 0, Q0 / _Sbase);
		_Ub.set(id + _nAgent, 0, qLim2 * (qLim2 > 0));
		_Lb.set(id + _nAgent, 0, qLim1 * (qLim1 < 0));

		_a.set(id, 0, cost1);
		_b.set(id, 0, cost2);
		_a.set(id + _nAgent, 0, costQ1);
		_b.set(id + _nAgent, 0, costQ2);

		_Pmin.set(id, 0, pLim1);
		_Pmax.set(id, 0, pLim2);
		_Pmin.set(id + _nAgent, 0, qLim1);
		_Pmax.set(id + _nAgent, 0, qLim2);
		_nVoisin.set(id, 0, nVoisin);
		_nVoisin.set(id + _nAgent, 0, _nAgent - 1);
	}

}

StudyCaseAgent::StudyCaseAgent(int nAgent, float P, float dP, float a, float da, float b, float db, float Gamma, float dGamma, float propCons, float propGenNFle, float propPro)
{
	srand(time(nullptr));
	if (propCons < 0 || propPro < 0 || propGenNFle < 0 || (1 - propCons - propPro - propGenNFle) < 0) {
		throw std::invalid_argument("propCons or propPro or propGenNFle are proportion <1 and >0");
	} if (dP > P || db > b || dGamma > Gamma) {
		throw std::invalid_argument("variation must be smaller than average");
	}
	if (P < 0) {
		P = -P;
	}
	_AC = false;
	_nAgent = nAgent ;
	_nCons = nAgent * propCons;
	_nPro = nAgent * propPro;
	_nGenNFle = nAgent * propGenNFle;
	_nGen = _nAgent - _nCons - _nPro;
	_nGenFle = _nGen - _nGenNFle;



	float pLim1, pLim2, cost1, cost2, P0, a0;

	initMat();

	genBetaRandom(Gamma, dGamma);

	int nVoisin;

	for (int id = 0; id < _nAgent; id++)
	{
		if (id < _nCons) { // consumer
			P0 = P + dP * 2 * (rand1() - 0.5);
			a0 = a + da * 2 * (rand1() - 0.5);
			pLim1 = -1.1 * P0 / _Sbase;
			pLim2 = -0.9 * P0 / _Sbase;
			cost1 = a0 * (_Sbase * _Sbase);
			cost2 = P0 * a0 * _Sbase;
			nVoisin = _nGen + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);
			_Pobj.set(id, 0, -P0 / _Sbase);
			_PobjD.set(id, 0, -P0 / _Sbase);
			_Ub.set(id, 0, 0);
			_Lb.set(id, 0, pLim1);
		}
		else if (id < _nCons + _nGenNFle) {
			P0 = P + dP * 2 * (rand1() - 0.5);
			a0 = a + da * 2 * (rand1() - 0.5);
			pLim1 = 0.9 * P0 / _Sbase;
			pLim2 = 1.1 * P0 / _Sbase;
			cost1 = a0 * (_Sbase * _Sbase);
			cost2 = -P0 * _Sbase * a0;
			nVoisin = _nCons + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);
			_Pobj.set(id, 0, P0 / _Sbase);
			_PobjD.set(id, 0, P0 / _Sbase);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, 0);
		}
		else if (id < _nCons + _nGen)// generator
		{
			P0 = P + dP * 2 * (rand1() - 0.5);
			a0 = a + da * 2 * (rand1() - 0.5);
			pLim1 = 0;
			pLim2 = P0 / _Sbase;

			cost1 = a0 * (_Sbase * _Sbase);
			cost2 = b * _Sbase + db * 2 * (rand1() - 0.5) * _Sbase;
			nVoisin = _nCons + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 2);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, 0);
			_Pobj.set(id, 0, P0 / _Sbase);
			_PobjD.set(id, 0, P0 / _Sbase);
		}
		else { // prosumer
			P0 = dP * 2 * (rand1() - 0.5); // positive or negative
			a0 = a + da * 2 * (rand1() - 0.5);
			
			pLim1 = -dP;
			pLim2 = dP;
	
			cost1 = a0 * (_Sbase * _Sbase);
			cost2 = -P0 * _Sbase * a0;
			
			nVoisin = _nGen + _nCons + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);
			_Pobj.set(id, 0, P0 / _Sbase);
			_PobjD.set(id, 0, P0 / _Sbase);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, pLim1);
		}

	

		_a.set(id, 0, cost1);
		_b.set(id, 0, cost2);
	
		_Pmin.set(id, 0, pLim1);
		_Pmax.set(id, 0, pLim2);
		
		_nVoisin.set(id, 0, nVoisin);
	}

}


StudyCaseAgent::StudyCaseAgent(int nAgent, float P0, float dP, float b, float db, float propCons)
{
	clock_t t = clock();
	srand(time(nullptr));
	if (propCons < 0 || (1 - propCons) < 0) {
		throw std::invalid_argument("propCons is a proportion <1 and >0");
	}
	_nAgent = nAgent;
	_nCons = nAgent * propCons;
	_nPro = 0;
	_nGen = nAgent - _nCons - _nPro;

	float pLim1, pLim2, cost1, cost2;
	initMat();
	genBetaUniforme(0);

	int nVoisin;

	bool impossible = true;

	while (impossible)
	{
		for (int id = 0; id < nAgent; id++)
		{
			if (id < _nCons) { // consumer
				pLim1 = -P0 - dP * (rand1() + 0.01);
				pLim2 = -P0 + dP * (rand1() + 0.01);
				cost1 = 1;
				cost2 = P0;
				nVoisin = _nGen + _nPro;
				_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, nAgent, 1);
				_Ub.set(id, 0, 0);
				_Lb.set(id, 0, pLim1);
			}
			else { // generator
				pLim1 = 0;
				pLim2 = P0 + dP * 2 * (rand1() - 0.5);
				cost1 = 0.1;
				cost2 = b + db * 2 * (rand1() - 0.5);
				nVoisin = _nCons + _nPro;
				_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, nAgent, 2);
				_Ub.set(id, 0, pLim2);
				_Lb.set(id, 0, 0);
			}
			_a.set(id, 0, cost1);
			_b.set(id, 0, cost2);

			_Pmin.set(id, 0, pLim1);
			_Pmax.set(id, 0, pLim2);
			_nVoisin.set(id, 0, nVoisin);
		}
		impossible = false;
	}
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}

StudyCaseAgent::StudyCaseAgent(const StudyCaseAgent& s)
{
	clock_t t = clock();
	_nAgent = s._nAgent;
	_nPro = s._nPro;
	_nGen = s._nGen;
	_nCons = s._nCons;
	_nGenNFle = s._nGenNFle;
	_Sbase = s._Sbase;
	temporalStep = s.temporalStep;
	_AC = s._AC;
	_name = s._name;

	_agents = new Agent[_nAgent];
	for (int i = 0; i < _nAgent; i++) {
		_agents[i] = s._agents[i];
	}

	_a = s._a;
	_b = s._b;
	_Ub = s._Ub;
	_connect =s._connect;
	_Lb = s._Lb;
	_Pmin = s._Pmin;
	_Pmax = s._Pmax;
	_nVoisin = s._nVoisin;
	_BETA = s._BETA;
	
	_Pobj = s._Pobj;
	_PobjD = s._PobjD;
	_PobjTemp = s._PobjTemp;
	_factor = s._factor;
	_PF = s._PF;
	
	

	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}

StudyCaseAgent::StudyCaseAgent(std::string fileName)
{
	clock_t t = clock();
	std::ifstream myFile(fileName, std::ios::in);
	std::cout << fileName << std::endl;
	if (myFile)
	{
		myFile >> _nAgent >> _nCons >> _nGen >> _nPro;
		int nbLineToRead = 7 + 2 * _nAgent;
		_agents = new Agent[_nAgent];
		_a = MatrixCPU(_nAgent, 1);
		_b = MatrixCPU(_nAgent, 1);
		_Lb = MatrixCPU(_nAgent, 1);
		_Ub = MatrixCPU(_nAgent, 1);
		_Pmin = MatrixCPU(_nAgent, 1);
		_Pmax = MatrixCPU(_nAgent, 1);
		_nVoisin = MatrixCPU(_nAgent, 1);
		_BETA = MatrixCPU(_nAgent, _nAgent);
		_connect = MatrixCPU(_nAgent, _nAgent);
		float c = 0;

		for (int i = 0; i < nbLineToRead; i++)
		{
			switch (i)
			{
			case 0:
				for (int j = 0; j < _nAgent; j++) {
					myFile >> c;
					_a.set(j, 0, c);
				}
				break;
			case 1:
				for (int j = 0; j < _nAgent; j++) {
					myFile >> c;
					_b.set(j, 0, c);
				}
				break;
			case 2:
				for (int j = 0; j < _nAgent; j++) {
					myFile >> c;
					_Lb.set(j, 0, c);
				}
				break;
			case 3:
				for (int j = 0; j < _nAgent; j++) {
					myFile >> c;
					_Ub.set(j, 0, c);
				}
				break;
			case 4:
				for (int j = 0; j < _nAgent; j++) {
					myFile >> c;
					_Pmin.set(j, 0, c);
				}
				break;
			case 5:
				for (int j = 0; j < _nAgent; j++) {
					myFile >> c;
					_Pmax.set(j, 0, c);
				}
				break;
			case 6:
				for (int j = 0; j < _nAgent; j++) {
					myFile >> c;
					_nVoisin.set(j, 0, c);
				}
				break;
			default:
				if (i < (6 + _nAgent)) {
					for (int j = 0; j < _nAgent; j++) {
						myFile >> c;
						_BETA.set(i, j, c);
					}
				}
				else {
					for (int j = 0; j < _nAgent; j++) {
						myFile >> c;
						_connect.set(i, j, c);
					}
				}
				break;
			}

		}
		myFile.close();

		for (int id = 0; id < _nAgent; id++)
		{
			if (id < _nCons) { // consumer
				_agents[id].setAgent(id, _Pmin.get(id, 0), _Pmax.get(id, 0), _a.get(id, 0), _b.get(id, 0), _nVoisin.get(id, 0), &_connect, _nAgent, 1);

			}
			else if (id < (_nCons + _nGen)) { // generator
				_agents[id].setAgent(id, _Pmin.get(id, 0), _Pmax.get(id, 0), _a.get(id, 0), _b.get(id, 0), _nVoisin.get(id, 0), &_connect, _nAgent, 2);
			}
			else { // prosumer
				_agents[id].setAgent(id, _Pmin.get(id, 0), _Pmax.get(id, 0), _a.get(id, 0), _b.get(id, 0), _nVoisin.get(id, 0), &_connect, _nAgent, 3);
			}
		}
	}
	else {
		throw std::invalid_argument("can't open this file");
	}
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;

}

StudyCaseAgent& StudyCaseAgent::operator= (const StudyCaseAgent& s)
{
	clock_t t = clock();
	_nAgent = s._nAgent;
	_nPro = s._nPro;
	_nGen = s._nGen;

	_nCons = s._nCons;
	
	_name = s._name;
	_nGenNFle = s._nGenNFle; 
	_nGenFle = s._nGenFle;

	temporalStep = s.temporalStep;
	_AC = s._AC;


	DELETEA(_agents);
	_agents = new Agent[_nAgent];
	for (int i = 0; i < _nAgent; i++) {
		_agents[i] = s._agents[i];
	}

	_a = s._a;
	_b = s._b;
	_Ub = s._Ub;
	_connect = s._connect;
	_Lb = s._Lb;
	_Pmin = s._Pmin;
	_Pmax = s._Pmax;
	_nVoisin = s._nVoisin;
	_BETA = s._BETA;

	_Pobj = s._Pobj;
	_PobjD = s._PobjD;
	_PobjTemp = s._PobjTemp;
	_factor = s._factor;
	_PF = s._PF;
	_Sbase = s._Sbase;

	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
	return *this;
}



void StudyCaseAgent::genAgents(int nAgent, float propCons, float Pconso, float dPconso, float bProd, float dbProd, float Pprod, float dPprod, float Gamma, float dGamma)
{
	clock_t t = clock();
	srand(time(nullptr));
	if (propCons < 0 || (1 - propCons) < 0) {
		throw std::invalid_argument("propCons is a proportion <1 and >0");
	} if (dPconso > Pconso || dbProd > bProd || dPprod > Pprod || dGamma > Gamma) {
		throw std::invalid_argument("variation must be smaller than average");
	}
	if (Pconso < 0) {
		Pconso = -Pconso;
	}
	_nAgent = nAgent;
	_nCons = nAgent * propCons;
	_nPro = 0;
	_nGen = nAgent - _nCons - _nPro;

	_agents = new Agent[nAgent];
	float pLim1;
	float pLim2;
	float cost1;
	float cost2;

	_a = MatrixCPU(nAgent, 1);
	_b = MatrixCPU(nAgent, 1);
	_Ub = MatrixCPU(nAgent, 1);
	_connect = MatrixCPU(nAgent, nAgent);
	_Lb = MatrixCPU(nAgent, 1);
	_Pmin = MatrixCPU(nAgent, 1);
	_Pmax = MatrixCPU(nAgent, 1);
	_nVoisin = MatrixCPU(nAgent, 1);
	_BETA = MatrixCPU(_nAgent, _nAgent);

	genConnec(&_connect);
	float gamma = Gamma + dGamma * 2 * (rand1() - 0.5);
	genBetaUniforme(gamma);

	int nVoisin;

	bool impossible = true;

	while (impossible)
	{
		for (int id = 0; id < nAgent; id++)
		{
			if (id < _nCons) { // consumer
				float P0 = Pconso + dPconso * 2 * (rand1() - 0.5);
				pLim1 = -1.1 * P0;
				pLim2 = -0.9 * P0;
				cost1 = 1;
				cost2 = P0;
				nVoisin = _nGen + _nPro;
				_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, nAgent, 1);
				_Ub.set(id, 0, 0);
				_Lb.set(id, 0, pLim1);
			}
			else { // generator
				float P0 = Pprod + dPprod * 2 * (rand1() - 0.5);
				pLim1 = 0;
				pLim2 = P0;
				cost1 = 0.1;
				cost2 = bProd + dbProd * 2 * (rand1() - 0.5);
				nVoisin = _nCons + _nPro;
				_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, nAgent, 2);
				_Ub.set(id, 0, pLim2);
				_Lb.set(id, 0, 0);
			}
			_a.set(id, 0, cost1);
			_b.set(id, 0, cost2);

			_Pmin.set(id, 0, pLim1);
			_Pmax.set(id, 0, pLim2);
			_nVoisin.set(id, 0, nVoisin);
		}
		impossible = false;
	}
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;



}

void StudyCaseAgent::genAgentsAC(int nAgent, float propCons, float propGenNFle, float Pconso, float dPconso, float bProd, float dQconso, float dbProd, float Pprod, float dPprod, float Gamma, float dGamma)
{

	
	if (propCons < 0 || propGenNFle < 0 || (1 - propCons - propGenNFle) < 0) {
		throw std::invalid_argument("propCons or propGenNFle are proportion <1 and >0");
	} if (dPconso > Pconso || dbProd > bProd || dPprod > Pprod || dGamma > Gamma) {
		throw std::invalid_argument("variation must be smaller than average");
	}
	if (Pconso < 0) {
		Pconso = -Pconso;
	}
	
	_nAgent = nAgent + 1;
	_nCons = nAgent * propCons + 1;
	_nPro = 0;
	_nGenNFle = nAgent * propGenNFle;
	_nGen = _nAgent - _nCons - _nPro;
	_nGenFle = _nGen - _nGenNFle;

	//std::cout << _nAgent << " " << _nCons << " " << _nGen << std::endl;
	
	float pLim1, pLim2, cost1, cost2, costQ1, costQ2, qLim1, qLim2, Q0;
		
	initMatAC();
	float gamma = Gamma + dGamma * 2 * (rand1() - 0.5);
	genBetaUniforme(gamma);

	int nVoisin;

	bool impossible = true;
	
	(_agents[0]).setAgent(0, 0, 0, 0, 0, _nGen + _nPro, &_connect, _nAgent, 1);

	_Lb.set(0, 0, -10000); // pour ne pas avoir besoin de le modifier
	_Ub.set(_nAgent, 0, 10000); // idem
	_Lb.set(_nAgent, 0, -10000);
	_nVoisin.set(0, 0, _nGen);
	_nVoisin.set(_nAgent, 0, _nAgent - 1);

	//std::cout << "set carac agent" << std::endl;
	for (int id = 1; id < _nAgent; id++)
	{
		if (id < _nCons) { // consumer
			float P0 = Pconso + dPconso * 2 * (rand1() - 0.5);
			Q0 = dQconso * 2 * (rand1() - 0.5);
			pLim1 = -1.1 * P0 / _Sbase;
			pLim2 = -0.9 * P0 / _Sbase;
			qLim1 = Q0 / _Sbase * (1.05 - 0.1 * (Q0 > 0));
			qLim2 = Q0 / _Sbase * (0.95 + 0.1 * (Q0 > 0));
			cost1 = 1 * (_Sbase * _Sbase);
			cost2 = P0 * cost1 /_Sbase;
			costQ1 = 0.1 * (_Sbase * _Sbase);
			costQ2 = -Q0 * costQ1 / _Sbase;
			nVoisin = _nGen + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);
			_Pobj.set(id, 0, -P0 / _Sbase);
			_PobjD.set(id, 0, -P0 / _Sbase);
			_Ub.set(id, 0, 0);
			_Lb.set(id, 0, pLim1);

		}
		else if (id < _nCons + _nGenNFle) {
			float P0 = Pconso + dPconso * 2 * (rand1() - 0.5);
			Q0 = dQconso * 2 * (rand1() - 0.5);
			pLim1 = 0.9 * P0 / _Sbase;
			pLim2 = 1.1 * P0 / _Sbase;
			qLim1 = - 2 * dQconso;
			qLim2 =   2 * dQconso;
			cost1 = 0.1 * (_Sbase * _Sbase);
			cost2 = - P0 * cost1 / _Sbase;
			costQ1 = 0.1 * (_Sbase * _Sbase);
			costQ2 = -Q0 * costQ1 / _Sbase;
			nVoisin = _nCons + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);
			_Pobj.set(id, 0, P0 / _Sbase);
			_PobjD.set(id, 0, P0 / _Sbase);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, 0);
		}
		else { // generator
			float P0 = Pconso + dPconso * 2 * (rand1() - 0.5);
			float P02 = Pprod + dPprod * 2 * (rand1() - 0.5);
			if (P0 < P02) {
				P0 = P02;
			}
			Q0 = dQconso * 2 * (rand1() - 0.5);
			pLim1 = 0;
			pLim2 = P0 / _Sbase;
			qLim1 = Q0 / _Sbase * (1.05 - 0.1 * (Q0 > 0));
			qLim2 = Q0 / _Sbase * (0.95 + 0.1 * (Q0 > 0));
			cost1 = 0.1 * (_Sbase * _Sbase);
			cost2 = bProd * _Sbase + dbProd * 2 * (rand1() - 0.5) * _Sbase;
			costQ1 = 0.1 * (_Sbase * _Sbase);
			costQ2 = -Q0 * costQ1 / _Sbase;
			nVoisin = _nCons + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 2);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, 0);
			_Pobj.set(id, 0, P0 / _Sbase);
			_PobjD.set(id, 0, P0 / _Sbase);
		}
		_PobjD.set(id + _nAgent, 0, Q0 / _Sbase);
		_Pobj.set(id + _nAgent, 0, Q0 / _Sbase);
		_Ub.set(id + _nAgent, 0, qLim2 * (qLim2 > 0));
		_Lb.set(id + _nAgent, 0, qLim1 * (qLim1 < 0));

		_a.set(id, 0, cost1);
		_b.set(id, 0, cost2);
		_a.set(id + _nAgent, 0, costQ1);
		_b.set(id + _nAgent, 0, costQ2);

		_Pmin.set(id, 0, pLim1);
		_Pmax.set(id, 0, pLim2);
		_Pmin.set(id + _nAgent, 0, qLim1);
		_Pmax.set(id + _nAgent, 0, qLim2);
		_nVoisin.set(id, 0, nVoisin);
		_nVoisin.set(id + _nAgent, 0, _nAgent - 1);
	}
	//std::cout << "fin carac agent" << std::endl;
}

void StudyCaseAgent::genAgentsFullRandom(int nAgent, float aMin, float aMax, float P0Min, float P0Max, float gammaMin, float gammaMax, float propConsoMin, float propConsoMax, float borneMin, float borneMax)
{
	clock_t t = clock();
	srand(time(nullptr));
	if (nAgent <= 0) {
		throw std::invalid_argument("nAgent must be positive");
	}
	if (aMin <= 0 || aMax < aMin) {
		throw std::invalid_argument("invalide value for aMin or aMax");
	}
	if (P0Max < P0Min || gammaMax < gammaMin || propConsoMax < propConsoMin || borneMax < borneMin || P0Min * P0Max < 0 || gammaMax * gammaMin < 0 || borneMin * borneMax < 0) {
		throw std::invalid_argument("yMax must be greater than yMin and same signe");
	}
	if (propConsoMin <= 0 || propConsoMax >= 1) {
		throw std::invalid_argument("proportion must be between 0 and 1 excluded");
	}
	if (borneMin < 0 ) {
		throw std::invalid_argument("borneMin must be positive");
	}
	if (P0Min < 0) {
		P0Min = -P0Min;
		P0Max = -P0Max;
	}

	_nAgent = nAgent;
	float propConso = randabfloat(propConsoMin, propConsoMax);
	_nCons = Mymin(Mymax(propConso, 1), _nAgent - 1);
	_nPro = 0;
	_nGen = _nAgent - _nCons;
	_nGenFle = _nGen;
	_nGenNFle = 0;

	initMat();

	float gammaMoy = (gammaMax + gammaMin) / 2;
	float dgamma = gammaMax - gammaMoy;
	
	genBetaRandom(gammaMoy, dgamma);

	int nVoisin;
	float pLim1, pLim2, cost1, cost2;

	
	for (int id = 0; id < nAgent; id++)
	{
		float P0 = randabfloat(P0Min, P0Max); // objective
		float borne = randabfloat(borneMin, borneMax);
		float cost1 = randabfloat(aMin, aMax);
		if (id < _nCons) { // consumer
				
			pLim1 =  - (1 + borne) * P0;
			pLim2 = Mymin((1 - borne) * (-P0), 0);
			cost2 =  P0 / cost1;
			nVoisin = _nGen + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, nAgent, 1);
			_Ub.set(id, 0, 0);
			_Lb.set(id, 0, pLim1);
		}
		else { // generator
			pLim1 = Mymax((1 - borne) * (P0), 0);
			pLim2 = (1 + borne) * P0;
			cost2 = -P0 / cost1;
			nVoisin = _nCons + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, nAgent, 2);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, 0);
		}
		_a.set(id, 0, cost1);
		_b.set(id, 0, cost2);

		_Pmin.set(id, 0, pLim1);
		_Pmax.set(id, 0, pLim2);
		_nVoisin.set(id, 0, nVoisin);
	}
	
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;

}





void StudyCaseAgent::nextStepPobj()
{
	temporalStep++;
	initCaseFromPobj();
}




void StudyCaseAgent::genConnec(MatrixCPU* connec)
{
	
#ifdef DEBUG_CONSTRUCTOR
	std::cout << " nPro =" << nPro << " nGen =" << nGen << " nCons =" << nCons << " nAgent =" << nAgent << std::endl;
	std::cout << " row =" << connec->getNLin() << " column =" << connec->getNCol() << std::endl;
#endif
	
	for (int i = 0; i < _nCons; i++) { // conso
		for (int j = _nCons; j < _nAgent;j++) {
			connec->set(i, j, 1); // gen and pro
		}
	}
	for (int i = _nCons; i < _nCons + _nGen; i++) { // gen
		for (int j = 0; j < _nCons;j++) { // conso
			connec->set(i, j, 1);
		}
		for (int j = _nCons + _nGen; j < _nAgent;j++) { //pro
			connec->set(i, j, 1); 
		}
	}
	for (int i = _nCons + _nGen; i < _nAgent;i++) { //pro
		for (int j = 0; j < _nCons + _nGen; j++) { // pro echange pas entre pro 
			if (i != j) {
				connec->set(i, j, 1);
			}
		}
	}
}

void StudyCaseAgent::genBetaUniforme(float beta)
{
	// Cons de 0 � nCons exclus
	// Prod de nCons � nCons+nGen exclus
	// Pro de nCons+nGen � nAgent exclus
	int offset = 0;
	if (_AC) {
		offset = 1;
	}

	for (int i = offset; i < _nCons; i++) { // conso
		for (int j = _nCons; j < _nAgent; j++) { // le reste
			_BETA.set(i, j, -beta *_Sbase);
		}
	}
	for (int i = _nCons; i < _nCons + _nGen; i++) { // gen
		for (int j = offset; j < _nCons; j++) { // conso
			_BETA.set(i, j, beta * _Sbase);
		}
		for (int j = _nCons + _nGen; j < _nAgent; j++) { // prosumer
			_BETA.set(i, j, beta * _Sbase);
		}
	}
	for (int i = _nCons + _nGen; i < _nAgent; i++) { // prosumer
		for (int j = offset; j < _nCons; j++) { // conso
			_BETA.set(i, j, beta * _Sbase);
		}
		for (int j = _nCons; j < _nCons + _nGen; j++) { // gen
			_BETA.set(i, j, -beta * _Sbase);
		}
	}

}

void StudyCaseAgent::genBetaDistance(float s, MatrixCPU* Distance)
{
	for (int i = 0; i < _nCons; i++) { // conso
		for (int j = _nCons; j < _nAgent; j++) { // le reste
			_BETA.set(i, j, -s * Distance->get(i, j) / 2);
		}
	}
	for (int i = _nCons; i < _nCons + _nGen; i++) { // gen
		for (int j = 0; j < _nCons; j++) { // conso
			_BETA.set(i, j, s * Distance->get(i, j) / 2);
		}
		for (int j = _nCons + _nGen; j < _nAgent; j++) { // prosumer
			_BETA.set(i, j, s * Distance->get(i, j) / 2);
		}
	}
	for (int i = _nCons + _nGen; i < _nAgent; i++) { // prosumer
		for (int j = 0; j < _nCons; j++) { // conso
			_BETA.set(i, j, s * Distance->get(i, j) / 2);
		}
		for (int j = _nCons; j < _nCons + _nGen; j++) { // gen
			_BETA.set(i, j, -s * Distance->get(i, j) / 2);
		}
	}
	_BETA.multiply(_Sbase);
}

void StudyCaseAgent::genBetaZone(MatrixCPU* Mats, MatrixCPU* Distance, MatrixCPU zones)
{
	std::cout << "work in progress" << std::endl;
	throw std::invalid_argument("WIP");
	/*int maxZone = _zoneBus.max2();
	if (maxZone >= Mats->getNCol() || maxZone >= Mats->getNLin()) {
		throw std::invalid_argument("Mats has not the a good size");
	}

	for (int i = 0; i < _nCons; i++) { // conso
		for (int j = _nCons; j < _nAgent; j++) { // le reste
			int zonei = _zoneBus.get(i, 0);
			int zonej = _zoneBus.get(j, 0);
			float s = Mats->get(zonei, zonej);
			_BETA.set(i, j, -s * _Distance.get(i, j) / 2);
		}
	}
	for (int i = _nCons; i < _nCons + _nGen; i++) { // gen
		for (int j = 0; j < _nCons; j++) { // conso
			int zonei = _zoneBus.get(i, 0);
			int zonej = _zoneBus.get(j, 0);
			float s = Mats->get(zonei, zonej);
			_BETA.set(i, j, s * _Distance.get(i, j) / 2);
		}
		for (int j = _nCons + _nGen; j < _nAgent; j++) { // prosumer
			int zonei = _zoneBus.get(i, 0);
			int zonej = _zoneBus.get(j, 0);
			float s = Mats->get(zonei, zonej);
			_BETA.set(i, j, s * _Distance.get(i, j) / 2);
		}
	}
	for (int i = _nCons + _nGen; i < _nAgent; i++) { // prosumer
		for (int j = 0; j < _nCons; j++) { // conso
			int zonei = _zoneBus.get(i, 0);
			int zonej = _zoneBus.get(j, 0);
			float s = Mats->get(zonei, zonej);
			_BETA.set(i, j, s * _Distance.get(i, j) / 2);
		}
		for (int j = _nCons; j < _nCons + _nGen; j++) { // gen
			int zonei = _zoneBus.get(i, 0);
			int zonej = _zoneBus.get(j, 0);
			float s = Mats->get(zonei, zonej);
			_BETA.set(i, j, -s * _Distance.get(i, j) / 2);
		}
	}*/
}

void StudyCaseAgent::genBetaRandom(float beta, float dbeta)
{
	for (int i = 0; i < _nCons; i++) { // conso
		for (int j = _nCons; j < _nAgent; j++) { // le reste
			float betaR = beta + dbeta * 2 * (rand1() - 0.5);
			_BETA.set(i, j, -betaR * _Sbase);
		}
	}
	for (int i = _nCons; i < _nCons + _nGen; i++) { // gen
		for (int j = 0; j < _nCons; j++) { // conso
			float betaR = beta + dbeta * 2 * (rand1() - 0.5);
			_BETA.set(i, j, betaR * _Sbase);
		}
		for (int j = _nCons + _nGen; j < _nAgent; j++) { // prosumer
			float betaR = beta + dbeta * 2 * (rand1() - 0.5);
			_BETA.set(i, j, betaR * _Sbase);
		}
	}
	for (int i = _nCons + _nGen; i < _nAgent; i++) { // prosumer
		for (int j = 0; j < _nCons; j++) { // conso
			float betaR = beta + dbeta * 2 * (rand1() - 0.5);
			_BETA.set(i, j, betaR * _Sbase);
		}
		for (int j = _nCons; j < _nCons + _nGen; j++) { // gen
			float betaR = beta + dbeta * 2 * (rand1() - 0.5);
			_BETA.set(i, j, -betaR * _Sbase);
		}
	}

}


MatrixCPU StudyCaseAgent::Set29node(bool AC)
{
	clock_t t = clock();
	_Sbase = 1; // MW car matlab
	int offset = 0;
	if (AC) {
		offset = 1;
	}

	_nAgent = 31 + offset;
	_nPro = 0;
	_nGen = 10;
	_nCons = 21 + offset;
	
	float Plim1[31] = { -146.4, -483, -750, -350.7, -783, -9.8, -12.8, -480, -493.5, -237, -1020, -411, -371.3, -462.9, -336, -208.5, -421.5, -309, -425.3,-13.8, -1656 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	float Plim2[31] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1040, 646, 725, 652, 508, 687, 580, 564, 865, 1100 };
	float Cost1[31] = { 67, 47, 47, 53, 82, 52, 87, 57, 50, 52, 71, 64, 57, 82, 69, 69, 86, 54, 78, 81, 59, 89, 67, 55, 82, 88, 76, 84, 77, 51, 87 };
	float Cost2[31] = { 64, 79, 71, 62, 65, 83, 63, 81, 73, 69, 62, 79, 60, 80, 78, 70, 62, 70, 66, 70, 71, 18, 21, 37, 25, 17, 38, 28, 36, 38, 19 };
	float Qobj[31] = { -14, -7, -14, 12, 9.7, 0.58, -0.54, 13, -14, 1.8, 3.5, 8, 8.9, -9.4, -0.3, -1.6, 4.4, 6.2, 7.6, -6.7, 5.4, 4.6, -10, -11, -0.05, 13.8, -4.7, 2.5, 8.3, 7.54, -7.4 };
	int BusAgent[31] = { 1, 3, 4, 7, 8, 9, 12, 15, 16, 18, 20, 21, 23, 24, 25, 26, 27, 28, 29, 31, 39, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39 };
	
	float gamma = 0;
	
	
	MatrixCPU coresBusAgentLin = MatrixCPU(_nAgent, 1);

	for (int i = offset; i < _nAgent; i++) {
		coresBusAgentLin.set(i, 0, BusAgent[i - offset] - 1);
	}

	for (int i = 0; i < _nAgent; i++) {
		Plim1[i] = Plim1[i] / _Sbase;
		Plim2[i] = Plim2[i] / _Sbase;
		Cost1[i] = Cost1[i] * _Sbase * _Sbase;
		Cost2[i] = Cost2[i] * _Sbase;
		Qobj[i] = Qobj[i] / _Sbase;
	}
	gamma = gamma * _Sbase;

	float pLim1, pLim2, cost1, cost2, costQ1, costQ2, qLim1, qLim2;
	
	if (AC) {
		initMatAC();
	}
	else {
		initMat();
	}
	genBetaUniforme(gamma);
	int nVoisin;
	
	if (AC) { // agent des pertes
		pLim1 = 0; // sera modifi� par l'algo selon le calcul des pertes
		pLim2 = 0; // idem
		cost1 = 0; // cela ne cause aucune division par 0 et la non convexit� ne joue pas comme on a qu'un seul point d�fini
		cost2 = 0; // idem
		nVoisin = _nGen + _nPro;

		qLim1 = 0; // idem que puissance active
		qLim2 = 0;

		costQ1 = 0;
		costQ2 = 0;

		(_agents[0]).setAgent(0, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);

		_Ub.set(0, 0, 0);
		_Lb.set(0, 0, -FLT_MAX); // pour ne pas avoir besoin de le modifier

		_Ub.set(0 + _nAgent, 0, FLT_MAX); // idem
		_Lb.set(0 + _nAgent, 0, -FLT_MAX);
	}

	int AgentIndice = 0;
	for (int id = 0; id < _nAgent-offset; id++)
	{
		AgentIndice = id + offset;
		if (id < _nCons - offset) { // consumer
			
			pLim1 = Plim1[id];
			pLim2 = Plim2[id];
			cost1 = Cost1[id] / 1000;
			cost2 = Cost2[id] ;
			nVoisin = _nGen + _nPro;
			
			qLim1 = Qobj[id] * (0.95 + 0.1 * (Qobj[id] < 0));
			qLim2 = Qobj[id] * (1.05 - 0.2 * (Qobj[id] < 0));
			if (qLim1 == 0) {
				qLim1 = -0.1;
			} if (qLim2 == 0) {
				qLim2 = 0.1;
			}
			costQ1 = 0.1;
			costQ2 = -costQ1 * Qobj[id];

			(_agents[AgentIndice]).setAgent(AgentIndice, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);
	
			_Ub.set(AgentIndice, 0, 0);
			_Lb.set(AgentIndice, 0, pLim1);
			if (AC) {
				_Ub.set(AgentIndice + _nAgent, 0, qLim2 * (qLim2 > 0));
				_Lb.set(AgentIndice + _nAgent, 0, qLim1 * (qLim1 < 0));
			}
			
			
		}
		else if (id < (_nCons + _nGen - offset)) { // generator
			pLim1 = Plim1[id];
			pLim2 = Plim2[id];
			cost1 = Cost1[id] / 1000;
			cost2 = Cost2[id] ;
			nVoisin = _nCons + _nPro;

			qLim1 = Qobj[id] * 2 * (Qobj[id] < 0);
			qLim2 = Qobj[id] * 2 * (Qobj[id] > 0);

			costQ1 = 0.1;
			costQ2 = -costQ1 * Qobj[id];
			if (qLim1 == 0) {
				qLim1 = -0.1;
			} if (qLim2 == 0) {
				qLim2 = 0.1;
			}
			_agents[AgentIndice].setAgent(AgentIndice, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 2);
			_Ub.set(AgentIndice, 0, pLim2);
			_Lb.set(AgentIndice, 0, 0);
			if (AC) {
				_Ub.set(AgentIndice + _nAgent, 0, qLim2);
				_Lb.set(AgentIndice + _nAgent, 0, qLim1);
			}
			
		}
		else { // prosumer
			pLim1 = Plim1[id];
			pLim2 = Plim2[id];
			cost1 = Cost1[id] / 1000;
			cost2 = Cost2[id];
			nVoisin = _nGen + _nCons;

			qLim1 = Qobj[id] * 2 * (Qobj[id] < 0);
			qLim2 = Qobj[id] * 2 * (Qobj[id] > 0);

			costQ1 = 0.1;
			costQ2 = -costQ1 * Qobj[id];
			if (qLim1 == 0) {
				qLim1 = -0.1;
			} if (qLim2 == 0) {
				qLim2 = 0.1;
			}
			_agents[AgentIndice].setAgent(AgentIndice, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 3);
			_Ub.set(AgentIndice, 0, pLim2);
			_Lb.set(AgentIndice, 0, pLim1);

			if(AC){
				_Ub.set(AgentIndice + _nAgent, 0, qLim2);
				_Lb.set(AgentIndice + _nAgent, 0, qLim1);
			}
		}
		_a.set(AgentIndice, 0, cost1);
		_b.set(AgentIndice, 0, cost2);

		_Pmin.set(AgentIndice, 0, pLim1);
		_Pmax.set(AgentIndice, 0, pLim2);
		_nVoisin.set(AgentIndice, 0, nVoisin);

		if (AC) {
			_a.set(AgentIndice + _nAgent, 0, costQ1);
			_b.set(AgentIndice + _nAgent, 0, costQ2);
			_Pmin.set(AgentIndice + _nAgent, 0, qLim1);
			_Pmax.set(AgentIndice + _nAgent, 0, qLim2);
			_nVoisin.set(AgentIndice + _nAgent, 0, _nAgent - 1);
		}
	}
	
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
	return coresBusAgentLin;
}
MatrixCPU StudyCaseAgent::Set3BusOld(bool AC) {
	
	if (AC) {
		throw std::invalid_argument("WIP not done yet");
	}
	clock_t t = clock();
	int offset = 0;
	if (AC) {
		offset = 1;
	}
	_nAgent = 3;
	_nPro = 0;
	_nGen = 2;
	_nCons = 1;
	DELETEA(_agents);
	_agents = new Agent[_nAgent];
	float Plim1[3] = { -1.51, 0, 0};
	float Plim2[3] = { -1.50, 1.00, 2.00};
	float Cost1[3] = { 1000, 0.01, 0.01 }; // on aimerait a = 0 mais on va �viter
	float Cost2[3] = { 1.50, 60, 120};
	
	int BusAgent[3] = { 0, 1, 2 };
	MatrixCPU coresBusAgentLin = MatrixCPU(_nAgent, 1);

	for (int i = offset; i < _nAgent; i++) {
		coresBusAgentLin.set(i, 0, BusAgent[i - offset]);
	}
	float pLim1;
	float pLim2;
	float cost1;
	float cost2;
	_connect = MatrixCPU(_nAgent, _nAgent);
	_a = MatrixCPU(_nAgent, 1);
	_b = MatrixCPU(_nAgent, 1);
	_Ub = MatrixCPU(_nAgent, 1);
	_Lb = MatrixCPU(_nAgent, 1);
	_Pmin = MatrixCPU(_nAgent, 1);
	_Pmax = MatrixCPU(_nAgent, 1);
	_nVoisin = MatrixCPU(_nAgent, 1);
	_BETA = MatrixCPU(_nAgent, _nAgent);
	genConnec(&_connect);
	genBetaUniforme(0);
	int nVoisin;

	for (int id = 0; id < _nAgent; id++)
	{
		if (id < _nCons) { // consumer
			pLim1 = Plim1[id];
			pLim2 = Plim2[id];
			cost1 = Cost1[id] / 1000;
			cost2 = Cost2[id];
			nVoisin = _nGen + _nPro;

			(_agents[id]).setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);

			_Ub.set(id, 0, 0);
			_Lb.set(id, 0, pLim1);

		}
		else if (id < (_nCons + _nGen)) { // generator
			pLim1 = Plim1[id];
			pLim2 = Plim2[id];
			cost1 = Cost1[id] / 1000;
			cost2 = Cost2[id];
			nVoisin = _nCons + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 2);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, 0);
		}
		else { // prosumer
			pLim1 = Plim1[id];
			pLim2 = Plim2[id];
			cost1 = Cost1[id] / 1000;
			cost2 = Cost2[id];
			nVoisin = _nGen + _nCons;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 3);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, pLim1);
		}
		_a.set(id, 0, cost1);
		_b.set(id, 0, cost2);

		_Pmin.set(id, 0, pLim1);
		_Pmax.set(id, 0, pLim2);
		_nVoisin.set(id, 0, nVoisin);
	}
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
	return coresBusAgentLin;
}

MatrixCPU StudyCaseAgent::Set3Bus(bool AC)
{
	clock_t t = clock();
	_Sbase = 1; // MW car matlab
	int offset = 0;
	if (AC) {
		offset = 1;
	}
	_nAgent = 3 + offset; // agent des pertes
	_nPro = 0;
	_nGen = 2;
	_nCons = 1 + offset; // agent des pertes
	DELETEA(_agents);
	_agents = new Agent[_nAgent];
	float Pobj[3] = { -2.00, 1.30,  0.70 };
	float Qobj[3] = { -1.20, 0.4137,  0.4701 };
	int BusAgent[3] = { 1, 0, 1 }; // avoir 1 bus sans agent et 1 bus avec 2 agents
	float pLim1, pLim2, cost1, cost2, costQ1, costQ2, qLim1, qLim2;
	MatrixCPU coresBusAgentLin = MatrixCPU(_nAgent, 1);


	for (int i = offset; i < _nAgent; i++) {
		Pobj[i - offset] = Pobj[i - offset] / _Sbase;
		Qobj[i - offset] = Qobj[i - offset] / _Sbase;
	}
	for (int i = offset; i < _nAgent; i++) {
		coresBusAgentLin.set(i, 0, BusAgent[i - offset]);
	}

	if (AC) {
		initMatAC();
	}
	else {
		initMat();
	}

	int nVoisin;

	if (AC) { // agent des pertes
		pLim1 = 0; // sera modifi� par l'algo selon le calcul des pertes
		pLim2 = 0; // idem
		cost1 = 0; // cela ne cause aucune division par 0 et la non convexit� ne joue pas comme on a qu'un seul point d�fini
		cost2 = 0; // idem
		nVoisin = _nGen + _nPro;

		qLim1 = 0; // idem que puissance active
		qLim2 = 0;

		costQ1 = 0;
		costQ2 = 0;

		(_agents[0]).setAgent(0, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);

		_Ub.set(0, 0, 0);
		_Lb.set(0, 0, -FLT_MAX); // pour ne pas avoir besoin de le modifier

		_Ub.set(0 + _nAgent, 0, FLT_MAX); // idem
		_Lb.set(0 + _nAgent, 0, -FLT_MAX);
		/*_a.set(0, 0, cost1);
		_b.set(0, 0, cost2);
		_a.set(id + _nAgent, 0, costQ1);
		_b.set(id + _nAgent, 0, costQ2);
		_Pmin.set(id + _nAgent, 0, qLim1);
		_Pmax.set(id + _nAgent, 0, qLim2);

		_Pmin.set(0, 0, pLim1);
		_Pmax.set(0, 0, pLim2);*/

		_nVoisin.set(0, 0, nVoisin);
		_nVoisin.set(_nAgent, 0, _nAgent - 1);
	}

	for (int id = offset; id < _nAgent; id++)
	{
		if (id < _nCons) { // consumer
			pLim1 = Pobj[id - offset] * (0.95 + 0.1 * (Pobj[id - offset] < 0));
			pLim2 = Pobj[id - offset] * (1.05 - 0.1 * (Pobj[id - offset] < 0));
			cost1 = 0.1 * _Sbase * _Sbase;
			cost2 = -cost1 * Pobj[id - offset];
			nVoisin = _nGen + _nPro;

			
			qLim1 = Qobj[id - offset] * (0.9 + 0.2 * (Qobj[id - offset] < 0));
			qLim2 = Qobj[id - offset] * (1.1 - 0.2 * (Qobj[id - offset] < 0));

			costQ1 = 0.1 * _Sbase * _Sbase;
			costQ2 = -costQ1 * Qobj[id - offset];

			(_agents[id]).setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);

			_Ub.set(id, 0, 0);
			_Lb.set(id, 0, pLim1);

			if (AC) {
				_Ub.set(id + _nAgent, 0, qLim2 * (qLim2 > 0));
				_Lb.set(id + _nAgent, 0, qLim1 * (qLim1 < 0));
			}
		}
		else if (id < (_nCons + _nGen)) { // generator

			pLim1 = 0;
			pLim2 = 2 * Pobj[id - offset];

			cost1 = 0.1 * _Sbase * _Sbase;
			cost2 = 0.5 * _Sbase;
			nVoisin = _nCons + _nPro;

			qLim1 = Qobj[id - offset] * 2 * (Qobj[id - offset] < 0);
			qLim2 = Qobj[id - offset] * 2 * (Qobj[id - offset] > 0);

			costQ1 = 0.1 * _Sbase * _Sbase;
			costQ2 = -costQ1 * Qobj[id - offset];

			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 2);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, 0);

			if (AC) {
				_Ub.set(id + _nAgent, 0, qLim2);
				_Lb.set(id + _nAgent, 0, qLim1);
			}

		}
		_a.set(id, 0, cost1);
		_b.set(id, 0, cost2);
		

		_Pmin.set(id, 0, pLim1);
		_Pmax.set(id, 0, pLim2);
		
		_nVoisin.set(id, 0, nVoisin);
		
		_Pobj.set(id, 0, Pobj[id - offset]);
		if (AC) {
			_a.set(id + _nAgent, 0, costQ1);
			_b.set(id + _nAgent, 0, costQ2);
			_Pmin.set(id + _nAgent, 0, qLim1);
			_Pmax.set(id + _nAgent, 0, qLim2);
			_nVoisin.set(id + _nAgent, 0, _nAgent - 1);
			_Pobj.set(id + _nAgent, 0, Qobj[id - offset]);
		}
		
	}
	_Pobj.toMatCPUD(_PobjD);

	return coresBusAgentLin;
}
MatrixCPU StudyCaseAgent::Set4Agent()
{
	clock_t t = clock();
	_nAgent = 4;
	_nPro = 1;
	_nGen = 1;
	_nCons = 2;
	DELETEA(_agents);
	_agents = new Agent[_nAgent];
	int BusAgent[4] = { 0, 0, 1, 1 };
	MatrixCPU coresBusAgentLin = MatrixCPU(_nAgent, 1);

	for (int i = 0; i < _nAgent; i++) {
		coresBusAgentLin.set(i, 0, BusAgent[i - 0]);
	}
	float Plim1[4] = {-30, -30, 0, -20 };
	float Plim2[4] = { -1, -2, 60, 30 };
	float Cost1[4] = { 1, 1, 0.7, 0.7 };
	float Cost2[4] = { 70, 70, 40, 60 };
	float pLim1;
	float pLim2;
	float cost1;
	float cost2;
	initMat();
	genBetaUniforme(1);
	int nVoisin;

	for (int id = 0; id < _nAgent; id++)
	{
		if (id < _nCons) { // consumer
			pLim1 = Plim1[id];
			pLim2 = Plim2[id];
			cost1 = Cost1[id];
			cost2 = Cost2[id];
			nVoisin = _nGen + _nPro;

			(_agents[id]).setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);

			_Ub.set(id, 0, 0);
			_Lb.set(id, 0, pLim1);

		}
		else if (id < (_nCons + _nGen)) { // generator
			pLim1 = Plim1[id];
			pLim2 = Plim2[id];
			cost1 = Cost1[id];
			cost2 = Cost2[id];
			nVoisin = _nCons + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 2);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, 0);
		}
		else { // prosumer
			pLim1 = Plim1[id];
			pLim2 = Plim2[id];
			cost1 = Cost1[id];
			cost2 = Cost2[id];
			nVoisin = _nGen + _nCons;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 3);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, pLim1);
		}
		_a.set(id, 0, cost1);
		_b.set(id, 0, cost2);

		_Pmin.set(id, 0, pLim1);
		_Pmax.set(id, 0, pLim2);
		_nVoisin.set(id, 0, nVoisin);
	}
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
	return coresBusAgentLin;
}


MatrixCPU StudyCaseAgent::Set2node(bool AC)
{
	clock_t t = clock();
	int offset = 0;
	if (AC) {
		offset = 1;
	}
	_AC = AC;
	_nAgent = 2 + offset;
	_nPro = 0;
	_nGen = 1;
	_nCons = 1 + offset;
	int BusAgent[2] = { 0, 1 };
	MatrixCPU coresBusAgentLin = MatrixCPU(_nAgent, 1);

	for (int i = offset; i < _nAgent; i++) {
		coresBusAgentLin.set(i, 0, BusAgent[i - offset]);
	}


	float Plim1[2] = { -30, 0 };
	float Plim2[2] = { 0, 60 };
	float Cost1[2] = { 1, 1 };
	float Cost2[2] = { 8, 4 };
	float Qobj[2] = { -1, 0 };
	float pLim1, pLim2, cost1, cost2, qLim1, qLim2, costQ1, costQ2;
	if (AC) {
		initMatAC();
	}
	else {
		initMat();
	}
	
	genBetaUniforme(0);
	int nVoisin;

	if (AC) { // agent des pertes
		pLim1 = 0; // sera modifi� par l'algo selon le calcul des pertes
		pLim2 = 0; // idem
		cost1 = 1; // Doit �tre convee
		cost2 = 0; // idem
		nVoisin = _nGen + _nPro;

		qLim1 = 0; // idem que puissance active
		qLim2 = 0;

		costQ1 = 0;
		costQ2 = 0;

		(_agents[0]).setAgent(0, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);
		_Ub.set(0, 0, 0);
		_Lb.set(0, 0, -FLT_MAX); // pour ne pas avoir besoin de le modifier

		_Ub.set(0 + _nAgent, 0, 0); // Pour �viter de la sp�culation
		_Lb.set(0 + _nAgent, 0, -0);
		_nVoisin.set(0, 0, nVoisin);
		_nVoisin.set(_nAgent, 0, _nAgent - 1);
	}

	int AgentIndice = 0;
	for (int id = 0; id < _nAgent - offset; id++)
	{
		AgentIndice = id + offset;
		if (id < _nCons-offset) { // consumer

			pLim1 = Plim1[id];
			pLim2 = Plim2[id];
			cost1 = Cost1[id];
			cost2 = Cost2[id];
			nVoisin = _nGen + _nPro;

			qLim1 = Qobj[id] * (0.9 + 0.2 * (Qobj[id] < 0));
			qLim2 = Qobj[id] * (1.1 - 0.2 * (Qobj[id] < 0));

			if (qLim1 == 0) {
				qLim1 = -1;
			} if (qLim2 == 0) {
				qLim2 = 1;
			}

			costQ1 = 0.1;
			costQ2 = -costQ1 * Qobj[id];

			(_agents[AgentIndice]).setAgent(AgentIndice, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);
			
			_Ub.set(AgentIndice, 0, 0);
			_Lb.set(AgentIndice, 0, pLim1);
			if (AC) {
				_Ub.set(AgentIndice + _nAgent, 0, qLim2 * (qLim2 > 0));
				_Lb.set(AgentIndice + _nAgent, 0, qLim1 * (qLim1 < 0));
			}


		}
		else if (id < (_nCons + _nGen - offset)) { // generator
			pLim1 = Plim1[id];
			pLim2 = Plim2[id];
			cost1 = Cost1[id];
			cost2 = Cost2[id];
			nVoisin = _nCons + _nPro;

			qLim1 = Qobj[id] * 2 * (Qobj[id] < 0);
			qLim2 = Qobj[id] * 2 * (Qobj[id] > 0);
			if (qLim1 == 0) {
				qLim1 = -1;
			} if (qLim2 == 0) {
				qLim2 = 1;
			}
			costQ1 = 0.1;
			costQ2 = -costQ1 * Qobj[id];
			
			_agents[AgentIndice].setAgent(AgentIndice, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 2);
			_Ub.set(AgentIndice, 0, pLim2);
			_Lb.set(AgentIndice, 0, 0);
			if (AC) {
				_Ub.set(AgentIndice + _nAgent, 0, qLim2);
				_Lb.set(AgentIndice + _nAgent, 0, qLim1);
			}

		}
		else { // prosumer
			pLim1 = Plim1[id];
			pLim2 = Plim2[id];
			cost1 = Cost1[id];
			cost2 = Cost2[id];
			nVoisin = _nGen + _nCons;

			qLim1 = Qobj[id] * 2 * (Qobj[id] < 0);
			qLim2 = Qobj[id] * 2 * (Qobj[id] > 0);
			if (qLim1 == 0) {
				qLim1 = -0.1;
			} if (qLim2 == 0) {
				qLim2 = 0.1;
			}
			costQ1 = 0.1;
			costQ2 = -costQ1 * Qobj[id];
			
			_agents[AgentIndice].setAgent(AgentIndice, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 3);
			_Ub.set(AgentIndice, 0, pLim2);
			_Lb.set(AgentIndice, 0, pLim1);

			if (AC) {
				_Ub.set(AgentIndice + _nAgent, 0, qLim2);
				_Lb.set(AgentIndice + _nAgent, 0, qLim1);
			}
		}
		
		_a.set(AgentIndice, 0, cost1);
		_b.set(AgentIndice, 0, cost2);

		_Pmin.set(AgentIndice, 0, pLim1);
		_Pmax.set(AgentIndice, 0, pLim2);
		_nVoisin.set(AgentIndice, 0, nVoisin);

		if (AC) {
			
			_a.set(AgentIndice + _nAgent, 0, costQ1);
			_b.set(AgentIndice + _nAgent, 0, costQ2);
			_Pmin.set(AgentIndice + _nAgent, 0, qLim1);
			_Pmax.set(AgentIndice + _nAgent, 0, qLim2);
			_nVoisin.set(AgentIndice + _nAgent, 0, _nAgent - 1);
			
		}
	}
	
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
	
	return coresBusAgentLin;
}


MatrixCPU StudyCaseAgent::Set4nodeBis(bool AC)
{
	// cas d'�tude pour simuler le cas d'EVA pendant son stage
	if (AC) {
		throw std::invalid_argument("WIP not done yet");
	}
	clock_t t = clock();
	_nAgent = 4;
	_nPro = 0;
	_nGen = 2;
	_nCons = 2;
	DELETEA(_agents);
	_agents = new Agent[_nAgent];
	float Plim1[4] = { -0.8001, -2.001, 0.4 , 0 };
	float Plim2[4] = { -0.8, -2, 10, 7 };
	float Cost1[4] = { 1, 1, 0.0001, 0.0001 };
	float Cost2[4] = { 0.8, 2, 7, 3 };
	int BusAgent[4] = { 1, 2, 0, 1 };
	MatrixCPU coresBusAgentLin = MatrixCPU(_nAgent, 1);

	for (int i = 0; i < _nAgent; i++) {
		coresBusAgentLin.set(i, 0, BusAgent[i - 0]);
	}
	float pLim1;
	float pLim2;
	float cost1;
	float cost2;
	_connect = MatrixCPU(_nAgent, _nAgent);
	_a = MatrixCPU(_nAgent, 1);
	_b = MatrixCPU(_nAgent, 1);
	_Ub = MatrixCPU(_nAgent, 1);
	_Lb = MatrixCPU(_nAgent, 1);
	_Pmin = MatrixCPU(_nAgent, 1);
	_Pmax = MatrixCPU(_nAgent, 1);
	_nVoisin = MatrixCPU(_nAgent, 1);
	_BETA = MatrixCPU(_nAgent, _nAgent);
	genConnec(&_connect);
	genBetaUniforme(1);
	int nVoisin;

	for (int id = 0; id < _nAgent; id++)
	{
		if (id < _nCons) { // consumer
			pLim1 = Plim1[id];
			pLim2 = Plim2[id];
			cost1 = Cost1[id];
			cost2 = Cost2[id];
			nVoisin = _nGen + _nPro;

			(_agents[id]).setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);

			_Ub.set(id, 0, 0);
			_Lb.set(id, 0, pLim1);

		}
		else if (id < (_nCons + _nGen)) { // generator
			pLim1 = Plim1[id];
			pLim2 = Plim2[id];
			cost1 = Cost1[id];
			cost2 = Cost2[id];
			nVoisin = _nCons + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 2);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, 0);
		}
		else { // prosumer
			pLim1 = Plim1[id];
			pLim2 = Plim2[id];
			cost1 = Cost1[id];
			cost2 = Cost2[id];
			nVoisin = _nGen + _nCons;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 3);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, pLim1);
		}
		_a.set(id, 0, cost1);
		_b.set(id, 0, cost2);

		_Pmin.set(id, 0, pLim1);
		_Pmax.set(id, 0, pLim2);
		_nVoisin.set(id, 0, nVoisin);
	}
	
	
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
	return coresBusAgentLin;

}


void StudyCaseAgent::SetEuropeP0(const std::string& path, MatrixCPU* P0)
{
	clock_t t = clock();
	_nGen = 969;
	_nCons = 1494;
	_nAgent = _nGen + _nCons;
	_nPro = 0;
	_name = "Europe";
	float dP = 0.1; // P = P0 +/- dP * P0 for the consumers

	initMat();
	
	int nVoisin;
	float pLim1, pLim2, cost1, cost2;

	MatrixCPU Pgen(_nGen, 1);
	MatrixCPU Cost(_nGen, 1);
	GenBus = MatrixCPU(_nGen, 1);
	std::string pathGen = path + _name + "/genCarac.txt";
	setGenFromFile(pathGen, &Pgen, &Cost, &GenBus);
	


	for (int id = 0; id < _nAgent; id++)
	{
		if (id < _nCons) { // consumer

			pLim1 = -(1 + dP) * P0->get(id, 0);
			pLim2 = -(1 - dP) * P0->get(id, 0);
			cost1 = 1;
			cost2 = P0->get(id, 0) * cost1;

			nVoisin = _nGen + _nPro;
			(_agents[id]).setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);
			_Ub.set(id, 0, 0);
			_Lb.set(id, 0, pLim1);

		}
		else { // generator

			pLim1 = 0;
			pLim2 = Pgen.get(id - _nCons, 0);
			cost1 = 0.1;
			cost2 = Cost.get(id - _nCons, 0);

			nVoisin = _nCons + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 2);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, 0);
		}
		_Pmin.set(id, 0, pLim1);
		_Pmax.set(id, 0, pLim2);
		_a.set(id, 0, cost1);
		_b.set(id, 0, cost2);
		_nVoisin.set(id, 0, nVoisin);

	}

	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}


void StudyCaseAgent::SetStudyCaseAgent(std::string path, std::string name, MatrixCPU* P0)
{
	// recuperation des tailles des fichiers
	std::string pathGen = path + "genCarac" + name + ".txt";

	_name = name;
	
	_nGen = getNFileline(pathGen);
	_nCons = P0->getNLin();
	
	clock_t t = clock();

	_nAgent = _nGen + _nCons;
	_nPro = 0;

	float dP = 0.1; // P = P0 +/- dP * P0 for the consumers
	
	initMat();

	int nVoisin;
	float pLim1, pLim2, cost1, cost2;
	genBetaUniforme(0);


	MatrixCPU Pgen(_nGen, 1);
	MatrixCPU Cost(_nGen, 1);
	GenBus = MatrixCPU(_nGen, 1);
	
	setGenFromFile(pathGen, &Pgen, &Cost, &GenBus);

	for (int id = 0; id < _nAgent; id++)
	{
		if (id < _nCons) { // consumer

			pLim1 = -(1 + dP) * P0->get(id, 0);
			pLim2 = -(1 - dP) * P0->get(id, 0);
			cost1 = 1;
			cost2 = P0->get(id, 0) * cost1;

			nVoisin = _nGen + _nPro;
			(_agents[id]).setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);
			_Ub.set(id, 0, 0);
			_Lb.set(id, 0, pLim1);

		}
		else { // generator

			pLim1 = 0;
			pLim2 = Pgen.get(id - _nCons, 0);
			cost1 = 0.1;
			cost2 = Cost.get(id - _nCons, 0);

			nVoisin = _nCons + _nPro;
			_agents[id].setAgent(id, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 2);
			_Ub.set(id, 0, pLim2);
			_Lb.set(id, 0, 0);
		}
		_Pmin.set(id, 0, pLim1);
		_Pmax.set(id, 0, pLim2);
		_a.set(id, 0, cost1);
		_b.set(id, 0, cost2);
		_nVoisin.set(id, 0, nVoisin);

	}
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}

MatrixCPU StudyCaseAgent::SetACFromFile(std::string name, std::string path)
{
	std::string fileName1 = path + "Case" + name + ".txt";
	std::string fileName2 = path + "Agent" + name + ".txt";

	MatrixCPUD Info(1, 9); // Sbase, Vbase, nAgent, nCons, nGenSup, nBus, nLine, V0, theta0
	Info.setFromFile(fileName1);
	_Sbase = Info.get(0, 0);
	_AC = true;

	_nAgent = Info.get(0, 2) + 1; // + the loss agent
	_nCons = Info.get(0, 3) + 1; // + the loss agent
	_nGen = _nAgent - _nCons;
	_nPro = 0;

	initMatAC();
	MatrixCPUD Mat(_nAgent - 1, 10); // bus, a, b, P, Pmin, Pmax, Qobj, Qmin, Qmax, zone
	Mat.setFromFile(fileName2);
	double pLim1, pLim2, cost1, cost2, qLim1, qLim2, costQ1, costQ2, Qobj;
	int nVoisin;
	int offsetbus = 1;

	(_agents[0]).setAgent(0, 0, 0, 0, 0, _nGen, &_connect, _nAgent, 1);
	_Lb.set(0, 0, -10000); // pour ne pas avoir besoin de le modifier
	_Ub.set(_nAgent, 0, 10000); // idem
	_Lb.set(_nAgent, 0, -10000);
	//_CoresBusAgent.set(0, 0, 1);
	_nVoisin.set(0, 0, _nGen);
	_nVoisin.set(_nAgent, 0, _nAgent - 1);

	MatrixCPU coresBusAgentLin = MatrixCPU(_nAgent, 1);

#ifdef DEBUG
	std::cout << "set Agent" << std::endl;
#endif // DEBUG

	for (int i = 1; i < _nAgent; i++) {
		int bus = Mat.get(i - 1, 0);
		coresBusAgentLin.set(i, 0, bus - offsetbus);
		cost1 = Mat.get(i - 1, 1) * (_Sbase * _Sbase);
		cost2 = Mat.get(i - 1, 2) * _Sbase;
		_Pobj.set(i, 0, Mat.get(i - 1, 3) / _Sbase);
		_PobjD.set(i, 0, Mat.get(i - 1, 3) / _Sbase);
		pLim1 = Mat.get(i - 1, 4) / _Sbase;
		pLim2 = Mat.get(i - 1, 5) / _Sbase;
		Qobj = Mat.get(i - 1, 6) / _Sbase;
		_Pobj.set(i + _nAgent, 0, Qobj);
		_PobjD.set(i + _nAgent, 0, Qobj);
		qLim1 = Mat.get(i - 1, 7) / _Sbase;
		qLim2 = Mat.get(i - 1, 8) / _Sbase;
		
		costQ1 = 0.1 * (_Sbase * _Sbase);
		costQ2 = -costQ1 * Qobj;
		
		if (i < _nCons) {
			nVoisin = _nGen;
			_Ub.set(i, 0, 0);
			_Lb.set(i, 0, pLim1);
		}
		else {
			nVoisin = _nCons;
			_Ub.set(i, 0, pLim2);
			_Lb.set(i, 0, 0);
		}
		_Pmin.set(i, 0, pLim1);
		_Pmax.set(i, 0, pLim2);

		(_agents[i]).setAgent(i, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);

		_a.set(i, 0, cost1);
		_b.set(i, 0, cost2);
		_a.set(i + _nAgent, 0, costQ1);
		_b.set(i + _nAgent, 0, costQ2);

		
		_Pmin.set(i + _nAgent, 0, qLim1);
		_Pmax.set(i + _nAgent, 0, qLim2);
		_nVoisin.set(i, 0, nVoisin);
		_nVoisin.set(i + _nAgent, 0, _nAgent - 1);
		_Ub.set(i + _nAgent, 0, qLim2 * (qLim2 > 0));
		_Lb.set(i + _nAgent, 0, qLim1 * (qLim1 < 0));
	}

	return coresBusAgentLin;
}


MatrixCPU StudyCaseAgent::SetACFromFileSimplify(std::string name, std::string path)
{
	std::string fileName1 = path + "Case" + name + ".txt";
	std::string fileName2 = path + "Agent" + name + ".txt";

	MatrixCPUD Info(1, 9); // Sbase, Vbase, nAgent, nCons, nGenSup, nBus, nLine, V0, theta0
	Info.setFromFile(fileName1);
	_Sbase = Info.get(0, 0);
	_AC = true;

	_nAgent = Info.get(0, 2) + 1; // + the loss agent
	_nCons = Info.get(0, 3) + 1; // + the loss agent
	_nGen = _nAgent - _nCons;
	_nPro = 0;

	initMatAC();
	MatrixCPUD Mat(_nAgent - 1, 10); // bus, a, b, P, Pmin, Pmax, Qobj, Qmin, Qmax, zone
	Mat.setFromFile(fileName2);
	double pLim1, pLim2, cost1, cost2, qLim1, qLim2, costQ1, costQ2, Qobj;
	int nVoisin;
	int offsetbus = 1;

	(_agents[0]).setAgent(0, 0, 0, 0, 0, _nGen, &_connect, _nAgent, 1);
	_Lb.set(0, 0, -10000); // pour ne pas avoir besoin de le modifier
	_Ub.set(_nAgent, 0, 10000); // idem
	_Lb.set(_nAgent, 0, -10000);
	//_CoresBusAgent.set(0, 0, 1);
	_nVoisin.set(0, 0, _nGen);
	_nVoisin.set(_nAgent, 0, _nAgent - 1);

	MatrixCPU coresBusAgentLin = MatrixCPU(_nAgent, 1);

#ifdef DEBUG
	std::cout << "set Agent" << std::endl;
#endif // DEBUG

	for (int i = 1; i < _nAgent; i++) {
		int bus = Mat.get(i - 1, 0);
		coresBusAgentLin.set(i, 0, bus - offsetbus);
		cost1 = Mat.get(i - 1, 1) * (_Sbase * _Sbase);
		cost2 = Mat.get(i - 1, 2) * _Sbase;
		_Pobj.set(i, 0, Mat.get(i - 1, 3) / _Sbase);
		_PobjD.set(i, 0, Mat.get(i - 1, 3) / _Sbase);
		pLim1 = Mat.get(i - 1, 4) / _Sbase;
		pLim2 = Mat.get(i - 1, 5) / _Sbase;
		Qobj = Mat.get(i - 1, 6) / _Sbase;
		_Pobj.set(i + _nAgent, 0, Qobj);
		_PobjD.set(i + _nAgent, 0, Qobj);
		qLim1 = Mat.get(i - 1, 7) / _Sbase;
		qLim2 = Mat.get(i - 1, 8) / _Sbase;

		costQ1 = 0.1 * (_Sbase * _Sbase);
		costQ2 = -costQ1 * Qobj;

		if (i < _nCons) {
			nVoisin = _nGen;
			_Ub.set(i, 0, 0);
			_Lb.set(i, 0, pLim1);
			pLim1 = pLim2;
			qLim1 = qLim2;/**/
		}
		else {
			nVoisin = _nCons;
			_Ub.set(i, 0, pLim2);
			_Lb.set(i, 0, 0);
			cost2 = 10;
			pLim1 = 0;
		}
		_Pmin.set(i, 0, pLim1);
		_Pmax.set(i, 0, pLim2);

		(_agents[i]).setAgent(i, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);

		_a.set(i, 0, cost1);
		_b.set(i, 0, cost2);
		_a.set(i + _nAgent, 0, costQ1);
		_b.set(i + _nAgent, 0, costQ2);


		_Pmin.set(i + _nAgent, 0, qLim1);
		_Pmax.set(i + _nAgent, 0, qLim2);
		_nVoisin.set(i, 0, nVoisin);
		_nVoisin.set(i + _nAgent, 0, _nAgent - 1);
		_Ub.set(i + _nAgent, 0, qLim2 * (qLim2 > 0));
		_Lb.set(i + _nAgent, 0, qLim1 * (qLim1 < 0));
	}

	return coresBusAgentLin;
}

MatrixCPU StudyCaseAgent::SetFromInterface(StudyCaseInterface* interface, bool DC)
{
    MatrixCPU Info = interface->getInfoCase();
	// Sbase, Vbase, nAgent, nCons, nGenSup, nBus, nLine, V0, theta0

	_Sbase = Info.get(0, Sbase_ind);
	_AC = !DC;
	int offset = 1 * (_AC);

	_nAgent = Info.get(0, nAgent_ind) + offset; // + the loss agent
	_nCons  = Info.get(0, nCons_ind)  + offset; // + the loss agent
	if(_nCons < offset+1){
		interface->checkCase();
		Info = interface->getInfoCase();
		_nCons = Info.get(0, nCons_ind) + offset;
	}

	_nGen = Info.get(0, nGen_ind);
	_nPro = _nAgent -_nCons - _nGen;

	if(_AC){
		initMatAC();
	} else{
		initMat();
	}

	
	if(interface->isConnexionDefined()){
		MatrixCPU connect_temp = interface->getConnexion(); // taille ((N-offset)) * (N-offset)
		_connect = MatrixCPU(_nAgent, _nAgent);
		if(_AC){
			for(int i=_nCons; i<_nAgent; i++){ // P
				_connect.set(0, i, 1);
				_connect.set(i, 0, 1);
			}
			/*for(int i=1; i<_nAgent; i++){ // Q
				_connect.set(_nAgent, i, 1);
				_connect.set(_nAgent + i, _nAgent, 1);
			}*/
		}
		
		for(int i=offset; i<_nAgent;i++){
			for(int j=offset; j<_nAgent; j++){
				_connect.set(i, j, connect_temp.get(i - offset, j - offset));
				/*if(_AC){
					_connect.set(_nAgent + i, j, connect_temp.get(_nAgent + i - 2, j - 1));
				}*/
			}
		}
	}

	if(interface->isBetaDefined()){
		MatrixCPU Beta_temp = interface->getBeta(); // taille N*N (without loss agent)
		_BETA = MatrixCPU(_nAgent, _nAgent);
		
		for(int i=offset; i<_nAgent;i++){
			for(int j=offset; j<_nAgent; j++){
				_BETA.set(i, j, Beta_temp.get(i - offset, j - offset));
			}
		}
	}


	MatrixCPU Mat = interface->getAgentCase();
	// bus, a, b, P, Pmin, Pmax, Qobj, Qmin, Qmax, zone
	
	double pLim1, pLim2, cost1, cost2, qLim1, qLim2, costQ1, costQ2, Qobj;
	int nVoisin;
	
	if(_AC){
		(_agents[0]).setAgent(0, 0, 0, 0, 0, _nAgent -_nCons, &_connect, _nAgent, 1);
		_Lb.set(0, 0, -10000); // pour ne pas avoir besoin de le modifier
		_Ub.set(_nAgent, 0, 10000); // idem
		_Lb.set(_nAgent, 0, -10000);
		//_CoresBusAgent.set(0, 0, 1);
		_nVoisin.set(0, 0,  _nAgent -_nCons);
		_nVoisin.set(_nAgent, 0, _nAgent - 1);
	}
	

	MatrixCPU coresBusAgentLin = MatrixCPU(_nAgent, 1);

#ifdef DEBUG
	std::cout << "set Agent" << std::endl;
#endif // DEBUG

	for (int i = offset; i < _nAgent; i++) {
		int bus = Mat.get(i - offset, PosBus_ind);
		coresBusAgentLin.set(i, 0, bus);
		
		cost1  = Mat.get(i - offset, a_ind) * (_Sbase * _Sbase);
		cost2  = Mat.get(i - offset, b_ind) * _Sbase;
		costQ1 = Mat.get(i - offset, aq_ind) * (_Sbase * _Sbase);
		costQ2 = Mat.get(i - offset, bq_ind) * _Sbase;

		_Pobj.set(i,  0, Mat.get(i - offset, Pobj_ind) / _Sbase);
		_PobjD.set(i, 0, Mat.get(i - offset, Pobj_ind) / _Sbase);

		pLim1 = Mat.get(i - offset, Pmin_ind) / _Sbase;
		pLim2 = Mat.get(i - offset, Pmax_ind) / _Sbase;
		Qobj  = Mat.get(i - offset, Qobj_ind) / _Sbase;

		qLim1 = Mat.get(i - offset, Qmin_ind) / _Sbase;
		qLim2 = Mat.get(i - offset, Qmax_ind) / _Sbase;
		
		
		if (i < _nCons) {
			nVoisin = _nGen;
			_Ub.set(i, 0, 0);
			_Lb.set(i, 0, pLim1);
		}
		else {
			nVoisin = _nCons;
			_Ub.set(i, 0, pLim2);
			_Lb.set(i, 0, 0);
		}
		
		_Pmin.set(i, 0, pLim1);
		_Pmax.set(i, 0, pLim2);

		(_agents[i]).setAgent(i, pLim1, pLim2, cost1, cost2, nVoisin, &_connect, _nAgent, 1);

		_a.set(i, 0, cost1);
		_b.set(i, 0, cost2);
		_nVoisin.set(i, 0, nVoisin);

		if(_AC){
			_a.set(i + _nAgent, 0, costQ1);
			_b.set(i + _nAgent, 0, costQ2);
			_Pmin.set(i + _nAgent, 0, qLim1);
			_Pmax.set(i + _nAgent, 0, qLim2);
			_Pobj.set(i + _nAgent, 0, Qobj);
			_PobjD.set(i + _nAgent, 0, Qobj);
			_nVoisin.set(i + _nAgent, 0, _nAgent - 1);
			_Ub.set(i + _nAgent, 0, qLim2 * (qLim2 > 0));
			_Lb.set(i + _nAgent, 0, qLim1 * (qLim1 < 0));
		}
	}


	if(interface->isTradeBoundDefined()){
		_Lbmat = interface->getLbMat();
		_Ubmat = interface->getUbMat();
		isBoundTradeDifferent = true;
		//throw std::invalid_argument("WIP : special tradeBound not yet implemented in solver");
	}
	//display(1);

	return coresBusAgentLin;
}
MatrixCPU StudyCaseAgent::SetEuropeTestFeeder(std::string path, int beggining)
{
	std::string fileName1 = path + "CaseTestFeeder.txt";
	std::string fileName2 = path + "AgentTestFeeder.txt";
	std::string fileName3 = path + "AgentConsumptionTestFeeder.txt";

	MatrixCPUD Info(1, 8); // Sbase, Vbase, Zbase, nAgent, nBus, nLine, V0, theta0
	Info.setFromFile(fileName1);
	_Sbase = Info.get(0, 0);

	_nAgent = Info.get(0, 3) + 2; // + the loss agent + the grid
	_nCons = Info.get(0, 3) + 1; // + the loss agent
	_nGen = 0;
	_nPro = 1;
	int nMaxShape = 100;
	int nMinute = 24 * 60;

	_PobjTemp = MatrixCPU(nMinute, nMaxShape);
	_PobjTemp.setFromFile(fileName3);


	MatrixCPU MatAgent(_nAgent - 2, 3); // bus, Pmult, facteur puissance
	MatAgent.setFromFile(fileName2);

	_AC = true;
	_factor = MatrixCPU(_nAgent - 2, 1);
	_PF = MatrixCPU(_nAgent - 2, 1);
	MatrixCPU CoresAgentBuslin(_nAgent, 1);

	for (int i = 0; i < _nAgent - 2; i++) {
		int bus = MatAgent.get(i, 0);
		_factor.set(i, 0, MatAgent.get(i, 1));
		_PF.set(i, 0, MatAgent.get(i, 2));
		CoresAgentBuslin.set(i + 1, 0, bus);
	}
	
	initMatAC();
	temporalStep = beggining;
	initCaseFromPobj();

	
	//std::cout << " agent des pertes" << std::endl;
	_Lb.set(0, 0, -10000); // pour ne pas avoir besoin de le modifier
	_Ub.set(_nAgent, 0, 10000); // idem
	_Lb.set(_nAgent, 0, -10000);
	_nVoisin.set(0, 0, _nGen + _nPro);
	_nVoisin.set(_nAgent, 0, _nAgent - 1);
	(_agents[0]).setAgent(0, 0, 0, 0, 0, _nGen + _nPro, &_connect, _nAgent, 1);

	// rajouter le prod
	// le pri spot pour la cost fonction ?
	std::cout << "reseau " << std::endl;
	_a.set(_nCons, 0, 0.01);
	_b.set(_nCons, 0, 0.03); // 30� / MWh	
	_Ub.set(_nCons, 0, 10000); // pour ne pas avoir besoin de le modifier
	_Lb.set(_nCons, 0, -10000); // pour ne pas avoir besoin de le modifier
	_Pmin.set(_nCons, 0, -10000);
	_Pmax.set(_nCons, 0, 10000);
	_Ub.set(_nCons + _nAgent, 0, 10000); // idem
	_Lb.set(_nCons + _nAgent, 0, -10000);
	_Pmax.set(_nCons + _nAgent, 0, 10000); // idem
	_Pmin.set(_nCons + _nAgent, 0, -10000);
	
	_nVoisin.set(_nCons, 0, _nCons);
	_nVoisin.set(_nAgent + _nCons, 0, _nAgent - 1);

	(_agents[_nCons]).setAgent(_nCons, -10000, 10000, 0.01, 0.03, _nCons, &_connect, _nAgent, 3);

	std::cout << "fin agent " << std::endl;
	return CoresAgentBuslin;
}



void StudyCaseAgent::UpdateP0(MatrixCPU* P0) 
{

	float dP = 0.1; // P = P0 +/- dP * P0 for the consumers
	
	if (P0->getNLin() != _nCons) {
		throw std::invalid_argument("P0 hasn't the good number of column");
	}
	for (int id = 0; id < _nCons; id++)
	{
		_agents[id].updateP0(P0->get(id, 0), dP);
		_Pmin.set(id, 0, _agents[id].getPmin());
		_Pmax.set(id, 0, _agents[id].getPmax());
		_a.set(id, 0, _agents[id].getA());
		_b.set(id, 0, _agents[id].getB());
		_Ub.set(id, 0, 0);
		_Lb.set(id, 0, _agents[id].getLb());
	}
	
}




///////////////////////////////////////////////////////////////////////////////
// Getter
///////////////////////////////////////////////////////////////////////////////

MatrixCPU StudyCaseAgent::getBeta() const
{
	return _BETA;
}

MatrixCPU StudyCaseAgent::getC() const
{
	return _connect;
}

MatrixCPU StudyCaseAgent::geta() const
{
	return _a;
}

MatrixCPU StudyCaseAgent::getb() const
{
	return _b;
}

MatrixCPU StudyCaseAgent::getUb() const
{
	if(isBoundTradeDifferent){
		return _Ubmat;
	} else {
		return _Ub;
	}
	
}

MatrixCPU StudyCaseAgent::getLb() const
{
	if(isBoundTradeDifferent){
		return _Lbmat;
	} else {
		return _Lb;
	}
}

MatrixCPU StudyCaseAgent::getPmin() const
{
	return _Pmin;
}

MatrixCPU StudyCaseAgent::getPmax() const
{
	return _Pmax;
}

MatrixCPU StudyCaseAgent::getNvoi() const
{
	return _nVoisin;
}

MatrixCPU StudyCaseAgent::getPobj()
{
	
	if (_Pobj.getNCol() == 0 || _Pobj.max2() == 0) {
		_Pobj = _Pmin;
		_Pobj.add(&_Pmax);
		_Pobj.multiply(0.5);
	}
	
	return _Pobj;
}

MatrixCPUD StudyCaseAgent::getPobjD()
{
	if (_PobjD.getNCol() == 0) {
		if (_Pobj.getNCol() == 0) {
			_Pobj = _Pmin;
			_Pobj.add(&_Pmax);
			_Pobj.multiply(0.5);
			_Pobj.toMatCPUD(_PobjD);
		}
		else {
			_Pobj.toMatCPUD(_PobjD);
		}
	}
	return _PobjD;
}

MatrixCPU StudyCaseAgent::getGenBus() const
{
	return GenBus;
}


float StudyCaseAgent::getTimeInit() const
{
	return _timeInit;
}

int StudyCaseAgent::getNagent() const
{
	return _nAgent;
}

int StudyCaseAgent::getNCons() const
{
	return _nCons;
}


MatrixCPU StudyCaseAgent::getVoisin(int agent) const
{
	if (agent > _nAgent) {
		throw std::invalid_argument("this agent doesn't exists");
	}
	return _agents[agent].getVoisin();
}
Agent StudyCaseAgent::getAgent(int agent) const
{
	return _agents[agent];
}

std::string StudyCaseAgent::getName() const
{
	return _name;
}


///////////////////////////////////////////////////////////////////////////////
// Study Case Modifier
///////////////////////////////////////////////////////////////////////////////



void StudyCaseAgent::removeLink(int i, int j)
{
	
	if (i>=_nAgent && j>=_nAgent) {
		throw std::invalid_argument("indice out of range");
	}
	if (_connect.get(i, j) == 0) {
		
		throw std::invalid_argument("agent already not linked"); 
	}
	_connect.set(i, j, 0);
	_connect.set(j, i, 0);


	_nVoisin.set(i, 0, _nVoisin.get(i, 0) - 1);
	_nVoisin.set(j, 0, _nVoisin.get(j, 0) - 1);

	int type = 3;

	if (i < _nCons) {
		type = 1;
	} 
	else if (i >= (_nCons + _nGen)) {
		type = 2;
	}
	_agents[i].setAgent(i, _Pmin.get(i, 0), _Pmax.get(i, 0), _a.get(i, 0), _b.get(i, 0), _nVoisin.get(i, 0), &_connect, _nAgent, type);
	type = 3;
	
	if (j < _nCons) {
		type = 1;
	}
	else if (j >= (_nCons + _nGen)) {
		type = 2;
	}
	_agents[j].setAgent(j, _Pmin.get(j, 0), _Pmax.get(j, 0), _a.get(j, 0), _b.get(j, 0), _nVoisin.get(j, 0), &_connect, _nAgent, type);

}


void StudyCaseAgent::addLink(int i, int j)
{
	//std::cout << "link " << i << " " << j << std::endl;
	if (i >= _nAgent && j >= _nAgent) {
		throw std::invalid_argument("indice out of range");
	}
	if (i == j) {
		throw std::invalid_argument("agent must not be link to itself");
	}
	if (_connect.get(i, j) == 1) {
		throw std::invalid_argument("agent already linked");
	}
	int type[2] = { 3 , 3 };
	if (i < _nCons) {
		type[0] = 1;
	}
	else if (i >= (_nCons + _nGen)) {
		type[0] = 2;
	}
	if (j < _nCons) {
		type[1] = 1;
	}
	else if (j >= (_nCons + _nGen)) {
		type[1] = 2;
	}
	if (type[0] == type[1]) {
		throw std::invalid_argument("agent must not be the same type");
	}

	_connect.set(i, j, 1);
	_connect.set(j, i, 1);


	_nVoisin.set(i, 0, _nVoisin.get(i, 0) + 1);
	_nVoisin.set(j, 0, _nVoisin.get(j, 0) + 1);

	_agents[i].setAgent(i, _Pmin.get(i, 0), _Pmax.get(i, 0), _a.get(i, 0), _b.get(i, 0), _nVoisin.get(i, 0), &_connect, _nAgent, type[0]);
	_agents[j].setAgent(j, _Pmin.get(j, 0), _Pmax.get(j, 0), _a.get(j, 0), _b.get(j, 0), _nVoisin.get(j, 0), &_connect, _nAgent, type[1]);

}


Agent StudyCaseAgent::removeAgent(int agent)
{
	if (agent > _nAgent || agent<0) {
		throw std::invalid_argument("this agent doesn't exist");
	}
	Agent agentRemove(_agents[agent]);
	MatrixCPU omega(getVoisin(agent));
	int Nvoisinmax = _nVoisin.get(agent, 0);
	for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
		int idVoisin = omega.get(voisin, 0);
		removeLink(agent, idVoisin);
	}

	_Pmin.set(agent, 0, 0);
	_Pmax.set(agent, 0, 0);
	int type =  3;
	if (agent < _nCons) {
		type = 1;
	}
	else if (agent >= (_nCons + _nGen)) {
		type = 2;
	}
	_agents[agent].setAgent(agent, _Pmin.get(agent, 0), _Pmax.get(agent, 0), _a.get(agent, 0), _b.get(agent, 0), _nVoisin.get(agent, 0), &_connect, _nAgent, type);

	return agentRemove;

}

void StudyCaseAgent::restoreAgent(Agent& agent, bool all) {


	int id = agent.getId();
	float Pmin = agent.getPmin();
	float Pmax = agent.getPmax();
	int type = agent.getType();

	_Pmin.set(id, 0, Pmin);
	_Pmax.set(id, 0, Pmax);
	if (!all) {
		std::cout << "Beware restore agent without giving it neighbor can lead to a non solvable problem" << std::endl;
	}
	else {
		MatrixCPU omega(agent.getVoisin());
		int nVoisin = agent.getNVoisin();
		for (int i = 0;i < nVoisin;i++) {
			int idVoisin = omega.get(i, 0);
			addLink(id, idVoisin);
		}
	}
	_agents[id].setAgent(id, _Pmin.get(id, 0), _Pmax.get(id, 0), _a.get(id, 0), _b.get(id, 0), _nVoisin.get(id, 0), &_connect, _nAgent, type);
}

void StudyCaseAgent::saveCSV(const std::string& fileName, bool all)
{

	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	MatrixCPU nombre(1, 7);
	nombre.set(0, 0, _nAgent);
	nombre.set(0, 1, _nCons);
	nombre.set(0, 2, _nGen);
	nombre.set(0, 3, _nPro);
	nombre.saveCSV(fileName, mode);


	MatrixCPU temp(1, _nAgent);
	MatrixCPU zero(1, _nAgent);
	temp.addTrans(&_a);
	temp.saveCSV(fileName, mode);
	temp.set(&zero);
	temp.addTrans(&_b);
	temp.saveCSV(fileName, mode);
	temp.set(&zero);
	temp.addTrans(&_Lb);
	temp.saveCSV(fileName, mode);
	temp.set(&zero);
	temp.addTrans(&_Ub);
	temp.saveCSV(fileName, mode);
	temp.set(&zero);
	temp.addTrans(&_Pmin);
	temp.saveCSV(fileName, mode);
	temp.set(&zero);
	temp.addTrans(&_Pmax);
	temp.saveCSV(fileName, mode);
	temp.set(&zero);
	temp.addTrans(&_nVoisin);
	temp.saveCSV(fileName, mode);
	temp.set(&zero);

	if (all) {
		_BETA.saveCSV(fileName, mode);
		_connect.saveCSV(fileName, mode);
	}


}

void StudyCaseAgent::display(int type) 
{
	
	std::cout << "Study Case : " << _nAgent << " agents and modelisation is " << (_AC ? "AC":"DC") << std::endl;
	std::cout << " a :" << std::endl;
	_a.display();
	std::cout << " b :" << std::endl;
	_b.display();
	std::cout << " lower bound: " << std::endl;
	_Lb.display();
	std::cout << " uper bound: " << std::endl;
	_Ub.display();
	std::cout << " Pmin :" << std::endl;
	_Pmin.display();
	std::cout << " Pmax :" << std::endl;
	_Pmax.display();
	std::cout << " Peers count :" << std::endl;
	_nVoisin.display();
	std::cout << "Pobj" << std::endl;
	(getPobj()).display();
	if (type == 1) {
		std::cout << " connection :" << std::endl;
		_connect.display();
		std::cout << " beta :" << std::endl;
		_BETA.display();
	
	}
}

StudyCaseAgent::~StudyCaseAgent()
{
#ifdef DEBUG_DESTRUCTOR
	std::cout << "case destructor" << std::endl;
#endif

	DELETEA(_agents);

}


