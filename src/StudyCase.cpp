#include "../head/StudyCase.h"





float StudyCase::rand1() const
{
	float a = (float)(rand()) / ((float)(RAND_MAX));
	return a;
}

int StudyCase::randab(int a, int b) const
{
	return a + (rand() % (b - a));
}

void StudyCase::genCoresBusAgent(bool all)
{

	if (all) {
		_CoresBusAgent = MatrixCPU(_nBus, _nAgent);
		for (int n = 0; n < _nAgent; n++) {
			int bus = _CoresBusAgentLin.get(n, 0);
			_CoresBusAgent.set(bus, n, 1);
		}
	}

	_nAgentByBus = MatrixCPU(_nBus, 1);
	_CoresAgentBusLinBegin = MatrixCPU(_nBus, 1);
	_CoresAgentBusLin = MatrixCPU(_nAgent, 1);

	for (int i = 0; i < _nAgent; i++) {
		int bus = _CoresBusAgentLin.get(i, 0);
		_nAgentByBus.increment(bus, 0, 1);
	}


	int debut = 0;
	int* decompteAgent = new int[_nBus];
	for (int b = 0; b < _nBus; b++) {
		_CoresAgentBusLinBegin.set(b, 0, debut);
		debut += _nAgentByBus.get(b, 0);
		decompteAgent[b] = 0;
	}

	for (int n = 0; n < _nAgent; n++) {
		int bus = _CoresBusAgentLin.get(n, 0);
		int indice = _CoresAgentBusLinBegin.get(bus, 0) + decompteAgent[bus];
		decompteAgent[bus]++;
		_CoresAgentBusLin.set(indice, 0, n);
	}
	DELETEA(decompteAgent);
}

void StudyCase::initMat()
{
	_nBus = SCGrid->getNBus();
	_nAgent = SCAg.getNagent();
	_nLine = SCGrid->getNLine();
	_SensiPower = MatrixCPU(_nLine, _nAgent); // G = A*I
	_CoresBusAgent = MatrixCPU(_nBus, _nAgent); // I
	_CoresBusAgentLin = MatrixCPU(_nAgent, 1);
	for (int id = 0; id < _nAgent / 2; id++) { // par defaut
		_CoresBusAgentLin.set(id, 0, 0);
	}
	for (int id = _nAgent / 2; id < _nAgent; id++) {
		_CoresBusAgentLin.set(id, 0, 1);
	}

}

void StudyCase::createGrid(bool _DC)
{
	
	DC = _DC;
	DELETEB(SCGrid);
	if (_DC) {
		SCGrid = new StudyCaseDCGrid;
	}
	else {
		SCGrid = new StudyCaseACGrid;
	}
	SCGrid->toReduce = toReduce;
}


void StudyCase::setReduce(bool toReduce1)
{
	
	SCGrid->toReduce = toReduce1;	
	toReduce = toReduce1;
}



StudyCase::StudyCase()
{
	 _timeInit = 0;
	 createGrid(true);
	 initMat();
	 srand(time(nullptr));
}
StudyCase::StudyCase(int nAgent, float P, float dP, float a, float da, float b, float db, float propCons, float propPro)
{
	clock_t t = clock();
	
	SCAg = StudyCaseAgent(nAgent, P, dP, a, da, b, db, propCons, propPro);
	createGrid(true);
	int nLineConstraint = SCGrid->getNLineConstraint();
	//std::cout << "nLineConstraint = " << nLineConstraint << std::endl;
	initMat();
	genCoresBusAgent(true);
	
	_SensiPower.multiply(&SCGrid->getPowerSensiBus(true), &_CoresBusAgent);
	_SensiPowerReduce = MatrixCPU(nLineConstraint, _nAgent); // Gred
	_SensiPowerReduce.multiply(&SCGrid->getPowerSensiBusReduce(), &_CoresBusAgent);

	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;

}
StudyCase::StudyCase(int nAgent, float P0, float dP, float b, float db, float propCons)
{
	clock_t t = clock();
	SCAg = StudyCaseAgent(nAgent, P0, dP, b, db, propCons);
	createGrid(true);

	int nLineConstraint = SCGrid->getNLineConstraint();
	initMat();
	genCoresBusAgent(true);
	_SensiPower.multiply(&SCGrid->getPowerSensiBus(true), &_CoresBusAgent);
	_SensiPowerReduce = MatrixCPU(nLineConstraint, _nAgent); // Gred
	_SensiPowerReduce.multiply(&SCGrid->getPowerSensiBusReduce(), &_CoresBusAgent);

	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}
StudyCase::StudyCase(int nAgent, float P, float dP, float a, float da, float b, float db, float Gamma, float dGamma, float propCons, float propGenNFle, float propPro) {

	clock_t t = clock();


	SCAg = StudyCaseAgent(nAgent, P, dP, a, da, b, db, Gamma, dGamma, propCons, propGenNFle, propPro);

	createGrid(true);
	_nAgent = SCAg.getNagent();
	_CoresBusAgentLin = MatrixCPU(_nAgent, 1);
	for (int id = 0; id < _nAgent / 2; id++) {
		_CoresBusAgentLin.set(id, 0, 0);
	}
	for (int id = _nAgent / 2; id < _nAgent; id++) {
		_CoresBusAgentLin.set(id, 0, 1);
	}

	_nBus = SCGrid->getNBus();


	genCoresBusAgent();


	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;

}
StudyCase::StudyCase(int nAgent, float P, float dP, float Q, float dQ, float a, float da, float aQ, float daQ, float b, float db, float Gamma, float dGamma, float propCons, float propGenNFle, float propPro){

	clock_t t = clock();

	
	SCAg = StudyCaseAgent(nAgent, P, dP, Q, dQ, a, da, a, daQ, b, db,Gamma, dGamma, propCons, propGenNFle, propPro);
	
	createGrid(false);
	_nAgent = SCAg.getNagent();
	_CoresBusAgentLin = MatrixCPU(_nAgent, 1);
	for (int id = 0; id < _nAgent / 2; id++) {
		_CoresBusAgentLin.set(id, 0, 0);
	}
	for (int id = _nAgent / 2; id < _nAgent; id++) {
		_CoresBusAgentLin.set(id, 0, 1);
	}
	
	_nBus = SCGrid->getNBus();
	

	genCoresBusAgent();
	

	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;

}


StudyCase::StudyCase(const StudyCase& s)
{
	clock_t t = clock();
	SCAg = s.SCAg;

	DELETEB(SCGrid);

	createGrid(s.DC);
	

	if (s.DC) {
		*SCGrid = *s.SCGrid;
	}
	else {
		((StudyCaseACGrid*)SCGrid)->copy((StudyCaseACGrid*)s.SCGrid); //dynamic_cast<StudyCaseACGrid*>
	}

	

	
	_name = s._name;
	_Sbase = s._Sbase;
	toReduce = s.toReduce;
	_nAgent = s._nAgent;
	_nBus = s._nBus;
	_nLine = s._nLine;

	
	_CoresBusAgent = s._CoresBusAgent; // I
	_SensiPower = s._SensiPower; // G
	_SensiPowerReduce = s._SensiPowerReduce;
	_Distance = s._Distance;

	_CoresBusAgentLin = s._CoresBusAgentLin;
	_CoresAgentBusLin = s._CoresAgentBusLin;
	_CoresAgentBusLinBegin = s._CoresAgentBusLinBegin;
	_nAgentByBus = s._nAgentByBus;
	

	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}
StudyCase::StudyCase(std::string fileName)
{
	clock_t t = clock();
	SCAg = StudyCaseAgent(fileName);
	createGrid(true);
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;

}
StudyCase& StudyCase::operator= (const StudyCase& s)
{
	clock_t t = clock();
	SCAg = s.SCAg;
	DELETEB(SCGrid);


	createGrid(s.DC);

	if (s.DC) {
		* SCGrid = *s.SCGrid;
	}
	else {
		((StudyCaseACGrid *)SCGrid)->copy((StudyCaseACGrid*) s.SCGrid); //dynamic_cast<StudyCaseACGrid*>
	}
	

	_name = s._name;
	_Sbase = s._Sbase;
	toReduce = s.toReduce;

	_nAgent = s._nAgent;
	_nBus = s._nBus;
	_nLine = s._nLine;


	_CoresBusAgent = s._CoresBusAgent; // I
	_SensiPower = s._SensiPower; // G
	_SensiPowerReduce = s._SensiPowerReduce;
	_Distance = s._Distance;

	_CoresBusAgentLin = s._CoresBusAgentLin;
	_CoresAgentBusLin = s._CoresAgentBusLin;
	_CoresAgentBusLinBegin = s._CoresAgentBusLinBegin;
	_nAgentByBus = s._nAgentByBus;

	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
	return *this;
}

void StudyCase::UpdateP0(MatrixCPU* P0)
{
	SCAg.UpdateP0(P0);
}
void StudyCase::genBetaUniforme(float beta)
{
	SCAg.genBetaUniforme(beta);
}

void StudyCase::genBetaDistance(float s)
{
	int nAgent = SCAg.getNagent();
	if (_Distance.getNLin() != nAgent) {
		setDistance();
	}

	SCAg.genBetaDistance(s, &_Distance);
}

void StudyCase::genBetaDistanceByZone(MatrixCPU* Mats)
{
	int nAgent = SCAg.getNagent();
	if (_Distance.getNLin() != nAgent) {
		setDistance();
	}
	// verification que c'est l�gal : 
	std::cout << "work in progress" << std::endl;
	throw std::invalid_argument("WIP");
	// differencier DC ou AC
	// 
	
	SCAg.genBetaZone(Mats, &_Distance, SCGrid->getZones());
	
	
}

void StudyCase::setDistance(bool alreadyDefine, std::string name )
{
	int nAgent = SCAg.getNagent();
	int nLine = SCGrid->getNLine(true);
	_Distance = MatrixCPU(nAgent, nAgent);
	if (alreadyDefine) {
		_Distance.setFromFile(name);
	}
	else {
		std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
		for (int i = 0; i < nAgent; i++) {
			for (int j = i + 1; j < nAgent; j++) {
				float sum = 0;
				for (int k = 0; k < nLine; k++) {
					float p = _SensiPower.get(k, i) - _SensiPower.get(k, j);
					sum = sum + fabs(p);
				}
				_Distance.set(i, j, sum);
				_Distance.set(j, i, sum);
			}
		}
		_Distance.display();
		_Distance.saveCSV(name, mode, 0, " ");
	}
}

bool StudyCase::isAC() const
{
	return !DC;
}

bool StudyCase::isRadial() const
{
	if (DC) {
		return false;
	}
	return ((StudyCaseACGrid*)SCGrid)->radial;
}

bool StudyCase::isCurrentLimit() const
{
	if (DC) {
		return false;
	}
	return ((StudyCaseACGrid*)SCGrid)->isCurrentLimit();
}



void StudyCase::setLineLimitMin(float min)
{
	SCGrid->setLineLimitMin(min);
}

void StudyCase::setLineLimitRelaxation(float eps)
{
	SCGrid->setLineLimitRelaxation(eps);
}

///////////////////////////////////////////////////////////////////////////////
// Set StudyCase
///////////////////////////////////////////////////////////////////////////////

void StudyCase::Set29node()
{
	
	clock_t t = clock();
	SCAg.Set29node();
	createGrid(true);
		
	_nBus = SCGrid->getNBus();
	_nLine = SCGrid->getNLine(true);
	_nAgent = SCAg.getNagent();

	int nLineConstraint = SCGrid->getNLineConstraint();
	
	_CoresBusAgentLin = MatrixCPU(_nAgent, 1);
	for (int id = 0; id < _nAgent/2; id++) {
		_CoresBusAgentLin.set(id, 0, 0);
	}
	for (int id = _nAgent / 2; id < _nAgent; id++) {
		_CoresBusAgentLin.set(id, 0, 1);
	}/**/
	_SensiPower = MatrixCPU(_nLine, _nAgent); // G
	
	genCoresBusAgent(true);
	_SensiPower.multiply(&SCGrid->getPowerSensiBus(true), &_CoresBusAgent);


	_SensiPowerReduce = MatrixCPU(nLineConstraint, _nAgent); // Gred
	_SensiPowerReduce.multiply(&SCGrid->getPowerSensiBusReduce(), &_CoresBusAgent);


	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}


void StudyCase::Set39Bus(std::string path, bool alreadyDefine)
{
	clock_t t = clock();
	_CoresBusAgentLin = SCAg.Set29node();
	createGrid(true);
	SCGrid->Set39Bus(path, alreadyDefine);
	_nBus = SCGrid->getNBus();
	_nAgent = SCAg.getNagent();
	_nLine = SCGrid->getNLine(true);
	int nLineConstraint = SCGrid->getNLineConstraint();

	
	_SensiPower = MatrixCPU(_nLine, _nAgent); // G

	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	std::string fileName2 = path + "SensiPower39node.txt";
	std::string fileName3 = path + "SensiPowerReduce39node.txt";

	genCoresBusAgent(true);

	if (alreadyDefine) {
		_SensiPowerReduce = MatrixCPU(nLineConstraint, _nAgent);
		_SensiPower.setFromFile(fileName2);
		_SensiPowerReduce.setFromFile(fileName3);
	}
	else {

		_SensiPower.multiply(&SCGrid->getPowerSensiBus(true), &_CoresBusAgent);
		_SensiPowerReduce = MatrixCPU(nLineConstraint, _nAgent); // Gred
		_SensiPowerReduce.multiply(&SCGrid->getPowerSensiBusReduce(), &_CoresBusAgent);
		_SensiPower.saveCSV(fileName2, mode, 0, " ");
		_SensiPowerReduce.saveCSV(fileName3, mode, 0, " ");

	}
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
	
}
void StudyCase::SetAC39Bus(std::string path, bool alreadyDefine)
{
	clock_t t = clock();

	_CoresBusAgentLin = SCAg.Set29node(true);
	createGrid(false);
	((StudyCaseACGrid*) SCGrid)->SetAC39Bus(path, alreadyDefine);
	
	_nBus = SCGrid->getNBus();
	_nAgent = SCAg.getNagent();
	_nLine = SCGrid->getNLine(true);

	genCoresBusAgent();


	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}

void StudyCase::SetAC3Bus(std::string path)
{
	//case3.m
	_Sbase = 100; // MW car matlab
	clock_t t = clock();
	
	_CoresBusAgentLin = SCAg.Set3Bus(true);
	createGrid(false);
	
	((StudyCaseACGrid*)SCGrid)->SetAC3Bus();

	
	_nBus = SCGrid->getNBus(); // 3
	_nLine = SCGrid->getNLine(true); // 3
	_nAgent = SCAg.getNagent(); // 4

	genCoresBusAgent();


	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}

void StudyCase::SetAC2node() {
	_Sbase = 1; // MW 
	clock_t t = clock();
	
	_CoresBusAgentLin = SCAg.Set2node(true);
	createGrid(false);

	_nBus = SCGrid->getNBus(); 
	_nLine = SCGrid->getNLine(true);
	_nAgent = SCAg.getNagent(); 


	genCoresBusAgent();

	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}

void StudyCase::SetACFromFile(std::string name, std::string path)
{
	_CoresBusAgentLin = SCAg.SetACFromFile(name, path);
	createGrid(false);
	((StudyCaseACGrid*)SCGrid)->SetACFromFile(name, path);
	_nBus = SCGrid->getNBus();
	_nAgent = SCAg.getNagent();
	_nLine = SCGrid->getNLine(true);
	
	genCoresBusAgent();
	
}

void StudyCase::SetACFromFileSimplify(std::string name, std::string path)
{
	_CoresBusAgentLin = SCAg.SetACFromFileSimplify(name, path);
	createGrid(false);
	((StudyCaseACGrid*)SCGrid)->SetACFromFile(name, path);
	_nBus = SCGrid->getNBus();
	_nAgent = SCAg.getNagent();
	_nLine = SCGrid->getNLine(true);

	genCoresBusAgent();

}

void StudyCase::SetStudyCase(std::string path, std::string name, MatrixCPU* P0, bool alreadyDefine)
{
	clock_t t = clock();
	std::string fileName = path + "SensiPower" + name + ".txt";
	std::string fileName3 = path + "SensiPowerReduce" + name + ".txt";
	_name = name;
	//std::cout << "agent" << std::endl;
	SCAg.SetStudyCaseAgent(path, name, P0);
	_nBus = SCAg.getNCons();
	_nAgent = SCAg.getNagent();
	//std::cout << "grid" << std::endl;
	createGrid(true);

	SCGrid->SetStudyCaseDCGrid(path, name, _nBus, alreadyDefine);
	
	_nLine = SCGrid->getNLine(true);

	//std::cout << "link" << std::endl;
	_SensiPower = MatrixCPU(_nLine, _nAgent); // G
	_CoresBusAgent = MatrixCPU(_nBus, _nAgent); // I
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;


	if (alreadyDefine) {
		int _nLineConstraint = SCGrid->getNLineConstraint();
		_SensiPowerReduce = MatrixCPU(_nLineConstraint, _nAgent);
		_SensiPower.setFromFile(fileName);
		_SensiPowerReduce.setFromFile(fileName3);
	}
	else {
		MatrixCPU fileCoresBus = SCGrid->getfileCoresBus();
		MatrixCPU GenBus = SCAg.getGenBus();

		int idBusMax = fileCoresBus.max2();
		MatrixCPU fileBusAgent(idBusMax + 1, 1, -1); // si reste � -1, le bus n'existe pas

		//std::cout << "agent" << std::endl;
		for (int i = 0; i < _nBus; i++) {
			int bus = fileCoresBus.get(i, 0);
			fileBusAgent.set(bus, 0, i);
		}

		//std::cout << _nLineConstraint << std::endl;
		//std::cout << "agent" << std::endl;
		for (int i = 0; i < _nAgent; i++) {
			if (i < _nBus) { // le bus correspond directement pour les conso
				int bus = i;
				_CoresBusAgent.set(bus, i, 1);
			}
			else {
				int idGen = i - _nBus;
				int bus = fileBusAgent.get(GenBus.get(idGen, 0), 0);
				_CoresBusAgent.set(bus, i, 1);
			}
		}

		MatrixCPU SensiBusLine = SCGrid->getPowerSensiBus(true);
		MatrixCPU SensiBusLineReduce = SCGrid->getPowerSensiBusReduce();
		_SensiPower.multiply(&SensiBusLine, &_CoresBusAgent);

		_SensiPower.saveCSV(fileName, mode);
		_SensiPowerReduce.saveCSV(fileName3, mode);
	}
	//_SensiPower.display();
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}

void StudyCase::SetACStudyCaseFromInterface(StudyCaseInterface* interface){
	//std::cout << "creation agent" << std::endl;
	_CoresBusAgentLin = SCAg.SetFromInterface(interface, false);
	createGrid(false);
	if(interface->getL()>= 1){
		//std::cout << "creation grid" << std::endl;
		SCGrid->setFromInterface(interface);
	}
	_nBus = SCGrid->getNBus();
	_nAgent = SCAg.getNagent();
	_nLine = SCGrid->getNLine(true);

	genCoresBusAgent();
	//display(1);
	//std::cout << "Fin creation cas" << std::endl;
}

void StudyCase::SetDCStudyCaseFromInterface(StudyCaseInterface* interface){
	//std::cout << "creation agent" << std::endl;
	_CoresBusAgentLin = SCAg.SetFromInterface(interface, true);
	createGrid(true);
	if(interface->getL()>= 1){
		//std::cout << "creation grid" << std::endl;
		(SCGrid)->setFromInterface(interface);
	}
	_nBus = SCGrid->getNBus();
	_nAgent = SCAg.getNagent();
	_nLine = SCGrid->getNLine(true);
	
	genCoresBusAgent(true);
	computeSensiPower();
	//display(0);
	//std::cout << "Fin creation cas" << std::endl;
}


void StudyCase::SetEuropeTestFeeder(std::string path, int typeOfAgentGen, int beggining)
{
	switch (typeOfAgentGen)
	{
	case 0: // 55 consumers, 1 productor as in the File
		std::cout << " agent " << std::endl;
		_CoresBusAgentLin = SCAg.SetEuropeTestFeeder(path, beggining);
		break;
	default:
		throw std::invalid_argument("WIP not implemented");
		break;
	}
	
	createGrid(false);
	((StudyCaseACGrid*)SCGrid)->SetEuropeTestFeeder(path);
	_nBus = SCGrid->getNBus();
	_nLine = SCGrid->getNLine(true);
	_nAgent = SCAg.getNagent();


	genCoresBusAgent();
	computeSensiPower();

}

void StudyCase::Set3Bus(std::string path) {
	clock_t t = clock();
	std::cout << "set agent" << std::endl;
	_CoresBusAgentLin = SCAg.Set3Bus();
	
	createGrid(true);
	std::cout << "set grid" << std::endl;
	SCGrid->Set3Bus();


	_nBus = SCGrid->getNBus();
	_nLine = SCGrid->getNLine(true);
	_nAgent = SCAg.getNagent();

	int nLineConstraint = SCGrid->getNLineConstraint();


	_SensiPower = MatrixCPU(_nLine, _nAgent); // G
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	

	genCoresBusAgent(true);
	
	
	_SensiPower.multiply(&SCGrid->getPowerSensiBus(true), &_CoresBusAgent);
	_SensiPowerReduce = MatrixCPU(nLineConstraint, _nAgent); // Gred
	_SensiPowerReduce.multiply(&SCGrid->getPowerSensiBusReduce(), &_CoresBusAgent);
	
	
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}

void StudyCase::Set4Agent()
{
	clock_t t = clock();
	
	SCAg.Set4Agent();
	createGrid(true);
	_nBus = SCGrid->getNBus();
	_nLine = SCGrid->getNLine(true);
	_nAgent = SCAg.getNagent();

	int nLineConstraint = SCGrid->getNLineConstraint();

	_CoresBusAgentLin = MatrixCPU(_nAgent, 1);
	for (int id = 0; id < _nAgent / 2; id++) {
		_CoresBusAgentLin.set(id, 0, 0);
	}
	for (int id = _nAgent / 2; id < _nAgent; id++) {
		_CoresBusAgentLin.set(id, 0, 1);
	}
	_SensiPower = MatrixCPU(_nLine, _nAgent); // G

	genCoresBusAgent(true);

	_SensiPower.multiply(&SCGrid->getPowerSensiBus(true), &_CoresBusAgent);


	_SensiPowerReduce = MatrixCPU(nLineConstraint, _nAgent); // Gred
	_SensiPowerReduce.multiply(&SCGrid->getPowerSensiBusReduce(), &_CoresBusAgent);

	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}

void StudyCase::Set2node() 
{
	clock_t t = clock();
	
	SCAg.Set2node();
	createGrid(true);
	_nBus = SCGrid->getNBus();
	_nLine = SCGrid->getNLine(true);
	_nAgent = SCAg.getNagent();

	int nLineConstraint = SCGrid->getNLineConstraint();

	_CoresBusAgentLin = MatrixCPU(_nAgent, 1);
	for (int id = 0; id < _nAgent / 2; id++) {
		_CoresBusAgentLin.set(id, 0, 0);
	}
	for (int id = _nAgent / 2; id < _nAgent; id++) {
		_CoresBusAgentLin.set(id, 0, 1);
	}
	_SensiPower = MatrixCPU(_nLine, _nAgent); // G

	genCoresBusAgent(true);

	_SensiPower.multiply(&SCGrid->getPowerSensiBus(true), &_CoresBusAgent);

	_SensiPower.display();

	_SensiPowerReduce = MatrixCPU(nLineConstraint, _nAgent); // Gred
	_SensiPowerReduce.multiply(&SCGrid->getPowerSensiBusReduce(), &_CoresBusAgent);


	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}

void StudyCase::Set4nodeBis(std::string path)
{
	// cas d'�tude pour simuler le cas d'EVA pendant son stage
	clock_t t = clock();
	
	
	
	_CoresBusAgentLin = SCAg.Set4nodeBis();
	createGrid(true);
	SCGrid->Set4nodeBis(path);
	_nBus = SCGrid->getNBus();
	_nLine = SCGrid->getNLine(true);
	_nAgent = SCAg.getNagent();

	int nLineConstraint = SCGrid->getNLineConstraint();
	
	
	_SensiPower = MatrixCPU(_nLine, _nAgent); // G
	_SensiPowerReduce = MatrixCPU(nLineConstraint, _nAgent); // Gred
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;



	genCoresBusAgent(true);
	
	_SensiPower.multiply(&SCGrid->getPowerSensiBus(true), &_CoresBusAgent);
	_SensiPowerReduce.multiply(&SCGrid->getPowerSensiBusReduce(), &_CoresBusAgent);

	
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}

void StudyCase::Set2nodeConstraint(float lim)
{
	clock_t t = clock();

	SCAg.Set2node();
	createGrid(true);
	SCGrid->Set2nodeConstraint(lim);
	_nBus = SCGrid->getNBus();
	_nLine = SCGrid->getNLine(true);
	_nAgent = SCAg.getNagent();

	int nLineConstraint = SCGrid->getNLineConstraint();
	
	_SensiPower = MatrixCPU(_nLine, _nAgent); // G
	_CoresBusAgentLin = MatrixCPU(_nAgent, 1);
	for (int id = 0; id < _nAgent / 2; id++) {
		_CoresBusAgentLin.set(id, 0, 0);
	}
	for (int id = _nAgent / 2; id < _nAgent; id++) {
		_CoresBusAgentLin.set(id, 0, 1);
	}
	genCoresBusAgent(true);
	

	_SensiPower.multiply(&SCGrid->getPowerSensiBus(true), &_CoresBusAgent);
	_SensiPowerReduce = MatrixCPU(nLineConstraint, _nAgent); // Gred
	_SensiPowerReduce.multiply(&SCGrid->getPowerSensiBusReduce(), &_CoresBusAgent);

}

void StudyCase::SetEuropeP0(const std::string& path, MatrixCPU* P0, bool alreadyDefine)
{
	clock_t t = clock();
	
	_name = "Europe";
	
	SCAg.SetEuropeP0(path, P0);
	createGrid(true);
	
	SCGrid->SetEuropeP0(path, alreadyDefine);
	_nBus = SCGrid->getNBus();
	_nLine = SCGrid->getNLine(true);
	_nAgent = SCAg.getNagent();

	int nLineConstraint = SCGrid->getNLineConstraint();
	int nCons = SCAg.getNCons();
	//genBetaUniforme(0);

	_SensiPower = MatrixCPU(_nLine, _nAgent); // G
	_SensiPowerReduce = MatrixCPU(nLineConstraint, _nAgent);

	_CoresBusAgentLin = MatrixCPU(_nAgent, 1); // Ilin
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	std::string fileName = path + "SensiPowerEurope.txt";
	std::string fileName3 = path + "SensiPowerReduceEurope.txt";
	
	if (alreadyDefine) {
		_SensiPower.setFromFile(fileName);
		_SensiPowerReduce.setFromFile(fileName3);
	}
	else {
		MatrixCPU fileCoresBus = SCGrid->getfileCoresBus();
		MatrixCPU GenBus = SCAg.getGenBus();

		int idBusMax = fileCoresBus.max2();
		MatrixCPU fileBusAgent(idBusMax + 1, 1, -1); // si reste � -1, le bus n'existe pas
		

		for (int i = 0; i < _nBus; i++) {
			int bus = fileCoresBus.get(i, 0);
			fileBusAgent.set(bus, 0, i);
		}
	
		//std::cout << _nLineConstraint << std::endl;
		for (int i = 0; i < _nAgent; i++) {
			if (i < nCons) { // le bus correspond directement pour les conso
				int bus = i;
				_CoresBusAgentLin.set(i, 0, bus);
			}
			else {
				int idGen = i - nCons;
				int bus = fileBusAgent.get(GenBus.get(idGen, 0), 0);
				_CoresBusAgentLin.set(i, 0, bus);
			}
		}

		genCoresBusAgent(true);
		
		_SensiPower.multiply(&SCGrid->getPowerSensiBus(true), &_CoresBusAgent);
		_SensiPowerReduce.multiply(&SCGrid->getPowerSensiBusReduce(), &_CoresBusAgent);
		
		_SensiPower.saveCSV(fileName, mode);
		_SensiPowerReduce.saveCSV(fileName3, mode);
	}
	//_SensiPower.display();
	t = clock() - t;
	//_timeInit = (float)t / CLOCKS_PER_SEC;
}

void StudyCase::SetEuropeP0WithoutConstraint(const std::string& path, MatrixCPU* P0)
{
	clock_t t = clock();

	_name = "Europe";

	SCAg.SetEuropeP0(path, P0);
	createGrid(true);
	_nBus = SCGrid->getNBus();
	_nLine = SCGrid->getNLine(true);
	_nAgent = SCAg.getNagent();

	int nLineConstraint = SCGrid->getNLineConstraint();

	_CoresBusAgentLin = MatrixCPU(_nAgent, 1);
	for (int id = 0; id < _nAgent / 2; id++) {
		_CoresBusAgentLin.set(id, 0, 0);
	}
	for (int id = _nAgent / 2; id < _nAgent; id++) {
		_CoresBusAgentLin.set(id, 0, 1);
	}
	_SensiPower = MatrixCPU(_nLine, _nAgent); // G

	genCoresBusAgent(true);

	_SensiPower.multiply(&SCGrid->getPowerSensiBus(true), &_CoresBusAgent);


	_SensiPowerReduce = MatrixCPU(nLineConstraint, _nAgent); // Gred
	_SensiPowerReduce.multiply(&SCGrid->getPowerSensiBusReduce(), &_CoresBusAgent);


	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}


///////////////////////////////////////////////////////////////////////////////
// Generator
///////////////////////////////////////////////////////////////////////////////


void StudyCase::genGrid(int _nBus, int _nMajorLine, int _minorLine, float ReacMajor, float DeltaReacMajor, float ReacMinor, float DeltaReacMinor, float LlimitMajor, float dLlimitMajor, float LlimitMinor, float dLlimitMinor)
{
	std::cout << "work in progress" << std::endl;
	throw std::invalid_argument("WIP");
}

void StudyCase::genGridBT(int nBus, int Nbranch, int Ndeep, float length, float dlength)
{
	_timeInit = clock();
	createGrid(false);
	((StudyCaseACGrid*)SCGrid)->genGridBT(nBus, Nbranch, Ndeep, length, dlength);
	_nBus = nBus;
	_nLine = nBus - 1;
	_name = "gridBT";

	_timeInit = (float) (clock() - _timeInit)/ CLOCKS_PER_SEC;
	
	
	/*
	// melange grid/agent

	MatrixCPU _SensiPower; // G = A*I
	MatrixCPU _CoresBusAgent; // I
	MatrixCPU _SensiPowerReduce; // Gred
	MatrixCPU _Distance; // Dij = sum |Pl|

	// Permettre de calculer W0
	MatrixCPU _CoresBusAgentLin;
	MatrixCPU _CoresAgentBusLin;
	MatrixCPU _CoresAgentBusLinBegin;
	MatrixCPU _nAgentByBus;*/
}

void StudyCase::genGridBTSpecial(int nBus, int Nbranch, int Ndeep, float length, float dlength, RadialType type)
{
	_timeInit = clock();
	createGrid(false);
	((StudyCaseACGrid*)SCGrid)->genGridBTSpecial(nBus, Nbranch, Ndeep, length, dlength, type);
	_nBus = nBus;
	_nLine = nBus - 1;
	_name = "gridBT";

	_timeInit = (float)(clock() - _timeInit) / CLOCKS_PER_SEC;

	/*
	// melange grid/agent

	MatrixCPU _SensiPower; // G = A*I
	MatrixCPU _CoresBusAgent; // I
	MatrixCPU _SensiPowerReduce; // Gred
	MatrixCPU _Distance; // Dij = sum |Pl|

	// Permettre de calculer W0
	MatrixCPU _CoresBusAgentLin;
	MatrixCPU _CoresAgentBusLin;
	MatrixCPU _CoresAgentBusLinBegin;
	MatrixCPU _nAgentByBus;*/
}

void StudyCase::genGridHTB(int nBus, int nLine, int dnLine, float length, float dlength)
{
	_timeInit = clock();
	createGrid(false);
	((StudyCaseACGrid*)SCGrid)->genGridHTB(nBus, nLine, dnLine, length, dlength);
	_nBus = nBus;
	_nLine = nLine;
	_name = "gridHTB";

	_timeInit = (float)(clock() - _timeInit) / CLOCKS_PER_SEC;
}


void StudyCase::genGridFromFile(std::string path, bool alreadyDefine)
{
	_timeInit = clock();
	createGrid(true);
	SCGrid->genGridFromFile(path, alreadyDefine);
	_nBus = SCGrid->getNBus();
	_nLine = SCGrid->getNLine();
	_timeInit = (float)(clock() - _timeInit) / CLOCKS_PER_SEC;
}

void StudyCase::genAgents(int nAgent, float propCons, float Pconso, float dPconso, float bProd, float dbProd, float Pprod, float dPprod, float Gamma, float dGamma)
{
	clock_t t = clock();
	SCAg.genAgents(nAgent, propCons, Pconso, dPconso, bProd, dbProd, Pprod, dPprod, Gamma, dGamma);
	_nAgent = nAgent;
	t = clock() - t;
	_timeInit += (float)t / CLOCKS_PER_SEC;
}

void StudyCase::genAgentsAC(int nAgent, float propCons, float propGenNFle, float Pconso, float dPconso, float bProd, float dbProd, float dQconso, float Pprod, float dPprod, float Gamma, float dGamma)
{
	clock_t t = clock();
	SCAg.genAgentsAC(nAgent, propCons, propGenNFle, Pconso, dPconso, bProd, dbProd, dQconso, Pprod, dPprod, Gamma, dGamma);
	_nAgent = nAgent;
	t = clock() - t;
	_timeInit += (float)t / CLOCKS_PER_SEC;
}

void StudyCase::genAgentsFullRandom(int nAgent, float aMin, float aMax, float P0Min, float P0Max, float gammaMin, float gammaMax, float propConsoMin, float propConsoMax, float borneMin, float borneMax)
{
	clock_t t = clock();
	SCAg.genAgentsFullRandom(nAgent, aMin, aMax, P0Min, P0Max, gammaMin, gammaMax, propConsoMin, propConsoMax, borneMin, borneMax);
	
	_nAgent = nAgent;
	t = clock() - t;
	_timeInit += (float)t / CLOCKS_PER_SEC;


}

void StudyCase::genLinkGridAgent()
{
	//std::cout << "gen link between agent" << std::endl;
	_nAgent = SCAg.getNagent();
	int nLineConstraint;
	
	_nBus = SCGrid->getNBus();
	_nLine = SCGrid->getNLine(true);

	if (DC) {
		nLineConstraint = SCGrid->getNLineConstraint();
		_CoresBusAgent = MatrixCPU(_nBus, _nAgent);
	}
	else {
		nLineConstraint = _nLine;
		_nAgentByBus = MatrixCPU(_nBus, 1);
		_nAgentByBus.set(0, 0, 1); // loss agent
		_CoresAgentBusLinBegin = MatrixCPU(_nBus, 1);
		_CoresBusAgentLin = MatrixCPU(_nAgent, 1);
		_CoresAgentBusLin = MatrixCPU(_nAgent, 1);
	}
	
	for (int n = 1; n < _nAgent; n++) {
		int bus =  rand() % _nBus;
		
		if (DC) {
			_CoresBusAgent.set(bus, n, 1);
		}
		else {
			_CoresBusAgentLin.set(n, 0, bus);
			_nAgentByBus.increment(bus, 0, 1);
		}
	}
	
	if (DC) {
		_SensiPower = MatrixCPU(_nLine, _nAgent); // G
		_SensiPowerReduce = MatrixCPU(nLineConstraint, _nAgent); // Gred
		_SensiPower.multiply(&SCGrid->getPowerSensiBus(true), &_CoresBusAgent);
		_SensiPowerReduce.multiply(&SCGrid->getPowerSensiBusReduce(), &_CoresBusAgent);
	}
	else {
		int debut = 0;
		int* decompteAgent = new int[_nBus];
		for (int b = 0; b < _nBus; b++) {
			_CoresAgentBusLinBegin.set(b, 0, debut);
			debut += _nAgentByBus.get(b, 0);
			decompteAgent[b] = 0;
		}
		for (int n = 0; n < _nAgent; n++) {
			int bus = _CoresBusAgentLin.get(n, 0);
			int indice = _CoresAgentBusLinBegin.get(bus, 0) + decompteAgent[bus];
			decompteAgent[bus]++;
			_CoresAgentBusLin.set(indice, 0, n);
		}
		DELETEA(decompteAgent);

	}
	//std::cout << "_nLine constraint " <<  nLineConstraint << std::endl;
	//std::cout << "Fin gen link between agent" << std::endl;
}

void StudyCase::computeSensiPower()
{
	if (!DC) {
		//std::cout << "Warning : this Study Case is AC" << std::endl;
		
		((StudyCaseACGrid*)SCGrid)->genDCGridFromAC();
	}
	
	
	int nLineConstraint = SCGrid->getNLineConstraint();
	//std::cout << "nLine Contraint" << nLineConstraint << std::endl;
	_SensiPower = MatrixCPU(_nLine, _nAgent); // G
	_SensiPowerReduce = MatrixCPU(nLineConstraint, _nAgent); // Gred

	

	_SensiPower.multiply(&SCGrid->getPowerSensiBus(true), &_CoresBusAgent);
	_SensiPowerReduce.multiply(&SCGrid->getPowerSensiBusReduce(), &_CoresBusAgent);
	
}

void StudyCase::genLineLimit(int nLine, float limit, float dlLimit)
{
	DC = true;
	SCGrid->genLineLimit(nLine, limit, dlLimit);
}

void StudyCase::genDCGridFromAC()
{
	if (DC) {
		throw std::invalid_argument("genDCFromAC : must be AC grid");
	}
	else {
		//std::cout << _nBus << " " << _nAgent << std::endl;
		_CoresBusAgent = MatrixCPU(_nBus, _nAgent);

		for (int n = 0; n < _nAgent; n++) {
			int bus = _CoresBusAgentLin.get(n, 0);
			_CoresBusAgent.set(bus, n, 1);
		}
		//_CoresBusAgent.display();
		computeSensiPower();
	}
}

///////////////////////////////////////////////////////////////////////////////
// Getter
///////////////////////////////////////////////////////////////////////////////
float StudyCase::getSbase() const
{
	return _Sbase;
}

int StudyCase::getLastBus() const
{
	if (DC) {
		throw std::invalid_argument("getLastBus : Must be AC grid");
		return 0;
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getLastBus();
	}
}


MatrixCPU StudyCase::getBeta() const
{
	return SCAg.getBeta();
}

MatrixCPU StudyCase::getC() const
{
	return SCAg.getC();
}

MatrixCPU StudyCase::geta() const
{
	return SCAg.geta();
}

MatrixCPU StudyCase::getb() const
{
	return SCAg.getb();
}

MatrixCPU StudyCase::getUb() const
{
	return SCAg.getUb();
}

MatrixCPU StudyCase::getLb() const
{
	return SCAg.getLb();
}

MatrixCPU StudyCase::getPmin() const
{
	return SCAg.getPmin();
}

MatrixCPU StudyCase::getPmax() const
{
	return SCAg.getPmax();
}

MatrixCPU StudyCase::getNvoi() const
{
	return SCAg.getNvoi();
}


MatrixCPU StudyCase::getPowerSensi() const
{
	if (toReduce) {
		return _SensiPowerReduce;
	}
	else {
		return _SensiPower;
	}
	
}

MatrixCPU StudyCase::getLineLimit() const
{
	return SCGrid->getLineLimit();
	
}

MatrixCPU StudyCase::getCurrentLimit() const
{
	if (DC) {
		throw std::invalid_argument("getCurrentLimit : Must be AC grid");
		return MatrixCPU();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getCurrentLimit();
	}
}

MatrixCPU StudyCase::getCoresLineBus(bool force) const
{
	return SCGrid->getCoresLineBus(force);
}

MatrixCPU StudyCase::getCoresBusAgent() const
{
	return _CoresBusAgent;
}

MatrixCPU StudyCase::getCoresBusAgentLin() const
{
	return _CoresBusAgentLin;
}

MatrixCPU StudyCase::getCoresAgentBusLin() const
{
	return _CoresAgentBusLin;
}

MatrixCPU StudyCase::getCoresAgentBusLinBegin() const
{
	return _CoresAgentBusLinBegin;
}

MatrixCPU StudyCase::getNagentByBus() const
{
	return _nAgentByBus;
}

MatrixCPU StudyCase::getLineSuceptance() const
{
	if (DC) {
		throw std::invalid_argument("getLineSuceptance : Must be AC grid");
		return MatrixCPU();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getLineSuceptance();
	}

	
}

MatrixCPU StudyCase::getLineReactance() const
{
	if (DC) {
		throw std::invalid_argument("getLineReactance : Must be AC grid");
		return MatrixCPU();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getLineReactance();
	}
}

MatrixCPUD StudyCase::getLineSuceptanceD() const
{
	if (DC) {
		throw std::invalid_argument("getLineSuceptanceD : Must be AC grid");
		return MatrixCPUD();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getLineSuceptanceD();
	}
}

MatrixCPUD StudyCase::getLineReactanceD() const
{
	if (DC) {
		throw std::invalid_argument("getLineReactanceD :Must be AC grid");
		return MatrixCPUD();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getLineReactanceD();
	}
}

MatrixCPU StudyCase::getUpperBound() const
{
	if (DC) {
		throw std::invalid_argument("getUpperBound : Must be AC grid");
		return MatrixCPU();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getUpperBound();
	}
}

MatrixCPU StudyCase::getLowerBound() const
{
	if (DC) {
		throw std::invalid_argument("getLowerBound : Must be AC grid");
		return MatrixCPU();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getLowerBound();
	}
}

MatrixCPU StudyCase::getPobj()
{
	return SCAg.getPobj();
}

MatrixCPUD StudyCase::getPobjD()
{

	return SCAg.getPobjD();
}


MatrixCPUD StudyCase::getSolPF() const
{
	if (DC) {
		throw std::invalid_argument("getSolPF : Must be AC grid");
		return MatrixCPUD();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getSolPF();
	}
}

double StudyCase::getV0() const
{
	if (DC) {
		return 1;
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getV0();
	}
}

double StudyCase::gettheta0() const
{
	if (DC) {
		return 0;
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->gettheta0();
	}
}

MatrixCPU StudyCase::getZsRe() const
{
	if (DC) {
		throw std::invalid_argument("getZsRe : Must be AC grid");
		return MatrixCPU();
		
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getZsRe();
	}
}

MatrixCPU StudyCase::getZsImag() const
{
	if (DC) {
		throw std::invalid_argument("getZsImag : Must be AC grid");
		return MatrixCPU();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getZsImag();
	}
}

MatrixCPU StudyCase::getYd() const
{
	if (DC) {
		throw std::invalid_argument("getYd : Must be AC grid");
		return MatrixCPU();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getYd();
	}
}

MatrixCPU StudyCase::getGlin() const
{
	if (DC) {
		throw std::invalid_argument("getGlin: Must be AC grid");
		return MatrixCPU();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getGlin();
	}
}

MatrixCPU StudyCase::getBlin() const
{
	if (DC) {
		throw std::invalid_argument("getBlin : Must be AC grid");
		return MatrixCPU();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getBlin();
	}
}

MatrixCPU StudyCase::getGlin2() const
{
	if (DC) {
		throw std::invalid_argument("getGlin2 : Must be AC grid");
		return MatrixCPU();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getGlin2();
	}
}

MatrixCPU StudyCase::getBlin2() const
{
	if (DC) {
		throw std::invalid_argument("getBlin2 : Must be AC grid");
		return MatrixCPU();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getBlin2();
	}
}

MatrixCPU StudyCase::getVoltageInit() const
{
	if (DC) {
		throw std::invalid_argument("getVoltageInit : Must be AC grid");
		return MatrixCPU();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getVoltageInit();
	}
}

MatrixCPUD StudyCase::getVoltageInitD() const
{
	if (DC) {
		throw std::invalid_argument("getVoltageInitD :Must be AC grid");
		return MatrixCPUD();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getVoltageInitD();
	}
}

MatrixCPUD StudyCase::getGlinD() const
{
	if (DC) {
		throw std::invalid_argument("getGlinD :Must be AC grid");
		return MatrixCPUD();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getGlinD();
	}
}

MatrixCPUD StudyCase::getBlinD() const
{
	if (DC) {
		throw std::invalid_argument("getBlinD : Must be AC grid");
		return MatrixCPUD();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getBlinD();
	}
}

MatrixCPU StudyCase::getNLines() const
{
	if (DC) {
		throw std::invalid_argument("getNLines : Must be AC grid");
		return MatrixCPU();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getNLines();
	}
}

MatrixCPU StudyCase::getNLinesBegin() const
{
	if (DC) {
		throw std::invalid_argument("getNLinesBegin : Must be AC grid");
		return MatrixCPU();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getNLinesBegin();
	}
}

MatrixCPU StudyCase::getCoresBusLin() const
{
	if (DC) {
		throw std::invalid_argument("getCoresBusLin :Must be AC grid");
		return MatrixCPU();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getCoresBusLin();
	}
}

MatrixCPU StudyCase::getCoresVoiLin() const
{
	if (DC) {
		throw std::invalid_argument("getCoresVoiLin : Must be AC grid");
		return MatrixCPU();
	}
	else {
		return ((StudyCaseACGrid*)SCGrid)->getCoresVoiLin();
	}
}

float StudyCase::getTimeInit() const
{
	return _timeInit;
}

int StudyCase::getNagent() const
{
	return SCAg.getNagent();
}

int StudyCase::getNCons() const
{
	return SCAg.getNCons();
}

int StudyCase::getNLine(bool force) const
{
	return SCGrid->getNLine(force);
}

int StudyCase::getNBus() const
{
	return SCGrid->getNBus();
}

MatrixCPU StudyCase::getVoisin(int agent) const
{
	return SCAg.getVoisin(agent);
}
Agent StudyCase::getAgent(int agent) const
{
	return SCAg.getAgent(agent);
}

std::string StudyCase::getName() const
{
	return _name;
}

void StudyCase::removeLink(int i, int j)
{
	SCAg.removeLink(i, j);
}


void StudyCase::addLink(int i, int j)
{
	SCAg.addLink(i, j);
}

void StudyCase::setLineLimit(int line, float limit)
{
	SCGrid->setLineLimit(line, limit);
	if (DC) {
		_SensiPowerReduce = MatrixCPU(SCGrid->getNLineConstraint(), SCAg.getNagent()); // Gred
		_SensiPowerReduce.multiply(&SCGrid->getPowerSensiBusReduce(), &_CoresBusAgent);
	}
	
}

Agent StudyCase::removeAgent(int agent)
{
	
	return SCAg.removeAgent(agent);

}

void StudyCase::restoreAgent(Agent& agent, bool all) {

	SCAg.restoreAgent(agent, all);
}

void StudyCase::saveCSV(const std::string& fileName, bool all)
{
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;

	SCAg.saveCSV(fileName, all);
	SCGrid->saveCSV(fileName);
	_SensiPowerReduce.saveCSV(fileName, mode);

}

void StudyCase::nextStepPobj()
{
	SCAg.nextStepPobj();
}


void StudyCase::display(int type) const
{
	if (type == 0) {
		SCAg.display();
		std::cout << "Agent to bus corespondance: " << std::endl;
		_CoresBusAgentLin.display();
		std::cout << "Bus to agent corespondance: " << std::endl;
		_CoresAgentBusLin.display();
	}
	else if (type == 1) {
		SCGrid->display();
		
	}
	else if (type == 2) {
		std::cout << "Grid Sensibility : " << std::endl;
		if (toReduce) {
			_SensiPowerReduce.display();
		}
		else {
			_SensiPower.display();
		}
	}
	
}

void StudyCase::displayLineCores(MatrixCPU* g, bool all)
{
	
	
	SCGrid->displayLineCores(g, all);
}

StudyCase::~StudyCase()
{
#ifdef DEBUG_DESTRUCTOR
	std::cout << "case destructor" << std::endl;
#endif
	DELETEB(SCGrid);

}


