#include "../head/StudyCaseDCGrid.h"
#include "../head/StudyCaseDCGrid.cuh"


void StudyCaseDCGrid::setGridFromFile(const std::string& path, MatrixCPU* fileCoresBus)
{
	MatrixCPU matFile(_nLine, 4);
	matFile.setFromFile(path);
	_nLineConstraint = 0;
	for (int i = 0; i < _nLine; i++) {
		int nodeFromFile = matFile.get(i, 0);
		int nodeToFile = matFile.get(i, 1);
		float react = matFile.get(i, 2); 
		if (react == 100000) { // cas pas de donn�e dans le r�seau europ�en
			react = 100;// que faire de ces "non" donn�es ?
		}
		float limit = matFile.get(i, 3) / _Sbase;

		int nodeFrom = fileCoresBus->get(nodeFromFile, 0);
		int nodeTo = fileCoresBus->get(nodeToFile, 0);

		//std::cout << " Ligne numero " << i << " entre bus " << nodeFromFile << " et " << nodeToFile << " dans le fichier mais en vrai c'est entre " << nodeFrom << " et " << nodeTo << " limite " << limit << " react " << react << std::endl;
		//std::cout << " Ligne numero " << i << " entre bus " << nodeFromFile << " et " << nodeToFile << " dans le fichier mais en vrai c'est entre " << nodeFrom << " et " << nodeTo << std::endl;
		
		_LineImpedance.set(i, i, react);
		_CoresBusLine.set(nodeFrom, i, 1);
		_CoresBusLine.set(nodeTo, i, -1);
		_CoresLineBus.set(i, 0, nodeFrom);
		_CoresLineBus.set(i, 1, nodeTo);
		if (limit > 0) {
			_nLineConstraint++;
			_lineLimits.set(i, 0, limit);
			_indiceLineConstraint.push_back(i);
		}
		else {
			_lineLimits.set(i, 0, LINELIMITMAX);
			_indiceLineNonConstraint.push_back(i);
		}
	}
}
void StudyCaseDCGrid::setBusFromFile(const std::string& path, MatrixCPU* fileCoresBus)
{
	int zone = 0;
	std::ifstream myfile(path, std::ios::in);
	bool found = false;
	int indice = zone;
	if (myfile)
	{
		for (int i = 0; i < _nBus; i++) {		
			int idAgent;
			int idBus;
			std::string country;
			myfile >> idAgent;
			myfile >> idBus;
			myfile >> country;
			fileCoresBus->set(idAgent, 0, idBus);
			
			found = false;
			indice = zone;
			for (int j = 0; j < _nameZone.size(); j++) {
				std::string value = _nameZone[j];
				if (!value.compare(country)) {
					found = true;
					indice = j;
					break;
				}
			}
			if (!found) {
				_nameZone.push_back(country);
				zone++;
			}
			_zoneBus.set(idAgent, 0, indice);

		}
		myfile.close();
	}
	else {
		throw std::invalid_argument("can't open this file");
	}
}
void StudyCaseDCGrid::CalcGridSensi()
{
	
	MatrixCPU temp1(_nLine, _nBus); // BC^T 
	MatrixCPU temp2(_nBus, _nBus); // CBC^T
	MatrixCPU temp3(_nBus, _nBus); // (CBC^T)^-1 avec mise � 0 de la ligne et colonne du noued de ref
	MatrixCPU temp33(_nBus - 1, _nBus - 1); //(CBC ^ T) ^ -1 sans la ligne et colonne du noued de ref
	MatrixCPU temp22(_nBus - 1, _nBus - 1); // on enl�ve la ligne du noeud de ref�rence
	MatrixCPU result(_nBus - 1, _nBus - 1);
	
	//MatrixCPU identity(_nBus - 1, _nBus - 1);
	//identity.setEyes(1);
	//std::cout << "nBus " << _nBus << " nLine " << _nLine <<  std::endl;



	temp1.multiplyTrans(&_LineImpedance, &_CoresBusLine);
	temp2.multiply(&_CoresBusLine, &temp1);
	
	temp2.getBloc(&temp22, 1, _nBus, 1, _nBus);

	//temp22.display();
	

	if (_invertMethod==1) {
		
		MatrixGPU temp33GPU(_nBus-1, _nBus-1, 0, 1);
		
		MatrixGPU temp22GPU(temp22,1);	
		temp33GPU.invertGaussJordan(&temp22GPU);
		temp33GPU.toMatCPU(temp33);


		//result.multiply(&temp33,&temp22);
		//float err = result.distance2(&identity);

		//std::cout << "err GPU " << err << std::endl;


	}
	else if(_invertMethod==2)
	{
		temp33.invertGaussJordan(&temp22);
		//result.multiply(&temp33, &temp22);
		//float err = result.distance2(&identity);

		//std::cout << "err CPU " << err << std::endl;
	}
	else {
	
		temp33.invertGaussJordan(&temp22);
		//result.multiply(&temp33, &temp22);
		//float err = result.distance2(&identity);

		//std::cout << "err Eigen " << err << std::endl;
	}

	temp3.setBloc(1, _nBus, 1, _nBus, &temp33);
	_SensiBusLine.multiply(&temp1, &temp3);

	

}
void StudyCaseDCGrid::ReduceSensi()
{
	
	_lineLimitsReduce = MatrixCPU(_nLineConstraint, 1);
	_SensiBusLineReduce = MatrixCPU(_nLineConstraint, _nBus); // Ared
	_CoresLineBusReduce = MatrixCPU(_nLineConstraint, 2);
	_indiceLineConstraint.clear();
	_indiceLineNonConstraint.clear();
	int indice = 0;
	for (int i = 0; i < _nLine; i++) {
		float lim = _lineLimits.get(i, 0);
		if (lim != 0 && lim != LINELIMITMAX) {
			_indiceLineConstraint.push_back(i);
			_lineLimitsReduce.set(indice, 0, lim);
			_CoresLineBusReduce.set(indice, 0, _CoresLineBus.get(i, 0));
			_CoresLineBusReduce.set(indice, 1, _CoresLineBus.get(i, 1));
			for (int j = 0; j < _nBus; j++) {
				_SensiBusLineReduce.set(indice, j, _SensiBusLine.get(i, j));
			}
			indice++;
		}
		else {
			_indiceLineNonConstraint.push_back(i);
			_lineLimits.set(i, 0, LINELIMITMAX);
		}
		
	}
}

float StudyCaseDCGrid::rand1() const
{
	float a = (float)(rand()) / ((float)(RAND_MAX));
	return a;
}
int StudyCaseDCGrid::randab(int a, int b) const
{
	return a + (rand() % (b - a));
}
int StudyCaseDCGrid::getNFileline(std::string nameFile)
{
	int number_of_lines = 0;
	std::string line;
	std::ifstream myfile(nameFile);

	while (std::getline(myfile, line))
		++number_of_lines;
	return number_of_lines;
}


StudyCaseDCGrid::StudyCaseDCGrid()
{
	
	_nBus = 2;
	_nLine = 1;
	_nLineConstraint = 0;
	float x = -0.01;
	
	_LineImpedance = MatrixCPU(_nLine, _nLine, -1 / x);
	_CoresLineBus = MatrixCPU(_nLine, 2);
	_CoresLineBus.set(0, 0, 0);
	_CoresLineBus.set(0, 1, 1);
	_CoresBusLine = MatrixCPU(_nBus, _nLine);
	_CoresBusLine.set(0, 0, 1);
	_CoresBusLine.set(1, 0, -1);
	_SensiBusLine = MatrixCPU(_nLine, _nBus); 
	_lineLimits = MatrixCPU(_nLine, 1, 0);
	
	//std::cout << "gridSensi " << std::endl;
	CalcGridSensi();
	//std::cout << "reduce " << std::endl;
	ReduceSensi();
	 
}
StudyCaseDCGrid::StudyCaseDCGrid(const StudyCaseDCGrid& s)
{
	clock_t t = clock();

	_nLine = s._nLine;
	_nBus = s._nBus;
	_nLineConstraint = s._nLineConstraint;
	_name = s._name;
	toReduce = s.toReduce;
	_Zbase = s._Zbase;
	_Sbase = s._Sbase;

	_lineLimitsReduce = s._lineLimitsReduce;
	_SensiBusLine = s._SensiBusLine;
	_SensiBusLineReduce = s._SensiBusLineReduce;



	_LineImpedance = s._LineImpedance; // B
	_CoresBusLine = s._CoresBusLine; // C
	_lineLimits = s._lineLimits; // l

	// min
	_lineLimitsChange = s._lineLimitsChange;
	lineMin = s.lineMin;
	lineoffset = s.lineoffset;

	_CoresLineBus = s._CoresLineBus;



	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}
StudyCaseDCGrid& StudyCaseDCGrid::operator= (const StudyCaseDCGrid& s)
{
	clock_t t = clock();
	//std::cout << " Copie egalite  DC " << std::endl;
	_nLine = s._nLine;
	_nBus = s._nBus;
	_nLineConstraint = s._nLineConstraint;
	_name = s._name;
	toReduce = s.toReduce;
	_Zbase = s._Zbase;
	_Sbase = s._Sbase;

	_lineLimitsReduce = s._lineLimitsReduce;
	_SensiBusLine = s._SensiBusLine;
	_SensiBusLineReduce = s._SensiBusLineReduce;



	_LineImpedance = s._LineImpedance; // B
	_CoresBusLine = s._CoresBusLine; // C
	_lineLimits = s._lineLimits; // l

	// min
	_lineLimitsChange = s._lineLimitsChange;
	lineMin = s.lineMin;
	lineoffset = s.lineoffset;

	_CoresLineBus = s._CoresLineBus;



	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;

	return *this;
}


void StudyCaseDCGrid::genGridFromFile(std::string path, bool alreadyDefine)
{
	// grid 
	_nBus = 1494;
	_nLine = 2156;
	//std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	std::string fileName = path + "SensiBusLineEurope.txt";
	std::string fileName2 = path + "lineLimitEurope.txt";
	std::string fileName3 = path + "SensiBusLineReduceEurope.txt";
	std::string fileName4 = path + "lineLimitReduceEurope.txt";
	std::string pathGrid = path + "Network.txt";
	std::string pathBus = path + "BusAgent.txt"; // corespondance entre les "idBus" du fichier et celui du code (exemple commence � 0 ou � 1)

	
	_LineImpedance = MatrixCPU(_nLine, _nLine); // B
	_CoresBusLine = MatrixCPU(_nBus, _nLine); // C
	_CoresLineBus = MatrixCPU(_nLine, 2); // from, to
	_lineLimits = MatrixCPU(_nLine, 1);
	_zoneBus = MatrixCPU(_nBus, 1);
	_SensiBusLine = MatrixCPU(_nLine, _nBus);
	MatrixCPU fileCoresBus(_nBus, 1);

	if (alreadyDefine) {
		setBusFromFile(pathBus, &fileCoresBus);
		int idBusMax = fileCoresBus.max2();
		MatrixCPU fileBusAgent(idBusMax + 1, 1, -1); // si reste � -1, le bus n'existe pas
		for (int i = 0; i < _nBus; i++) {
			int bus = fileCoresBus.get(i, 0);
			fileBusAgent.set(bus, 0, i);
		}
		setGridFromFile(pathGrid, &fileBusAgent);


		_nLineConstraint = getNFileline(fileName4);
		_SensiBusLineReduce = MatrixCPU(_nLineConstraint, _nBus);
		_lineLimitsReduce = MatrixCPU(_nLineConstraint, 1);
		_SensiBusLine.setFromFile(fileName);
		_lineLimits.setFromFile(fileName2);
		_SensiBusLineReduce.setFromFile(fileName3);
		_lineLimitsReduce.setFromFile(fileName4);
	}
	else {

		std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
		std::string pathGrid = path + "Network.txt";
		std::string pathBus = path + "BusAgent.txt"; // corespondance entre les "idBus" du fichier et celui du code (exemple commence � 0 ou � 1)


		setBusFromFile(pathBus, &fileCoresBus);

		int idBusMax = fileCoresBus.max2();
		MatrixCPU fileBusAgent(idBusMax + 1, 1, -1); // si reste � -1, le bus n'existe pas


		for (int i = 0; i < _nBus; i++) {
			int bus = fileCoresBus.get(i, 0);
			fileBusAgent.set(bus, 0, i);
		}


		setGridFromFile(pathGrid, &fileBusAgent);



		CalcGridSensi();
		ReduceSensi();
		_SensiBusLine.saveCSV(fileName, mode);
		_lineLimits.saveCSV(fileName2, mode);
		_SensiBusLineReduce.saveCSV(fileName3, mode);
		_lineLimitsReduce.saveCSV(fileName4, mode);
	}
}
void StudyCaseDCGrid::genLineLimit(int nLine, float limit, float dlLimit)
{
	//std::cout << "genLineLimit" << std::endl;
	if (nLine > _nLine) {
		throw std::invalid_argument("nLine is too big");
	}
	if (_nLine < 0) {
		throw std::invalid_argument("nLine must be positive");
	}
	if (_indiceLineNonConstraint.size() + _indiceLineConstraint.size() != _nLine) {
		_indiceLineNonConstraint.clear();
		_indiceLineConstraint.clear();
		//std::cout << _indiceLineNonConstraint.size() << " " << _indiceLineConstraint.size() << " " << _nLine << std::endl;
		//_lineLimits.display();
		for (int i = 0; i < _nLine; i++) {
			float lim = _lineLimits.get(i, 0);
			if (lim != 0 && lim != LINELIMITMAX) {
				_indiceLineConstraint.push_back(i);
			}
			else {
				_indiceLineNonConstraint.push_back(i);
			}

		}
	}
	
	if (nLine > _nLineConstraint) // doit augmenter le nombre de ligne contrainte
	{
		int dl = nLine - _nLineConstraint;
		for (int i = 0; i < dl; i++) {
			int indice = rand() % (_indiceLineNonConstraint.size());
			int j = _indiceLineNonConstraint[indice];
			_indiceLineNonConstraint.erase(_indiceLineNonConstraint.begin() + indice);
			_indiceLineConstraint.push_back(j);
			float l = limit + 2 * dlLimit * (rand1() - 0.5);
			_lineLimits.set(j, 0, l);
		}

	}
	else { // doit diminuer le nombre de ligne contrainte
		if (nLine == 0) {
			int dl = _nLineConstraint;
			for (int i = 0; i < dl; i++) {
				int j = _indiceLineConstraint[dl - i - 1];
				_indiceLineConstraint.pop_back();
				_indiceLineNonConstraint.push_back(j);
				_lineLimits.set(j, 0, 0);
			}
		}
		else {
			int dl = _nLineConstraint - nLine;
			for (int i = 0; i < dl; i++) {
				//std::cout << "-";
				int indice = rand() % _indiceLineConstraint.size();
				int j = _indiceLineConstraint[indice];

				_indiceLineConstraint.erase(_indiceLineConstraint.begin() + indice);
				_indiceLineNonConstraint.push_back(j);
				_lineLimits.set(j, 0, 0);
			}
		}

	}
	//std::cout << std::endl;
	_nLineConstraint = _indiceLineConstraint.size();
	ReduceSensi();

}



void StudyCaseDCGrid::Set39Bus(std::string path, bool alreadyDefine)
{
	clock_t t = clock();

	// grid 
	_nBus = 39;
	_nLine = 46;
	_Sbase = 1; //MW
	float zoneBus[39] = { 0, 0, 0, 2, 2, 2, 2, 2, 2, 2,
						  2, 2, 2, 2, 3, 3, 3, 0, 3, 3,
						  3, 3, 3, 3, 0, 1, 1, 1, 1, 0,
						  2, 2, 3, 3, 3, 3, 0, 1, 0 };
	
	std::string filename = path + "Network39.txt";
	_LineImpedance = MatrixCPU(_nLine, _nLine); // B
	_CoresBusLine = MatrixCPU(_nBus, _nLine); // C
	_CoresLineBus = MatrixCPU(_nLine, 2); // from, to
	_lineLimits = MatrixCPU(_nLine, 1);
	_SensiBusLine = MatrixCPU(_nLine, _nBus);
	
	_zoneBus = MatrixCPU(_nBus, 1);
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;

	MatrixCPU fileCoresBus(_nBus + 1, 1);
	for (int i = 0; i < _nBus + 1; i++) {
		fileCoresBus.set(i, 0, i - 1);
	}
	setGridFromFile(filename, &fileCoresBus);
	for (int id = 0; id < _nBus; id++) {
			_zoneBus.set(id, 0, zoneBus[id]);
	}

	std::string fileName = path + "SensiBusLine39node.txt";
	std::string fileName3 = path + "SensiBusLineReduce39node.txt";
	std::string fileName4 = path + "lineLimitReduce39node.txt";

	if (alreadyDefine) {
		_SensiBusLineReduce = MatrixCPU(_nLineConstraint, _nBus);
		_lineLimitsReduce = MatrixCPU(_nLineConstraint, 1);
		_SensiBusLine.setFromFile(fileName);
		_SensiBusLineReduce.setFromFile(fileName3);
		_lineLimitsReduce.setFromFile(fileName4);
	}
	else {

		CalcGridSensi();
		//_SensiPower.multiply(&_SensiBusLine, &_CoresBusAgent);
		ReduceSensi();
		_SensiBusLine.saveCSV(fileName, mode, 0, " ");
		_SensiBusLineReduce.saveCSV(fileName3, mode, 0, " ");
		_lineLimitsReduce.saveCSV(fileName4, mode, 0, " ");
	}
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}
void StudyCaseDCGrid::Set3Bus(std::string path) {
	clock_t t = clock();
	// grid 
	_nBus = 3;
	_nLine = 3;
	std::string filename = path + "Network3.txt";
	
	_LineImpedance = MatrixCPU(_nLine, _nLine); // B
	_CoresBusLine = MatrixCPU(_nBus, _nLine); // C
	_CoresLineBus = MatrixCPU(_nLine, 2); // from, to
	_lineLimits = MatrixCPU(_nLine, 1);
	_SensiBusLine = MatrixCPU(_nLine, _nBus); // A
	//std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	
	MatrixCPU fileCoresBus(_nBus, 1);
	for (int i = 0; i < _nBus; i++) {
		fileCoresBus.set(i, 0, i);
	}
	setGridFromFile(filename, &fileCoresBus);
	CalcGridSensi();

	/*for (int id = 0; id < _nAgent; id++) {
		_CoresBusAgent.set(fileCoresBus.get(BusAgent[id], 0), id, 1);
	}
	
	_SensiPower.multiply(&_SensiBusLine, &_CoresBusAgent);*/

	ReduceSensi();
	
	
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}
void StudyCaseDCGrid::Set4nodeBis(std::string path)
{
	// cas d'�tude pour simuler le cas d'EVA pendant son stage
	clock_t t = clock();
	
	// grid 
	_nBus = 4;
	_nLine = 4;
	std::string filename = path + "Network4.txt";
	
	_LineImpedance = MatrixCPU(_nLine, _nLine); // B
	_CoresBusLine = MatrixCPU(_nBus, _nLine); // C
	_CoresLineBus = MatrixCPU(_nLine, 2); // from, to
	_lineLimits = MatrixCPU(_nLine, 1);
	//_CoresBusAgent = MatrixCPU(_nBus, _nAgent);
	_SensiBusLine = MatrixCPU(_nLine, _nBus); // A
	//std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;



	MatrixCPU fileCoresBus(_nBus, 1);
	for (int i = 0; i < _nBus; i++) {
		fileCoresBus.set(i, 0, i);
	}
	setGridFromFile(filename, &fileCoresBus);


	/*for (int id = 0; id < _nAgent; id++) {
		_CoresBusAgent.set(fileCoresBus.get(BusAgent[id], 0), id, 1);
	}
	_SensiPower.multiply(&_SensiBusLine, &_CoresBusAgent);*/
	CalcGridSensi();
	
	ReduceSensi();

	

	
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;


}
void StudyCaseDCGrid::Set2nodeConstraint(float lim)
{
	clock_t t = clock();
	
	// grid 
	_nBus = 2;
	_nLine = 1;
	

	_LineImpedance = MatrixCPU(_nLine, _nLine); // B
	_CoresBusLine = MatrixCPU(_nBus, _nLine); // C
	_lineLimits = MatrixCPU(_nLine, 1);
	_SensiBusLine = MatrixCPU(_nLine, _nBus); // A
	_zoneBus = MatrixCPU(_nBus, 1);
	_CoresLineBus = MatrixCPU(_nLine, 2); // from, to
	

	/*_CoresBusAgent.set(0, 0, 1);
	_CoresBusAgent.set(1, 1, 1);*/
	
	
	_LineImpedance.set(0, 0, 1); //bii = 1; se simplifie
	

	_SensiBusLine.set(0, 1, -1);
	_lineLimits.set(0, 0, lim);

	
	_nLineConstraint = _nLine;
	_SensiBusLineReduce = MatrixCPU(_SensiBusLine);
	_lineLimitsReduce = MatrixCPU(_lineLimits);
	
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}

void StudyCaseDCGrid::SetEuropeP0(const std::string& path, bool alreadyDefine)
{
	clock_t t = clock();
	_nBus = 1494;
	_nLine = 2156;
	
	_LineImpedance = MatrixCPU(_nLine, _nLine); // B
	_CoresBusLine = MatrixCPU(_nBus, _nLine); // C
	_CoresLineBus = MatrixCPU(_nLine, 2); // from, to
	_lineLimits = MatrixCPU(_nLine, 1);
	_SensiBusLine = MatrixCPU(_nLine, _nBus); // A
	_zoneBus = MatrixCPU(_nBus, 1);
	
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	std::string fileName = path  + "SensiBusLineEurope.txt";
	std::string fileName2 = path + "lineLimitEurope.txt";
	std::string fileName3 = path + "SensiBusLineReduceEurope.txt";
	std::string fileName4 = path + "lineLimitReduceEurope.txt";

	if (alreadyDefine) {
		_nLineConstraint = getNFileline(fileName4);
		_SensiBusLineReduce = MatrixCPU(_nLineConstraint, _nBus);
		_lineLimitsReduce = MatrixCPU(_nLineConstraint, 1);
		_SensiBusLine.setFromFile(fileName);
		_lineLimits.setFromFile(fileName2);
		_SensiBusLineReduce.setFromFile(fileName3);
		_lineLimitsReduce.setFromFile(fileName4);
	}
	else {

		std::string pathGrid = path + "Network.txt";
		fileCoresBus = MatrixCPU(_nBus, 1);
		std::string pathBus = path + "BusAgent.txt"; // corespondance entre les "idBus" du fichier et celui du code (exemple commence � 0 ou � 1)
		setBusFromFile(pathBus, &fileCoresBus);
		int idBusMax = fileCoresBus.max2();
		MatrixCPU fileBusAgent(idBusMax + 1, 1, -1); // si reste � -1, le bus n'existe pas
		//std::cout << idBusMax << std::endl;
		for (int i = 0; i < _nBus; i++) {
			int bus = fileCoresBus.get(i, 0);
			fileBusAgent.set(bus, 0, i);
		}
		setGridFromFile(pathGrid, &fileBusAgent);
		CalcGridSensi();
		ReduceSensi();
		

		_SensiBusLine.saveCSV(fileName, mode);
		_lineLimits.saveCSV(fileName2, mode);
		_SensiBusLineReduce.saveCSV(fileName3, mode);
		_lineLimitsReduce.saveCSV(fileName4, mode);
	}
	//_SensiPower.display();
	t = clock() - t;
	//_timeInit = (float)t / CLOCKS_PER_SEC;
}
void StudyCaseDCGrid::SetStudyCaseDCGrid(std::string path, std::string name, int nBus, bool alreadyDefine)
{
	
	clock_t t = clock();

	std::string fileName = path + "SensiBusLine" + name + ".txt";
	std::string fileName2 = path + "lineLimit" + name + ".txt";
	std::string fileName3 = path + "SensiBusLineReduce" + name + ".txt";
	std::string fileName4 = path + "lineLimitReduce" + name + ".txt";
	
	_name = name;
	
	
	// grid 
	_nBus = nBus; //_nBus = _nCons;
	
	std::string pathGrid = path + "Network" + name + ".txt";
	_nLine = getNFileline(pathGrid);
	//std::cout << "nb de ligne " <<  _nLine << std::endl;
	_LineImpedance = MatrixCPU(_nLine, _nLine); // B
	_CoresBusLine = MatrixCPU(_nBus, _nLine); // C
	_lineLimits = MatrixCPU(_nLine, 1);
	_SensiBusLine = MatrixCPU(_nLine, _nBus); // A
	_zoneBus = MatrixCPU(_nBus, 1);
	_CoresLineBus = MatrixCPU(_nLine, 2);
	
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	

	if (alreadyDefine) {
		_nLineConstraint = getNFileline(fileName4);
		_SensiBusLineReduce = MatrixCPU(_nLineConstraint, _nBus);
		_lineLimitsReduce = MatrixCPU(_nLineConstraint, 1);
		_SensiBusLine.setFromFile(fileName);
		_lineLimits.setFromFile(fileName2);
		_SensiBusLineReduce.setFromFile(fileName3);
		_lineLimitsReduce.setFromFile(fileName4);
	}
	else {

		std::string pathGrid = path + "Network" + name + ".txt";
		fileCoresBus = MatrixCPU(_nBus, 1);
		std::string pathBus = path + "BusAgent" + name + ".txt"; // corespondance entre les "idBus" du fichier et celui du code (exemple commence � 0 ou � 1)
		
		setBusFromFile(pathBus, &fileCoresBus);
		int idBusMax = fileCoresBus.max2();
		MatrixCPU fileBusAgent(idBusMax + 1, 1, -1); // si reste � -1, le bus n'existe pas
		//std::cout << idBusMax << std::endl;
		
		for (int i = 0; i < _nBus; i++) {
			int bus = fileCoresBus.get(i, 0);
			fileBusAgent.set(bus, 0, i);
		}
		
		setGridFromFile(pathGrid, &fileBusAgent);
	
		CalcGridSensi();

		ReduceSensi();

		/*std::cout << _nLineConstraint << std::endl;
		for (int i = 0; i < _nAgent; i++) {
			if (i < _nCons) { // le bus correspond directement pour les conso
				int bus = i;
				_CoresBusAgent.set(bus, i, 1);
			}
			else {
				int idGen = i - _nCons;
				int bus = fileBusAgent.get(GenBus.get(idGen, 0), 0);
				_CoresBusAgent.set(bus, i, 1);
			}
		}
		_SensiPower.multiply(&_SensiBusLine, &_CoresBusAgent);
		*/
		
		_SensiBusLine.saveCSV(fileName, mode);
		_lineLimits.saveCSV(fileName2, mode);
		_SensiBusLineReduce.saveCSV(fileName3, mode);
		_lineLimitsReduce.saveCSV(fileName4, mode);
	}
	//_SensiPower.display();
	t = clock() - t;
	//_timeInit = (float)t / CLOCKS_PER_SEC;
}

void StudyCaseDCGrid::setFromInterface(StudyCaseInterface* interface){
	clock_t t = clock();

	_name = interface->getName();
	
	MatrixCPU infoCase = interface->getInfoCase();
	_Sbase = infoCase.get(0, Sbase_ind);

	// grid 
	_nBus = interface->getB(); //_nBus = _nCons;
	_nLine = interface->getL();
	//std::cout << "nb de ligne " <<  _nLine << std::endl;
	_LineImpedance = MatrixCPU(_nLine, _nLine); // B
	_CoresBusLine = MatrixCPU(_nBus, _nLine); // C
	_lineLimits = MatrixCPU(_nLine, 1);
	_SensiBusLine = MatrixCPU(_nLine, _nBus); // A
	_zoneBus = MatrixCPU(_nBus, 1);
	_CoresLineBus = MatrixCPU(_nLine, 2);
	
	//std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	
	MatrixCPU branchCase = interface->getBranchCase();
	MatrixCPU busCase    = interface->getBusCase();

	_nLineConstraint = 0;
	for (int i = 0; i < _nLine; i++) {
		int nodeFrom = branchCase.get(i, From_ind);
		int nodeTo   = branchCase.get(i, To_ind);
		float react = branchCase.get(i, ZsIm_ind); 
		
		float limit = branchCase.get(i, lim_ind) / _Sbase;

		_LineImpedance.set(i, i, -react);
		_CoresBusLine.set(nodeFrom, i, 1);
		_CoresBusLine.set(nodeTo, i, -1);
		_CoresLineBus.set(i, 0, nodeFrom);
		_CoresLineBus.set(i, 1, nodeTo);
		if (limit > 0) {
			_nLineConstraint++;
			_lineLimits.set(i, 0, limit);
			_indiceLineConstraint.push_back(i);
		}
		else {
			_lineLimits.set(i, 0, LINELIMITMAX);
			_indiceLineNonConstraint.push_back(i);
		}
	}

	CalcGridSensi();

	ReduceSensi();

	/*std::cout << _nLineConstraint << std::endl;
	for (int i = 0; i < _nAgent; i++) {
		if (i < _nCons) { // le bus correspond directement pour les conso
			int bus = i;
			_CoresBusAgent.set(bus, i, 1);
		}
		else {
			int idGen = i - _nCons;
			int bus = fileBusAgent.get(GenBus.get(idGen, 0), 0);
			_CoresBusAgent.set(bus, i, 1);
		}
	}
	_SensiPower.multiply(&_SensiBusLine, &_CoresBusAgent);
	*/
		
	
	//_SensiPower.display();
	t = clock() - t;
	//_timeInit = (float)t / CLOCKS_PER_SEC;
	
}




void StudyCaseDCGrid::setLineLimitMin(float min)
{
	if (lineoffset) {
		lineoffset = 0;
	}
	_lineLimitsChange = getLineLimit();
	lineMin = min;
	for (int l = 0; l < getNLine(); l++) {
		if (_lineLimitsChange.get(l, 0) < lineMin) {
			_lineLimitsChange.set(l, 0, lineMin);
		}
	}
}
void StudyCaseDCGrid::setLineLimitRelaxation(float eps)
{
	if (lineMin) {
		lineMin = 0;
	}
	_lineLimitsChange = getLineLimit();
	lineoffset = eps;
	for (int l = 0; l < getNLine(); l++) {
		
		_lineLimitsChange.increment(l, 0, -eps);
		
	}
}
void StudyCaseDCGrid::setLineLimit(int line, float limit)
{
	if (line > _nLine) {
		throw std::invalid_argument("this line doesn't exist");
	}
	else {
		float oldLimit = _lineLimits.get(line, 0);
		if (oldLimit == LINELIMITMAX) {
			_nLineConstraint++;
		}
		_lineLimits.set(line, 0, limit/_Sbase);
		ReduceSensi();
	}
}


MatrixCPU StudyCaseDCGrid::getPowerSensiBus(bool force) const
{
	if (force) {
		return _SensiBusLine;
	}
	if (toReduce) {
		return _SensiBusLineReduce;
	}
	else {
		return _SensiBusLine;
	}
	
}
MatrixCPU StudyCaseDCGrid::getPowerSensiBusReduce() const
{
	return _SensiBusLineReduce;
}
MatrixCPU StudyCaseDCGrid::getLineLimit() const
{
	if (lineMin || lineoffset) {
		return _lineLimitsChange;
	} 
	else if (toReduce) {
		return _lineLimitsReduce;
	}
	else {
		return _lineLimits;
	}
	
}
MatrixCPU StudyCaseDCGrid::getCoresLineBus(bool force) const
{
	if (toReduce && !force) {
		return _CoresLineBusReduce;
	}
	else {
		return _CoresLineBus;
	}
}
MatrixCPU StudyCaseDCGrid::getfileCoresBus() const
{
	return fileCoresBus;
}
MatrixCPU StudyCaseDCGrid::getZones() const
{
	return _zoneBus;
}
float StudyCaseDCGrid::getTimeInit() const
{
	return _timeInit;
}
int StudyCaseDCGrid::getNLine(bool force) const
{
	if (force) {
		return _nLine;
	}
	if (toReduce) {
		return _nLineConstraint;
	}
	else {
		return _nLine;
	}
	
}
int StudyCaseDCGrid::getNLineConstraint() const
{
	return _nLineConstraint;
}
int StudyCaseDCGrid::getNBus() const
{
	return _nBus;
}
std::string StudyCaseDCGrid::getName() const
{
	return _name;
}



void StudyCaseDCGrid::saveCSV(const std::string& fileName)
{

	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	MatrixCPU nombre(1, 3);
	nombre.set(0, 0, _nLine);
	nombre.set(0, 1, _nLineConstraint);
	nombre.set(0, 2, _nBus);
	nombre.saveCSV(fileName, mode);


	MatrixCPU temp2(1, _nLineConstraint);
	temp2.addTrans(&_lineLimitsReduce);
	temp2.saveCSV(fileName, mode);
	_SensiBusLineReduce.saveCSV(fileName, mode);

}

void StudyCaseDCGrid::display(int type) const
{
	std::cout << "Study Case : " << _nBus << " bus and " << _nLine << " lines " << std::endl;
	std::cout << "and " << _nLineConstraint << " and reduced " << toReduce << std::endl;

	if (_nLine < 100 && _nBus < 100) {
		std::cout << " B :" << std::endl;
		_LineImpedance.display();
		std::cout << " C :" << std::endl;
		_CoresBusLine.display();
		if (toReduce) {
			std::cout << " Line limit :" << std::endl;
			_lineLimitsReduce.display();
			std::cout << " Sensibility :" << std::endl;
			_SensiBusLineReduce.display();
		}
		else {
			std::cout << " Line limit :" << std::endl;
			_lineLimits.display();
			std::cout << " Sensibility :" << std::endl;
			_SensiBusLine.display();
		}
		
	}
	
	
	if(type ==1){

		MatrixCPU temp1(_nLine, _nBus); // BC^T 
		MatrixCPU temp2(_nBus, _nBus); // CBC^T
		MatrixCPU temp33(_nBus - 1, _nBus - 1); //(CBC ^ T) ^ -1 sans la ligne et colonne du noued de ref
		MatrixCPU temp22(_nBus - 1, _nBus - 1); // on enl�ve la ligne du noeud de ref�rence

		MatrixCPU LineImpedance(_LineImpedance);
		MatrixCPU CoresBusLine(_CoresBusLine);

		temp1.multiplyTrans(&LineImpedance, &CoresBusLine);
		temp2.multiply(&CoresBusLine, &temp1);
		temp2.getBloc(&temp22, 1, _nBus, 1, _nBus);
		temp33.invertGaussJordan(&temp22);
		std::cout << " BC^T :" << std::endl;
		temp1.display();
		std::cout << " CBC^T :" << std::endl;
		temp2.display();
		std::cout << " (CBC^T)^-1 :" << std::endl;
		temp33.display();
	
	}
	
}

void StudyCaseDCGrid::displayLineCores(MatrixCPU* g, bool all)
{
	std::cout << "Line correspendance : " << std::endl;
	
	MatrixCPU Cores(getCoresLineBus());
	//Cores.display();
	
	MatrixCPU Limit(getLineLimit());
	//Limit.display();
	if (all) {
		if (Cores.getNLin() == 0) {
			for (int l = 0; l < getNLine(); l++) {
				if (fabs(g->get(l, 0)) > Limit.get(l, 0)) {
					std::cout << "******* OverFlow : Line n " << l << " line limit " << Limit.get(l, 0)
						<< " flow " << g->get(l, 0) << "********" << std::endl;
				}
				else if (Limit.get(l, 0) - fabs(g->get(l, 0)) < 0.1) {
					std::cout << "+++ Close bounds : Line n " << l << " line limit " << Limit.get(l, 0)
						<< " flow " << g->get(l, 0) << "+++" << std::endl;
				}
				else {
					std::cout << "Line n " << l << " line limit " << Limit.get(l, 0)
						<< " flow " << g->get(l, 0) << std::endl;
				}
				
			}
		}
		else {
			for (int l = 0; l < getNLine(); l++) {
				if (fabs(g->get(l, 0)) > Limit.get(l, 0)) {
					std::cout << "********* OverFlow :Line n " << l << " from node " << Cores.get(l, 0)
						<< " to node " << Cores.get(l, 1) << " line limit " << Limit.get(l, 0)
						<< " flow " << g->get(l, 0) << "********" << std::endl;
				}
				else if (Limit.get(l, 0) - fabs(g->get(l, 0)) < 0.1) {
					std::cout << "+++ Close bounds :Line n " << l << " from node " << Cores.get(l, 0)
						<< " to node " << Cores.get(l, 1) << " line limit " << Limit.get(l, 0)
						<< " flow " << g->get(l, 0) << std::endl;
				}
				else {
					std::cout << "Line n " << l << " from node " << Cores.get(l, 0)
						<< " to node " << Cores.get(l, 1) << " line limit " << Limit.get(l, 0)
						<< " flow " << g->get(l, 0) << std::endl;
				}
				
			}
		}
	}
	else {
		if (Cores.getNLin() == 0) {
			for (int l = 0; l < getNLine(); l++) {
				if (fabs(g->get(l, 0)) > Limit.get(l, 0)) {
					std::cout << "********* OverFlow : Line n " << l << " line limit " << Limit.get(l, 0)
						<< " flow " << g->get(l, 0) << "*********"<< std::endl;
				}
				else if (Limit.get(l, 0) - fabs(g->get(l, 0)) < 0.1) {
					std::cout << "+++ Close bounds  : Line n " << l << " line limit " << Limit.get(l, 0)
						<< " flow " << g->get(l, 0) << "+++" << std::endl;
				}
				
			}
		}
		else {
			for (int l = 0; l < getNLine(); l++) {
				if (fabs(g->get(l, 0)) > Limit.get(l, 0)) {
					std::cout << "********* OverFlow :Line n " << l << " from node " << Cores.get(l, 0)
						<< " to node " << Cores.get(l, 1) << " line limit " << Limit.get(l, 0)
						<< " flow " << g->get(l, 0) << "*********" << std::endl;
				}
				else if (Limit.get(l, 0) - fabs(g->get(l, 0)) < 0.1) {
					std::cout << "+++ Close bounds :Line n " << l << " from node " << Cores.get(l, 0)
						<< " to node " << Cores.get(l, 1) << " line limit " << Limit.get(l, 0)
						<< " flow " << g->get(l, 0) << std::endl;
				}
				
			}
		}
	}
	
	
}

StudyCaseDCGrid::~StudyCaseDCGrid()
{
#ifdef DEBUG_DESTRUCTOR
	std::cout << "case destructor" << std::endl;
#endif

}


