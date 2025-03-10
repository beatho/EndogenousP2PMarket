#include "../head/StudyCaseACGrid.h"
float  BTLINE[] = { 0.306, 0.29, 13.2, 0.35 };// 94-AL1/15-ST1A 0.4 : r x b Imax (ohm, ohm, nF) /km kA
float HTALINE[] = { 0.5939, 0.372, 9.5, 0.21 };// 48-AL1/8-ST1A 20 : r x b Imax (ohm, ohm, nF) /km  kA
float HTBLINE[] = { 0.059, 0.253, 11, 0.96 };// 490-AL1/64-ST1A 380.0 : r x b Imax (ohm, ohm, nF) /km kA


void StudyCaseACGrid::genGridLine(int nBus, float length, float dlength)
{
	
	_Sbase = 1; // 1MW
	_Vbase = 0.4; // 400V
	_Zbase = _Vbase * _Vbase / _Sbase;
	int billion = 1000000000;

	_nBus = nBus;
	_nLine = nBus - 1;
	_nConstraint = _nLine + 2 * _nBus;
	_nLineConstraint = 0;
	_V0 = 1;
	_theta0 = 0;

	initMat();


	for (int l = 0; l < _nLine; l++) {
		int j = l + 1;
		int i = l; // tous en ligne
		creatLine(l, i, j, length, dlength);
	}
	hasCurrentLimit = true;



	setDefaultConstraint();
	LinearizeImp();
}

void StudyCaseACGrid::genGridOneStep(int nBus, float length, float dlength)
{

	_Sbase = 1; // 1MW
	_Vbase = 0.4; // 400V
	_Zbase = _Vbase * _Vbase / _Sbase;
	int billion = 1000000000;

	_nBus = nBus;
	_nLine = nBus - 1;
	_nConstraint = _nLine + 2 * _nBus;
	_nLineConstraint = 0;
	_V0 = 1;
	_theta0 = 0;

	initMat();


	for (int l = 0; l < _nLine; l++) {
		int j = l + 1;
		int i = 0; // tous lie a 0
		creatLine(l, i, j, length, dlength);
	}
	hasCurrentLimit = true;


	
	setDefaultConstraint();
	LinearizeImp();
}


void StudyCaseACGrid::genGridBalance(int nBus, float length, float dlength)
{

	_Sbase = 1; // 1MW
	_Vbase = 0.4; // 400V
	_Zbase = _Vbase * _Vbase / _Sbase;
	int billion = 1000000000;

	_nBus = nBus;
	_nLine = nBus - 1;
	_nConstraint = _nLine + 2 * _nBus;
	_nLineConstraint = 0;
	_V0 = 1;
	_theta0 = 0;

	initMat();

	int Nstep = getBalanceNChild(nBus);
	//std::cout << "nombre de bus " << nBus << " nombre d'�tage et d'enfant par bus " << Nstep << std::endl;
	int line = 0;
	int from = 0;
	int nBusOnStep = 1;
	int i = 0;
	hasCurrentLimit = true;
	for (int step = 1; step < Nstep; step++) {
		//std::cout << "�tage " << step << " nombre de bus de letage precedent " << nBusOnStep << std::endl;
		int futurNBusOnStep = 0;
		for (int bus = 0; bus < nBusOnStep; bus++) {
			i = from + bus;
			for (int child = 0; child < Nstep; child++) {
				if (line < _nLine) {
					int j = line + 1;
					//std::cout << "create line " << line << "between bus " << i << " " << j << std::endl;
					creatLine(line, i, j, length, dlength);
					line++;
					futurNBusOnStep++;
				}
				else {
					
					setDefaultConstraint();
					LinearizeImp();
					return;

				}
			}
		}
		nBusOnStep = futurNBusOnStep;
		from = i + 1;
	}

	
}


int StudyCaseACGrid::getBalanceNChild(int nBus)
{
	for (int i = 0; i < nBus; i++) {
		int accumulateur = 0;
		int puissance = 1;
		for (int j = 0; j < i; j++) {
			accumulateur += puissance;
			puissance *= i;
		}
		if (accumulateur >= nBus) {
			return i;
		}

	}

	return 0;
}

void StudyCaseACGrid::creatLine(int line, int from, int to, float length, float dlength)
{
	int l = line;
	int i = from;
	int j = to;
	int billion = 1000000000;


	_CoresLineBus.set(l, 0, i);
	_CoresLineBus.set(l, 1, j);

	double L = length + 2 * (rand1() - 0.5) * dlength;
	double ZsRe = L * BTLINE[0] / _Zbase;
	double ZsIm = L * BTLINE[1] / _Zbase;
	double limitLine = 2 * BTLINE[3] * _Vbase / _Sbase; // limite puissance active
	double limitI = limitLine * limitLine; // limite courrant
	double YlsRe = ZsRe / (ZsRe * ZsRe + ZsIm * ZsIm);
	double YlsIm = -ZsIm / (ZsRe * ZsRe + ZsIm * ZsIm);
	double Ylp = 314.15 * L * BTLINE[2] / (2 * billion * _Zbase); // b/2


	//std::cout << " Ligne numero " << l << " entre bus " << i << " et " << j << " re " << YlsRe << " im " << YlsIm << " b " << Ylp << " limit "<< limitLine << std::endl;



	_lineReactance.increment(i, j, -YlsRe);
	_lineReactance.increment(j, i, -YlsRe);
	_lineSuceptance.increment(j, i, -YlsIm);
	_lineSuceptance.increment(i, j, -YlsIm);

	_lineReactance.increment(j, j, YlsRe);
	_lineSuceptance.increment(j, j, YlsIm + Ylp);
	_lineReactance.increment(i, i, YlsRe);
	_lineSuceptance.increment(i, i, YlsIm + Ylp);


	_lineReactanceD.increment(i, j, -YlsRe);
	_lineReactanceD.increment(j, i, -YlsRe);
	_lineSuceptanceD.increment(j, i, -YlsIm);
	_lineSuceptanceD.increment(i, j, -YlsIm);

	_lineReactanceD.increment(j, j, YlsRe);
	_lineSuceptanceD.increment(j, j, YlsIm + Ylp);
	_lineReactanceD.increment(i, i, YlsRe);
	_lineSuceptanceD.increment(i, i, YlsIm + Ylp);


	_busSuceptance.increment(i, 0, Ylp);
	_busSuceptance.increment(j, 0, Ylp);


	_lineImpedanceReal.increment(l, 0, ZsRe);
	_lineImpedanceImag.increment(l, 0, ZsIm);
	if (limitLine > 0) {
		_nLineConstraint++;
		_lineLimits.set(l, 0, limitLine); // dans le cas DC : S = P = UI avec U = 1 pu, donc P = I
		_currentLimit.set(l, 0, limitI);
	}
	else {
		_lineLimits.set(l, 0, LINELIMITMAX);
	}
}

void StudyCaseACGrid::setDefaultConstraint()
{
	for (int i = 0; i < _nBus; i++) { // bound on voltage angle rad
		_upperBound.set(i, 0, 3);
		_lowerBound.set(i, 0, -3);
		_upperBound.set(i + _nBus, 0, 1.1 * _V0);
		_lowerBound.set(i + _nBus, 0, 0.9 * _V0);
		_VoltageInit.set(i, 0, 0);
		_VoltageInit.set(i + _nBus, 0, _V0);
		_VoltageInitD.set(i, 0, 0);
		_VoltageInitD.set(i + _nBus, 0, _V0);

	}


	int indice = 0;
	for (int i = 2 * _nBus; i < _nConstraint; i++) { // bound on power flow
		_upperBound.set(i, 0, _lineLimits.get(indice, 0));
		_lowerBound.set(i, 0, -_lineLimits.get(indice, 0));
		indice++;
	}
}

void StudyCaseACGrid::setGridACFromFile(const std::string& path, MatrixCPU* fileBusAgent)
{
	MatrixCPU matFile(_nLine, 6);
	matFile.setFromFile(path);
	
	_nLineConstraint = 0;
	for (int i = 0; i < _nLine; i++) {
		int nodeFromFile = matFile.get(i, 0);
		int nodeToFile = matFile.get(i, 1);
		float YlsRe = matFile.get(i, 2);
		float YlsIm = matFile.get(i, 3);
		float Ylp = matFile.get(i, 4);
		float limit = matFile.get(i, 5);

		int nodeFrom = fileBusAgent->get(nodeFromFile, 0);
		int nodeTo = fileBusAgent->get(nodeToFile, 0);

		
		_lineReactance.set(nodeFrom, nodeTo, -YlsRe);
		_lineReactance.set(nodeTo, nodeFrom, -YlsRe);
		_lineSuceptance.set(nodeFrom, nodeTo, -YlsIm);
		_lineSuceptance.set(nodeTo, nodeFrom, -YlsIm);
		
		_lineReactance.increment(nodeTo, nodeTo, YlsRe);
		_lineSuceptance.increment(nodeTo, nodeTo, YlsIm + Ylp);
		_lineReactance.increment(nodeFrom, nodeFrom, YlsRe);
		_lineSuceptance.increment(nodeFrom, nodeFrom, YlsIm + Ylp);
		
		_CoresLineBus.set(i, 0, nodeFrom);
		_CoresLineBus.set(i, 1, nodeTo);
		
		
		if (limit > 0) {
			_nLineConstraint++;
			_lineLimits.set(i, 0, limit);
		}
		else {
			_lineLimits.set(i, 0, LINELIMITMAX);
		}
		
	}

	std::cout << "nLine Contraint" << _nLineConstraint << std::endl;

}
void StudyCaseACGrid::setBusFromFile(const std::string& path, MatrixCPU* fileCoresBus)
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


void StudyCaseACGrid::initMat()
{
	_lineLimits = MatrixCPU(_nLine, 1);
	//_CoresBusAgent = MatrixCPU(_nBus, _nAgent);

	_lineReactance = MatrixCPU(_nBus, _nBus);
	_lineSuceptance = MatrixCPU(_nBus, _nBus);
	_lineReactanceD = MatrixCPUD(_nBus, _nBus);
	_lineSuceptanceD = MatrixCPUD(_nBus, _nBus);

	
	_CoresLineBus = MatrixCPU(_nLine, 2); // from, to
	_upperBound = MatrixCPU(_nConstraint, 1);
	_lowerBound = MatrixCPU(_nConstraint, 1);
	_currentLimit = MatrixCPU(_nLine, 1, 10000);

	_zoneBus = MatrixCPU(_nBus, 1);

	_busSuceptance = MatrixCPU(_nBus, 1);

	_lineImpedanceReal = MatrixCPU(_nLine, 1);
	_lineImpedanceImag = MatrixCPU(_nLine, 1);

	_VoltageInit = MatrixCPU(2 * _nBus, 1);
	_VoltageInitD = MatrixCPUD(2 * _nBus, 1);

	for (int i = 0; i < _nBus; i++) {
		_VoltageInit.set(i + _nBus, 0, _V0);
		_VoltageInitD.set(i + _nBus, 0, _V0);
	}

	_lineSuceptanceLin = MatrixCPU(_nBus + 2 * _nLine, 1);
	_lineReactanceLin = MatrixCPU(_nBus + 2 * _nLine, 1);
	_lineSuceptanceLinD = MatrixCPUD(_nBus + 2 * _nLine, 1);
	_lineReactanceLinD = MatrixCPUD(_nBus + 2 * _nLine, 1);
	_CoresVoiLin = MatrixCPU(_nBus + 2 * _nLine, 1);
	_CoresBusLin = MatrixCPU(_nBus, 1);
	_nLines = MatrixCPU(_nBus, 1);
	_nLinesBegin = MatrixCPU(_nBus, 1);

	_lineSuceptanceLin2 = MatrixCPU(_nLine, 1);
	_lineReactanceLin2 = MatrixCPU(_nLine, 1);

}

void StudyCaseACGrid::copy(StudyCaseACGrid* s)
{
	clock_t t = clock();

	_nLine = s->_nLine;
	_nBus = s->_nBus;
	_nLineConstraint = s->_nLineConstraint;
	_lineLimits = s->_lineLimits;
	_CoresLineBus = s->_CoresLineBus;

	_name = s->_name;
	_currentLimit = s->_currentLimit;
	hasCurrentLimit = s->hasCurrentLimit;
	toReduce = s->toReduce;

	// AC
	_Zbase = s->_Zbase;
	_Sbase = s->_Sbase;
	_Vbase = s->_Vbase;

	_lineSuceptance = s->_lineSuceptance;
	_lineReactance = s->_lineReactance;
	_lineSuceptanceD = s->_lineSuceptanceD;
	_lineReactanceD = s->_lineReactanceD;

	_lineSuceptanceLin = s->_lineSuceptanceLin;
	_lineReactanceLin = s->_lineReactanceLin;
	_lineSuceptanceLinD = s->_lineSuceptanceLinD;
	_lineReactanceLinD = s->_lineReactanceLinD;

	_upperBound = s->_upperBound;
	_lowerBound = s->_lowerBound;

	_CoresVoiLin = s->_CoresVoiLin;
	_CoresBusLin = s->_CoresBusLin;
	_nLines = s->_nLines;



	//distribution network
	_busSuceptance = s->_busSuceptance;
	_lineImpedanceImag = s->_lineImpedanceImag;
	_lineImpedanceReal = s->_lineImpedanceReal;


	//Sol
	_SolutionPF = s->_SolutionPF;
	_VoltageInit = s->_VoltageInit;
	_VoltageInitD = s->_VoltageInitD;

	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;

	
}

void StudyCaseACGrid::LinearizeImp()
{
	_lineSuceptanceLin = MatrixCPU(_nBus + 2 * _nLine, 1);
	_lineReactanceLin = MatrixCPU(_nBus + 2 * _nLine, 1);
	_lineSuceptanceLinD = MatrixCPUD(_nBus + 2 * _nLine, 1);
	_lineReactanceLinD = MatrixCPUD(_nBus + 2 * _nLine, 1);
	_CoresVoiLin = MatrixCPU(_nBus + 2 * _nLine, 1);
	_CoresBusLin = MatrixCPU(_nBus, 1);
	_nLines = MatrixCPU(_nBus, 1);
	_nLinesBegin = MatrixCPU(_nBus, 1);

	_lineSuceptanceLin2 = MatrixCPU(_nLine, 1);
	_lineReactanceLin2 = MatrixCPU(_nLine, 1);

	int indice = 0;
	int line = 0;
	for (int i = 0; i < _nBus; i++) {
		_CoresBusLin.set(i, 0, indice);
		_lineSuceptanceLin.set(indice, 0, _lineSuceptance.get(i, i));
		_lineReactanceLin.set(indice, 0, _lineReactance.get(i, i));
		_lineSuceptanceLinD.set(indice, 0, _lineSuceptanceD.get(i, i));
		_lineReactanceLinD.set(indice, 0, _lineReactanceD.get(i, i));
		_CoresVoiLin.set(indice, 0, i);
		_nLines.increment(i, 0, 1);
		indice++;
		for (int j = 0; j < _nBus; j++) {
			if (i != j) {
				if (abs(_lineReactance.get(i, j)) > 0 || abs(_lineSuceptance.get(i, j)) > 0) {
					_nLines.increment(i, 0, 1);
					_lineSuceptanceLin.set(indice, 0, _lineSuceptance.get(i, j));
					_lineReactanceLin.set(indice, 0, _lineReactance.get(i, j));
					_lineSuceptanceLinD.set(indice, 0, _lineSuceptanceD.get(i, j));
					_lineReactanceLinD.set(indice, 0, _lineReactanceD.get(i, j));
					_CoresVoiLin.set(indice, 0, j);
					indice++;
					if (i < j) {
						_lineSuceptanceLin2.set(line, 0, _lineSuceptance.get(i, j));
						_lineReactanceLin2.set(line, 0, _lineReactance.get(i, j));
						line++;
					}
				}
			}
		}
		_nLinesBegin.set(i, 0, indice - (i + 1));
	}
}

void StudyCaseACGrid::genDCGridFromAC()
{
	
	_LineImpedance = MatrixCPU(_nLine, _nLine); // B
	_CoresBusLine = MatrixCPU(_nBus, _nLine); // C
	_SensiBusLine = MatrixCPU(_nLine, _nBus);
	_nLineConstraint = 0;
	for (int i = 0; i < _nLine; i++) {
		float limit = _lineLimits.get(i, 0);
		
		int nodeFrom = _CoresLineBus.get(i, 0);
		int nodeTo = _CoresLineBus.get(i, 1);

		float react = _lineSuceptance.get(nodeFrom, nodeTo);
		_LineImpedance.set(i, i, react);
		_CoresBusLine.set(nodeFrom, i, 1);
		_CoresBusLine.set(nodeTo, i, -1);
		
		if (limit > 0 && limit != LINELIMITMAX) {
			_nLineConstraint++;
			_indiceLineConstraint.push_back(i);
		}
		else {
			_indiceLineNonConstraint.push_back(i);
		}
	}
	
	//std::cout << "gridSensi " << std::endl;
	CalcGridSensi();
	
	ReduceSensi();

}


StudyCaseACGrid::StudyCaseACGrid() : StudyCaseDCGrid()
{
	_Sbase = 1;
	_Vbase = 1;
	_Zbase = _Vbase * _Vbase / _Sbase;

	// AC et DC
	_nBus = 2;
	_nLine = 1;
	_nLineConstraint = 1;
	_zoneBus = MatrixCPU(_nBus,1); // taille B*1 indique pour chaque agent la zone o� il est 
		
	_CoresLineBus = MatrixCPU(_nLine, 2); // Perso from, to
	_CoresLineBus.set(0, 1, 1);
	
	_lineLimits = MatrixCPU(_nLine, 1, 0.8); // l
	_currentLimit = MatrixCPU(_nLine, 1, 200);
	_nConstraint = 2*_nBus + _nLine;

	hasCurrentLimit = false;
	radial = true;

	_busSuceptance = MatrixCPU(_nBus, 1, 0); // Yd for distribution network
	_SolutionPF = MatrixCPUD(_nBus, 4);
	
	_SolutionPF.set(0, 0, 1);
	_SolutionPF.set(1, 0, 1);
	_SolutionPF.set(0, 2, 1);
	_SolutionPF.set(1, 2, -1);
	_SolutionPF.set(0, 3, 0.9);
	_SolutionPF.set(1, 3, -0.9);
	
	_VoltageInitD = MatrixCPUD(2 * _nBus, 1);
	_VoltageInitD.set(2, 0, 1);
	_VoltageInitD.set(3, 0, 1);
	_VoltageInit = _VoltageInitD;
	
	float x =  0.01;
	float r =  0.005;

	float YlsRe =   r / (x * x + r * r);
	float YlsIm = - x / (x * x + r * r);
	float Ylp = 0;
    
	_lineImpedanceReal = MatrixCPU(_nLine, 1, r); // real(zs) for distribution network
	_lineImpedanceImag = MatrixCPU(_nLine, 1, x); // imag(zs) for distribution network
	_lineReactance = MatrixCPU(_nBus, _nBus);
	_lineSuceptance = MatrixCPU(_nBus, _nBus);


	_lineReactance.increment(0, 1, -YlsRe);
	_lineReactance.increment(1, 0, -YlsRe);
	_lineSuceptance.increment(1, 0, -YlsIm);
	_lineSuceptance.increment(0, 1, -YlsIm);

	_lineReactance.increment(1, 1, YlsRe);
	_lineSuceptance.increment(1, 1, YlsIm + Ylp);
	_lineReactance.increment(0, 0, YlsRe);
	_lineSuceptance.increment(0, 0, YlsIm + Ylp);

	_lineReactance.toMatCPUD(_lineReactanceD);
	_lineSuceptance.toMatCPUD(_lineSuceptanceD);
	LinearizeImp();
	_upperBound = MatrixCPU(_nConstraint, 1); // overline(Y) : angle, voltage, powerFlow
	_lowerBound = MatrixCPU(_nConstraint, 1);; // underline(Y)  : angle, voltage, powerFlow

	
	for (int i = 1; i < _nBus; i++) { // bound on voltage angle rad
		_upperBound.set(i, 0, 3);
		_lowerBound.set(i, 0, -3);
	}
	_upperBound.set(_nBus, 0, 1);
	_lowerBound.set(_nBus, 0, 1);
	for (int i = _nBus + 1; i < 2 * _nBus; i++) { // bound on voltage 
		_upperBound.set(i, 0, 1.1);
		_lowerBound.set(i, 0, 0.9);

	}
	for (int i = 0; i < _nLine; i++) {
		_upperBound.set(2 * _nBus + i, 0, _lineLimits.get(i, 0));
		_lowerBound.set(2 * _nBus + i, 0, -_lineLimits.get(i, 0));
	}
	 
}

void StudyCaseACGrid::SetAC39Bus(std::string path, bool alreadyDefine)
{
	clock_t t = clock();
	
	// grid 
	_nBus = 39;
	_nLine = 46;
	_nConstraint = _nLine + 2 * _nBus;
	float zoneBus[39] = { 0, 0, 0, 2, 2, 2, 2, 2, 2, 2,
						  2, 2, 2, 2, 3, 3, 3, 0, 3, 3,
						  3, 3, 3, 3, 0, 1, 1, 1, 1, 0,
						  2, 2, 3, 3, 3, 3, 0, 1, 0 };
	std::string filename = path + "NetworkAC39.txt";

	initMat();

	
	_zoneBus = MatrixCPU(_nBus, 1);
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	std::string fileName = path + "lineReactance39.txt";
	std::string fileName2 = path + "lineSuceptance39.txt";
	std::string fileName3 = path + "upperBound39node.txt";
	std::string fileName4 = path + "lowerBound39node.txt";

	if (alreadyDefine) {
		
		_lineReactance.setFromFile(fileName);
		_lineSuceptance.setFromFile(fileName2);
		_upperBound.setFromFile(fileName3);
		_lowerBound.setFromFile(fileName4);
	}
	else {
		MatrixCPU fileCoresBus(_nBus + 1, 1);
		for (int i = 0; i < _nBus + 1; i++) {
			fileCoresBus.set(i, 0, i - 1);
		}
		setGridACFromFile(filename, &fileCoresBus);


		/*for (int id = 0; id < _nAgent; id++) {
			_CoresBusAgent.set(fileCoresBus.get(BusAgent[id], 0), id, 1);
			_zoneBus.set(BusAgent[id]-1, 0, zoneBus[BusAgent[id] - 1]);
		}*/

		//_lineLimitsReduce = MatrixCPU(_nLineConstraint, 1);
		//_CoresLineBusReduce = MatrixCPU(_nLineConstraint, 2);
		
		for (int i = 0; i < _nLine; i++) {
			int lim = _lineLimits.get(i, 0);
			if (lim == 0) {
				_lineLimits.set(i, 0, LINELIMITMAX);
			}
			else {
				_lineLimits.set(i, 0, lim / _Sbase);
			}
		}
		for (int i = 0; i < _nBus; i++) { // bound on voltage angle rad
			_upperBound.set(i, 0, 3);
			_lowerBound.set(i, 0, -3);
		}
		for (int i = _nBus; i < 2*_nBus; i++) { // bound on voltage 
			_upperBound.set(i, 0, 1.05);
			_lowerBound.set(i, 0, 0.95);
			
		}
		int indice = 0;
		for (int i = 2 * _nBus; i < _nConstraint; i++) { // bound on power flow
			_upperBound.set(i, 0, _lineLimits.get(indice,0));
			_lowerBound.set(i, 0, -_lineLimits.get(indice, 0));
			indice++;
		}

		_lineReactance.saveCSV(fileName, mode);
		_lineSuceptance.saveCSV(fileName2, mode);
		_upperBound.saveCSV(fileName3, mode);
		_lowerBound.saveCSV(fileName4, mode);
	}
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}

void StudyCaseACGrid::SetAC3Bus(std::string path)
{
	//case3.m
	clock_t t = clock();
	_Vbase = 230; // kV car matlab
	_Sbase = 100; // MW car matlab
	_Zbase = _Vbase * _Vbase / _Sbase; // kV*kV/MW les puissances de 10 s'en vont
	
	radial = false;
	
	// grid 
	_nBus = 3;
	_nLine = 3;
	_nConstraint = _nLine + 2 * _nBus;
	float zoneBus[3] = { 0, 0, 0};
	std::string filename = path + "NetworkAC3.txt";

	
	initMat();
	

	_lineReactanceD = MatrixCPUD(_nBus, _nBus);
	_lineSuceptanceD = MatrixCPUD(_nBus, _nBus);
	
	


	
	
	MatrixCPU fileCoresBus(_nBus + 1, 1);
	for (int i = 0; i < _nBus + 1; i++) {
		fileCoresBus.set(i, 0, i - 1);
	}
	
	setGridACFromFile(filename, &fileCoresBus);

	_lineReactance.toMatCPUD(_lineReactanceD);
	_lineSuceptance.toMatCPUD(_lineSuceptanceD);
	
	LinearizeImp();
	
	
	
	for (int i = 0; i < _nLine; i++) {
		int lim = _lineLimits.get(i, 0);
		if (lim == 0 || lim == LINELIMITMAX) {
			_lineLimits.set(i, 0, LINELIMITMAX);
		}
		else {
			_lineLimits.set(i, 0, lim / _Sbase);
		}
	}
	for (int i = 0; i < _nBus; i++) { // bound on voltage angle rad
		_upperBound.set(i, 0, 3);
		_lowerBound.set(i, 0, -3);
	}
	for (int i = _nBus; i < 2 * _nBus; i++) { // bound on voltage 
		_upperBound.set(i, 0, 1.05);
		_lowerBound.set(i, 0, 0.95);

	}
	int indice = 0;
	for (int i = 2 * _nBus; i < _nConstraint; i++) { // bound on power flow
		_upperBound.set(i, 0, _lineLimits.get(indice, 0));
		_lowerBound.set(i, 0, -_lineLimits.get(indice, 0));
		indice++;
	}
	
		
	
	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}

void StudyCaseACGrid::SetACFromFile(std::string name, std::string path)
{
	std::string fileName1 = path + "Case"	+ name + ".txt";
	std::string fileName3 = path + "Branch" + name + ".txt";
	std::string fileName4 = path + "Bus"	+ name + ".txt";
	std::string fileName5 = path + "Sol"	+ name + ".txt";
	std::string fileName6 = path + "Bgrid"  + name + ".txt";
	std::string fileName7 = path + "Ggrid"  + name + ".txt";

	MatrixCPUD Info(1, 9); // Sbase, Vbase, nAgent, nCons, nGenSup, nBus, nLine, V0, theta0
	Info.setFromFile(fileName1);
	//Info.display();
	
	
	_Sbase = Info.get(0, 0);
	_Vbase = Info.get(0, 1);
	_Zbase = _Vbase * _Vbase / _Sbase;

	_nBus = Info.get(0, 5);
	_nLine = Info.get(0, 6);
	_nConstraint = _nLine + 2 * _nBus;
	_V0 = Info.get(0, 7);
	_theta0 = Info.get(0, 8);

	std::cout << "V0 " << _V0 << " theta0 " << _theta0 << std::endl;
	_zoneBus = MatrixCPU(_nBus, 1);
	

	MatrixCPUD MatLine(_nLine, 10); // from, to, Ys Real, Ys Im, Yp , tau, theta, Limit=0, zs Real, Zs imag
	MatLine.setFromFile(fileName3);
	

	MatrixCPU inversionLine(_nLine, 1, 0);
	if (_nBus == _nLine + 1) {

		for (int l = 0; l < _nLine; l++) {
			if (MatLine.get(l, 0) > MatLine.get(l, 1)) {
				if (MatLine.get(l, 5) > 0 && MatLine.get(l, 5) != 1) {
					std::cout << "la presence de transformateur est peut �tre mal prise en compte ?" << std::endl;
					inversionLine.set(l, 0, 1);
				}
				int temp = MatLine.get(l, 0);
				MatLine.set(l, 0, MatLine.get(l, 1));
				MatLine.set(l, 1, temp);
				MatLine.set(l, 6, -MatLine.get(l, 6));	// ????
			}
		}
		// il faut ordonner pour que Matline(k,1) = k+1;
		std::cout << "changement de l'ordre" << std::endl;
		radial = true;
		for (int l = 0; l < _nLine; l++) {
			while (MatLine.get(l, 1) != l + 2) {
				if (MatLine.get(l, 1) == l) {
					throw std::invalid_argument("problem with branch for distribution network");
				}
				if (MatLine.get(l, 1) - 2 < 0) {
					std::cout << "arret du tri des branchs" << std::endl;
					radial = false;
					break;
				}
				else {
					MatLine.swapLine(l, MatLine.get(l, 1) - 2);
				}
				
			}
		}
	}
	else {
		radial = false;
	}
	
	int offsetbus = 1; 

	//Mat.display();
	
	
	MatrixCPUD MatBus(_nBus, 6); // Gs, Bs, min, max, V0, theta0, zone � rajouter ici ?
	MatBus.setFromFile(fileName4);
	
	
	initMat();
	
	
	bool impedanceToBeDefined = false;
	std::cout << "set Line" << std::endl;
	try
	{
		_lineReactanceD.setFromFile(fileName7);
		_lineSuceptanceD.setFromFile(fileName6);
		_lineReactance = _lineReactanceD;
		_lineSuceptance = _lineSuceptanceD;
	}
	catch (const std::exception&)
	{
		std::cout << " echec du chargement des impedances" << std::endl;
		impedanceToBeDefined = true;
	}

	for (int l = 0; l < _nLine; l++) {
		int i = MatLine.get(l, 0) - offsetbus;
		int j = MatLine.get(l, 1) - offsetbus;
		
		_CoresLineBus.set(l, 0, i);
		_CoresLineBus.set(l, 1, j);
		double limit = MatLine.get(l, 7);
		double ZsRe = MatLine.get(l, 8);
		double ZsIm = MatLine.get(l, 9);
		if (impedanceToBeDefined) {
			double YlsRe = MatLine.get(l, 2); // re(1/Z)
			double YlsIm = MatLine.get(l, 3); // re(1/Z)
			double Ylp = MatLine.get(l, 4); // b/2
			double tau = MatLine.get(l, 5);
			double theta = MatLine.get(l, 6);

			//std::cout << " Ligne numero " << l << " entre bus " << i << " et " << j << " re " << YlsRe << " im " << YlsIm << " b " << Ylp << std::endl;
			if (tau > 0) {
				if (inversionLine.get(l, 0)) {
					int temp = i;
					i = j;
					j = temp;
				}
				double YijRe = (YlsRe * cos(theta) - YlsIm * sin(theta)) / tau;
				double YijImag = (YlsIm * cos(theta) + YlsRe * sin(theta)) / tau;

				double YjiRe = (YlsRe * cos(theta) + YlsIm * sin(theta)) / tau;
				double YjiImag = (YlsIm * cos(theta) - YlsRe * sin(theta)) / tau;
				_lineReactance.increment(i, j, -YijRe);  // Yft
				_lineSuceptance.increment(i, j, -YijImag);

				_lineReactance.increment(j, i, -YjiRe); // Ytf
				_lineSuceptance.increment(j, i, -YjiImag);

				_lineReactance.increment(j, j, YlsRe); // Ytt
				_lineSuceptance.increment(j, j, YlsIm + Ylp);

				_lineReactance.increment(i, i, YlsRe / (tau * tau)); // Yff
				_lineSuceptance.increment(i, i, (YlsIm + Ylp) / (tau * tau));

				_lineReactanceD.increment(i, j, -YijRe);  // Yft
				_lineSuceptanceD.increment(i, j, -YijImag);

				_lineReactanceD.increment(j, i, -YjiRe); // Ytf
				_lineSuceptanceD.increment(j, i, -YjiImag);

				_lineReactanceD.increment(j, j, YlsRe); // Ytt
				_lineSuceptanceD.increment(j, j, YlsIm + Ylp);

				_lineReactanceD.increment(i, i, YlsRe / (tau * tau)); // Yff
				_lineSuceptanceD.increment(i, i, (YlsIm + Ylp) / (tau * tau));

			}
			else {
				_lineReactance.increment(i, j, -YlsRe);
				_lineReactance.increment(j, i, -YlsRe);
				_lineSuceptance.increment(j, i, -YlsIm);
				_lineSuceptance.increment(i, j, -YlsIm);

				_lineReactance.increment(j, j, YlsRe);
				_lineSuceptance.increment(j, j, YlsIm + Ylp);
				_lineReactance.increment(i, i, YlsRe);
				_lineSuceptance.increment(i, i, YlsIm + Ylp);


				_lineReactanceD.increment(i, j, -YlsRe);
				_lineReactanceD.increment(j, i, -YlsRe);
				_lineSuceptanceD.increment(j, i, -YlsIm);
				_lineSuceptanceD.increment(i, j, -YlsIm);

				_lineReactanceD.increment(j, j, YlsRe);
				_lineSuceptanceD.increment(j, j, YlsIm + Ylp);
				_lineReactanceD.increment(i, i, YlsRe);
				_lineSuceptanceD.increment(i, i, YlsIm + Ylp);


				_busSuceptance.increment(i, 0, Ylp);
				_busSuceptance.increment(j, 0, Ylp);
			}
		}
		_lineImpedanceReal.increment(l, 0, ZsRe);
		_lineImpedanceImag.increment(l, 0, ZsIm);
		if (limit > 0) {
			_nLineConstraint++;
			_lineLimits.set(l, 0, limit);
		}
		else {
			_lineLimits.set(l, 0, LINELIMITMAX);
			//_lineLimits.set(i, 0, FLT_MAX);
		}
	}
	
	
	std::cout << "set bus" <<std::endl;

	for (int i = 0; i < _nBus; i++) { // bound on voltage angle rad
		_upperBound.set(i, 0, 3);
		_lowerBound.set(i, 0, -3);
		_VoltageInit.set(i, 0, MatBus.get(i, 5));
		_VoltageInit.set(i + _nBus, 0, MatBus.get(i, 4));
		_VoltageInitD.set(i, 0, MatBus.get(i, 5));
		_VoltageInitD.set(i + _nBus, 0, MatBus.get(i, 4));
		if (impedanceToBeDefined) {
			_lineReactance.increment( i, i, MatBus.get(i, 0) / _Sbase);
			_lineSuceptance.increment(i, i, MatBus.get(i, 1) / _Sbase);

			_lineReactanceD.increment(i, i, MatBus.get(i, 0) / _Sbase);
			_lineSuceptanceD.increment(i, i, MatBus.get(i, 1) / _Sbase);
		}
		

		_busSuceptance.increment(i, 0, MatBus.get(i, 1) / _Sbase);
	}
	for (int i = _nBus; i < 2 * _nBus; i++) { // bound on voltage 
		_upperBound.set(i, 0, MatBus.get(i - _nBus, 3));
		_lowerBound.set(i, 0, MatBus.get(i - _nBus, 2));

	}
	
	int indice = 0;
	for (int i = 2 * _nBus; i < _nConstraint; i++) { // bound on power flow
		_upperBound.set(i, 0,  _lineLimits.get(indice, 0));
		_lowerBound.set(i, 0, -_lineLimits.get(indice, 0));
		indice++;
	}

	_SolutionPF = MatrixCPUD(_nBus, 4);
	_SolutionPF.setFromFile(fileName5);


	LinearizeImp();

}

void StudyCaseACGrid::SetEuropeTestFeeder(std::string path)
{
	std::string fileName1 = path + "CaseTestFeeder.txt";
	std::string fileName2 = path + "BranchTestFeeder.txt";
	

	MatrixCPUD Info(1, 8); // Sbase, Vbase, Zbase, nAgent, nBus, nLine, V0, theta0
	Info.setFromFile(fileName1);
	//Info.display();

	_Sbase = Info.get(0, 0);
	_Vbase = Info.get(0, 1);
	_Zbase = Info.get(0, 2);

	_nBus = Info.get(0, 4);
	_nLine = Info.get(0, 5);
	_nConstraint = _nLine + 2 * _nBus;
	_V0 = Info.get(0, 6);
	_theta0 = Info.get(0, 7);
	std::cout << _nBus << " " << _nLine << " " << _nConstraint << std::endl;

	MatrixCPUD MatLine(_nLine, 10); // from, to, Ys real, Ys Im, Yp, tau, thetha , limit, zs real, zs Imag
	MatLine.setFromFile(fileName2);
	MatrixCPU inversionLine(_nLine, 1, 0);

	for (int l = 0; l < _nLine; l++) {
		if (MatLine.get(l, 0) > MatLine.get(l, 1)) {
			int temp = MatLine.get(l, 0);
			MatLine.set(l, 0, MatLine.get(l, 1));
			MatLine.set(l, 1, temp);
		}
	}
	
	radial = (_nBus==(_nLine+1));
	if (!radial) {
		throw std::invalid_argument("problem with the study case, not radial");
	}
	for (int l = 0; l < _nLine; l++) { // s'il y a des bus non reli� et donc un resau pas radial cela ne marchera pas
		while (MatLine.get(l, 1) != l + 1) {
			MatLine.swapLine(l, MatLine.get(l, 1) - 1);
		}
	}

	initMat();
	
	
	
	for (int l = 0; l < _nLine; l++) {
		int i = MatLine.get(l, 0);
		int j = MatLine.get(l, 1);

		_CoresLineBus.set(l, 0, i);
		_CoresLineBus.set(l, 1, j);
		double limit = MatLine.get(l, 7);
		double ZsRe = MatLine.get(l, 8);
		double ZsIm = MatLine.get(l, 9);
		
		double YlsRe = MatLine.get(l, 2); // re(1/Z)
		double YlsIm = MatLine.get(l, 3); // re(1/Z)
		double Ylp = MatLine.get(l, 4); // b/2
		double tau = MatLine.get(l, 5);
		double theta = MatLine.get(l, 6);

		if (tau > 0 || theta > 0) {
			throw std::invalid_argument("WIP, transformers not taking into account for now in this case");
		}

		_lineReactance.increment(i, j, -YlsRe);
		_lineReactance.increment(j, i, -YlsRe);
		_lineSuceptance.increment(j, i, -YlsIm);
		_lineSuceptance.increment(i, j, -YlsIm);

		_lineReactance.increment(j, j, YlsRe);
		_lineSuceptance.increment(j, j, YlsIm + Ylp);
		_lineReactance.increment(i, i, YlsRe);
		_lineSuceptance.increment(i, i, YlsIm + Ylp);

		_lineReactanceD.increment(i, j, -YlsRe);
		_lineReactanceD.increment(j, i, -YlsRe);
		_lineSuceptanceD.increment(j, i, -YlsIm);
		_lineSuceptanceD.increment(i, j, -YlsIm);

		_lineReactanceD.increment(j, j, YlsRe);
		_lineSuceptanceD.increment(j, j, YlsIm + Ylp);
		_lineReactanceD.increment(i, i, YlsRe);
		_lineSuceptanceD.increment(i, i, YlsIm + Ylp);

		_busSuceptance.increment(i, 0, Ylp);
		_busSuceptance.increment(j, 0, Ylp);
			
		_lineImpedanceReal.increment(l, 0, ZsRe);
		_lineImpedanceImag.increment(l, 0, ZsIm);
		if (limit > 0) {
			_nLineConstraint++;
			_lineLimits.set(l, 0, limit);
		}
		else {
			_lineLimits.set(l, 0, LINELIMITMAX);
		}
	}


	std::cout << "set bus" << std::endl;

	for (int i = 0; i < _nBus; i++) { // bound on voltage angle rad
		_upperBound.set(i, 0, 3);
		_lowerBound.set(i, 0, -3);
		_VoltageInit.set(i, 0, _theta0);
		_VoltageInitD.set(i, 0, _theta0);
	}
	for (int i = _nBus; i < 2 * _nBus; i++) { // bound on voltage 
		_upperBound.set(i, 0, 1.1);
		_lowerBound.set(i, 0, 0.9);
		_VoltageInit.set(i, 0, _V0);
		_VoltageInitD.set(i, 0, _V0);
	}
	_upperBound.set(0, 0, 0);
	_lowerBound.set(0, 0, 0);
	_upperBound.set(_nBus, 0, 1);
	_lowerBound.set(_nBus, 0, 1);

	int indice = 0;
	for (int i = 2 * _nBus; i < _nConstraint; i++) { // bound on power flow
		_upperBound.set(i, 0, _lineLimits.get(indice, 0));
		_lowerBound.set(i, 0, -_lineLimits.get(indice, 0));
		indice++;
	}

	LinearizeImp();

}

void StudyCaseACGrid::setLineLimit(int line, float limit)
{
	if (line > _nLine) {
		throw std::invalid_argument("this line doesn't exist");
	}
	if(limit<=0) {
		throw std::invalid_argument("limit must be positive");
	}
	else {
		float oldLimit = _lineLimits.get(line, 0);
		if (oldLimit == LINELIMITMAX) {
			_nLineConstraint++;
		}
		_lineLimits.set(line, 0, limit);
		_upperBound.set(2 * _nBus + line, 0, limit);
		_lowerBound.set(2 * _nBus + line, 0, -limit);
	}
}

void StudyCaseACGrid::setFromInterface(StudyCaseInterface interface) {
	

	MatrixCPU Info = interface.getInfoCase();
	 // Sbase, Vbase, nAgent, nCons, nGenSup, nBus, nLine, V0, theta0
	
	//Info.display();
	
	
	_Sbase = Info.get(0, Sbase_ind);
	_Vbase = Info.get(0, Vbase_ind);
	_Zbase = _Vbase * _Vbase / _Sbase;

	_nBus = Info.get(0, nBus_ind);
	_nLine = Info.get(0, nLine_ind);
	_nConstraint = _nLine + 2 * _nBus;
	_V0 = Info.get(0, V0_ind);
	_theta0 = Info.get(0, theta0_ind);

	std::cout << "V0 " << _V0 << " theta0 " << _theta0 << std::endl;
	_zoneBus = MatrixCPU(_nBus, 1);
	

	MatrixCPU MatLine = interface.getBranchCase();
	// from, to, Ys Real, Ys Im, Yp , tau, theta, Limit=0, zs Real, Zs imag

	MatrixCPU inversionLine(_nLine, 1, 0);

	if (_nBus == _nLine + 1) {
		for (int l = 0; l < _nLine; l++) {
			if (MatLine.get(l, From_ind) > MatLine.get(l, To_ind)) {
				if (MatLine.get(l, theta_ind) > 0 && MatLine.get(l, tau_ind) != 1) {
					std::cout << "la presence de transformateur est peut etre mal prise en compte ?" << std::endl;
					inversionLine.set(l, 0, 1);
				}
				int temp = MatLine.get(l, From_ind);
				MatLine.set(l, From_ind, MatLine.get(l, To_ind));
				MatLine.set(l, To_ind, temp);
				MatLine.set(l, theta_ind, -MatLine.get(l, theta_ind));	// ????
			}
		}
		// il faut ordonner pour que Matline(k,1) = k+1;
		std::cout << "changement de l'ordre" << std::endl;
		radial = true;
		for (int l = 0; l < _nLine; l++) {
			while (MatLine.get(l, To_ind) != l + 2) {
				if (MatLine.get(l,To_ind) == l) {
					throw std::invalid_argument("problem with branch for distribution network");
				}
				if (MatLine.get(l, To_ind) - 2 < 0) {
					std::cout << "arret du tri des branchs" << std::endl;
					radial = false;
					break;
				}
				else {
					MatLine.swapLine(l, MatLine.get(l, To_ind) - 2);
				}
				
			}
		}
	}
	else {
		radial = false;
	}
	
	int offsetbus = 100; 
	
	for(int l =0; l<_nLine;l++){
		if(MatLine.get(l,From_ind)<offsetbus){
			offsetbus = MatLine.get(l,From_ind);
		}
	}

	//Mat.display();
	
	
	MatrixCPU MatBus = interface.getBusCase();
	 // Gs, Bs, min, max, V0, theta0, zone � rajouter ici ?
		
	initMat();
	
	
	bool impedanceToBeDefined = false;
	std::cout << "set Line" << std::endl;
	
	if(interface.isImpedanceDefined()){
		//_lineReactanceD.setFromFile(fileName7);
		//_lineSuceptanceD.setFromFile(fileName6);
		_lineReactance = interface.getGmat();
		_lineSuceptance = interface.getBmat();
		_lineReactance.toMatCPUD(_lineReactanceD);
		_lineSuceptance.toMatCPUD(_lineSuceptanceD);
	} else{
		impedanceToBeDefined = true;
	}
	

	for (int l = 0; l < _nLine; l++) {
		int i = MatLine.get(l, From_ind) - offsetbus;
		int j = MatLine.get(l, To_ind) - offsetbus;
		
		_CoresLineBus.set(l, 0, i);
		_CoresLineBus.set(l, 1, j);
		double limit = MatLine.get(l, lim_ind);
		double ZsRe = MatLine.get(l, ZsRe_ind);
		double ZsIm = MatLine.get(l, ZsIm_ind);
		if (impedanceToBeDefined) {
			double YlsRe = MatLine.get(l, YsRe_ind); // re(1/Z)
			double YlsIm = MatLine.get(l, YsIm_ind); // re(1/Z)
			double Ylp = MatLine.get(l, Yp_ind); // b/2
			double tau = MatLine.get(l, tau_ind);
			double theta = MatLine.get(l, theta_ind);

			//std::cout << " Ligne numero " << l << " entre bus " << i << " et " << j << " re " << YlsRe << " im " << YlsIm << " b " << Ylp << std::endl;
			if (tau > 0) {
				if (inversionLine.get(l, 0)) {
					int temp = i;
					i = j;
					j = temp;
				}
				double YijRe = (YlsRe * cos(theta) - YlsIm * sin(theta)) / tau;
				double YijImag = (YlsIm * cos(theta) + YlsRe * sin(theta)) / tau;

				double YjiRe = (YlsRe * cos(theta) + YlsIm * sin(theta)) / tau;
				double YjiImag = (YlsIm * cos(theta) - YlsRe * sin(theta)) / tau;
				_lineReactance.increment(i, j, -YijRe);  // Yft
				_lineSuceptance.increment(i, j, -YijImag);

				_lineReactance.increment(j, i, -YjiRe); // Ytf
				_lineSuceptance.increment(j, i, -YjiImag);

				_lineReactance.increment(j, j, YlsRe); // Ytt
				_lineSuceptance.increment(j, j, YlsIm + Ylp);

				_lineReactance.increment(i, i, YlsRe / (tau * tau)); // Yff
				_lineSuceptance.increment(i, i, (YlsIm + Ylp) / (tau * tau));

				_lineReactanceD.increment(i, j, -YijRe);  // Yft
				_lineSuceptanceD.increment(i, j, -YijImag);

				_lineReactanceD.increment(j, i, -YjiRe); // Ytf
				_lineSuceptanceD.increment(j, i, -YjiImag);

				_lineReactanceD.increment(j, j, YlsRe); // Ytt
				_lineSuceptanceD.increment(j, j, YlsIm + Ylp);

				_lineReactanceD.increment(i, i, YlsRe / (tau * tau)); // Yff
				_lineSuceptanceD.increment(i, i, (YlsIm + Ylp) / (tau * tau));

			}
			else {
				_lineReactance.increment(i, j, -YlsRe);
				_lineReactance.increment(j, i, -YlsRe);
				_lineSuceptance.increment(j, i, -YlsIm);
				_lineSuceptance.increment(i, j, -YlsIm);

				_lineReactance.increment(j, j, YlsRe);
				_lineSuceptance.increment(j, j, YlsIm + Ylp);
				_lineReactance.increment(i, i, YlsRe);
				_lineSuceptance.increment(i, i, YlsIm + Ylp);


				_lineReactanceD.increment(i, j, -YlsRe);
				_lineReactanceD.increment(j, i, -YlsRe);
				_lineSuceptanceD.increment(j, i, -YlsIm);
				_lineSuceptanceD.increment(i, j, -YlsIm);

				_lineReactanceD.increment(j, j, YlsRe);
				_lineSuceptanceD.increment(j, j, YlsIm + Ylp);
				_lineReactanceD.increment(i, i, YlsRe);
				_lineSuceptanceD.increment(i, i, YlsIm + Ylp);


				_busSuceptance.increment(i, 0, Ylp);
				_busSuceptance.increment(j, 0, Ylp);
			}
		}
		_lineImpedanceReal.increment(l, 0, ZsRe);
		_lineImpedanceImag.increment(l, 0, ZsIm);
		if (limit > 0) {
			_nLineConstraint++;
			_lineLimits.set(l, 0, limit);
		}
		else {
			_lineLimits.set(l, 0, LINELIMITMAX);
			//_lineLimits.set(i, 0, FLT_MAX);
		}
	}
	
	
	std::cout << "set bus" <<std::endl;

	for (int i = 0; i < _nBus; i++) { // bound on voltage angle rad
		_upperBound.set(i, 0, MatBus.get(i,thetamax_ind));
		_lowerBound.set(i, 0, MatBus.get(i,thetamin_ind)); 
		_VoltageInit.set(i, 0, MatBus.get(i, thetainit_ind));
		_VoltageInit.set(i + _nBus, 0, MatBus.get(i, Vinit_ind));
		_VoltageInitD.set(i, 0, MatBus.get(i, thetainit_ind));
		_VoltageInitD.set(i + _nBus, 0, MatBus.get(i, Vinit_ind));
		if (impedanceToBeDefined) {
			_lineReactance.increment( i, i, MatBus.get(i, Gs_ind) / _Sbase);
			_lineSuceptance.increment(i, i, MatBus.get(i, Bs_ind) / _Sbase);

			_lineReactanceD.increment(i,  i, MatBus.get(i, Gs_ind) / _Sbase);
			_lineSuceptanceD.increment(i, i, MatBus.get(i, Bs_ind) / _Sbase);
		}
		

		_busSuceptance.increment(i, 0, MatBus.get(i, Bs_ind) / _Sbase);
	}
	for (int i = _nBus; i < 2 * _nBus; i++) { // bound on voltage 
		_upperBound.set(i, 0, MatBus.get(i - _nBus, Vmax_ind));
		_lowerBound.set(i, 0, MatBus.get(i - _nBus, Vmin_ind));

	}
	
	int indice = 0;
	for (int i = 2 * _nBus; i < _nConstraint; i++) { // bound on power flow
		_upperBound.set(i, 0,  _lineLimits.get(indice, 0));
		_lowerBound.set(i, 0, -_lineLimits.get(indice, 0));
		indice++;
	}

	_SolutionPF = MatrixCPUD(_nBus, 4);
	
	LinearizeImp();
}



void StudyCaseACGrid::genGrid(int _nBus, int _nMajorLine, int _minorLine, float ReacMajor, float DeltaReacMajor, float ReacMinor, float DeltaReacMinor, float LlimitMajor, float dLlimitMajor, float LlimitMinor, float dLlimitMinor)
{
	std::cout << "work in progress" << std::endl;
	throw std::invalid_argument("WIP");
}

void StudyCaseACGrid::genGridBT(int nBus, int Nbranch, int Ndeep, float length, float dlength)
{
	if (nBus > Nbranch * Ndeep) {
		std::cout << nBus << " " << Nbranch << " " << Ndeep << std::endl;
		throw std::invalid_argument("nBus > Nbranch * Ndeep, impossible to build a grid");
	}
	float rapContr = Ndeep / Nbranch;
	float proba = rapContr / (1 + rapContr);


	_Sbase = 1; // 1MW
	_Vbase = 0.4; // 400V
	_Zbase = _Vbase * _Vbase / _Sbase;
	int billion = 1000000000;
	



	_nBus = nBus;
	_nLine = nBus - 1;
	_nConstraint = _nLine + 2 * _nBus;
	_nLineConstraint = 0;
	_V0 = 1;
	_theta0 = 0;

	initMat();

	std::vector<int> branch;
	std::vector<int> distZero(nBus,0);
	int sizeVector = 1;
	branch.push_back(0);
	
	for (int l = 0; l < _nLine; l++) {
		int j = l + 1;
		int i = 0; // random dans vecteur !!!
		float nRandom = rand1();
		int dist = distZero[j];
		sizeVector = branch.size();
		if (sizeVector > Nbranch || nRandom < proba) {
			int indice = 0;
			do
			{
				if (sizeVector > Nbranch) {
					indice = rand1() * (branch.size() - 2) + 1;
				}
				else {
					indice = rand1() * (branch.size() - 1);
				}
				i = branch[indice];
				dist = distZero[i] + 1;
			} while (dist>Ndeep);
			distZero[j] = dist;
			
			if (indice) {
				branch[indice] = j;
			}
			else {
				branch.push_back(j);
			}
		}
		else 
		{
			do
			{
				i = rand1() * l;
				dist = distZero[i] + 1;
			} while (dist > Ndeep);
			distZero[j] = dist;
			branch.push_back(j);
		}
		
		creatLine(l, i, j, length, dlength);
	}
	hasCurrentLimit = true;
	/* Rajouter ici impedance shunt
	_lineReactance.increment(i, i, MatBus.get(i, 0) / _Sbase);
		_lineSuceptance.increment(i, i, MatBus.get(i, 1) / _Sbase);

		_lineReactanceD.increment(i, i, MatBus.get(i, 0) / _Sbase);
		_lineSuceptanceD.increment(i, i, MatBus.get(i, 1) / _Sbase);
	
		_busSuceptance.increment(i, 0, MatBus.get(i, 1) / _Sbase);*/
	setDefaultConstraint();
	
	LinearizeImp();

}


void StudyCaseACGrid::genGridBTSpecial(int nBus, int Nbranch, int Ndeep, float length, float dlength, RadialType type) {

	switch (type)
	{
	case Normal:
		genGridBT(nBus, Nbranch, Ndeep, length, dlength);
		break;
	case Line:
		genGridLine(nBus, length, dlength);
		break;
	case Balance:
		genGridBalance(nBus, length, dlength);
		break;
	case OneStep:
		genGridOneStep(nBus, length, dlength);
		break;
	default:
		throw std::invalid_argument("unknown radialType");
		break;
	}

	 
}


MatrixCPU StudyCaseACGrid::getCurrentLimit() const
{
	return _currentLimit;
}

void StudyCaseACGrid::genGridHTB(int nBus, int nLine, int dnLine, float length, float dlength)
{
	//std::cout << "work in progress" << std::endl;
	//throw std::invalid_argument("WIP");
	float rapLineBus = 2.0 * nLine / nBus; // il vaut mieux que cela tombe juste sinon les arrondis vont tous casser

	if (dnLine >= rapLineBus || nLine < nBus) {
		std::cout << dnLine << " " << nLine << " " << nBus << std::endl;
		throw std::invalid_argument("dnLine > rapLineBus or nLine< nBus, impossible to build a grid");
	}
	if (nBus * (nBus - 1) / 2 < nLine) {
		std::cout << nBus << " " << nLine << " " << nBus * nBus / 2 << std::endl;
		throw std::invalid_argument("too many lines, impossible to build a grid");
	}
	
	
	
	_Sbase = 1; // 1MW
	_Vbase = 380; // 380kV
	_Zbase = _Vbase * _Vbase / _Sbase;
	int billion = 1000000000;

	_nBus = nBus;
	_nLine = nLine;
	_nConstraint = _nLine + 2 * _nBus;
	_V0 = 1;
	_theta0 = 0;

	initMat();

	
	std::vector<std::vector<int>> alreadyLinked;

	for (int b = 0; b < nBus; b++) {
		std::vector<int> v;
		v.push_back(b);
		alreadyLinked.push_back(v);
	}
	
		
	int indiceLine = 0;
	for (int b = 0; b < nBus-1; b++) {
		//std::cout << " bus " << b << " il y a " << indiceLine << " ligne et il reste " << (_nLine - indiceLine) << " ligne a faire ";
		int nLineBus = randab(rapLineBus - dnLine, rapLineBus + dnLine + 1);
		//std::cout << nLineBus << " ,";
		int oldnLineBus = _nLines.get(b, 0); // nombre de ligne qu'il a d�j�
		int nLineToDo = 0;
		if (nLineBus > oldnLineBus) { // s'il y a besoin d'ajouter de ligne
			nLineToDo = nLineBus - oldnLineBus;
		}
		//std::cout << nLineToDo << " ,";
		nLineToDo = Mymin(nLineToDo, nBus - b - 1);
		//std::cout << nLineToDo << " ,";
		nLineToDo = Mymin(nLineToDo, _nLine - indiceLine - ( nBus-b ) );

		//std::cout << nLineToDo << std::endl;


		int i = b; // bus from
		_nLines.increment(i, 0, nLineToDo); 
		
		for (int l = 0; l < nLineToDo; l++) {
			int j = 0;
			bool mustRerand = false;
			do
			{
				mustRerand = false;
				j = randab(b + 1, nBus); // bus to
				//std::cout << j << std::endl;
				for (int k = 0; k < alreadyLinked[b].size(); k++) {
					if (alreadyLinked[b][k] == j) {
						mustRerand = true;
						break;
					}
				}

			} while (mustRerand);
			alreadyLinked[b].push_back(j);
			alreadyLinked[j].push_back(b);
			
			_nLines.increment(j, 0, 1);
			

			_CoresLineBus.set(indiceLine, 0, b);
			_CoresLineBus.set(indiceLine, 1, j);

			double limit = 0;
			double L = length + 2 * (rand1() - 0.5) * dlength;
			double ZsRe = L * HTBLINE[0] / _Zbase;
			double ZsIm = L * HTBLINE[1] / _Zbase;

			double YlsRe = ZsRe / (ZsRe * ZsRe + ZsIm * ZsIm);
			double YlsIm = -ZsIm / (ZsRe * ZsRe + ZsIm * ZsIm);
			double Ylp = 314.15 * L * HTBLINE[2] / (2 * billion * _Zbase); // b/2

			//std::cout << " Ligne numero " << l << " entre bus " << i << " et " << j << " re " << YlsRe << " im " << YlsIm << " b " << Ylp << std::endl;



			_lineReactance.increment(i, j, -YlsRe);
			_lineReactance.increment(j, i, -YlsRe);
			_lineSuceptance.increment(j, i, -YlsIm);
			_lineSuceptance.increment(i, j, -YlsIm);

			_lineReactance.increment(j, j, YlsRe);
			_lineSuceptance.increment(j, j, YlsIm + Ylp);
			_lineReactance.increment(i, i, YlsRe);
			_lineSuceptance.increment(i, i, YlsIm + Ylp);


			_lineReactanceD.increment(i, j, -YlsRe);
			_lineReactanceD.increment(j, i, -YlsRe);
			_lineSuceptanceD.increment(j, i, -YlsIm);
			_lineSuceptanceD.increment(i, j, -YlsIm);

			_lineReactanceD.increment(j, j, YlsRe);
			_lineSuceptanceD.increment(j, j, YlsIm + Ylp);
			_lineReactanceD.increment(i, i, YlsRe);
			_lineSuceptanceD.increment(i, i, YlsIm + Ylp);


			_busSuceptance.increment(i, 0, Ylp);
			_busSuceptance.increment(j, 0, Ylp);


			_lineImpedanceReal.increment(indiceLine, 0, ZsRe);
			_lineImpedanceImag.increment(indiceLine, 0, ZsIm);
			if (limit > 0) {
				_nLineConstraint++;
				_lineLimits.set(i, 0, limit);
			}
			else {
				_lineLimits.set(i, 0, LINELIMITMAX);
			}
			indiceLine++;

		}
	}

	for (int l = indiceLine; l < _nLine; l++) {
		int i = 0;
		int j = 0;

		bool mustRerand = false;
		do
		{
			mustRerand = false;
			i = randab(0, nBus);
			if (_nLines.get(i, 0) == _nBus - 1) {
				mustRerand = true;
			}
		} while (mustRerand);


		do
		{
			mustRerand = false;
			j = randab(0, nBus); // bus to

			for (int k = 0; k < alreadyLinked[i].size(); k++) {
				if (alreadyLinked[i][k] == j) {
					mustRerand = true;
					break;
				}
			}

		} while (mustRerand);
		alreadyLinked[i].push_back(j);
		alreadyLinked[j].push_back(i);

		_nLines.increment(j, 0, 1);
		_nLines.increment(i, 0, 1);

		_CoresLineBus.set(l, 0, i);
		_CoresLineBus.set(l, 1, j);

		double limit = 0;
		double L = length + 2 * (rand1() - 0.5) * dlength;
		double ZsRe = L * HTBLINE[0] / _Zbase;
		double ZsIm = L * HTBLINE[1] / _Zbase;

		double YlsRe = ZsRe / (ZsRe * ZsRe + ZsIm * ZsIm);
		double YlsIm = -ZsIm / (ZsRe * ZsRe + ZsIm * ZsIm);
		double Ylp = 314.15 * L * HTBLINE[2] / (2 * billion * _Zbase); // b/2

		//std::cout << " Ligne numero " << l << " entre bus " << i << " et " << j << " re " << YlsRe << " im " << YlsIm << " b " << Ylp << std::endl;



		_lineReactance.increment(i, j, -YlsRe);
		_lineReactance.increment(j, i, -YlsRe);
		_lineSuceptance.increment(j, i, -YlsIm);
		_lineSuceptance.increment(i, j, -YlsIm);

		_lineReactance.increment(j, j, YlsRe);
		_lineSuceptance.increment(j, j, YlsIm + Ylp);
		_lineReactance.increment(i, i, YlsRe);
		_lineSuceptance.increment(i, i, YlsIm + Ylp);


		_lineReactanceD.increment(i, j, -YlsRe);
		_lineReactanceD.increment(j, i, -YlsRe);
		_lineSuceptanceD.increment(j, i, -YlsIm);
		_lineSuceptanceD.increment(i, j, -YlsIm);

		_lineReactanceD.increment(j, j, YlsRe);
		_lineSuceptanceD.increment(j, j, YlsIm + Ylp);
		_lineReactanceD.increment(i, i, YlsRe);
		_lineSuceptanceD.increment(i, i, YlsIm + Ylp);


		_busSuceptance.increment(i, 0, Ylp);
		_busSuceptance.increment(j, 0, Ylp);


		_lineImpedanceReal.increment(l, 0, ZsRe);
		_lineImpedanceImag.increment(l, 0, ZsIm);
		if (limit > 0) {
			_nLineConstraint++;
			_lineLimits.set(i, 0, limit);
		}
		else {
			_lineLimits.set(i, 0, LINELIMITMAX);
		}

	}


	//std::cout << _nLines.sum() << " " << _nLine << std::endl;


	for (int i = 0; i < _nBus; i++) { // bound on voltage angle rad
		_upperBound.set(i, 0, 3);
		_lowerBound.set(i, 0, -3);
		_upperBound.set(i + _nBus, 0, 1.1 * _V0);
		_lowerBound.set(i + _nBus, 0, 0.9 * _V0);
		_VoltageInit.set(i, 0, 0);
		_VoltageInit.set(i + _nBus, 0, _V0);
		_VoltageInitD.set(i, 0, 0);
		_VoltageInitD.set(i + _nBus, 0, _V0);

		/* Rajouter ici impedance shunt
		_lineReactance.increment(i, i, MatBus.get(i, 0) / _Sbase);
		_lineSuceptance.increment(i, i, MatBus.get(i, 1) / _Sbase);

		_lineReactanceD.increment(i, i, MatBus.get(i, 0) / _Sbase);
		_lineSuceptanceD.increment(i, i, MatBus.get(i, 1) / _Sbase);

		_busSuceptance.increment(i, 0, MatBus.get(i, 1) / _Sbase);*/
	}


	int indice = 0;
	for (int i = 2 * _nBus; i < _nConstraint; i++) { // bound on power flow
		_upperBound.set(i, 0, _lineLimits.get(indice, 0));
		_lowerBound.set(i, 0, -_lineLimits.get(indice, 0));
		indice++;
	}


	LinearizeImp();


}



StudyCaseACGrid::StudyCaseACGrid(const StudyCaseACGrid& s)
{
	clock_t t = clock();
	
	_nLine = s._nLine;
	_nBus = s._nBus;
	_nLineConstraint = s._nLineConstraint;
	_lineLimits = s._lineLimits;
	_CoresLineBus = s._CoresLineBus;

	_name = s._name;
	_currentLimit = s._currentLimit;
	hasCurrentLimit = s.hasCurrentLimit;
	// AC
	_Zbase = s._Zbase;
	_Sbase = s._Sbase;
	_Vbase = s._Vbase;

	_lineSuceptance = s._lineSuceptance; 
	_lineReactance = s._lineReactance;
	_lineSuceptanceD = s._lineSuceptanceD;
	_lineReactanceD = s._lineReactanceD;

	_lineSuceptanceLin = s._lineSuceptanceLin;
	_lineReactanceLin = s._lineReactanceLin;
	_lineSuceptanceLinD = s._lineSuceptanceLinD;
	_lineReactanceLinD = s._lineReactanceLinD;

	_upperBound = s._upperBound;
	_lowerBound = s._lowerBound;

	_CoresVoiLin = s._CoresVoiLin;
	_CoresBusLin = s._CoresBusLin;
	_nLines = s._nLines;

	

	//distribution network
	_busSuceptance = s._busSuceptance;
	_lineImpedanceImag = s._lineImpedanceImag;
	_lineImpedanceReal = s._lineImpedanceReal;


	//Sol
	_SolutionPF = s._SolutionPF;
	_VoltageInit = s._VoltageInit;
	_VoltageInitD = s._VoltageInitD;

	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
}


StudyCaseACGrid& StudyCaseACGrid::operator= (const StudyCaseACGrid& s) 
{
	clock_t t = clock();
	std::cout << " Copie �galit� AC " << std::endl;
	_nLine = s._nLine;
	_nBus = s._nBus;
	_nLineConstraint = s._nLineConstraint;
	_lineLimits = s._lineLimits;
	_CoresLineBus = s._CoresLineBus;

	_name = s._name;
	_currentLimit = s._currentLimit;
	hasCurrentLimit = s.hasCurrentLimit;

	// AC
	_Zbase = s._Zbase;
	_Sbase = s._Sbase;
	_Vbase = s._Vbase;
	_V0 = s._V0;
	_theta0 = s._V0;
	radial = s.radial;

	_lineSuceptance = s._lineSuceptance;
	_lineReactance = s._lineReactance;
	_lineSuceptanceD = s._lineSuceptanceD;
	_lineReactanceD = s._lineReactanceD;

	_lineSuceptanceLin = s._lineSuceptanceLin;
	_lineReactanceLin = s._lineReactanceLin;
	_lineSuceptanceLinD = s._lineSuceptanceLinD;
	_lineReactanceLinD = s._lineReactanceLinD;

	
	_upperBound = s._upperBound;
	_lowerBound = s._lowerBound;

	_CoresVoiLin = s._CoresVoiLin;
	_CoresBusLin = s._CoresBusLin;
	_nLines = s._nLines;



	//distribution network
	_busSuceptance = s._busSuceptance;
	_lineImpedanceImag = s._lineImpedanceImag;
	_lineImpedanceReal = s._lineImpedanceReal;


	//Sol
	_SolutionPF = s._SolutionPF;
	_VoltageInit = s._VoltageInit;
	_VoltageInitD = s._VoltageInitD;

	t = clock() - t;
	_timeInit = (float)t / CLOCKS_PER_SEC;
	return *this;
}





MatrixCPU StudyCaseACGrid::getLineLimit() const
{
	return _lineLimits;
}

MatrixCPU StudyCaseACGrid::getCoresLineBus() const
{
	return _CoresLineBus;
}


MatrixCPU StudyCaseACGrid::getLineSuceptance() const
{
	return _lineSuceptance;
}

MatrixCPU StudyCaseACGrid::getLineReactance() const
{
	return _lineReactance;
}

MatrixCPUD StudyCaseACGrid::getLineSuceptanceD() const
{
	return _lineSuceptanceD;
}

MatrixCPUD StudyCaseACGrid::getLineReactanceD() const
{
	return _lineReactanceD;
}

MatrixCPU StudyCaseACGrid::getUpperBound() const
{
	return _upperBound;
}

MatrixCPU StudyCaseACGrid::getLowerBound() const
{
	return _lowerBound;
}



MatrixCPUD StudyCaseACGrid::getSolPF() const
{
	return _SolutionPF;
}


double StudyCaseACGrid::getV0() const
{
	return _V0;
}

double StudyCaseACGrid::gettheta0() const
{
	return _theta0;
}

MatrixCPU StudyCaseACGrid::getZsRe() const
{
	return _lineImpedanceReal;
}

MatrixCPU StudyCaseACGrid::getZsImag() const
{
	return _lineImpedanceImag;
}

MatrixCPU StudyCaseACGrid::getYd() const
{
	return _busSuceptance;
}

MatrixCPU StudyCaseACGrid::getGlin() const
{
	return _lineReactanceLin;
}

MatrixCPU StudyCaseACGrid::getBlin() const
{
	return _lineSuceptanceLin;
}

bool StudyCaseACGrid::isCurrentLimit() const
{
	return hasCurrentLimit;
}

MatrixCPU StudyCaseACGrid::getGlin2() const
{
	return _lineReactanceLin2;
}

MatrixCPU StudyCaseACGrid::getBlin2() const
{
	return _lineSuceptanceLin2;
}

MatrixCPU StudyCaseACGrid::getVoltageInit() const
{
	return _VoltageInit;
}

MatrixCPUD StudyCaseACGrid::getVoltageInitD() const
{
	return _VoltageInitD;
}

MatrixCPUD StudyCaseACGrid::getGlinD() const
{
	return _lineReactanceLinD;
}

MatrixCPUD StudyCaseACGrid::getBlinD() const
{
	return _lineSuceptanceLinD;
}

MatrixCPU StudyCaseACGrid::getNLines() const
{
	return _nLines;
}

MatrixCPU StudyCaseACGrid::getCoresBusLin() const
{
	return _CoresBusLin;
}

MatrixCPU StudyCaseACGrid::getCoresVoiLin() const
{
	return _CoresVoiLin;
}

MatrixCPU StudyCaseACGrid::getZones() const
{
	return _zoneBus;
}

MatrixCPU StudyCaseACGrid::getNLinesBegin() const
{
	return _nLinesBegin;
}

float StudyCaseACGrid::getTimeInit() const
{
	return _timeInit;
}

int StudyCaseACGrid::getLastBus() const
{
	if (!LastBus) {

		int* tabTemp = new int[_nBus];

		for (int i = 0; i < _nBus; i++)
		{
			tabTemp[i] = 0;
		}

		int deepMa = 0;
		int bus = 0;
		for (int lold = 0; lold < _nLine; lold++) {
			int busTo = _CoresLineBus.get(lold, 1);
			int busFrom = _CoresLineBus.get(lold, 0);
			int deep = tabTemp[busFrom] + 1;
			if (deep > deepMa) {
				deepMa = deep;
				bus = busTo;
			}
			tabTemp[busTo] = deep;
		}
		// LastBus = bus; peut pas car const
		DELETEA(tabTemp);
		return bus;
	}
	return LastBus;
}



int StudyCaseACGrid::getNLine() const
{
	return _nLine;
}

int StudyCaseACGrid::getNBus() const
{
	return _nBus;
}



std::string StudyCaseACGrid::getName() const
{
	return _name;
}


void StudyCaseACGrid::saveCSV(const std::string& fileName, bool all)
{

	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	MatrixCPU nombre(1, 3);
	nombre.set(0, 0, _nLine);
	nombre.set(0, 1, _nLineConstraint);
	nombre.set(0, 2, _nBus);
	nombre.saveCSV(fileName, mode);


	std::cout << "work in progress" << std::endl;
	throw std::invalid_argument("WIP");
	


}

void StudyCaseACGrid::display(int type) 
{
	std::cout << "work in progress" << std::endl;
	throw std::invalid_argument("WIP");
	
}

void StudyCaseACGrid::displayLineCores(MatrixCPU* g, bool all)
{
	
	MatrixCPU Cores(getCoresLineBus());
	//Cores.display();
	
	MatrixCPU Limit(getLineLimit());
	//Limit.display();
	if (all) {
		if (Cores.getNLin() == 0) {
			for (int l = 0; l < getNLine(); l++) {
				if (fabs(g->get(l, 0)) > Limit.get(l, 0)) {
					std::cout << "*******Limites depasse : Line n " << l << " line limit " << Limit.get(l, 0)
						<< " flow " << g->get(l, 0) << "********" << std::endl;
				}
				else if (Limit.get(l, 0) - fabs(g->get(l, 0)) < 0.1) {
					std::cout << "+++Limites proche : Line n " << l << " line limit " << Limit.get(l, 0)
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
					std::cout << "*********Limites depasse :Line n " << l << " from node " << Cores.get(l, 0)
						<< " to node " << Cores.get(l, 1) << " line limit " << Limit.get(l, 0)
						<< " flow " << g->get(l, 0) << "********" << std::endl;
				}
				else if (Limit.get(l, 0) - fabs(g->get(l, 0)) < 0.1) {
					std::cout << "+++Limites proche :Line n " << l << " from node " << Cores.get(l, 0)
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
					std::cout << "*********Limites depasse : Line n " << l << " line limit " << Limit.get(l, 0)
						<< " flow " << g->get(l, 0) << "*********"<< std::endl;
				}
				else if (Limit.get(l, 0) - fabs(g->get(l, 0)) < 0.1) {
					std::cout << "+++Limites proche : Line n " << l << " line limit " << Limit.get(l, 0)
						<< " flow " << g->get(l, 0) << "+++" << std::endl;
				}
				
			}
		}
		else {
			for (int l = 0; l < getNLine(); l++) {
				if (fabs(g->get(l, 0)) > Limit.get(l, 0)) {
					std::cout << "*********Limites depass� :Line n " << l << " from node " << Cores.get(l, 0)
						<< " to node " << Cores.get(l, 1) << " line limit " << Limit.get(l, 0)
						<< " flow " << g->get(l, 0) << "*********" << std::endl;
				}
				else if (Limit.get(l, 0) - fabs(g->get(l, 0)) < 0.1) {
					std::cout << "+++Limites proche :Line n " << l << " from node " << Cores.get(l, 0)
						<< " to node " << Cores.get(l, 1) << " line limit " << Limit.get(l, 0)
						<< " flow " << g->get(l, 0) << std::endl;
				}
				
			}
		}
	}
	
	
}

StudyCaseACGrid::~StudyCaseACGrid()
{
#ifdef DEBUG_DESTRUCTOR
	std::cout << "case destructor" << std::endl;
#endif

}


