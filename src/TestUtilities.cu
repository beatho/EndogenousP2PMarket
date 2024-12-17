

#include "../head/TestUtilities.cuh"




int testUtilities() {
	int n = 1;
	if (!testcoefPolynome3From4to2coef1()) return n;
	n++;
	if (!testcoefPolynome3From4to2coef2()) return n;
	n++;
	if (!testcoefPolynome3From4to2coef3()) return n;
	n++;
	if (!testresolveRealPolynome3without2term1()) return n;
	n++;
	if (!testresolveRealPolynome3without2term2()) return n;
	n++;
	if (!testresolveRealPolynome3without2term3()) return n;
	n++;
	if (!testresolveRealPolynome4without2term()) return n;
	n++;
	if (!testresolveRealPolynome4without2term2()) return n;
	n++;
	if (!testresolveRealPolynome4without2termLagrange()) return n;
	n++;
	if (!testresolveRealPolynome4without2term2Lagrange()) return n;
		n++;
	std::cout << "---------- GPU ------------------" << std::endl;
	if (!testresolveRealPolynome3without2termGPU()) return n;
	n++; 
	if (!testresolveRealPolynome4without2termGPU()) return n;
	n++;
	if (!testresolveRealPolynome4without2termGPULagrange()) return n;
	n++;
	std::cout << "---------- Eigen------------------" << std::endl;
	if (!testPolyEigen3()) return n;
	n++;
	if (!testPolyEigen4()) return n;
	n++;
	if (!testresolveRealPolynome3without2termEigen()) return n;
	n++;
	if (!testresolveRealPolynome4without2termEigen()) return n;
	n++;
	std::cout << "---------- NEWTON------------------" << std::endl;
	if (!testresolveRealPolynome3Newton1()) return n;
	n++;
	if (!testresolveRealPolynome3Newton2()) return n;
	n++;
	if (!testresolveRealPolynome3Newton3()) return n;
	n++;
	if (!testresolveRealPolynome4Newton1()) return n;
	n++;
	if (!testresolveRealPolynome4Newton2()) return n;
	n++; 
	std::cout << "---------- Halley ------------------" << std::endl;
	if (!testresolveRealPolynome3Halley1()) return n;
	n++;
	if (!testresolveRealPolynome3Halley2()) return n;
	n++;
	if (!testresolveRealPolynome3Halley3()) return n;
	n++;
	if (!testresolveRealPolynome4Halley1()) return n;
	n++;
	if (!testresolveRealPolynome4Halley2()) return n;
	n++; // 20
	std::cout << "---------- GPU 2------------------" << std::endl;
	if (!testresolveRealPolynome3GPU()) return n;
	n++;  
	if (!testresolveRealPolynome4GPU()) return n;
	n++;
	std::cout << "---------- Laguerre------------------" << std::endl;
	if (!testresolveRealPolynome3Laguerre1()) return n;
	n++;
	if (!testresolveRealPolynome3Laguerre2()) return n;
	n++;
	if (!testresolveRealPolynome3Laguerre3()) return n;
	n++;

	return 0;
}

void compareCPUGPU()
{
	std::string fileName = "rootCPUGPU2.csv";


	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;
	


	// faire varier p et q
	double qmax = 100;
	double qmin = 0;
	double pmax = 100;
	double pmin = -100;
	int nQ = 100;
	int nP = 200;

	int blocksize = 512;

	int N = nQ * nP; // nombre total de polynome à résoudre

	int Numblock = ceil((N + blocksize - 1) / blocksize);

	MatrixCPU Param(1, 10);
	Param.set(0, 0, nQ);
	Param.set(0, 1, nP);
	Param.set(0, 2, qmin);
	Param.set(0, 3, qmax);
	Param.set(0, 4, pmin);
	Param.set(0, 5, pmax);
	Param.set(0, 6, blocksize);
	Param.set(0, 7, Numblock);


	// calcul des pas
	double dP = (pmax - pmin) / nP;
	double dQ = (qmax - qmin) / nQ;

	MatrixGPUD coef(4, N);
	MatrixGPUD rootsGPU(3, N, nanf(""), 1);
	MatrixCPU  rootsCPU(3, N, nanf(""));
	MatrixCPU nRootCPU(1, N);
	MatrixGPUD nRootGPU(1, N, 0, 1);

	int poly = 0;
	for (int i = 0; i < nP; i++) {
		float p = pmin + i * dP;
		for (int j = 0; j < nQ; j++) {
			float q = qmin + j * dQ;
			coef.set(0, poly, 1);
			coef.set(1, poly, 0);
			coef.set(2, poly, p);
			coef.set(3, poly, q);
			poly++;
		}
	}
	
	// resolution CPU
	double coef2[2];
	double root[3];
	t1 = std::chrono::high_resolution_clock::now();
	for (int n = 0; n < N; n++) {
		coef2[0] = coef.get(2, n);
		coef2[1] = coef.get(3, n);
		int nRoot = resolveRealPolynome3without2term(root, coef2);
		nRootCPU.set(0, n, nRoot);
		for (int i = 0; i < nRoot; i++) {
			rootsCPU.set(i, n, root[i]);
		}
	}
	t2 = std::chrono::high_resolution_clock::now();
	
	Param.set(0, 8, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()/ BILLION);

	//resolution GPU

	t1 = std::chrono::high_resolution_clock::now();
	coef.transferGPU();
	resolveSeveralRealPolynome3termGPU << <Numblock, blocksize >> > (nRootGPU._matrixGPU, rootsGPU._matrixGPU, coef._matrixGPU, N);
	cudaDeviceSynchronize();
	t2 = std::chrono::high_resolution_clock::now();


	Param.set(0, 9, (float)std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / BILLION);
	coef.transferCPU();
	rootsGPU.transferCPU();
	nRootGPU.transferCPU();
	//save

	Param.saveCSV(fileName);
	coef.saveCSV(fileName);
	nRootCPU.saveCSV(fileName);
	rootsCPU.saveCSV(fileName);
	nRootGPU.saveCSV(fileName);
	rootsGPU.saveCSV(fileName);




}

bool testcoefPolynome3From4to2coef1()
{
	double coef4[4] = { 1.5, 2.2, -3.5, -4 };
	double coef2[2];

	double pSol = -3.0504;
	double qSol = -1.2922;
	coefPolynome3From4to2coef(coef4, coef2);

	if (abs(pSol - coef2[0])>0.001 || abs(qSol - coef2[1])>0.001) {
		std::cout << "p " << pSol << " " << coef2[0] << " q " << qSol << " " << coef2[1] << std::endl;
		return false;
	}

	return true;
}

bool testcoefPolynome3From4to2coef2()
{
	double coef4[4] = { 0, 2.2, -3.5, -4 };
	double coef2[2];
	try
	{
		coefPolynome3From4to2coef(coef4, coef2);
	}
	catch (const std::exception&)
	{
		return true;
	}

	return false;
}

bool testcoefPolynome3From4to2coef3()
{
	double coef4[4] = { 1, 0, -3.5, -4 };
	double coef2[2];

	double pSol = -3.5;
	double qSol = -4;
	coefPolynome3From4to2coef(coef4, coef2);

	if (abs(pSol - coef2[0]) > 0.001 || abs(qSol - coef2[1]) > 0.001) {
		std::cout << "p " << pSol << " " << coef2[0] << " q " << qSol << " " << coef2[1] << std::endl;
		return false;
	}

	return true;
}

bool testresolveRealPolynome3without2term1() {

	double a = 1;
	double b = -5;
	
	double coef4[4] = { 1, -5, 3, 1 };
	double coef2[2];
	double root1 = 2 - sqrt(5);
	double root2 = 1;
	double root3 = 2 + sqrt(5);

	coefPolynome3From4to2coef(coef4, coef2);
	double root[3];

	int nRoot = resolveRealPolynome3without2term(root, coef2);
	if (nRoot != 3) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}
	for (int k = 0; k < 3; k++) {
		root[k] += -b / (3 * a);
	}

	bool find[3] = { false, false, false };
	for (int k = 0; k < 3; k++) {
		
		if (abs(root[k] - root1)<0.001) {
			find[0] = true;
		}
		else if (abs(root[k] - root2) < 0.001) {
			find[1] = true;
		}
		else if (abs(root[k] - root3) < 0.001) {
			find[2] = true;
		}
		else {
			std::cout << "wrong root " << root[k] << " against " << root1 << " " << root2 << " " << root3 << std::endl;
			return false;
		}
	}

	for (int k = 0; k < 3; k++) {
		if (!find[k]) {
			std::cout << "wrong root " << root[0] << " " << root[1] << " " << root[2] << " against " << root1 << " " << root2 << " " << root3 << std::endl;
			return false;
		}
	}


	return true;
}
bool testresolveRealPolynome3without2term2() {
	double a = 1;
	double b = 0;
	double c = 3;
	double d = 1;

	double coef4[4] = { 1, 0, 3, 1 };
	double coef2[2];
	double root1 = cbrt((-1 + sqrt(5)) / 2) + cbrt((-1 - sqrt(5)) / 2);
	
	coefPolynome3From4to2coef(coef4, coef2);
	double root[1];

	int nRoot = resolveRealPolynome3without2term(root, coef2);
	if (nRoot != 1) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}
	
	root[0] += -b / (3 * a);
	
	
	if (abs(root[0] - root1) > 0.001) {
		std::cout << "wrong root " << root[0] << " against " << root1 << std::endl;
		return false;
	}
		

	return true;
}
bool testresolveRealPolynome3without2term3() {
	double a = 1;
	double b = 2;
	double c = -12.75;
	double d = 11.25;

	double coef4[4] = { 1, 2,  -12.75, 11.25 };
	double coef2[2];
	double root1 = 1.5;
	double root2 = -5;
	

	coefPolynome3From4to2coef(coef4, coef2);
	
	std::cout << "poly 3 " << coef2[0] << " " << coef2[1] << std::endl;
	
	double root[3];

	int nRoot = resolveRealPolynome3without2term(root, coef2);
	if (nRoot == 1) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}
	for (int k = 0; k < 2; k++) {
		root[k] += -b / (3 * a);
	}

	bool find[2] = { false, false };
	for (int k = 0; k < 2; k++) {

		if (abs(root[k] - root1) < 0.001) {
			find[0] = true;
		}
		else if (abs(root[k] - root2) < 0.001) {
			find[1] = true;
		}
		else {
			std::cout << "wrong root " << root[k] << " against " << root1 << " " << root2 << " " << std::endl;
			return false;
		}
	}

	for (int k = 0; k < 2; k++) {
		if (!find[k]) {
			std::cout << "wrong root " << root[0] << " " << root[1] << " against " << root1 << " " << root2 << std::endl;
			return false;
		}
	}
	return true;
}


bool testresolveRealPolynome4without2term()
{
	double rootbis[3];
	double coef4[4] = { 1, 7, 7, -6 };
	double coef3[3] = { 6, -13, 6 };
	double coef2[2];
	coefPolynome3From4to2coef(coef4, coef2);
	int nroot = resolveRealPolynome3without2term(rootbis, coef2);

	for (int k = 0; k < nroot; k++) {
		rootbis[k] += -coef4[1] / (3 * coef4[0]);
	}

	double root1 = rootbis[0];
	double root2 = rootbis[1];
	double root3 = rootbis[2];
	double root4 = 1;

	double root[4];

	int nRoot = resvolveRealPolynome4without2term(root, coef3);
	if (nRoot != 4) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}
	
	bool find[4] = { false, false, false, false};
	for (int k = 0; k < 4; k++) { // si racine multiple ne renvoie pas d'erreur si on ne trouve pas la bonne multiplicité mais les bonnes racines -> pas grave dans notre cas.

		if (abs(root[k] - root1) < 0.001) {
			find[0] = true;
		}
		else if (abs(root[k] - root2) < 0.001) {
			find[1] = true;
		}
		else if (abs(root[k] - root3) < 0.001) {
			find[2] = true;
		}
		else if (abs(root[k] - root4) < 0.001) {
			find[3] = true;
		}
		else {
			std::cout << "wrong root " << root[k] << " against " << root1 << " " << root2 << " " << root3 << " " << root4 << std::endl;
			return false;
		}
	}

	for (int k = 0; k < 4; k++) {
		if (!find[k]) {
			std::cout << "wrong root " << root[0] << " " << root[1] << " " << root[2] << " " << root[3] << " against " << root1 << " " << root2 << " " << root3 << " " << root4 << std::endl;
			return false;
		}
	}


	return true;
}
bool testresolveRealPolynome4without2term2()
{
	double root[4];
	double coef3[3] = { -109.778, -4260.6, -3051.76 };
	
	
	int nRoot = resvolveRealPolynome4without2term(root, coef3);
	if (nRoot != 2) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}
	double rootBis[2] = { -0.707106, 110.132 };

	for (int i = 0; i < nRoot; i++) {
		double r = root[i];
		double poly = coef3[2] + coef3[1] * r + coef3[0] * r * r * r + r * r * r * r;
		if (abs(poly) > 0.000001) {
			std::cout << "wrong root poly= "<< poly << " " << root[0] << " " << root[1] << " against " << rootBis[0] << " " << rootBis[1] << std::endl;
			return false;
		}
	}

	return true;
}


bool testresolveRealPolynome4without2termLagrange()
{
	double rootbis[3];
	double coef4[4] = { 1, 7, 7, -6 };
	double coef3[3] = { 6, -13, 6 };
	double coef2[2];
	coefPolynome3From4to2coef(coef4, coef2);
	int nroot = resolveRealPolynome3without2term(rootbis, coef2);

	for (int k = 0; k < nroot; k++) {
		rootbis[k] += -coef4[1] / (3 * coef4[0]);
	}

	double root1 = rootbis[0];
	double root2 = rootbis[1];
	double root3 = rootbis[2];
	double root4 = 1;

	double root[4];

	int nRoot = resvolveRealPolynome4without2termLagrange(root, coef3);
	if (nRoot != 4) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}

	bool find[4] = { false, false, false, false };
	for (int k = 0; k < 4; k++) { // si racine multiple ne renvoie pas d'erreur si on ne trouve pas la bonne multiplicité mais les bonnes racines -> pas grave dans notre cas.

		if (abs(root[k] - root1) < 0.001) {
			find[0] = true;
		}
		else if (abs(root[k] - root2) < 0.001) {
			find[1] = true;
		}
		else if (abs(root[k] - root3) < 0.001) {
			find[2] = true;
		}
		else if (abs(root[k] - root4) < 0.001) {
			find[3] = true;
		}
		else {
			std::cout << "wrong root " << root[k] << " against " << root1 << " " << root2 << " " << root3 << " " << root4 << std::endl;
			return false;
		}
	}

	for (int k = 0; k < 4; k++) {
		if (!find[k]) {
			std::cout << "wrong root " << root[0] << " " << root[1] << " " << root[2] << " " << root[3] << " against " << root1 << " " << root2 << " " << root3 << " " << root4 << std::endl;
			return false;
		}
	}


	return true;
}
bool testresolveRealPolynome4without2term2Lagrange()
{
	double root[4];
	double coef3[3] = { -109.778, -4260.6, -3051.76 };


	int nRoot = resvolveRealPolynome4without2termLagrange(root, coef3);
	if (nRoot != 2) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}
	double rootBis[2] = { -0.707106, 110.132 };

	for (int i = 0; i < nRoot; i++) {
		double r = root[i];
		double poly = coef3[2] + coef3[1] * r + coef3[0] * r * r * r + r * r * r * r;
		if (abs(poly) > 0.000001) {
			std::cout << "wrong root poly= " << poly << " " << root[0] << " " << root[1] << " against " << rootBis[0] << " " << rootBis[1] << std::endl;
			return false;
		}
	}

	return true;
}

bool testresolveRealPolynome3without2termGPU() {
	
	int nPoly = 3;
	MatrixGPUD coefs(4, nPoly);
	MatrixGPUD roots(3, nPoly, 0, 1);
	MatrixGPUD rootToFind(3, nPoly);
	MatrixGPUD nRoot(nPoly, 1, 0, 1);
	MatrixGPUD nRootToFind(nPoly, 1);

	int poly = 0;
	// --------poly 1-----------
	//double coef4[4] = { 1, -5, 3, 1 };
	coefs.set(0, poly, 1);
	coefs.set(1, poly, -5);
	coefs.set(2, poly, 3);
	coefs.set(3, poly, 1);
	// double root1 = 2 - sqrt(5); 	double root2 = 1; 	double root3 = 2 + sqrt(5);
	rootToFind.set(0, poly, 2 + sqrt(5));
	rootToFind.set(1, poly, 2 - sqrt(5));
	rootToFind.set(2, poly, 1);
	
	nRootToFind.set(poly, 0, 3);
	poly++;
	//--------poly 2-----------
	//double coef4[4] = { 1, 0, 3, 1 };
	coefs.set(0, poly, 1);
	coefs.set(1, poly, 0);
	coefs.set(2, poly, 3);
	coefs.set(3, poly, 1);
	// root1 = cbrt((-1 + sqrt(5)) / 2) + cbrt((-1 - sqrt(5)) / 2);
	rootToFind.set(0, poly, cbrt((-1 + sqrt(5)) / 2) + cbrt((-1 - sqrt(5)) / 2));
	nRootToFind.set(poly, 0, 1);
	poly++;


	// poly 3
	//double coef4[4] = { 1, 2,  -12.75, 11.25 };
	coefs.set(0, poly, 1);
	coefs.set(1, poly, 2);
	coefs.set(2, poly, -12.75);
	coefs.set(3, poly, 11.25);
	// double root1 = 1.5; 	double root2 = -5;
	rootToFind.set(0, poly, 1.5);
	rootToFind.set(1, poly, -5);
	nRootToFind.set(poly, 0, 2);
	poly++;

	coefs.transferGPU();

	resolveSeveralRealPolynome3termGPU << <1, 32 >> > (nRoot._matrixGPU, roots._matrixGPU, coefs._matrixGPU, nPoly);

	nRoot.transferCPU();
	roots.transferCPU();

	nRoot.display();
	nRootToFind.display();

	roots.display();
	rootToFind.display();

	return true;


}
bool testresolveRealPolynome4without2termGPU() {

	int nPoly = 2;
	MatrixGPUD coefs(4, nPoly);
	MatrixGPUD roots(4, nPoly, 0, 1);
	MatrixGPUD rootToFind(4, nPoly);
	MatrixGPUD nRoot(nPoly, 1, 0, 1);
	MatrixGPUD nRootToFind(nPoly, 1);

	int poly = 0;
	// --------poly 1-----------
	// double coef3[3] = { 6, -13, 6 };
		// determination des racines
	double rootbis[3];
	double coef4[4] = { 1, 7, 7, -6 };
	double coef2[2];
	coefPolynome3From4to2coef(coef4, coef2);
	int nroot = resolveRealPolynome3without2term(rootbis, coef2);

	coefs.set(0, poly, 6);
	coefs.set(2, poly, -13);
	coefs.set(3, poly, 6);

	for (int k = 0; k < nroot; k++) {
		rootbis[k] += -coef4[1] / (3 * coef4[0]);
	}
	rootToFind.set(0, poly, rootbis[2]);
	rootToFind.set(1, poly, rootbis[1]);
	rootToFind.set(2, poly, 1);
	rootToFind.set(3, poly, rootbis[0]);

	nRootToFind.set(poly, 0, 4);
	poly++;

	// --------poly 2-----------
// double coef3[3] = { -109.778, -4260.6, -3051.76 };

	coefs.set(0, poly, -109.778);
	coefs.set(2, poly, -4260.6);
	coefs.set(3, poly, -3051.76);	
	rootToFind.set(0, poly, 110.132);
	rootToFind.set(1, poly, -0.707106);


	nRootToFind.set(poly, 0, 2);
	coefs.transferGPU();

	resolveSeveralRealPolynome4WO2termGPU << <1, 32 >> > (nRoot._matrixGPU, roots._matrixGPU, coefs._matrixGPU, nPoly);

	nRoot.transferCPU();
	roots.transferCPU();

	nRoot.display();
	nRootToFind.display();

	roots.display();
	rootToFind.display();

	return true;
}

bool testresolveRealPolynome4without2termGPULagrange() {

	int nPoly = 2;
	MatrixGPUD coefs(4, nPoly);
	MatrixGPUD roots(4, nPoly, 0, 1);
	MatrixGPUD rootToFind(4, nPoly);
	MatrixGPUD nRoot(nPoly, 1, 0, 1);
	MatrixGPUD nRootToFind(nPoly, 1);

	int poly = 0;
	// --------poly 1-----------
	// double coef3[3] = { 6, -13, 6 };
		// determination des racines
	double rootbis[3];
	double coef4[4] = { 1, 7, 7, -6 };
	double coef2[2];
	coefPolynome3From4to2coef(coef4, coef2);
	int nroot = resolveRealPolynome3without2term(rootbis, coef2);

	coefs.set(0, poly, 6);
	coefs.set(2, poly, -13);
	coefs.set(3, poly, 6);

	for (int k = 0; k < nroot; k++) {
		rootbis[k] += -coef4[1] / (3 * coef4[0]);
	}
	rootToFind.set(0, poly, rootbis[2]);
	rootToFind.set(1, poly, rootbis[1]);
	rootToFind.set(2, poly, 1);
	rootToFind.set(3, poly, rootbis[0]);

	nRootToFind.set(poly, 0, 4);
	poly++;

	// --------poly 2-----------
// double coef3[3] = { -109.778, -4260.6, -3051.76 };

	coefs.set(0, poly, -109.778);
	coefs.set(2, poly, -4260.6);
	coefs.set(3, poly, -3051.76);
	rootToFind.set(0, poly, 110.132);
	rootToFind.set(1, poly, -0.707106);


	nRootToFind.set(poly, 0, 2);
	coefs.transferGPU();

	resolveSeveralRealPolynome4WO2termGPULagrange << <1, 32 >> > (nRoot._matrixGPU, roots._matrixGPU, coefs._matrixGPU, nPoly);

	nRoot.transferCPU();
	roots.transferCPU();

	nRoot.display();
	nRootToFind.display();

	roots.display();
	rootToFind.display();

	return true;
}

bool testPolyEigen3()
{
	Eigen::Vector4d coeff(1, 3, -5, 1); //double coef4[4] = { 1, -5, 3, 1 };
	Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
	solver.compute(coeff);
	const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType& r = solver.roots();

	std::cout << r << std::endl;
	
	double root1 = 2 - sqrt(5);
	double root2 = 1;
	double root3 = 2 + sqrt(5);
	double root[3];

	if (r(0).imag() || r(1).imag() || r(2).imag()) {
		std::cout << r(0).imag() << " " << r(1).imag() << " " << r(2).imag() << std::endl;
		return false;
	}
	bool find[3] = { false, false, false };
	
	for (int k = 0; k < 3; k++) {
		root[k] = r(k).real();
	}
	
	for (int k = 0; k < 3; k++) {

		if (abs(root[k] - root1) < 0.001) {
			find[0] = true;
		}
		else if (abs(root[k] - root2) < 0.001) {
			find[1] = true;
		}
		else if (abs(root[k] - root3) < 0.001) {
			find[2] = true;
		}
		else {
			std::cout << "wrong root " << root[k] << " against " << root1 << " " << root2 << " " << root3 << std::endl;
			return false;
		}
	}

	for (int k = 0; k < 3; k++) {
		if (!find[k]) {
			std::cout << "wrong root " << root[0] << " " << root[1] << " " << root[2] << " against " << root1 << " " << root2 << " " << root3 << std::endl;
			return false;
		}
	}
	std::cout << std::endl << std::endl;

	return true;
}

bool testPolyEigen4()
{
	Eigen::VectorXd coeff(5);
	coeff(0) = 6;
	coeff(1) = -13;
	coeff(2) = 0;
	coeff(3) = 6;
	coeff(4) = 1;
	
	Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
	solver.compute(coeff);
	const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType& r = solver.roots();

	std::cout << r << std::endl;

	double rootbis[3];
	double coef4[4] = { 1, 7, 7, -6 };
	double coef3[3] = { 6, -13, 6 };
	double coef2[2];
	coefPolynome3From4to2coef(coef4, coef2);
	int nroot = resolveRealPolynome3without2term(rootbis, coef2);

	for (int k = 0; k < nroot; k++) {
		rootbis[k] += -coef4[1] / (3 * coef4[0]);
	}

	double root1 = rootbis[0];
	double root2 = rootbis[1];
	double root3 = rootbis[2];
	double root4 = 1;

	double root[4];

	bool find[4] = { false, false, false, false };

	for (int k = 0; k < 4; k++) {
		if (r(k).imag()) {
			std::cout << r(k).imag() << std::endl;
			return false;
		}
		root[k] = r(k).real();
	}


	for (int k = 0; k < 4; k++) { // si racine multiple ne renvoie pas d'erreur si on ne trouve pas la bonne multiplicité mais les bonnes racines -> pas grave dans notre cas.

		if (abs(root[k] - root1) < 0.001) {
			find[0] = true;
		}
		else if (abs(root[k] - root2) < 0.001) {
			find[1] = true;
		}
		else if (abs(root[k] - root3) < 0.001) {
			find[2] = true;
		}
		else if (abs(root[k] - root4) < 0.001) {
			find[3] = true;
		}
		else {
			std::cout << "wrong root " << root[k] << " against " << root1 << " " << root2 << " " << root3 << " " << root4 << std::endl;
			return false;
		}
	}

	for (int k = 0; k < 4; k++) {
		if (!find[k]) {
			std::cout << "wrong root " << root[0] << " " << root[1] << " " << root[2] << " " << root[3] << " against " << root1 << " " << root2 << " " << root3 << " " << root4 << std::endl;
			return false;
		}
	}

	std::cout << std::endl << std::endl;

	return true;
}

bool testresolveRealPolynome3without2termEigen()
{
	{
		double a = 1;
		double b = -5;

		double coef4[4] = { 1, -5, 3, 1 };
		double coef2[2];
		double root1 = 2 - sqrt(5);
		double root2 = 1;
		double root3 = 2 + sqrt(5);

		coefPolynome3From4to2coef(coef4, coef2);

		
		double root[3];

		int nRoot = resolveRealPolynome3without2termEigen(root, coef2);
		if (nRoot != 3) {
			std::cout << "wrong number of root " << nRoot << std::endl;
			return false;
		}
		for (int k = 0; k < 3; k++) {
			root[k] += -b / (3 * a);
		}

		bool find[3] = { false, false, false };
		for (int k = 0; k < 3; k++) {

			if (abs(root[k] - root1) < 0.001) {
				find[0] = true;
			}
			else if (abs(root[k] - root2) < 0.001) {
				find[1] = true;
			}
			else if (abs(root[k] - root3) < 0.001) {
				find[2] = true;
			}
			else {
				std::cout << "wrong root " << root[k] << " against " << root1 << " " << root2 << " " << root3 << std::endl;
				return false;
			}
		}

		for (int k = 0; k < 3; k++) {
			if (!find[k]) {
				std::cout << "wrong root " << root[0] << " " << root[1] << " " << root[2] << " against " << root1 << " " << root2 << " " << root3 << std::endl;
				return false;
			}
		}
		std::cout << "root Eigen ";
		for (int k = 0; k < nRoot; k++) {

			std::cout << root[k] << " ";
		}
		std::cout << std::endl;
	}

	{
		double a = 1;
		double b = 0;
		double c = 3;
		double d = 1;

		double coef4[4] = { 1, 0, 3, 1 };
		double coef2[2];
		double root1 = cbrt((-1 + sqrt(5)) / 2) + cbrt((-1 - sqrt(5)) / 2);

		coefPolynome3From4to2coef(coef4, coef2);
		double root[1];

		int nRoot = resolveRealPolynome3without2term(root, coef2);
		if (nRoot != 1) {
			std::cout << "wrong number of root " << nRoot << std::endl;
			return false;
		}

		root[0] += -b / (3 * a);


		if (abs(root[0] - root1) > 0.001) {
			std::cout << "wrong root " << root[0] << " against " << root1 << std::endl;
			return false;
		}
		std::cout << "root Eigen ";
		for (int k = 0; k < nRoot; k++) {

			std::cout << root[k] << " ";
		}
		std::cout << std::endl;
	}


	{
		double a = 1;
		double b = 2;
		double c = -12.75;
		double d = 11.25;

		double coef4[4] = { 1, 2,  -12.75, 11.25 };
		double coef2[2];
		double root1 = 1.5;
		double root2 = -5;


		coefPolynome3From4to2coef(coef4, coef2);

		std::cout << "poly 3 " << coef2[0] << " " << coef2[1] << std::endl;

		double root[3];

		int nRoot = resolveRealPolynome3without2term(root, coef2);
		if (nRoot == 1) {
			std::cout << "wrong number of root " << nRoot << std::endl;
			return false;
		}
		for (int k = 0; k < 2; k++) {
			root[k] += -b / (3 * a);
		}

		bool find[2] = { false, false };
		for (int k = 0; k < 2; k++) {

			if (abs(root[k] - root1) < 0.001) {
				find[0] = true;
			}
			else if (abs(root[k] - root2) < 0.001) {
				find[1] = true;
			}
			else {
				std::cout << "wrong root " << root[k] << " against " << root1 << " " << root2 << " " << std::endl;
				return false;
			}
		}

		for (int k = 0; k < 2; k++) {
			if (!find[k]) {
				std::cout << "wrong root " << root[0] << " " << root[1] << " against " << root1 << " " << root2 << std::endl;
				return false;
			}
		}
		std::cout << "root Eigen ";
		for (int k = 0; k < nRoot; k++) {

			std::cout << root[k] << " ";
		}
		std::cout << std::endl;
	}
	return true;
}


bool testresolveRealPolynome4without2termEigen()
{
	{
		double rootbis[3];
		double coef4[4] = { 1, 7, 7, -6 };
		double coef3[3] = { 6, -13, 6 };
		double coef2[2];
		coefPolynome3From4to2coef(coef4, coef2);
		int nroot = resolveRealPolynome3without2termEigen(rootbis, coef2);

		for (int k = 0; k < nroot; k++) {
			rootbis[k] += -coef4[1] / (3 * coef4[0]);
		}

		double root1 = rootbis[0];
		double root2 = rootbis[1];
		double root3 = rootbis[2];
		double root4 = 1;

		double root[4];

		int nRoot = resvolveRealPolynome4without2termEigen(root, coef3);
		if (nRoot != 4) {
			std::cout << "wrong number of root " << nRoot << std::endl;
			return false;
		}

		bool find[4] = { false, false, false, false };
		for (int k = 0; k < 4; k++) { // si racine multiple ne renvoie pas d'erreur si on ne trouve pas la bonne multiplicité mais les bonnes racines -> pas grave dans notre cas.

			if (abs(root[k] - root1) < 0.001) {
				find[0] = true;
			}
			else if (abs(root[k] - root2) < 0.001) {
				find[1] = true;
			}
			else if (abs(root[k] - root3) < 0.001) {
				find[2] = true;
			}
			else if (abs(root[k] - root4) < 0.001) {
				find[3] = true;
			}
			else {
				std::cout << "wrong root " << root[k] << " against " << root1 << " " << root2 << " " << root3 << " " << root4 << std::endl;
				return false;
			}
		}

		for (int k = 0; k < 4; k++) {
			if (!find[k]) {
				std::cout << "wrong root " << root[0] << " " << root[1] << " " << root[2] << " " << root[3] << " against " << root1 << " " << root2 << " " << root3 << " " << root4 << std::endl;
				return false;
			}
		}
		std::cout << "root Eigen ";
		for (int k = 0; k < nRoot; k++) {
			
			std::cout << root[k] << " ";
		}
		std::cout << std::endl;
	}


	{
		double root[4];
		double coef3[3] = { -109.778, -4260.6, -3051.76 };


		int nRoot = resvolveRealPolynome4without2termEigen(root, coef3);
		if (nRoot != 2) {
			std::cout << "wrong number of root " << nRoot << std::endl;
			return false;
		}
		double rootBis[2] = { -0.707106, 110.132 };

		for (int i = 0; i < nRoot; i++) {
			double r = root[i];
			double poly = coef3[2] + coef3[1] * r + coef3[0] * r * r * r + r * r * r * r;
			if (abs(poly) > 0.000001) {
				std::cout << "wrong root poly= " << poly << " " << root[0] << " " << root[1] << " against " << rootBis[0] << " " << rootBis[1] << std::endl;
				return false;
			}
		}

		std::cout << "root Eigen ";
		for (int k = 0; k < nRoot; k++) {

			std::cout << root[k] << " ";
		}
		std::cout << std::endl;
	}
	return true;
}


bool testresolveRealPolynome3Newton1() {

	double a = 1;
	double b = -5;

	double coef3[3] = {-5, 3, 1 };
	double root1 = 2 - sqrt(5);
	double root2 = 1;
	double root3 = 2 + sqrt(5);

	double root[3];

	int nRoot = resolveRealPolynome3Newton(root, coef3);
	if (nRoot != 3) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}
	
	bool find[3] = { false, false, false };
	for (int k = 0; k < 3; k++) {

		if (abs(root[k] - root1) < 0.001) {
			find[0] = true;
		}
		else if (abs(root[k] - root2) < 0.001) {
			find[1] = true;
		}
		else if (abs(root[k] - root3) < 0.001) {
			find[2] = true;
		}
		else {
			std::cout << "wrong root " << root[k] << " against " << root1 << " " << root2 << " " << root3 << std::endl;
			return false;
		}
	}

	for (int k = 0; k < 3; k++) {
		if (!find[k]) {
			std::cout << "wrong root " << root[0] << " " << root[1] << " " << root[2] << " against " << root1 << " " << root2 << " " << root3 << std::endl;
			return false;
		}
	}


	return true;
}
bool testresolveRealPolynome3Newton2() {
	
	double a = 1;
	double b = 0;
	double c = 3;
	double d = 1;

	double coef3[3] = { 0, 3, 1 };
	double root1 = cbrt((-1 + sqrt(5)) / 2) + cbrt((-1 - sqrt(5)) / 2);


	double root[3];

	int nRoot = resolveRealPolynome3Newton(root, coef3);
	
	if (nRoot != 1) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}



	if (abs(root[0] - root1) > 0.001) {
		std::cout << "wrong root " << root[0] << " against " << root1 << std::endl;
		return false;
	}


	return true;
}
bool testresolveRealPolynome3Newton3() {
	double a = 1;
	double b = 2;
	double c = -12.75;
	double d = 11.25;

	double coef3[3] = { 2,  -12.75, 11.25 };
	
	double root1 = 1.5;
	double root2 = -5;



	double root[3];

	int nRoot = resolveRealPolynome3Newton(root, coef3);
	if (nRoot == 1) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}
	
	bool find[2] = { false, false };
	for (int k = 0; k < nRoot; k++) {

		if (abs(root[k] - root1) < 0.001) {
			find[0] = true;
		}
		else if (abs(root[k] - root2) < 0.001) {
			find[1] = true;
		}
		else {
			std::cout << "wrong root " << root[k] << " against " << root1 << " " << root2 << " " << std::endl;
			return false;
		}
	}

	for (int k = 0; k < 2; k++) {
		if (!find[k]) {
			std::cout << "wrong root " << root[0] << " " << root[1] << " against " << root1 << " " << root2 << std::endl;
			return false;
		}
	}
	return true;
}



bool testresolveRealPolynome4Newton1()
{
	double rootbis[3];
	double coef4[4] = { 1, 7, 7, -6 }; //poly 3 avec tous les coefs
	
	double coef4bis[4] = { 6, 0, -13, 6 }; // poly 4 avec les coef sans le premier unitaire

	// juste pour trouver les racines
	double coef2[2];
	coefPolynome3From4to2coef(coef4, coef2);
	int nroot = resolveRealPolynome3without2term(rootbis, coef2);

	for (int k = 0; k < nroot; k++) {
		rootbis[k] += -coef4[1] / (3 * coef4[0]);
	}

	double root1 = rootbis[0];
	double root2 = rootbis[1];
	double root3 = rootbis[2];
	double root4 = 1;


	// vrai test
	double root[4];

	int nRoot = resolveRealPolynome4Newton(root, coef4bis);
	if (nRoot != 4) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}

	bool find[4] = { false, false, false, false };
	for (int k = 0; k < 4; k++) { // si racine multiple ne renvoie pas d'erreur si on ne trouve pas la bonne multiplicité mais les bonnes racines -> pas grave dans notre cas.

		if (abs(root[k] - root1) < 0.001) {
			find[0] = true;
		}
		else if (abs(root[k] - root2) < 0.001) {
			find[1] = true;
		}
		else if (abs(root[k] - root3) < 0.001) {
			find[2] = true;
		}
		else if (abs(root[k] - root4) < 0.001) {
			find[3] = true;
		}
		else {
			std::cout << "wrong root " << root[k] << " against " << root1 << " " << root2 << " " << root3 << " " << root4 << std::endl;
			return false;
		}
	}

	for (int k = 0; k < 4; k++) {
		if (!find[k]) {
			std::cout << "wrong root " << root[0] << " " << root[1] << " " << root[2] << " " << root[3] << " against " << root1 << " " << root2 << " " << root3 << " " << root4 << std::endl;
			return false;
		}
	}


	return true;
}
bool testresolveRealPolynome4Newton2()
{
	double root[4];
	double coef4[4] = { -109.778, 0, -4260.6, -3051.76 };


	int nRoot = resolveRealPolynome4Newton(root, coef4);
	if (nRoot != 2) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}
	double rootBis[2] = { -0.707106, 110.132 };

	for (int i = 0; i < nRoot; i++) {
		double r = root[i];
		double poly = coef4[3] + coef4[2] * r + coef4[0] * r * r * r + r * r * r * r;
		if (abs(poly) > 0.000001) {
			std::cout << "wrong root poly= " << poly << " " << root[0] << " " << root[1] << " against " << rootBis[0] << " " << rootBis[1] << std::endl;
			return false;
		}
	}

	return true;
}



bool testresolveRealPolynome3Laguerre1() {

	double a = 1;
	double b = -5;

	double coef3[3] = { -5, 3, 1 };
	double root1 = 2 - sqrt(5);
	double root2 = 1;
	double root3 = 2 + sqrt(5);

	double root[3];

	int nRoot = resolveRealPolynome3Laguerre(root, coef3);
	if (nRoot != 3) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}

	bool find[3] = { false, false, false };
	for (int k = 0; k < 3; k++) {

		if (abs(root[k] - root1) < 0.001) {
			find[0] = true;
		}
		else if (abs(root[k] - root2) < 0.001) {
			find[1] = true;
		}
		else if (abs(root[k] - root3) < 0.001) {
			find[2] = true;
		}
		else {
			std::cout << "wrong root " << root[k] << " against " << root1 << " " << root2 << " " << root3 << std::endl;
			return false;
		}
	}

	for (int k = 0; k < 3; k++) {
		if (!find[k]) {
			std::cout << "wrong root " << root[0] << " " << root[1] << " " << root[2] << " against " << root1 << " " << root2 << " " << root3 << std::endl;
			return false;
		}
	}


	return true;
}
bool testresolveRealPolynome3Laguerre2() {

	double a = 1;
	double b = 0;
	double c = 3;
	double d = 1;

	double coef3[3] = { 0, 3, 1 };
	double root1 = cbrt((-1 + sqrt(5)) / 2) + cbrt((-1 - sqrt(5)) / 2);


	double root[3];

	int nRoot = resolveRealPolynome3Laguerre(root, coef3);

	if (nRoot != 1) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}



	if (abs(root[0] - root1) > 0.001) {
		std::cout << "wrong root " << root[0] << " against " << root1 << std::endl;
		return false;
	}


	return true;
}
bool testresolveRealPolynome3Laguerre3() {
	double a = 1;
	double b = 2;
	double c = -12.75;
	double d = 11.25;

	double coef3[3] = { 2,  -12.75, 11.25 };

	double root1 = 1.5;
	double root2 = -5;



	double root[3];

	int nRoot = resolveRealPolynome3Laguerre(root, coef3);
	if (nRoot == 1) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}

	bool find[2] = { false, false };
	for (int k = 0; k < nRoot; k++) {

		if (abs(root[k] - root1) < 0.001) {
			find[0] = true;
		}
		else if (abs(root[k] - root2) < 0.001) {
			find[1] = true;
		}
		else {
			std::cout << "wrong root " << root[k] << " against " << root1 << " " << root2 << " " << std::endl;
			return false;
		}
	}

	for (int k = 0; k < 2; k++) {
		if (!find[k]) {
			std::cout << "wrong root " << root[0] << " " << root[1] << " against " << root1 << " " << root2 << std::endl;
			return false;
		}
	}
	return true;
}


bool testresolveRealPolynome3Halley1() {

	double a = 1;
	double b = -5;

	double coef3[3] = { -5, 3, 1 };
	double root1 = 2 - sqrt(5);
	double root2 = 1;
	double root3 = 2 + sqrt(5);

	double root[3];

	int nRoot = resolveRealPolynome3Halley(root, coef3);
	if (nRoot != 3) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}

	bool find[3] = { false, false, false };
	for (int k = 0; k < 3; k++) {

		if (abs(root[k] - root1) < 0.001) {
			find[0] = true;
		}
		else if (abs(root[k] - root2) < 0.001) {
			find[1] = true;
		}
		else if (abs(root[k] - root3) < 0.001) {
			find[2] = true;
		}
		else {
			std::cout << "wrong root " << root[k] << " against " << root1 << " " << root2 << " " << root3 << std::endl;
			return false;
		}
	}

	for (int k = 0; k < 3; k++) {
		if (!find[k]) {
			std::cout << "wrong root " << root[0] << " " << root[1] << " " << root[2] << " against " << root1 << " " << root2 << " " << root3 << std::endl;
			return false;
		}
	}


	return true;
}
bool testresolveRealPolynome3Halley2() {

	double a = 1;
	double b = 0;
	double c = 3;
	double d = 1;

	double coef3[3] = { 0, 3, 1 };
	double root1 = cbrt((-1 + sqrt(5)) / 2) + cbrt((-1 - sqrt(5)) / 2);


	double root[3];

	int nRoot = resolveRealPolynome3Halley(root, coef3);

	if (nRoot != 1) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}



	if (abs(root[0] - root1) > 0.001) {
		std::cout << "wrong root " << root[0] << " against " << root1 << std::endl;
		return false;
	}


	return true;
}
bool testresolveRealPolynome3Halley3() {
	double a = 1;
	double b = 2;
	double c = -12.75;
	double d = 11.25;

	double coef3[3] = { 2,  -12.75, 11.25 };

	double root1 = 1.5;
	double root2 = -5;



	double root[3];

	int nRoot = resolveRealPolynome3Halley(root, coef3);
	if (nRoot == 1) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		std::cout << root[0] << std::endl;
		return false;
	}

	bool find[2] = { false, false };
	for (int k = 0; k < nRoot; k++) {

		if (abs(root[k] - root1) < 0.001) {
			find[0] = true;
		}
		else if (abs(root[k] - root2) < 0.001) {
			find[1] = true;
		}
		else {
			std::cout << "wrong root " << root[k] << " against " << root1 << " " << root2 << " " << std::endl;
			return false;
		}
	}

	for (int k = 0; k < 2; k++) {
		if (!find[k]) {
			std::cout << "wrong root " << root[0] << " " << root[1] << " against " << root1 << " " << root2 << std::endl;
			return false;
		}
	}
	return true;
}

bool testresolveRealPolynome4Halley1()
{
	double rootbis[3];
	double coef4[4] = { 1, 7, 7, -6 }; //poly 3 avec tous les coefs

	double coef4bis[4] = { 6, 0, -13, 6 }; // poly 4 avec les coef sans le premier unitaire

	// juste pour trouver les racines
	double coef2[2];
	coefPolynome3From4to2coef(coef4, coef2);
	int nroot = resolveRealPolynome3without2term(rootbis, coef2);

	for (int k = 0; k < nroot; k++) {
		rootbis[k] += -coef4[1] / (3 * coef4[0]);
	}

	double root1 = rootbis[0];
	double root2 = rootbis[1];
	double root3 = rootbis[2];
	double root4 = 1;


	// vrai test
	double root[4];

	int nRoot = resolveRealPolynome4Halley(root, coef4bis);
	if (nRoot != 4) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		return false;
	}

	bool find[4] = { false, false, false, false };
	for (int k = 0; k < 4; k++) { // si racine multiple ne renvoie pas d'erreur si on ne trouve pas la bonne multiplicité mais les bonnes racines -> pas grave dans notre cas.

		if (abs(root[k] - root1) < 0.001) {
			find[0] = true;
		}
		else if (abs(root[k] - root2) < 0.001) {
			find[1] = true;
		}
		else if (abs(root[k] - root3) < 0.001) {
			find[2] = true;
		}
		else if (abs(root[k] - root4) < 0.001) {
			find[3] = true;
		}
		else {
			std::cout << "wrong root " << root[k] << " against " << root1 << " " << root2 << " " << root3 << " " << root4 << std::endl;
			return false;
		}
	}

	for (int k = 0; k < 4; k++) {
		if (!find[k]) {
			std::cout << "wrong root " << root[0] << " " << root[1] << " " << root[2] << " " << root[3] << " against " << root1 << " " << root2 << " " << root3 << " " << root4 << std::endl;
			return false;
		}
	}


	return true;
}
bool testresolveRealPolynome4Halley2()
{
	double root[4];
	double coef4[4] = { -109.778, 0, -4260.6, -3051.76 };


	int nRoot = resolveRealPolynome4Halley(root, coef4);
	if (nRoot != 2) {
		std::cout << "wrong number of root " << nRoot << std::endl;
		for (int i = 0; i < nRoot; i++) {
			std::cout << root[i] << " ";
		}
		std::cout<<std::endl;
		return false;
	}
	double rootBis[2] = { -0.707106, 110.132 };

	for (int i = 0; i < nRoot; i++) {
		double r = root[i];
		double poly = coef4[3] + coef4[2] * r + coef4[0] * r * r * r + r * r * r * r;
		if (abs(poly) > 0.000001) {
			std::cout << "wrong root poly= " << poly << " " << root[0] << " " << root[1] << " against " << rootBis[0] << " " << rootBis[1] << std::endl;
			return false;
		}
	}

	return true;
}



bool testresolveRealPolynome3GPU() {

	int nPoly = 3;
	MatrixGPUD coefs(4, nPoly);
	MatrixGPUD roots(3, nPoly, 0, 1);
	MatrixGPUD rootToFind(3, nPoly);
	MatrixGPUD nRoot(nPoly, 1, 0, 1);
	MatrixGPUD nRootToFind(nPoly, 1);

	int poly = 0;
	// --------poly 1-----------
	//double coef4[4] = { 1, -5, 3, 1 };
	coefs.set(0, poly, 1);
	coefs.set(1, poly, -5);
	coefs.set(2, poly, 3);
	coefs.set(3, poly, 1);
	// double root1 = 2 - sqrt(5); 	double root2 = 1; 	double root3 = 2 + sqrt(5);
	rootToFind.set(0, poly, 2 + sqrt(5));
	rootToFind.set(1, poly, 2 - sqrt(5));
	rootToFind.set(2, poly, 1);

	nRootToFind.set(poly, 0, 3);
	poly++;
	//--------poly 2-----------
	//double coef4[4] = { 1, 0, 3, 1 };
	coefs.set(0, poly, 1);
	coefs.set(1, poly, 0);
	coefs.set(2, poly, 3);
	coefs.set(3, poly, 1);
	// root1 = cbrt((-1 + sqrt(5)) / 2) + cbrt((-1 - sqrt(5)) / 2);
	rootToFind.set(0, poly, cbrt((-1 + sqrt(5)) / 2) + cbrt((-1 - sqrt(5)) / 2));
	nRootToFind.set(poly, 0, 1);
	poly++;


	// poly 3
	//double coef4[4] = { 1, 2,  -12.75, 11.25 };
	coefs.set(0, poly, 1);
	coefs.set(1, poly, 2);
	coefs.set(2, poly, -12.75);
	coefs.set(3, poly, 11.25);
	// double root1 = 1.5; 	double root2 = -5;
	rootToFind.set(0, poly, 1.5);
	rootToFind.set(1, poly, -5);
	nRootToFind.set(poly, 0, 2);
	poly++;

	coefs.transferGPU();

	resolveSeveralRealPolynome3GPU << <1, 32 >> > (nRoot._matrixGPU, roots._matrixGPU, coefs._matrixGPU, nPoly);

	nRoot.transferCPU();
	roots.transferCPU();

	nRoot.display();
	nRootToFind.display();

	roots.display();
	rootToFind.display();

	return true;


}
bool testresolveRealPolynome4GPU() {

	int nPoly = 2;
	MatrixGPUD coefs(4, nPoly);
	MatrixGPUD roots(4, nPoly, 0, 1);
	MatrixGPUD rootToFind(4, nPoly);
	MatrixGPUD nRoot(nPoly, 1, 0, 1);
	MatrixGPUD nRootToFind(nPoly, 1);

	int poly = 0;
	// --------poly 1-----------
	// double coef3[3] = { 6, -13, 6 };
		// determination des racines
	double rootbis[3];
	double coef4[4] = { 1, 7, 7, -6 };
	double coef2[2];
	coefPolynome3From4to2coef(coef4, coef2);
	int nroot = resolveRealPolynome3without2term(rootbis, coef2);

	coefs.set(0, poly, 6);
	coefs.set(2, poly, -13);
	coefs.set(3, poly, 6);

	for (int k = 0; k < nroot; k++) {
		rootbis[k] += -coef4[1] / (3 * coef4[0]);
	}
	rootToFind.set(0, poly, rootbis[2]);
	rootToFind.set(1, poly, rootbis[1]);
	rootToFind.set(2, poly, 1);
	rootToFind.set(3, poly, rootbis[0]);

	nRootToFind.set(poly, 0, 4);
	poly++;

	// --------poly 2-----------
// double coef3[3] = { -109.778, -4260.6, -3051.76 };

	coefs.set(0, poly, -109.778);
	coefs.set(2, poly, -4260.6);
	coefs.set(3, poly, -3051.76);
	rootToFind.set(0, poly, 110.132);
	rootToFind.set(1, poly, -0.707106);


	nRootToFind.set(poly, 0, 2);
	coefs.transferGPU();

	resolveSeveralRealPolynome4GPU << <1, 32 >> > (nRoot._matrixGPU, roots._matrixGPU, coefs._matrixGPU, nPoly);

	nRoot.transferCPU();
	roots.transferCPU();

	nRoot.display();
	nRootToFind.display();

	roots.display();
	rootToFind.display();

	return true;
}
