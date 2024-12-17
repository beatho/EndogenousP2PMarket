#include "../head/TestADMMGPUConst3.cuh"
#define NSTEPLOCAL 5
#define NMAXPEERPERTRHREAD 5

int testADMMGPUConst3()
{
	int n = 1;

	if (!testADMMGPUConst3Contruct1()) return n;
	n++;
	if (!testADMMGPUConst3Contruct2()) return n;
	n++;
	if (!testADMMGPUConst3Contruct3()) return n;
	n++;
	if (!testADMMGPUConst3LAMBDA()) return n;
	n++;
	if (!testADMMGPUConst3Bt1()) return n;
	n++;
	if (!testADMMGPUConst3TradeP()) return n;
	n++;
	if (!testADMMGPUConst3UpdateRes()) return n;
	n++;
	if (!testADMMGPUConst3CalcRes()) return n;
	n++;
	if (!testADMMGPUConst3Solve1()) return n;
	n++;
	if (!testADMMGPUConst3Solve2()) return n;
	n++;


	return 0;
}

void testADMMGPUConst3Time(int test) {

	switch (test)
	{
	case 0:
		testADMMGPUConst3TimeLAMBDA();
		break;
	case 1:
		testADMMGPUConst3TimeBt1();
		break;
	case 2:
		testADMMGPUConst3TimeTradeP();
		break;
	case 3:
		testADMMGPUConst3TimeUpdateRes();
		break;
	case 4:
		testADMMGPUConst3TimeCalcRes();
		break;
	default:
		std::cout << "No valid input " << std::endl;
		break;
	}
	
}

bool testADMMGPUConst3Contruct1()
{
	std::cout << "contructeur par defaut" << std::endl;
	ADMMGPUConst3 a;
	return true;
}

bool testADMMGPUConst3Contruct2()
{
	float rho = 2;

	std::cout << "contructeur avec parametres" << std::endl;
	ADMMGPUConst3 a(rho);
	return true;
}
bool testADMMGPUConst3Contruct3()
{
	float rho = 2;

	std::cout << "contructeur en deux temps" << std::endl;
	ADMMGPUConst3 a;
	a = ADMMGPUConst3(rho);
	return true;
}

bool testADMMGPUConst3Solve1()
{
	//solve(Simparam* result, Simparam sim, StudyCase cas);
	StudyCase cas;
	cas.Set2node();
	//cas.display();
	int nAgent = cas.getNagent();
	Simparam param(nAgent);
	param.setRho(1);
	Simparam res(param);

	ADMMGPUConst3 a;

	a.solve(&res, param, cas);
	res.display();
	MatrixCPU Trade(nAgent, nAgent);
	Trade.set(0, 1, -1);
	Trade.set(1, 0, 1);
	MatrixCPU Res(res.getRes());
	Res.display();

	MatrixCPU trade = res.getTrade();
	trade.display();
	return trade.isEqual(&Trade, 0.001);

}
bool testADMMGPUConst3Solve2()
{
	//solve(Simparam* result, Simparam sim, StudyCase cas);
	StudyCase cas;
	cas.Set29node();
	//cas.display();
	int nAgent = cas.getNagent();

	Simparam param(nAgent, cas.getNLine());
	float epsG = 0.00002f;
	float epsL = 0.000002f;
	param.setRho(10000);
	param.setEpsL(epsL);
	param.setEpsG(epsG);
	param.setStep(1, 1);
	Simparam res(param);
	ADMMGPUConst3 a;
	a.solve(&res, param, cas);
	res.display();


	float Pn[31] = { -1.008853555,-4.62966156,-2.927534103,-0.8979898691,-0.9462603927,-0.09805059433,-0.127968356,-4.168303013,-3.151874542,-2.261414766,-0.670329392,-3.399893284,-0.4841034412,-2.775528431,-3.008597374,-1.849177122,-0.5534118414,-2.362840891,-1.122991204,-0.1379692554,-2.332088947,4.406820297,5.406073093,3.676487684,3.929354668,4.570535183,2.529039145,3.478654861,2.755935192,3.768760443,4.393183708, };

	MatrixCPU P(31, 1);
	for (int i = 0; i < 31; i++) {
		P.set(i, 0, Pn[i]);
	}

	MatrixCPU P2 = res.getPn();
	P2.display();
	MatrixCPU temp(P2);
	temp.subtract(&P);
	temp.display();
	return P2.isEqual(&P, 0.01);

}

bool testADMMGPUConst3LAMBDA()
{
	int nAgent = 3; // 2 conso et un prod
	int ntrade = 4;
	int blockSize = 256;
	int numBlocks = ceil((nAgent + blockSize - 1) / blockSize);
	float value1 = 2;
	float value2 = -8;
	float value3 = 1.5;
	float value4 = 4;
	MatrixGPU Bt1(ntrade, 1, 0);
	MatrixGPU Bt11(ntrade, 1, -value4 - value1 / value3);
	MatrixGPU LAMBDALin(ntrade, 1, value1);
	MatrixGPU trade(ntrade, 1, value2);
	MatrixGPU CoresLinTrans(ntrade, 1);
	
	MatrixGPU LAMBDALin2(ntrade, 1, value1 + 0.5 * value3 * (value2 + value4));
	float rho = value3;


	CoresLinTrans.set(0, 0, 2);
	CoresLinTrans.set(1, 0, 3);
	CoresLinTrans.set(2, 0, 0);
	CoresLinTrans.set(3, 0, 1);

	trade.set(2, 0, value4);
	trade.set(3, 0, value4);

	
	trade.transferGPU();
	LAMBDALin.transferGPU();
	CoresLinTrans.transferGPU();
	
	updateLAMBDAGPU << <numBlocks, blockSize >> > (LAMBDALin._matrixGPU, trade._matrixGPU, rho, CoresLinTrans._matrixGPU, ntrade);

	
	LAMBDALin.transferCPU();

	return (LAMBDALin.isEqual(&LAMBDALin2));
}

bool testADMMGPUConst3Bt1()
{
	int nAgent = 3; // 2 conso et un prod
	int ntrade = 4;
	int blockSize = 256;
	int numBlocks = ceil((nAgent + blockSize - 1) / blockSize);
	float value1 = 2;
	float value2 = -8;
	float value3 = 1.5;
	float value4 = 4;
	MatrixGPU Bt1(ntrade, 1, 0);
	MatrixGPU Bt11(ntrade, 1, -value4 - value1 / value3);
	MatrixGPU trade(ntrade, 1, value2);
	MatrixGPU CoresLinTrans(ntrade, 1);

	MatrixGPU LAMBDALin(ntrade, 1, value1 + 0.5 * value3 * (value2 + value4));
	float rho = value3;


	CoresLinTrans.set(0, 0, 2);
	CoresLinTrans.set(1, 0, 3);
	CoresLinTrans.set(2, 0, 0);
	CoresLinTrans.set(3, 0, 1);

	trade.set(2, 0, value4);
	trade.set(3, 0, value4);
	Bt11.set(2, 0, -value2 - value1 / value3);
	Bt11.set(3, 0, -value2 - value1 / value3);



	Bt1.transferGPU();
	trade.transferGPU();
	LAMBDALin.transferGPU();

	CoresLinTrans.transferGPU();



	updateBt1GPU << <numBlocks, blockSize >> > (Bt1._matrixGPU, trade._matrixGPU, rho, LAMBDALin._matrixGPU, CoresLinTrans._matrixGPU, ntrade);

	Bt1.transferCPU();
	LAMBDALin.transferCPU();

	return Bt1.isEqual(&Bt11);
}

bool testADMMGPUConst3TradeP()
{
	int nAgent = 4;
	int ntrade = 8;
	int blockSize = 512;
	//int numBlocks = ceil((nAgent + blockSize - 1) / blockSize);
	float value1 = 2;
	float value2 = 3;
	float value3 = 1;
	float value4 = 1;
	float value5 = -1;
	float value6 = -30;
	float value7 = 10;
	float value8 = 5;
	float value9 = value8 - value3 + value4 - value5;
	MatrixGPU Bt1(ntrade, 1, value1);
	float at1 = value3;
	float at2 = value4;
	MatrixGPU Ct(ntrade, 1, value5);
	MatrixGPU Lb(ntrade, 1, value6);
	MatrixGPU Ub(ntrade, 1, value7);
	MatrixGPU Tlocal(ntrade, 1);
	
	MatrixGPU Tlocal_pre(ntrade, 1, value8);
	//Tlocal_pre.set(0,0,value8+1);
	MatrixGPU Tlocal2(Tlocal_pre);
	MatrixGPU CoresLinAgent(ntrade, 1);
	CoresLinAgent.set(1, 0, 0);
	CoresLinAgent.set(2, 0, 1);
	CoresLinAgent.set(3, 0, 1);
	CoresLinAgent.set(4, 0, 2);
	CoresLinAgent.set(5, 0, 2);
	CoresLinAgent.set(6, 0, 3);
	CoresLinAgent.set(7, 0, 3);
	MatrixGPU Bp1(nAgent, 1);
	MatrixGPU Bp11(nAgent, 1);
	MatrixGPU Ap1(nAgent, 1, value2);
	MatrixGPU Ap2(nAgent, 1, value3);
	MatrixGPU Ap12(nAgent, 1, value2 + value3);
	MatrixGPU Cp(nAgent, 1, value5);
	MatrixGPU PLb(nAgent, 1, value6/2);
	MatrixGPU PUb(nAgent, 1, value7/2);

	MatrixGPU nVoisin(nAgent, 1, 2);
	MatrixGPU Tmoy(nAgent, 1, value3);
	MatrixGPU Tmoy2(Tmoy);
	MatrixGPU MU(nAgent, 1, value5);
	MatrixGPU MU2(MU);
	MatrixGPU P(nAgent, 1, value4);
	MatrixGPU P2(P);

	MatrixGPU CoresAgentLin(nAgent + 1, 1);
	CoresAgentLin.set(1, 0, 2);
	CoresAgentLin.set(2, 0, 4);
	CoresAgentLin.set(3, 0, 6);
	CoresAgentLin.set(4, 0, 8);

	for (int iter = 0; iter < NSTEPLOCAL; iter++) {
		for (int i = 0; i < nAgent; i += 1) // 1 bloc = 1 agent
		{
			float s = 0.0;
			for (int j = CoresAgentLin.get(i, 0); j < CoresAgentLin.get(i + 1, 0); j++) // on parcourt les trades de l'agent i
			{
				float m = Tlocal2.get(j, 0) - Tmoy2.get(i, 0) + P2.get(i, 0) - MU2.get(i, 0);
				float r = (Bt1.get(j, 0) * at1 + m * at2 - Ct.get(j, 0)) / (at1 + at2);
				float ub = Ub.get(j, 0);
				float lb = Lb.get(j, 0);
				float t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
				Tlocal2.set(j, 0, t);
				s += t;
			}
			float r = s / nVoisin.get(i, 0);
			Tmoy2.set(i, 0, r);
			Bp11.set(i, 0, r + MU2.get(i, 0));
			float p = (Ap1.get(i, 0) * Bp11.get(i, 0) - Cp.get(i, 0)) / (Ap12.get(i, 0));
			float ub = PUb.get(i, 0);
			float lb = PLb.get(i, 0);
			p = (ub - p) * (p > ub) + (lb - p) * (p < lb) + p;
			P2.set(i, 0, p);
			MU2.set(i, 0, MU2.get(i, 0) + r - P2.get(i, 0));
		}
	}
	


	Bt1.transferGPU();
	Tlocal.transferGPU();
	Ct.transferGPU();
	Lb.transferGPU();
	Ub.transferGPU();
	Tlocal_pre.transferGPU();
	Tmoy.transferGPU();
	P.transferGPU();
	MU.transferGPU();
	CoresLinAgent.transferGPU();
	Ap1.transferGPU();
	Ap2.transferGPU();
	Ap12.transferGPU();
	Bp1.transferGPU();
	Cp.transferGPU();
	PLb.transferGPU();
	PUb.transferGPU();
	nVoisin.transferGPU();
	CoresAgentLin.transferGPU();

	std::cout << " fin transfert " << std::endl;


	updateTradePGPUShared<256> <<<nAgent, blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, at1, at2, Bt1._matrixGPU, Ct._matrixGPU,
		Lb._matrixGPU, Ub._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, PLb._matrixGPU, PUb._matrixGPU, CoresAgentLin._matrixGPU);

	Tlocal.transferCPU();
	P.transferCPU();
	MU.transferCPU();
	Tmoy.transferCPU();
	
	Tlocal2.display();
	Tlocal.display();
	Tmoy2.display();
	Tmoy.display();
	P2.display();
	P.display();
	MU2.display();
	MU.display();

	Bp11.display();


	return ((Tlocal.isEqual(&Tlocal2)) && (P.isEqual(&P2)) && (MU.isEqual(&MU2)) && (Tmoy.isEqual(&Tmoy2)));
}

bool testADMMGPUConst3UpdateRes()
{	/*float ADMMGPU5::updateRes(MatrixCPU* res, MatrixGPU* Tlocal, MatrixGPU* trade, int iter, MatrixGPU* CoresLinAgent, MatrixGPU* CoresLinVoisin, MatrixGPU* CoresMatLin, MatrixGPU* tempNN)
{
	tempNN->subtract(Tlocal, trade);
	//cudaDeviceSynchronize();
	float resS = tempNN->distance2();

	updateDiffGPU <<<_numBlocks, _blockSize >>> (tempNN->_matrixGPU, Tlocal->_matrixGPU, CoresLinTrans->_matrixGPU, _N);
	//cudaDeviceSynchronize();
	float resR = tempNN->distance2();

	res->set(0, iter, resR);
	res->set(1, iter, resS);


	return resR * (resR > resS) + resS * (resR <= resS);*/
	int nAgent = 3;
	int ntrade = 4;
	int blockSize = 15;
	int numBlocks = ceil((ntrade + blockSize - 1) / blockSize);
	std::cout << "blockSize " << blockSize << " numBlocks " << numBlocks << std::endl;
	float value1 = 4;
	float value2 = 2.5;
	float value3 = -2;
	float value4 = value3 + value2 - value1;
	MatrixCPU res(2, 1);
	MatrixCPU res2(2, 1);
	MatrixGPU Tlocal(ntrade, 1, value1);
	MatrixGPU Tlocal_pre(ntrade, 1, value2);

	MatrixGPU CoresLinTrans(ntrade, 1);

	MatrixGPU tempN(numBlocks, 1, 0, 1);
	MatrixGPU tempN2(numBlocks, 1, 0, 1);
	
	

	CoresLinTrans.set(0, 0, 2);
	CoresLinTrans.set(1, 0, 3);
	CoresLinTrans.set(2, 0, 0);
	CoresLinTrans.set(3, 0, 1);

	Tlocal.set(2, 0, value3);
	Tlocal.set(3, 0, value3);
	Tlocal_pre.set(2, 0, value4);
	Tlocal_pre.set(3, 0, value4);



	res2.set(0, 0, sqrtf((value1 + value3) * (value1 + value3) ));
	res2.set(1, 0, sqrtf((value1 - value2) * (value1 - value2) ));
	int iter = 0;

	Tlocal.transferGPU();
	Tlocal_pre.transferGPU();
	CoresLinTrans.transferGPU();
	
	float resS = Tlocal.max2(&Tlocal_pre);
	updateDiffGPU << <numBlocks, blockSize >> > (tempN._matrixGPU, Tlocal._matrixGPU, CoresLinTrans._matrixGPU, ntrade);
	float resR = tempN.max2();
	
	res.set(0, 0, resR);
	res.set(1, 0, resS);

	return res2.isEqual(&res);
}
bool testADMMGPUConst3CalcRes() {
	/*
	* float ADMMGPU5::calcRes( MatrixGPU* Tlocal, MatrixGPU* P, MatrixGPU* tempN1, MatrixGPU* tempNN)
{
	 tempNN->subtract(Tlocal, &Tlocal_pre);
	 tempN1->subtract(&Tmoy, P);

	 float d1 = tempN1->max2();
	 float d2 = tempNN->max2();


	 return d1* (d1 > d2) + d2 * (d2 >= d1);
}*/

	int nAgent = 3;
	int ntrade = 4;
	int blockSize = 256;
	int numBlocks = ceil((ntrade + blockSize - 1) / blockSize);
	float value1 = 5;
	float value2 = 2;
	float value3 = -3;
	float value4 = -1;
	MatrixGPU Tlocal(ntrade, 1, value1);
	MatrixGPU Tlocal_pre(ntrade, 1, value2);
	MatrixGPU Tmoy(nAgent, 1, value3);
	MatrixGPU P(nAgent, 1, value4);


	Tlocal.transferGPU();
	Tlocal_pre.transferGPU();
	Tmoy.transferGPU();
	P.transferGPU();
	

	float d11 = Tlocal.max2(&Tlocal_pre);
	float d22 = P.max2(&Tmoy);
	float d = d11 * (d11 > d22) + d22 * (d22 >= d11);



	float d1 = fabs(value1 - value2);
	float d2 = fabs(value3 - value4);
	float df = d1 * (d1 > d2) + d2 * (d2 >= d1);



	return (df == d);
}


void testADMMGPUConst3TimeLAMBDA()
{
	std::string fileName = "TempsLAMBDA.csv";
	//cudaEvent_t start, stop;
	//float elapsedTime;
	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	unsigned int time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNAgent = 7;
	const int nSimu = 100;
	const int nRepet = 10;
	int nAgent[nNAgent] = { 10, 100, 500, 1000, 5000, 10000, 40000 }; // autant conso que de prod, la dernière veleur ne "marche" pas (trop rapide)
	int ntrade[nNAgent];
	int blockSize = 256;
	float values1[nSimu];
	float values2[nSimu];
	float rhos[nSimu];
	MatrixCPU temps(nNAgent, nSimu, 0);
	
	for (int j = 0; j < nSimu;j++) {
		values1[j] = (float) (rand()) / rand();
		values2[j] = (float)(rand()) / rand();
		rhos[j] = (float)(rand()%100) / rand();
	}

	for (int i = 0; i < nNAgent; i++) {
		
		ntrade[i] = nAgent[i] * nAgent[i] / 2;
		std::cout << "iteration " << i << " nAgent " << nAgent[i] << " ntrade " << ntrade[i] << std::endl;
		int numBlocks = ceil((ntrade[i] + blockSize - 1) / blockSize);
		
		MatrixGPU CoresLinAgent(ntrade[i], 1);
		MatrixGPU CoresLinVoisin(ntrade[i], 1);
		MatrixGPU CoresMatLin(nAgent[i], nAgent[i], -1);
		MatrixGPU CoresLinTrans(ntrade[i], 1);
		

		int indice = 0;
		int idVoisin = 0;

		for (int idAgent = 0; idAgent < nAgent[i]; idAgent++) {
			int Nvoisinmax = nAgent[i] / 2;
			if (idAgent < Nvoisinmax) {
				idVoisin = 0;
			}
			else {
				idVoisin = nAgent[i] / 2;
			}
			for (int voisin = idVoisin; voisin < Nvoisinmax; voisin++) {
				CoresLinAgent.set(indice, 0, idAgent);
				CoresLinVoisin.set(indice, 0, voisin);
				CoresMatLin.set(idAgent, voisin, indice);
				indice = indice + 1;
			}
		}
		for (int lin = 0;lin < ntrade[i];lin++) {
			int i = CoresLinAgent.get(lin, 0);
			int j = CoresLinVoisin.get(lin, 0);
			int k = CoresMatLin.get(j, i);
			CoresLinTrans.set(lin, 0, k);
		}
		CoresLinAgent.transferGPU();
		CoresLinVoisin.transferGPU();
		CoresMatLin.transferGPU();
		CoresLinTrans.transferGPU();
		clock_t t = clock();
		
		for (int simu = 0; simu < nSimu; simu++) {
			MatrixGPU LAMBDALin(ntrade[i], 1, values1[simu], 1);
			MatrixGPU trade(ntrade[i], 1, values2[simu], 1);
			float rho = rhos[simu];
			/*cudaEventCreate(&start);
			cudaEventRecord(start, 0);*/
			cudaDeviceSynchronize();
			time = 0;
			for (int repet = 0; repet < nRepet; repet++) {
				MatrixGPU LAMBDALinCopy(LAMBDALin);
				MatrixGPU tradeCopy(trade);
				MatrixGPU CoresLinTransCopy(CoresLinTrans);
				cudaDeviceSynchronize();
				a = std::chrono::high_resolution_clock::now();
				updateLAMBDAGPU <<<numBlocks, blockSize >>> (LAMBDALinCopy._matrixGPU, tradeCopy._matrixGPU, rho, CoresLinTransCopy._matrixGPU, ntrade[i]);
				cudaDeviceSynchronize();
				b = std::chrono::high_resolution_clock::now();
				time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
			}
			temps.set(i, simu, (float) time/nRepet);
		}
	}
	temps.saveCSV(fileName, mode);
}

void testADMMGPUConst3TimeBt1()
{
	std::string fileName = "TempsBt.csv";
	//cudaEvent_t start, stop;
	//float elapsedTime;
	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	unsigned int time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNAgent = 7;
	const int nSimu = 100;
	const int nRepet = 10;
	int nAgent[nNAgent] = { 10, 100, 500, 1000, 5000, 10000, 40000 }; // autant conso que de prod, la dernière veleur ne "marche" pas (trop rapide)
	int ntrade[nNAgent];
	int blockSize = 256;
	float values1[nSimu];
	float values2[nSimu];
	float rhos[nSimu];
	MatrixCPU temps(nNAgent, nSimu, 0);

	for (int j = 0; j < nSimu;j++) {
		values1[j] = (float)(rand()) / rand();
		values2[j] = (float)(rand()) / rand();
		rhos[j] = (float)(rand() % 100) / rand();
	}

	for (int i = 0; i < nNAgent; i++) {

		ntrade[i] = nAgent[i] * nAgent[i] / 2;
		std::cout << "iteration " << i << " nAgent " << nAgent[i] << " ntrade " << ntrade[i] << std::endl;
		int numBlocks = ceil((ntrade[i] + blockSize - 1) / blockSize);

		MatrixGPU Bt1(ntrade[i], 1, 0, 1);

		MatrixGPU CoresLinAgent(ntrade[i], 1);
		MatrixGPU CoresLinVoisin(ntrade[i], 1);
		MatrixGPU CoresMatLin(nAgent[i], nAgent[i], -1);
		MatrixGPU CoresLinTrans(ntrade[i], 1);


		int indice = 0;
		int idVoisin = 0;

		for (int idAgent = 0; idAgent < nAgent[i]; idAgent++) {
			int Nvoisinmax = nAgent[i] / 2;
			if (idAgent < Nvoisinmax) {
				idVoisin = 0;
			}
			else {
				idVoisin = nAgent[i] / 2;
			}
			for (int voisin = idVoisin; voisin < Nvoisinmax; voisin++) {
				CoresLinAgent.set(indice, 0, idAgent);
				CoresLinVoisin.set(indice, 0, voisin);
				CoresMatLin.set(idAgent, voisin, indice);
				indice = indice + 1;
			}
		}
		for (int lin = 0;lin < ntrade[i];lin++) {
			int i = CoresLinAgent.get(lin, 0);
			int j = CoresLinVoisin.get(lin, 0);
			int k = CoresMatLin.get(j, i);
			CoresLinTrans.set(lin, 0, k);
		}
		CoresLinAgent.transferGPU();
		CoresLinVoisin.transferGPU();
		CoresMatLin.transferGPU();
		CoresLinTrans.transferGPU();
		clock_t t = clock();

		for (int simu = 0; simu < nSimu; simu++) {
			MatrixGPU LAMBDALin(ntrade[i], 1, values1[simu], 1);
			MatrixGPU trade(ntrade[i], 1, values2[simu], 1);
			float rho = rhos[simu];
			/*cudaEventCreate(&start);
			cudaEventRecord(start, 0);*/
			cudaDeviceSynchronize();
			time = 0;
			for (int repet = 0; repet < nRepet; repet++) {
				MatrixGPU LAMBDALinCopy(LAMBDALin);
				MatrixGPU tradeCopy(trade);
				MatrixGPU CoresLinTransCopy(CoresLinTrans);
				MatrixGPU Bt1Copy(Bt1);
				cudaDeviceSynchronize();
				a = std::chrono::high_resolution_clock::now();
				updateBt1GPU << <numBlocks, blockSize >> > (Bt1Copy._matrixGPU, tradeCopy._matrixGPU, rho, LAMBDALinCopy._matrixGPU, CoresLinTransCopy._matrixGPU, ntrade[i]);
				cudaDeviceSynchronize();
				b = std::chrono::high_resolution_clock::now();
				time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
			}
			temps.set(i, simu, (float)time / nRepet);
		}
	}
	temps.saveCSV(fileName, mode);
}


void testADMMGPUConst3TimeTradeP() {
	std::string fileName = "TempsTradeP10_F5.csv";
	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	unsigned int time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNAgent = 6;
	const int nSimu = 100;
	const int nRepet = 10;
	const int nVar = 14;
	int nAgent[nNAgent] = { 10, 100, 500, 1000, 5000, 10000 }; // autant conso que de prod, 
	int ntrade[nNAgent];
	int blockSize = 256;
	float values[nVar][nSimu];
	float rhos[nSimu];
	MatrixCPU temps(nSimu, nNAgent, 0);
	MatrixCPU nAgentMat(1, nNAgent, 0);
	MatrixCPU nTradeMat(1, nNAgent, 0);
	MatrixCPU nPro(1, nNAgent, 0);

	for (int j = 0; j < nSimu; j++) {
		for (int var = 0; var < nVar; var++) {
			values[var][j] = (float)(rand()) / rand();
		}
		rhos[j] = (float)(rand() % 100) / rand();
	}

	for (int i = 0; i < nNAgent; i++) {

		ntrade[i] = nAgent[i] * nAgent[i] / 2;
		std::cout << "iteration " << i << " nAgent " << nAgent[i] << " ntrade " << ntrade[i] << std::endl;
		int numBlocks = nAgent[i];
		nAgentMat.set(0, i, nAgent[i]);
		nTradeMat.set(0, i, ntrade[i]);
		nPro.set(0, i, numBlocks * blockSize);
		MatrixGPU Tlocal(ntrade[i], 1, 0, 1);
		MatrixGPU CoresLinVoisin(ntrade[i], 1);
		MatrixGPU CoresAgentLin(nAgent[i] + 1, 1);

		int indice = 0;
		int idVoisin = 0;

		for (int idAgent = 0; idAgent < nAgent[i]; idAgent++) {
			int Nvoisinmax = nAgent[i] / 2;
			if (idAgent < Nvoisinmax) {
				idVoisin = 0;
			}
			else {
				idVoisin = nAgent[i] / 2;
			}
			for (int voisin = idVoisin; voisin < Nvoisinmax; voisin++) {
				CoresLinVoisin.set(indice, 0, voisin);

				indice = indice + 1;
			}
			CoresAgentLin.set(idAgent + 1, 0, indice);


		}
		CoresLinVoisin.transferGPU();
		CoresAgentLin.transferGPU();
		clock_t t = clock();

		for (int simu = 0; simu < nSimu; simu++) {

			MatrixGPU Tlocal_pre(ntrade[i], 1, values[0][simu], 1);
			MatrixGPU Bt1(ntrade[i], 1, values[1][simu], 1);
			MatrixGPU Ct(ntrade[i], 1, values[2][simu], 1);
			MatrixGPU Lb(ntrade[i], 1, values[3][simu], 1);
			MatrixGPU Ub(ntrade[i], 1, values[4][simu], 1);

			MatrixGPU Tmoy(nAgent[i], 1, values[5][simu], 1);
			MatrixGPU P(nAgent[i], 1, values[6][simu], 1);
			MatrixGPU MU(nAgent[i], 1, values[7][simu], 1);
			MatrixGPU nVoisin(nAgent[i], 1, nAgent[i] / 2, 1);
			MatrixGPU Ap1(nAgent[i], 1, values[9][simu], 1);
			MatrixGPU Ap12(nAgent[i], 1, values[10][simu], 1);
			MatrixGPU Cp(nAgent[i], 1, values[11][simu], 1);
			MatrixGPU PLb(nAgent[i], 1, values[12][simu], 1);
			MatrixGPU PUb(nAgent[i], 1, values[13][simu], 1);

			float rho = rhos[simu];
			float at1 = 2 * rho;
			float at2 = 3 * rho;


			time = 0;
			for (int repet = 0; repet < nRepet; repet++) {
				int _n = nAgent[i];
				MatrixGPU Tlocal_preCopy(Tlocal_pre);
				MatrixGPU Bt1Copy(Bt1);
				MatrixGPU CtCopy(Ct);
				MatrixGPU LbCopy(Lb);
				MatrixGPU UbCopy(Ub);

				MatrixGPU TmoyCopy(Tmoy);
				MatrixGPU PCopy(P);
				MatrixGPU MUCopy(MU);
				MatrixGPU nVoisinCopy(nVoisin);
				MatrixGPU Ap1Copy(Ap1);
				MatrixGPU Ap12Copy(Ap12);
				MatrixGPU CpCopy(Cp);
				MatrixGPU PLbCopy(PLb);
				MatrixGPU PUbCopy(PUb);

				MatrixGPU TlocalCopy(Tlocal);
				MatrixGPU CoresLinVoisinCopy(CoresLinVoisin);
				MatrixGPU CoresAgentLinCopy(CoresAgentLin);

				cudaDeviceSynchronize();
				a = std::chrono::high_resolution_clock::now();
				updateTradePGPUShared<256> << <_n, blockSize >> > (TlocalCopy._matrixGPU, Tlocal_preCopy._matrixGPU, TmoyCopy._matrixGPU, PCopy._matrixGPU,
					MUCopy._matrixGPU, nVoisinCopy._matrixGPU, at1, at2, Bt1Copy._matrixGPU, CtCopy._matrixGPU, LbCopy._matrixGPU, UbCopy._matrixGPU,
					Ap1Copy._matrixGPU, Ap12Copy._matrixGPU, CpCopy._matrixGPU, PLbCopy._matrixGPU, PUbCopy._matrixGPU, CoresAgentLinCopy._matrixGPU);
				cudaDeviceSynchronize();
				b = std::chrono::high_resolution_clock::now();
				time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
			}
			temps.set(simu, i, (float)time / nRepet);
		}
	}

	nAgentMat.saveCSV(fileName, mode);
	nTradeMat.saveCSV(fileName, mode);
	nPro.saveCSV(fileName, mode);
	temps.saveCSV(fileName, mode);
}



void testADMMGPUConst3TimeUpdateRes() {
	std::string fileName = "TempsResG.csv";
	
	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	unsigned int time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNAgent = 7;
	const int nSimu = 100;
	const int nRepet = 10;
	int nAgent[nNAgent] = { 10, 100, 500, 1000, 5000, 10000, 40000 }; // autant conso que de prod, la dernière veleur ne "marche" pas (trop rapide)
	int ntrade[nNAgent];
	int blockSize = 256;
	float values1[nSimu];
	float values2[nSimu];
	
	MatrixCPU temps(nNAgent, nSimu, 0);
	MatrixCPU res(2, 1, 0);

	for (int j = 0; j < nSimu;j++) {
		values1[j] = (float)(rand()) / rand();
		values2[j] = (float)(rand()) / rand();
	}

	for (int i = 0; i < nNAgent; i++) {

		ntrade[i] = nAgent[i] * nAgent[i] / 2;
		std::cout << "iteration " << i << " nAgent " << nAgent[i] << " ntrade " << ntrade[i] << std::endl;
		int numBlocks = ceil((ntrade[i] + blockSize - 1) / blockSize);

		MatrixGPU tempN(blockSize, 1, 0, 1);
		MatrixGPU tempN2(blockSize, 1, 0, 1);

		MatrixGPU CoresLinAgent(ntrade[i], 1);
		MatrixGPU CoresLinVoisin(ntrade[i], 1);
		MatrixGPU CoresMatLin(nAgent[i], nAgent[i], -1);
		MatrixGPU CoresLinTrans(ntrade[i], 1);


		int indice = 0;
		int idVoisin = 0;

		for (int idAgent = 0; idAgent < nAgent[i]; idAgent++) {
			int Nvoisinmax = nAgent[i] / 2;
			if (idAgent < Nvoisinmax) {
				idVoisin = 0;
			}
			else {
				idVoisin = nAgent[i] / 2;
			}
			for (int voisin = idVoisin; voisin < Nvoisinmax; voisin++) {
				CoresLinAgent.set(indice, 0, idAgent);
				CoresLinVoisin.set(indice, 0, voisin);
				CoresMatLin.set(idAgent, voisin, indice);
				indice = indice + 1;
			}
		}
		for (int lin = 0;lin < ntrade[i];lin++) {
			int i = CoresLinAgent.get(lin, 0);
			int j = CoresLinVoisin.get(lin, 0);
			int k = CoresMatLin.get(j, i);
			CoresLinTrans.set(lin, 0, k);
		}
		CoresLinAgent.transferGPU();
		CoresLinVoisin.transferGPU();
		CoresMatLin.transferGPU();
		CoresLinTrans.transferGPU();
		clock_t t = clock();

		for (int simu = 0; simu < nSimu; simu++) {
			MatrixGPU Tlocal_pre(ntrade[i], 1, values1[simu], 1);
			MatrixGPU Tlocal(ntrade[i], 1, values2[simu], 1);
			
			cudaDeviceSynchronize();
			time = 0;
			for (int repet = 0; repet < nRepet; repet++) {
				MatrixGPU Tlocal_preCopy(Tlocal_pre);
				MatrixGPU TlocalCopy(Tlocal);
				MatrixGPU CoresLinTransCopy(CoresLinTrans);
				MatrixGPU tempNCopy(tempN);
				MatrixGPU tempN2Copy(tempN2);
				cudaDeviceSynchronize();
				a = std::chrono::high_resolution_clock::now();
				float resS = TlocalCopy.max2(&Tlocal_preCopy);
				updateDiffGPU << <numBlocks, blockSize >> > (tempNCopy._matrixGPU, TlocalCopy._matrixGPU, CoresLinTransCopy._matrixGPU, ntrade[i]);
				float resR = tempNCopy.max2();
				cudaDeviceSynchronize();
				b = std::chrono::high_resolution_clock::now();
				time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
			}
			temps.set(i, simu, (float)time / nRepet);
		}
	}
	temps.saveCSV(fileName, mode);

}

void testADMMGPUConst3TimeCalcRes() {
	std::string fileName = "TempsResL.csv";

	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	unsigned int time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNAgent = 7;
	const int nSimu = 100;
	const int nRepet = 10;
	int nAgent[nNAgent] = { 10, 100, 500, 1000, 5000, 10000, 40000 }; // autant conso que de prod, la dernière veleur ne "marche" pas (trop rapide)
	int ntrade[nNAgent];
	int blockSize = 256;
	float values1[nSimu];
	float values2[nSimu];

	MatrixCPU temps(nNAgent, nSimu, 0);
	
	for (int j = 0; j < nSimu;j++) {
		values1[j] = (float)(rand()) / rand();
		values2[j] = (float)(rand()) / rand();
	}

	for (int i = 0; i < nNAgent; i++) {

		ntrade[i] = nAgent[i] * nAgent[i] / 2;
		std::cout << "iteration " << i << " nAgent " << nAgent[i] << " ntrade " << ntrade[i] << std::endl;
		int numBlocks = ceil((ntrade[i] + blockSize - 1) / blockSize);

		MatrixGPU tempN(blockSize, 1, 0, 1);
		MatrixGPU tempN2(blockSize, 1, 0, 1);

		clock_t t = clock();

		for (int simu = 0; simu < nSimu; simu++) {
			MatrixGPU Tlocal_pre(ntrade[i], 1, values1[simu], 1);
			MatrixGPU Tlocal(ntrade[i], 1, values2[simu], 1);
			MatrixGPU P(nAgent[i], 1, values1[simu], 1);
			MatrixGPU Tmoy(nAgent[i], 1, values2[simu], 1);

			cudaDeviceSynchronize();
			time = 0;
			for (int repet = 0; repet < nRepet; repet++) {
				MatrixGPU Tlocal_preCopy(Tlocal_pre);
				MatrixGPU TlocalCopy(Tlocal);
				MatrixGPU PCopy(P);
				MatrixGPU TmoyCopy(Tmoy);
				MatrixGPU tempNCopy(tempN);
				MatrixGPU tempN2Copy(tempN2);
				cudaDeviceSynchronize();
				a = std::chrono::high_resolution_clock::now();
				
				float d11 = Tlocal.max2(&Tlocal_pre);
				float d22 = P.max2(&Tmoy);
			
				cudaDeviceSynchronize();
				b = std::chrono::high_resolution_clock::now();
				time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
			}


			/*cudaEventCreate(&stop);
			cudaEventRecord(stop, 0);
			cudaEventSynchronize(stop);
			cudaEventElapsedTime(&elapsedTime, start, stop);*/
			temps.set(i, simu, (float)time / nRepet);
		}
	}
	temps.saveCSV(fileName, mode);

}



