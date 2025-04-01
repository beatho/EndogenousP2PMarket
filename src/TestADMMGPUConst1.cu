#include "../head/TestADMMGPUConst1.cuh"
#define NSTEPLOCAL 5
#define NMAXPEERPERTRHREAD 5

int testADMMGPUConst1()
{
	int n = 1;

	if (!testADMMGPUConst1Contruct1()) return n;
	n++;
	
	if (!testADMMGPUConst1Contruct2()) return n;
	n++;
	
	if (!testADMMGPUConst1Contruct3()) return n;
	n++;
	
	if (!testADMMGPUConst1LAMBDA()) return n;
	n++;
	if (!testADMMGPUConst1Kappa()) return n;
	n++;
	if (!testADMMGPUConst1Bt1()) return n;
	n++;
	if (!testADMMGPUConst1CP()) return n;
	n++; 
	if (!testADMMGPUConstCpb()) return n;
	n++;
	if (!testADMMGPUConst1TradeP()) return n;
	n++; //10
	if (!testADMMGPUConst1Q()) return n;
	n++; 
	if (!testADMMGPUConst1alpha()) return n;
	n++;
	//std::cout << n << std::endl;
	if (!testADMMGPUConst1UpdateRes()) return n;
	n++;
	if (!testADMMGPUConst1CalcRes()) return n;
	n++;
	if (!testADMMGPUConst1Solve1()) return n;
	n++; 
	if (!testADMMGPUConst1Solve2()) return n;
	n++;
	if (!testADMMGPUConst1Solve3()) return n;
	n++;
	return 0;
}

void testADMMGPUConst1Time(int test) {

	switch (test)
	{
	case 0:
		testADMMGPUConst1TimeLAMBDA();
		break;
	case 1:
		testADMMGPUConst1TimeBt1();
		break;
	case 2:
		testADMMGPUConst1TimeTradeP();
		break;
	case 3:
		testADMMGPUConst1TimeUpdateRes();
		break;
	case 4:
		testADMMGPUConst1TimeCalcRes();
		break;
	default:
		std::cout << "No valid input " << std::endl;
		break;
	}
	
}

bool testADMMGPUConst1Contruct1()
{
	std::cout << "contructeur par defaut" << std::endl;
	ADMMGPUConst1 a;
	return true;
}

bool testADMMGPUConst1Contruct2()
{
	float rho = 2;

	std::cout << "contructeur avec parametres" << std::endl;
	ADMMGPUConst1 a(rho);
	return true;
}
bool testADMMGPUConst1Contruct3()
{
	float rho = 2;

	std::cout << "contructeur en deux temps" << std::endl;
	ADMMGPUConst1 a;
	a = ADMMGPUConst1(rho);
	return true;
}

bool testADMMGPUConst1Solve1()
{
	//solve(Simparam* result, Simparam sim, StudyCase cas);
	std::cout << "-------------------------------------------------------- " << std::endl;
	StudyCase cas;
	cas.Set2node();
	//cas.display();
	int nAgent = cas.getNagent();
	Simparam param(nAgent, 1);
	param.setRho(1);
	Simparam res(param);

	ADMMGPUConst1 a;

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
bool testADMMGPUConst1Solve2()
{
	//solve(Simparam* result, Simparam sim, StudyCase cas);
	std::cout << "-------------------------------------------------------- " << std::endl;
	StudyCase cas;
	cas.Set29node();
	//cas.display();
	int nAgent = cas.getNagent();

	Simparam param(nAgent, cas.getNLine());
	float epsG = 0.00002f;
	float epsL = 0.000002f;
	param.setEpsL(epsL);
	param.setEpsG(epsG);
	param.setStep(1, 1);
	param.setRho(10000);
	Simparam res(param);
	ADMMGPUConst1 a;
	a.solve(&res, param, cas);
	res.display();
	MatrixCPU Trade = res.getTrade();
	MatrixCPU P2(29, 1, 0);
	P2.sum(&Trade);

	float Pn[31] = { -1.008853555,-4.62966156,-2.927534103,-0.8979898691,-0.9462603927,-0.09805059433,-0.127968356,-4.168303013,-3.151874542,-2.261414766,-0.670329392,-3.399893284,-0.4841034412,-2.775528431,-3.008597374,-1.849177122,-0.5534118414,-2.362840891,-1.122991204,-0.1379692554,-2.332088947,4.406820297,5.406073093,3.676487684,3.929354668,4.570535183,2.529039145,3.478654861,2.755935192,3.768760443,4.393183708, };

	MatrixCPU P(31, 1);
	for (int i = 0; i < 31; i++) {
		P.set(i, 0, Pn[i]);
	}
	MatrixCPU P22 = res.getPn(); 
	

	return (P2.isEqual(&P, 0.01) && P2.isEqual(&P22, 0.01));

}

bool testADMMGPUConst1Solve3()
{
	std::cout << "-------------------------------------------------------- " << std::endl;
	StudyCase cas;
	float lim = 0.8;
	cas.Set2nodeConstraint(lim);
	int nAgent = cas.getNagent();
	Simparam param(nAgent, 1);
	Simparam res(param);
	
	param.setRho1(50);
	float value = (1 - lim) * (lim > 1) + lim;

	ADMMGPUConst1 a;

	MatrixCPU Trade(nAgent, nAgent);
	Trade.set(0, 1, -value);
	Trade.set(1, 0, value);
	a.solve(&res, param, cas);

	MatrixCPU trade = res.getTrade();
	res.display();
	trade.display();
	return trade.isEqual(&Trade, 0.001);
}

bool testADMMGPUConst1LAMBDA()
{
	int nAgent = 3; // 2 conso et un prod
	int ntrade = 4;
	int blockSize = 256;
	int numBlocks = ceil((nAgent + blockSize - 1) / blockSize);
	float value1 = 2;
	float value2 = -8;
	float value3 = 1.5;
	float value4 = 4;
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

	//LAMBDALin.display();
	//std::cout << "--------------" << std::endl;
	//LAMBDALin2.display();


	return (LAMBDALin.isEqual(&LAMBDALin2));
}

bool testADMMGPUConst1Kappa()
{
	
	int _nLine = 3;
	
	float value1 = 2;
	float value2 = 3;
	float value3 = 1;
	float value4 = -2;
	
	int _blockSize = 256;
	int _numBlocksL = ceil((_nLine + _blockSize - 1) / _blockSize);

	
	
	
	MatrixCPU Qtot(_nLine, 1, value1);
	MatrixCPU Llimit(_nLine, 1, value2);
	MatrixCPU Kappa1(_nLine, 1, value3);
	MatrixCPU Kappa2(_nLine, 1, value4);
	
	Kappa1.projectNeg();
	Kappa1.add(&Llimit);
	Kappa1.subtract(&Qtot);
	
	Kappa2.projectNeg();
	Kappa2.add(&Llimit);
	Kappa2.add(&Qtot);


	MatrixGPU QtotGPU(_nLine, 1, value1, 1);
	MatrixGPU LlimitGPU(_nLine, 1, value2, 1);
	MatrixGPU Kappa1GPU(_nLine, 1, value3, 1);
	MatrixGPU Kappa2GPU(_nLine, 1, value4, 1);

	MatrixCPU Kappa1Result(_nLine, 1);
	MatrixCPU Kappa2Result(_nLine, 1);

	updateKappaGPU << <_numBlocksL, _blockSize >> > (Kappa1GPU._matrixGPU, Kappa2GPU._matrixGPU, LlimitGPU._matrixGPU, QtotGPU._matrixGPU, _nLine);

	Kappa1GPU.toMatCPU(Kappa1Result);
	Kappa2GPU.toMatCPU(Kappa2Result);

	

	return (Kappa1.isEqual(&Kappa1Result) && Kappa2.isEqual(&Kappa2Result));
}

bool testADMMGPUConst1Bt1()
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

bool testADMMGPUConst1TradeP()
{
	int nAgent = 4;
	int ntrade = 8;
	int blockSize = 256;
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


	updateTradePGPU<256> <<<nAgent, blockSize >> > (Tlocal._matrixGPU, Tlocal_pre._matrixGPU, Tmoy._matrixGPU, P._matrixGPU, MU._matrixGPU, nVoisin._matrixGPU, at1, at2, Bt1._matrixGPU, Ct._matrixGPU,
		Lb._matrixGPU, Ub._matrixGPU, Ap1._matrixGPU, Ap12._matrixGPU, Cp._matrixGPU, PLb._matrixGPU, PUb._matrixGPU, CoresAgentLin._matrixGPU, CoresLinAgent._matrixGPU, nAgent);

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

bool testADMMGPUConst1alpha()
{
	int _nAgent = 2;
	int _nLine = 3;

	float value1 = 2;
	float value2 = 3;
	
	int _blockSize = 256;
	int _numBlocksNL = ceil((_nLine*_nAgent + _blockSize - 1) / _blockSize);


	MatrixCPU Pn(_nAgent, 1, value1);
	MatrixCPU G(_nLine, _nAgent, value2);
	MatrixCPU alpha(_nLine, _nAgent);
	
	alpha.multiplyTVector(&G, &Pn, 0);
	
	MatrixGPU PnGPU(_nAgent, 1, value1, 1);
	MatrixGPU GGPU(_nLine, _nAgent, value2, 1);
	MatrixGPU alphaGPU(_nLine, _nAgent, 0, 1);

	updateAlpha << < _numBlocksNL, _blockSize >> > (alphaGPU._matrixGPU, GGPU._matrixGPU, PnGPU._matrixGPU, _nLine, _nAgent);

	MatrixCPU alphaResult(_nLine, _nAgent);
	
	alphaGPU.toMatCPU(alphaResult);
	


	return alpha.isEqual(&alphaResult);
}

bool testADMMGPUConst1Q()
{
	int _nAgent = 10;
	int _nLine = 7;

	float value1 = 2;
	
	int _blockSize = 256;
	int _numBlocksL = ceil((_nLine + _blockSize - 1) / _blockSize);

	MatrixCPU alpha(_nLine, _nAgent, value1);
	MatrixCPU Qpart(_nLine, _nAgent);
	MatrixCPU Qtot(_nLine, 1);
	

	for (int l = 0; l < _nLine; l++) {
		float qt = 0;
		for (int n = _nAgent - 1; n >= 0; n--) {
			qt += alpha.get(l, n);
			if (n > 0) {
				Qpart.set(l, n - 1, qt);
			}
		}
		Qtot.set(l, 0, qt);
	}

	MatrixGPU alphaGPU(_nLine, _nAgent, value1, 1);
	MatrixGPU QpartGPU(_nLine, _nAgent, 0, 1);
	MatrixGPU alphaGPUTrans(_nAgent, _nLine, value1, 1);
	MatrixGPU QpartGPUTrans(_nAgent, _nLine, 0, 1);
	MatrixGPU QtotGPU(_nLine, 1, 0, 1);
	MatrixGPU QtotGPU1(_nLine, 1, 0, 1);
	MatrixGPU QtotGPU2(_nLine, 1, 0, 1);

	
	updateQpart <<< _nLine, _blockSize, _nAgent * sizeof(float) >> > (QpartGPU._matrixGPU, alphaGPU._matrixGPU, _nAgent);
	calculQpartAgentBlocTrans <<< _nAgent, _blockSize >>> (QpartGPUTrans._matrixGPU, alphaGPUTrans._matrixGPU, _nLine, _nAgent);

	QtotGPU.sum(&alphaGPU);
	updateQtotTest <<< _numBlocksL, _blockSize >> > (QtotGPU1._matrixGPU, QpartGPU._matrixGPU, alphaGPU._matrixGPU, _nLine, _nAgent);
	updateQtotTestTrans << < _numBlocksL, _blockSize >> > (QtotGPU2._matrixGPU, QpartGPUTrans._matrixGPU, alphaGPU._matrixGPU, _nLine);

	MatrixCPU QpartResult(_nLine, _nAgent);
	MatrixCPU QpartResultTrans(_nAgent, _nLine);
	MatrixCPU QtotResult(_nLine, 1);
	MatrixCPU QtotResult1(_nLine, 1);
	MatrixCPU QtotResult2(_nLine, 1);
	
	QtotGPU.toMatCPU(QtotResult);
	QtotGPU1.toMatCPU(QtotResult1);
	QtotGPU2.toMatCPU(QtotResult2);
	QpartGPU.toMatCPU(QpartResult);
	QpartGPUTrans.toMatCPU(QpartResultTrans);

	std::cout << "-------------------------------------------------------- " << std::endl;
	alpha.display();
	QpartResult.display();
	QpartResultTrans.display();

	std::cout << "-------------------------------------------------------- " << std::endl;
	QtotResult.display();
	QtotResult1.display();
	QtotResult2.display();

	std::cout << (Qtot.isEqual(&QtotResult)) <<" " << (Qtot.isEqual(&QtotResult1)) << " " << (Qtot.isEqual(&QtotResult2)) << " " << (Qpart.isEqual(&QpartResult)) << std::endl;

	return (Qtot.isEqual(&QtotResult)) && (Qtot.isEqual(&QtotResult1)) && (Qtot.isEqual(&QtotResult2)) && (Qpart.isEqual(&QpartResult));
}

bool testADMMGPUConst1CP()
{
	int _nAgent = 2;
	int _nLine = 3;
	float _rho1 = 1.5;
	float value1 = 2;
	float value2 = 3;
	float value3 = 1;
	float value4 = -2;
	float value5 = -1;
	float value6 = -30;
	//float value7 = 10;
	//float value8 = 5;
	int _blockSize = 256;
	int numBlocks = _nAgent;
	int _numBlocksN = ceil((_nAgent + _blockSize - 1) / _blockSize);
	int _numBlocksL = ceil((_nLine + _blockSize - 1) / _blockSize);
	
	MatrixCPU Cp(_nAgent, 1);
	MatrixCPU tempN1(_nAgent, 1);
	MatrixCPU Cp1(_nAgent, 1, value1);
	MatrixCPU Cp2(_nAgent, 1, value2);
	MatrixCPU tempL1(_nLine, 1);
	MatrixCPU Kappa1(_nLine, 1, value3);
	MatrixCPU Kappa2(_nLine, 1, value4);
	MatrixCPU G(_nLine, _nAgent, value5);
	MatrixCPU Qpart(_nLine, _nAgent, value6);
	MatrixCPU nVoisin(_nAgent, 1, 1);


	ADMMGPUConst1 a;
	tempL1.subtractAbs(&Kappa1, &Kappa2);
	//Cp2->multiplyTrans(G, tempL1, 0);

	float r = 0;
	for (int i = 0; i < _nAgent; i++)
	{
		r = 0;
		for (int k = 0; k < _nLine; ++k)
		{
			r +=  G.get(k, i) * (tempL1.get(k, 0) + 2 * Qpart.get(k, i));
		}
		Cp2.set(i, 0, r);
	}

	Cp2.multiply(_rho1);
	Cp2.multiplyT(&nVoisin);

	Cp.add(&Cp1, &Cp2);
	

	MatrixGPU CpGPU(_nAgent, 1, 0, 1);
	MatrixGPU tempN1GPU(_nAgent, 1, 0, 1);
	MatrixGPU Cp1GPU(_nAgent, 1, value1, 1);
	MatrixGPU Cp2GPU(_nAgent, 1, value2, 1);
	MatrixGPU tempL1GPU(_nLine, 1, 0, 1);
	MatrixGPU Kappa1GPU(_nLine, 1, value3, 1);
	MatrixGPU Kappa2GPU(_nLine, 1, value4, 1);
	MatrixGPU GGPU(_nLine, _nAgent, value5, 1);
	MatrixGPU QpartGPU(_nLine, _nAgent, value6, 1);
	MatrixCPU CpResult(_nAgent, 1);
	MatrixGPU nVoisinGPU(_nAgent, 1, 1, 1);

	diffKappa << <_numBlocksL, _blockSize >> > (tempL1GPU._matrixGPU, Kappa1GPU._matrixGPU, Kappa2GPU._matrixGPU, _nLine);
	updateCp2a<256> << <numBlocks, _blockSize >> > (Cp2GPU._matrixGPU, tempL1GPU._matrixGPU, GGPU._matrixGPU, _nLine, _nAgent);
	updateCp2b<256> << <numBlocks, _blockSize >> > (tempN1GPU._matrixGPU, GGPU._matrixGPU, QpartGPU._matrixGPU, _nLine, _nAgent);
	updateCpOld << <_numBlocksN, _blockSize >> > (CpGPU._matrixGPU, Cp1GPU._matrixGPU, Cp2GPU._matrixGPU, tempN1GPU._matrixGPU, nVoisinGPU._matrixGPU, _rho1, _nAgent);

	CpGPU.toMatCPU(CpResult);

	Cp.display();
	CpResult.display();

	return Cp.isEqual(&CpResult);
}

bool testADMMGPUConstCpb()
{
	int _nAgent = 100;
	int _nLine = 200;

	float value1 = 2;
	float value2 = 3;
	int _blockSize = 256;
	int numBlocks = _nAgent;

	
	MatrixCPU Qpart(_nLine, _nAgent, value1);
	MatrixCPU G(_nLine, _nAgent, value2);
	MatrixCPU Cpb(_nAgent, 1);

	for (int n = 0; n < _nAgent; n++) {
		float sum = 0;
		for (int l = 0; l < _nLine; l++) {
			sum += G.get(l, n) * Qpart.get(l, n);
		}
		Cpb.set(n, 0, 2*sum);
	}
	

	MatrixGPU GGPU(_nLine, _nAgent, value2, 1);
	MatrixGPU GGPUTrans(_nAgent, _nLine, value2, 1);
	MatrixGPU QpartGPU(_nLine, _nAgent, value1, 1);
	MatrixGPU QpartGPUTrans(_nAgent, _nLine, value1, 1);

	MatrixGPU CpbGPU(_nAgent, 1, 0, 1);
	MatrixGPU CpbGPU1(_nAgent, 1, 0, 1);



	
	updateCp2bTest<256> << <numBlocks, _blockSize >> > (CpbGPU._matrixGPU, GGPU._matrixGPU, QpartGPU._matrixGPU, _nLine, _nAgent);
	updateCp2bTestTrans<256> << <numBlocks, _blockSize >> > (CpbGPU1._matrixGPU, GGPUTrans._matrixGPU, QpartGPUTrans._matrixGPU, _nLine, _nAgent);

	MatrixCPU CpbResult(_nAgent, 1);
	MatrixCPU CpbResult2(_nAgent, 1);


	CpbGPU.toMatCPU(CpbResult);
	CpbGPU1.toMatCPU(CpbResult2);
	

	
	return (Cpb.isEqual(&CpbResult)) && (Cpb.isEqual(&CpbResult2)) ;

}

bool testADMMGPUConst1UpdateRes()
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
bool testADMMGPUConst1CalcRes() {
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


void testADMMGPUConst1TimeLAMBDA()
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
	int nAgent[nNAgent] = { 10, 100, 500, 1000, 5000, 10000, 40000 }; // autant conso que de prod, la derni�re veleur ne "marche" pas (trop rapide)
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

void testADMMGPUConst1TimeBt1()
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
	int nAgent[nNAgent] = { 10, 100, 500, 1000, 5000, 10000, 40000 }; // autant conso que de prod, la derni�re veleur ne "marche" pas (trop rapide)
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


void testADMMGPUConst1TimeTradeP() {
	std::string fileName = "TempsTradeP8.csv";
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
	MatrixCPU temps(nNAgent, nSimu, 0);

	for (int j = 0; j < nSimu; j++) {
		for (int var = 0; var < nVar; var++) {
			values[var][j] = (float)(rand()) / rand();
		}
		rhos[j] = (float)(rand() % 100) / rand();
	}

	for (int i = 0; i < nNAgent; i++) {

		ntrade[i] = nAgent[i] * nAgent[i] / 2;
		std::cout << "iteration " << i << " nAgent " << nAgent[i] << " ntrade " << ntrade[i] << std::endl;
	
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
			MatrixGPU Ap1(nAgent[i], 1, values[8][simu], 1);
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
				updateTradePGPU<256> << <_n, blockSize >> > (TlocalCopy._matrixGPU, Tlocal_preCopy._matrixGPU, TmoyCopy._matrixGPU, PCopy._matrixGPU,
					MUCopy._matrixGPU, nVoisinCopy._matrixGPU, at1, at2, Bt1Copy._matrixGPU, CtCopy._matrixGPU, LbCopy._matrixGPU, UbCopy._matrixGPU,
					Ap1Copy._matrixGPU, Ap12Copy._matrixGPU, CpCopy._matrixGPU, PLbCopy._matrixGPU, PUbCopy._matrixGPU, CoresAgentLinCopy._matrixGPU, CoresLinVoisinCopy._matrixGPU, _n);
				cudaDeviceSynchronize();
				b = std::chrono::high_resolution_clock::now();
				time += std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
			}
			temps.set(i, simu, (float)time / nRepet);
		}
	}

	temps.saveCSV(fileName, mode);
}



void testADMMGPUConst1TimeUpdateRes() {
	std::string fileName = "TempsResG8.csv";
	
	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	unsigned int time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNAgent = 7;
	const int nSimu = 100;
	const int nRepet = 10;
	int nAgent[nNAgent] = { 10, 100, 500, 1000, 5000, 10000, 40000 }; // autant conso que de prod, la derni�re veleur ne "marche" pas (trop rapide)
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
				updateDiffGPU<< <numBlocks, blockSize >> > (tempNCopy._matrixGPU, TlocalCopy._matrixGPU, CoresLinTransCopy._matrixGPU, ntrade[i]);
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

void testADMMGPUConst1TimeCalcRes() {
	std::string fileName = "TempsResL8.csv";

	std::chrono::high_resolution_clock::time_point a;
	std::chrono::high_resolution_clock::time_point b;
	unsigned int time;
	std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
	const int nNAgent = 7;
	const int nSimu = 100;
	const int nRepet = 10;
	int nAgent[nNAgent] = { 10, 100, 500, 1000, 5000, 10000, 40000 }; // autant conso que de prod, la derni�re veleur ne "marche" pas (trop rapide)
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



