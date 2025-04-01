#include "../head/TestADMMGPUConst1T.cuh"
#define NSTEPLOCAL 5
#define NMAXPEERPERTRHREAD 5

int testADMMGPUConst1T()
{
	int n = 1;

	if (!testADMMGPUConst1TContruct1()) return n;
	n++;
	if (!testADMMGPUConst1TContruct2()) return n;
	n++;
	if (!testADMMGPUConst1TContruct3()) return n;
	n++;
	if (!testADMMGPUConst1TLAMBDA()) return n;
	n++;
	if (!testADMMGPUConst1TKappa()) return n;
	n++;
	if (!testADMMGPUConst1TBt1()) return n;
	n++;
	if (!testADMMGPUConst1TCP()) return n;
	n++; 
	if (!testADMMGPUConst1TCpb()) return n;
	n++;
	if (!testADMMGPUConst1TTradeP()) return n;
	n++; //10
	if (!testADMMGPUConst1TQ()) return n;
	n++; 
	if (!testADMMGPUConst1Talpha()) return n;
	n++;
	//std::cout << n << std::endl;
	if (!testADMMGPUConst1TUpdateRes()) return n;
	n++;
	if (!testADMMGPUConst1TCalcRes()) return n;
	n++;
	if (!testADMMGPUConst1TSolve1()) return n;
	n++; 
	if (!testADMMGPUConst1TSolve2()) return n;
	n++;
	if (!testADMMGPUConst1TSolve3()) return n;
	n++;
	return 0;
}

bool testADMMGPUConst1TContruct1()
{
	std::cout << "contructeur par defaut" << std::endl;
	ADMMGPUConst1T a;
	return true;
}

bool testADMMGPUConst1TContruct2()
{
	float rho = 2;

	std::cout << "contructeur avec parametres" << std::endl;
	ADMMGPUConst1T a(rho);
	return true;
}
bool testADMMGPUConst1TContruct3()
{
	float rho = 2;

	std::cout << "contructeur en deux temps" << std::endl;
	ADMMGPUConst1T a;
	a = ADMMGPUConst1T(rho);
	return true;
}

bool testADMMGPUConst1TSolve1()
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

	ADMMGPUConst1T a;

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
bool testADMMGPUConst1TSolve2()
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
	ADMMGPUConst1T a;
	a.solve(&res, param, cas);
	res.display();
	MatrixCPU Trade = res.getTrade();
	MatrixCPU P2(31, 1, 0);
	P2.sum(&Trade);

	float Pn[31] = { -1.008853555,-4.62966156,-2.927534103,-0.8979898691,-0.9462603927,-0.09805059433,-0.127968356,-4.168303013,-3.151874542,-2.261414766,-0.670329392,-3.399893284,-0.4841034412,-2.775528431,-3.008597374,-1.849177122,-0.5534118414,-2.362840891,-1.122991204,-0.1379692554,-2.332088947,4.406820297,5.406073093,3.676487684,3.929354668,4.570535183,2.529039145,3.478654861,2.755935192,3.768760443,4.393183708, };

	MatrixCPU P(31, 1);
	for (int i = 0; i < 31; i++) {
		P.set(i, 0, Pn[i]);
	}
	MatrixCPU P22 = res.getPn(); 
	

	return (P2.isEqual(&P, 0.01) && P2.isEqual(&P22, 0.01));

}

bool testADMMGPUConst1TSolve3()
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

	ADMMGPUConst1T a;

	MatrixCPU Trade(nAgent, nAgent);
	Trade.set(0, 1, -value);
	Trade.set(1, 0, value);
	a.solve(&res, param, cas);

	MatrixCPU trade = res.getTrade();
	res.display();
	trade.display();
	return trade.isEqual(&Trade, 0.001);
}

bool testADMMGPUConst1TLAMBDA()
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

	return (LAMBDALin.isEqual(&LAMBDALin2));
}

bool testADMMGPUConst1TKappa()
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

bool testADMMGPUConst1TBt1()
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

bool testADMMGPUConst1TTradeP()
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

bool testADMMGPUConst1Talpha()
{
	int _nAgent = 2;
	int _nLine = 3;

	float value1 = 2;
	float value2 = 3;
	
	int _blockSize = 256;
	int _numBlocksNL = ceil((_nLine*_nAgent + _blockSize - 1) / _blockSize);


	MatrixCPU Pn(_nAgent, 1, value1);
	MatrixCPU G(_nAgent, _nLine, value2);
	MatrixCPU alpha(_nAgent, _nLine);
	
	alpha.multiplyTVector(&G, &Pn, 1);
	
	MatrixGPU PnGPU(_nAgent, 1, value1, 1);
	MatrixGPU GGPU(_nAgent, _nLine, value2, 1);
	MatrixGPU alphaGPU(_nAgent, _nLine, 0, 1);

	updateAlphaTrans <<< _numBlocksNL, _blockSize >> > (alphaGPU._matrixGPU, GGPU._matrixGPU, PnGPU._matrixGPU, _nLine, _nAgent);


	MatrixCPU alphaResult(_nAgent, _nLine);
	
	alphaGPU.toMatCPU(alphaResult);
	
	std::cout << "------------------------------------ " << std::endl;
	Pn.display();
	G.display();
	alpha.display();
	

	std::cout << "------------------------------------ " << std::endl;
	alphaResult.display();


	return alpha.isEqual(&alphaResult);
}

bool testADMMGPUConst1TQ()
{
	int _nAgent = 10;
	int _nLine = 7;

	float value1 = 2;

	int _blockSize = 256;
	int _numBlocksL = ceil((_nLine + _blockSize - 1) / _blockSize);

	MatrixCPU alpha(_nAgent, _nLine, value1);
	MatrixCPU Qpart(_nAgent, _nLine);
	MatrixCPU Qtot(_nLine, 1);
	

	for (int l = 0; l < _nLine; l++) {
		float qt = 0;
		for (int n = _nAgent - 1; n >= 0; n--) {
			qt += alpha.get(n, l);
			if (n > 0) {
				Qpart.set(n - 1, l , qt);
			}
		}
		Qtot.set(l, 0, qt);
	}

	MatrixGPU alphaGPU(_nAgent, _nLine, value1, 1);
	MatrixGPU QpartGPU(_nAgent, _nLine, 0, 1);
	MatrixGPU QtotGPU(_nLine, 1, 0, 1);


	updateQpartTrans << < _nLine, _blockSize, _nAgent * sizeof(float) >> > (QpartGPU._matrixGPU, alphaGPU._matrixGPU, _nAgent, _nLine);
	updateQtotTrans << <_numBlocksL, _blockSize >> > (QtotGPU._matrixGPU, QpartGPU._matrixGPU, alphaGPU._matrixGPU, _nLine);


	MatrixCPU QpartResult(_nAgent, _nLine);
	MatrixCPU QtotResult(_nLine, 1);

	QtotGPU.toMatCPU(QtotResult);
	QpartGPU.toMatCPU(QpartResult);
	

	/*std::cout << "------------------------------------ " << std::endl;
	alpha.display();
	Qpart.display();
	Qtot.display();

	std::cout << "------------------------------------ " << std::endl;
	QpartResult.display();
	QtotResult.display();*/

	

	return (Qtot.isEqual(&QtotResult)) && (Qpart.isEqual(&QpartResult));
}

bool testADMMGPUConst1TCP()
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
	MatrixCPU G(_nAgent, _nLine, value5);
	MatrixCPU Qpart(_nAgent, _nLine, value6);
	MatrixCPU nVoisin(_nAgent, 1, 1);


	tempL1.subtractAbs(&Kappa1, &Kappa2);
	float r = 0;
	for (int i = 0; i < _nAgent; ++i)
	{
		r = 0;
		for (int k = 0; k < _nLine; ++k)
		{
			r += G.get(i, k) * (tempL1.get(k, 0) + 2 * Qpart.get(i, k));
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
	MatrixGPU GGPU(_nAgent, _nLine, value5, 1);
	MatrixGPU QpartGPU(_nAgent, _nLine, value6, 1);
	MatrixCPU CpResult(_nAgent, 1);
	MatrixGPU nVoisinGPU(_nAgent, 1, 1, 1);

	diffKappa << <_numBlocksL, _blockSize >> > (tempL1GPU._matrixGPU, Kappa1GPU._matrixGPU, Kappa2GPU._matrixGPU, _nLine);
	updateCp2aTrans<256> << <numBlocks, _blockSize >> > (Cp2GPU._matrixGPU, tempL1GPU._matrixGPU, GGPU._matrixGPU, _nLine, _nAgent);
	updateCp2bTrans<256> << <numBlocks, _blockSize >> > (tempN1GPU._matrixGPU, GGPU._matrixGPU, QpartGPU._matrixGPU, _nLine, _nAgent);
	updateCpOld << <_numBlocksN, _blockSize >> > (CpGPU._matrixGPU, Cp1GPU._matrixGPU, Cp2GPU._matrixGPU, tempN1GPU._matrixGPU, nVoisinGPU._matrixGPU, _rho1, _nAgent);

	CpGPU.toMatCPU(CpResult);

	Cp.display();
	CpResult.display();

	return Cp.isEqual(&CpResult);
}

bool testADMMGPUConst1TCpb()
{
	int _nAgent = 100;
	int _nLine = 200;

	float value1 = 2;
	float value2 = 3;

	int _blockSize = 256;
	int numBlocks = _nAgent;

	
	MatrixCPU Qpart(_nAgent, _nLine, value1);
	MatrixCPU G(_nAgent, _nLine, value2);
	MatrixCPU Cpb(_nAgent, 1);

	for (int n = 0; n < _nAgent; n++) {
		float sum = 0;
		for (int l = 0; l < _nLine; l++) {
			sum += G.get(n, l) * Qpart.get(n, l);
		}
		Cpb.set(n, 0, 2*sum);
	}
	

	MatrixGPU GGPU(_nAgent, _nLine, value2, 1);
	MatrixGPU QpartGPU(_nAgent, _nLine, value1, 1);
	MatrixGPU CpbGPU(_nAgent, 1, 0, 1);

	updateCp2bTrans<256> << <numBlocks, _blockSize >> > (CpbGPU._matrixGPU, GGPU._matrixGPU, QpartGPU._matrixGPU, _nLine, _nAgent);
	
	MatrixCPU CpbResult(_nAgent, 1);
	
	CpbGPU.toMatCPU(CpbResult);
	
	return Cpb.isEqual(&CpbResult) ;

}

bool testADMMGPUConst1TUpdateRes()
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
bool testADMMGPUConst1TCalcRes() {
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
