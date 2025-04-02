#include "../head/MethodP2PGPU.cuh"
 


MethodP2PGPU::MethodP2PGPU() : Method(){

}
MethodP2PGPU::~MethodP2PGPU(){

}


float MethodP2PGPU::updateRes(int iter)
{
	float resS = Tlocal.max2(&tradeLin);
	updateDiffGPU << <_numBlocksM, _blockSize >> > (tempNN._matrixGPU, tradeLin._matrixGPU, CoresLinTrans._matrixGPU, _nTrade);
	float resR = tempNN.max2();

	updateResX << <_numBlocksL, _blockSize >> > (tempL1._matrixGPU, Kappa1._matrixGPU, Kappa2._matrixGPU, Kappa1_pre._matrixGPU, Kappa2_pre._matrixGPU, _nLine);

	resF.set(0, iter, resR);
	resF.set(1, iter, resS);
	if (iter > 0 && _tau > 1) {
		if (resR > _mu * resS) {
			_rhog = _tau * _rhog;
			_at1 = _rhog;
			//std::cout << iter << ", rho augmente :" << _rhog << std::endl;
		}
		else if (resS > _mu * resR) {// rho = rho / tau_inc;
			_rhog = _rhog / _tau;
			_at1 = _rhog;
			//std::cout << iter << ", rho diminue :" << _rhog << std::endl;
		}
	}
	
	return MYMAX(resS, resR);
}
float MethodP2PGPU::updateResEndo(int iter)
{
	float resS = Tlocal.max2(&tradeLin);
	updateDiffGPU << <_numBlocksM, _blockSize >> > (tempNN._matrixGPU, tradeLin._matrixGPU, CoresLinTrans._matrixGPU, _nTrade);
	float resR = tempNN.max2();

	updateResX << <_numBlocksL, _blockSize >> > (tempL1._matrixGPU, Kappa1._matrixGPU, Kappa2._matrixGPU, Kappa1_pre._matrixGPU, Kappa2_pre._matrixGPU, _nLine);


	float resXf = _ratioEps * sqrt(tempL1.max2());
	

	resF.set(0, iter, resR);
	resF.set(1, iter, resS);
	resF.set(2, iter, resXf);
	
	return MYMAX(MYMAX(resXf, resS), resR);
}
float MethodP2PGPU::calcRes()
{
	float d1 = Tlocal.max2(&Tlocal_pre);
	float d2 = P.max2(&Tmoy);

	return d1* (d1 > d2) + d2 * (d2 >= d1);
}


void MethodP2PGPU::updateLAMBDA(MatrixGPU* LAMBDA, MatrixGPU* trade, float rho, MatrixGPU* tempNN)
{
	tempNN->set(trade);
	tempNN->addTrans(trade);
	tempNN->multiply(rho);
	tempNN->multiply(0.5);
	LAMBDA->add(LAMBDA, tempNN);
}
void MethodP2PGPU::updateKappa()
{
	Kappa1.projectNeg();
	Kappa1.add(&lLimit);
	Kappa1.subtract(&Qtot);
	Kappa2.projectNeg();
	Kappa2.add(&lLimit);
	Kappa2.add(&Qtot);
}
void MethodP2PGPU::updatePn()
{
	Pn.set(&Tmoy);
	Pn.multiplyT(&nVoisin);
}


void MethodP2PGPU::updateP0(const StudyCase& cas)
{
	_id = _id + 1;
#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
	// Change : Power Limits, cost function


	MatrixGPU Lb(cas.getLb());
	MatrixGPU Ub(cas.getUb());
	MatrixCPU BETA(cas.getBeta());

	matLb.transferCPU();
	matUb.transferCPU();
	Ct.transferCPU();
	CoresLinVoisin.transferCPU();

	if (cas.isAC() && !isAC) {
		MatrixGPU aT(cas.geta(), 1);
		MatrixGPU bT(cas.getb(), 1);
		MatrixGPU PminT(cas.getPmin(), 1);
		MatrixGPU PmaxT(cas.getPmax(), 1);

		a.setFromBloc(0, _nAgent, 0, 1, &aT);
		b.setFromBloc(0, _nAgent, 0, 1, &bT);
		Pmin.setFromBloc(0, _nAgent, 0, 1, &PminT);
		Pmax.setFromBloc(0, _nAgent, 0, 1, &PmaxT);
	}
	else if ((!cas.isAC()) && isAC){
		throw std::invalid_argument("updateP0 : Study Case is not AC, but this method require AC information");
	}
	else {
		a = cas.geta();
		b = cas.getb();
		Pmin = cas.getPmin();
		Pmax = cas.getPmax();
	}
	Cp1 = b;
	int indice = 0;

	// hypothese : ce sont les mêmes voisins !!!
	for (int idAgent = 0; idAgent < _nAgent; idAgent++) {
		int Nvoisinmax = (int) nVoisinCPU.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			int idVoisin = (int) CoresLinVoisin.get(indice, 0);
			if(Lb.getNCol()==1){
				matLb.set(indice, 0, Lb.get(idAgent, 0));
				matUb.set(indice, 0, Ub.get(idAgent, 0));
			} else {
				matLb.set(indice, 0, Lb.get(idAgent, idVoisin));
				matUb.set(indice, 0, Ub.get(idAgent, idVoisin));
			}
			Ct.set(indice, 0, BETA.get(idAgent, idVoisin));
			indice = indice + 1;
		}
	}
	for (int idAgent = _nAgentTrue; idAgent < _nAgent; idAgent++) {
		for (int voisin = 0; voisin < (_nAgent - 1); voisin++) {
			int idVoisin = (int) CoresLinVoisin.get(indice, 0);
			if(Lb.getNCol()==1){
				matLb.set(indice, 0, Lb.get(idAgent, 0));
				matUb.set(indice, 0, Ub.get(idAgent, 0));
			} else {
				matLb.set(indice, 0, Lb.get(idAgent, idVoisin));
				matUb.set(indice, 0, Ub.get(idAgent, idVoisin));
			}
			Ct.set(indice, 0, BETA.get(idAgent, idVoisin));
			indice = indice + 1;
		}
	}


	matLb.transferGPU();
	matUb.transferGPU();
	Ct.transferGPU();
	CoresLinVoisin.transferGPU();

	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);
	Cp1.multiplyT(&nVoisin);

	Ap2a = a;
	Ap2.add(&Ap2a, &Ap2b);
	Ap2.multiplyT(&nVoisin);
	Ap2.multiplyT(&nVoisin);
	Ap12.add(&Ap1, &Ap2);
	Ap123.add(&Ap12, &Ap3);
	Cp.add(&Cp1, &Cp2);

#ifdef INSTRUMENTATION
	cudaDeviceSynchronize();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
	occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION
	
}

void MethodP2PGPU::initLinForm( const Simparam& sim, const StudyCase& cas){

	MatrixCPU BETA(cas.getBeta());
	MatrixCPU Ub(cas.getUb());
	MatrixCPU Lb(cas.getLb());
	LAMBDA = sim.getLambda(); 
	trade = sim.getTrade();

	// Rem : si matrice deja existante, elles sont deja sur GPU donc bug pour les get

	if (Ct.getPos()) { // une copie en trop mais pour l'instant c'est ok...
		CoresMatLin.transferCPU();

		CoresLinAgent.transferCPU();
		CoresAgentLin.transferCPU();
		CoresLinVoisin.transferCPU();
		CoresLinTrans.transferCPU();

		Tlocal_pre.transferCPU();
		tradeLin.transferCPU();
		LAMBDALin.transferCPU();

		matLb.transferCPU();
		matUb.transferCPU();
		Ct.transferCPU();

	}

	CoresMatLin = MatrixGPU(_nAgent, _nAgentTrue, -1);
	CoresAgentLin = MatrixGPU(_nAgent + 1, 1);
	CoresLinAgent = MatrixGPU(_nTrade, 1);
	CoresLinVoisin = MatrixGPU(_nTrade, 1);
	CoresLinTrans = MatrixGPU(_nTrade, 1);

	Tlocal_pre = MatrixGPU(_nTrade, 1);
	tradeLin = MatrixGPU(_nTrade, 1);
	LAMBDALin = MatrixGPU(_nTrade, 1);

	matLb = MatrixGPU(_nTrade, 1);
	matUb = MatrixGPU(_nTrade, 1);
	Ct = MatrixGPU(_nTrade, 1);
	

	int indice = 0;
	//std::cout << " P " << std::endl;
	for (int idAgent = 0; idAgent < _nAgentTrue; idAgent++) { // P
		MatrixCPU omega(cas.getVoisin(idAgent));
		int Nvoisinmax = (int) nVoisinCPU.get(idAgent, 0);
		for (int voisin = 0; voisin < Nvoisinmax; voisin++) {
			int idVoisin = (int) omega.get(voisin, 0);
			if(Lb.getNCol()==1){
				matLb.set(indice, 0, Lb.get(idAgent, 0));
				matUb.set(indice, 0, Ub.get(idAgent, 0));
			} else {
				matLb.set(indice, 0, Lb.get(idAgent, idVoisin));
				matUb.set(indice, 0, Ub.get(idAgent, idVoisin));
			}
			Ct.set(indice, 0, BETA.get(idAgent, idVoisin));
			tradeLin.set(indice, 0, trade.get(idAgent, idVoisin));
			Tlocal_pre.set(indice, 0, trade.get(idAgent, idVoisin));
			LAMBDALin.set(indice, 0, LAMBDA.get(idAgent, idVoisin));
			CoresLinAgent.set(indice, 0, idAgent);
			CoresLinVoisin.set(indice, 0, idVoisin);
			CoresMatLin.set(idAgent, idVoisin, indice);
			indice = indice + 1;
		}
		CoresAgentLin.set(idAgent + 1, 0, indice);
	}
	//std::cout << " Q " << std::endl;
	for (int idAgent = _nAgentTrue; idAgent < _nAgent; idAgent++) { // Q
		for (int idVoisin = 0; idVoisin < _nAgentTrue; idVoisin++) {
			if (idVoisin != (idAgent - _nAgentTrue)) {
				if(Lb.getNCol()==1){
					matLb.set(indice, 0, Lb.get(idAgent, 0));
					matUb.set(indice, 0, Ub.get(idAgent, 0));
				} else {
					matLb.set(indice, 0, Lb.get(idAgent, idVoisin));
					matUb.set(indice, 0, Ub.get(idAgent, idVoisin));
				}
				tradeLin.set(indice, 0, trade.get(idAgent, idVoisin));
				Tlocal_pre.set(indice, 0, trade.get(idAgent, idVoisin));
				LAMBDALin.set(indice, 0, LAMBDA.get(idAgent, idVoisin));
				CoresLinAgent.set(indice, 0, idAgent);
				CoresLinVoisin.set(indice, 0, idVoisin + _nAgentTrue);
				CoresMatLin.set(idAgent, idVoisin, indice);
				indice = indice + 1;
			}
		}

		CoresAgentLin.set(idAgent + 1, 0, indice);
	}
	for (int lin = 0; lin < _nTrade; lin++) {
		int i = (int) CoresLinAgent.get(lin, 0);
		int j = (int) CoresLinVoisin.get(lin, 0);
		if (lin >= _nTradeP) {
			i -= _nAgentTrue;
		}

		int k = (int) CoresMatLin.get(j, i);
		CoresLinTrans.set(lin, 0, k);
	}

		
	// transfert des mises lineaire
	matUb.transferGPU();
	matLb.transferGPU();
	Ct.transferGPU();

	Tlocal_pre.transferGPU();
	tradeLin.transferGPU();
	LAMBDALin.transferGPU();

	CoresAgentLin.transferGPU();
	CoresLinAgent.transferGPU();
	CoresLinVoisin.transferGPU();
	CoresMatLin.transferGPU();
	CoresLinTrans.transferGPU();
}

void MethodP2PGPU::initSize(const StudyCase& cas){
	_nAgentTrue = cas.getNagent();
	_nAgent = _nAgentTrue + isAC * _nAgentTrue;
	if (cas.isAC() && !isAC) {
		MatrixCPU nVoisinT = cas.getNvoi();
		nVoisinCPU = MatrixCPU(_nAgent, 1);
		for (int n = 0; n < _nAgent; n++) {
			nVoisinCPU.set(n, 0, nVoisinT.get(n, 0));
		}
	}else if(!cas.isAC() && isAC){
		throw std::invalid_argument("initSize : Study Case is not AC, but this method require AC information");
	}
	else {
		nVoisinCPU = cas.getNvoi();
	}
	nVoisin = MatrixGPU(nVoisinCPU, 1);
	nVoisin.preallocateReduction();

	_nLine = cas.getNLine();
	_nBus = cas.getNBus();
	_nTrade = (int) nVoisin.sum();
	_nTradeP = 0;
	if(!isAC){
		_nTradeP = _nTrade;
		_nTradeQ = 0;
	} else{
		for (int n = 0; n < _nAgentTrue; n++) {
			_nTradeP += (int) nVoisin.get(n, 0);
		}
		_nTradeQ = _nTrade - _nTradeP;
	}
	_numBlocksN = ceil((_nAgent + _blockSize - 1) / _blockSize);
	_numBlocksM = ceil((_nTrade + _blockSize - 1) / _blockSize);
	_numBlocksL = ceil((_nLine + _blockSize - 1) / _blockSize);
	_numBlocksNL = ceil((_nAgent*_nLine + _blockSize - 1) / _blockSize);


}

void MethodP2PGPU::initSimParam(const Simparam& sim){
	
	_rhog = sim.getRho();
	_rho1 = sim.getRho1();
	_rhol = _rho;
	if (_rho == 0) {
		_rhol = _rhog;
	}

	_iterG = sim.getIterG();
	_iterL = sim.getIterL();
	_iterIntern = sim.getIterIntern();

	_stepG = sim.getStepG();
	_stepL = sim.getStepL();
	_stepIntern = sim.getStepIntern();

	_epsG = sim.getEpsG();
	_epsX = sim.getEpsGC();
	_epsIntern = sim.getEpsIntern();
	_epsL = sim.getEpsL();
	_ratioEps = _epsG / _epsX;

	resF = MatrixCPU(3, (_iterG / _stepG) + 1);
	resX = MatrixCPU(4, (_iterG / _stepG) + 1);

	tempNN = MatrixGPU(_nTrade, 1, 0, 1);
	tempN1 = MatrixGPU(_nAgent, 1, 0, 1); // plut�t que de re-allouer de la m�moire � chaque utilisation
	tempL1 = MatrixGPU(_nLine, 1, 0, 1);
	tempL2 = MatrixGPU(_nLine, 1, 0, 1);

	
	tempNN.preallocateReduction();
	tempL1.preallocateReduction();

}

void MethodP2PGPU::initDCEndoGrid(const StudyCase& cas){
	
	Kappa1 = MatrixGPU(_nLine, 1, 0, 1);
	Kappa2 = MatrixGPU(_nLine, 1, 0, 1);
	Kappa1_pre = MatrixGPU(_nLine, 1, 0, 1);
	Kappa2_pre = MatrixGPU(_nLine, 1, 0, 1);
	Qpart = MatrixGPU(_nAgent, _nLine, 0, 1);
	Qtot = MatrixGPU(_nLine, 1, 0, 1);
	alpha = MatrixGPU(_nAgent, _nLine, 0, 1);

	G = MatrixGPU(cas.getPowerSensi());

	lLimit = MatrixGPU(cas.getLineLimit(), 1);

	GTrans = MatrixGPU(_nAgent, _nLine);
	if (GTrans.getPos()) {
		GTrans.transferCPU();
		G.transferCPU();
	}


	GTrans.setTrans(&G);

	G.transferGPU();
	GTrans.transferGPU();


	G2 = GTrans;
	G2.multiplyT(&GTrans);

}

void MethodP2PGPU::initDCEndoMarket(){
	initP2PMarket();

	Ap2a = a;
	Ap2b = MatrixGPU(_nAgent, 1, 0, 1);
	Ap3 = MatrixGPU(_nAgent, 1, 0, 1); // not used by default but exists
	Ap123 = MatrixGPU(_nAgent, 1, 0, 1); // idem

	Cp2 = MatrixGPU(_nAgent, 1, 0, 1);
	Cp1 = b;

	Cp1.multiplyT(&nVoisin);
	

	Ap2b.sum(&G2);
	Ap2b.multiply(2 * _rho1);
	Ap2.add(&Ap2a, &Ap2b);

	Ap2.multiplyT(&nVoisin);
	Ap2.multiplyT(&nVoisin);
	Ap12.add(&Ap1, &Ap2);
	Cp = Cp1;
}
void MethodP2PGPU::initP2PMarket(){
	_at1 = _rhog; 
	_at2 = _rhol;
	Ap2 = a;
	Ap1 = nVoisin;
	Ap12 = MatrixGPU(_nAgent, 1, 0, 1);

	Bt1 = MatrixGPU(_nTrade, 1, 0, 1);
	Cp = b;

	
	Pmin.divideT(&nVoisin);
	Pmax.divideT(&nVoisin);
	Ap1.multiply(_rhol);
	Cp.multiplyT(&nVoisin);
	Tmoy.divideT(&nVoisin);
	
	Ap2.multiplyT(&nVoisin);
	Ap2.multiplyT(&nVoisin);
	Ap12.add(&Ap1, &Ap2);

}

void MethodP2PGPU::initCaseParam(const Simparam& sim,const StudyCase& cas){

	Tlocal = MatrixGPU(_nTrade, 1, 0, 1);
	P = MatrixGPU(_nAgent, 1, 0, 1); // moyenne des trades
	Pn = MatrixGPU(_nAgent, 1, 0, 1); // somme des trades

	// si cas AC, a, b , Nvoisin, Pmin, Pmax n'ont pas la bonne taille !!!
	if (cas.isAC() && !isAC) {
		MatrixGPU aT(cas.geta(), 1);
		MatrixGPU bT(cas.getb(), 1);
		MatrixGPU PminT(cas.getPmin(), 1);
		MatrixGPU PmaxT(cas.getPmax(), 1);
		MatrixGPU MUT(sim.getMU(), 1); // facteur reduit i.e lambda_l/_rho
		MatrixGPU TmoyT(sim.getPn(), 1);
		a = MatrixGPU(_nAgent, 1, 0, 1);
		b = MatrixGPU(_nAgent, 1, 0, 1);
		Pmin = MatrixGPU(_nAgent, 1, 0, 1);
		Pmax = MatrixGPU(_nAgent, 1, 0, 1);
		MU = MatrixGPU(_nAgent, 1, 0, 1);
		Tmoy = MatrixGPU(_nAgent, 1, 0, 1);

		
		a.setFromBloc(0, _nAgent, 0, 1, &aT);
		b.setFromBloc(0, _nAgent, 0, 1, &bT);
		Pmax.setFromBloc(0, _nAgent, 0, 1, &PmaxT);
		Pmin.setFromBloc(0, _nAgent, 0, 1, &Pmin);
		MU.setFromBloc(0, _nAgent, 0, 1, &MUT);
		Tmoy.setFromBloc(0, _nAgent, 0, 1, &TmoyT);
	}
	else if(!cas.isAC() && isAC){
		throw std::invalid_argument("initCaseParam : Study Case is not AC, but this method require AC information");
	}
	else {
		a = MatrixGPU(cas.geta(), 1);
		b = MatrixGPU(cas.getb(), 1);

		Pmin = MatrixGPU(cas.getPmin(), 1);
		Pmax = MatrixGPU(cas.getPmax(), 1);
		MU = MatrixGPU(sim.getMU(), 1); // facteur reduit i.e lambda_l/_rho
		Tmoy = MatrixGPU(sim.getPn(), 1);
	}
	Pn = Tmoy;

	Tlocal.preallocateReduction();
	P.preallocateReduction();
}



void MethodP2PGPU::solveWithMinPower(Simparam* result, const Simparam& sim, const StudyCase& cas)
{
	std::cout << "solveWithMinPower : should not be called" << std::endl;
}




float MethodP2PGPU::calcFc()
{
	
	tempN1.set(&a);
	
	tempN1.multiply(0.5);
	tempN1.multiplyT(&Pn);
	
	tempN1.add(&b);
	
	tempN1.multiplyT(&Pn);
	
	float fc = tempN1.sum();
	

	tempNN.set(&trade);
	
	tempNN.multiplyT(&Ct);
	
	fc = fc + tempNN.sum();



	//std::cout << "fc " << fc << std::endl;
	return fc;

}


void MethodP2PGPU::display(){
	std::cout << " resolution par la methode " << _name << std::endl;
}



__global__ void updateLAMBDAGPU(float* LAMBDALin, float* tradeLin, float rho, float* CoresLinTrans, int const N)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int l = index; l < N; l += step)
	{
		float m = LAMBDALin[l];
		int k = CoresLinTrans[l];
		LAMBDALin[l] = m + 0.5 * rho * (tradeLin[l] + tradeLin[k]);
	}
}
__global__ void updateBt1GPU(float* Bt1, float* tradeLin, float rho, float* LAMBDA, float* CoresLinTrans, int const N)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int l = index; l < N; l += step)
	{
		int k = CoresLinTrans[l];
		Bt1[l] = 0.5 * (tradeLin[l] - tradeLin[k]) - LAMBDA[l] / rho;
	}

}

__global__ void updateLAMBDABt1GPU(float* Bt1, float* LAMBDA, float* tradeLin, float rho, float* CoresLinTrans, int const N) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int l = index; l < N; l += step)
	{
		int k = CoresLinTrans[l];
		float m = LAMBDA[l] + 0.5 * rho * (tradeLin[l] + tradeLin[k]);
		Bt1[l] = 0.5 * (tradeLin[l] - tradeLin[k]) - m / rho;
		LAMBDA[l] = m;
	}
}


__global__ void updateDiffGPU(float* tempN, float* Tlocal, float* CoresLinTrans, int const N)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int l = index; l < N; l += step)
	{
		int k = CoresLinTrans[l];
		tempN[l] = (Tlocal[l] + Tlocal[k]);
	}
}

__global__ void updateResKappa(float* result, float* Kappa1, float* Kappa2, float* Kappapre1, float* Kappapre2 ,float ratio, int const L)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int l = index; l < L; l += step)
	{
		float kappaNeg1 = Kappa1[l] < 0 ? Kappa1[l] : 0;
		float kappaNeg2 = Kappa2[l] < 0 ? Kappa2[l] : 0;
		float kappaNegpre1 = Kappapre1[l] < 0 ? Kappapre1[l] : 0;
		float kappaNegpre2 = Kappapre2[l] < 0 ? Kappapre2[l] : 0;

		float res1 = kappaNeg1 - kappaNegpre1;
		res1 *= res1;
		float res2 = kappaNeg2 - kappaNegpre2;
		res2 *= res2;

		result[l] = ratio * sqrt(res1 + res2);
	}
}


__global__ void selectResidual(float* res, unsigned int id1, unsigned int id2, unsigned int id3, float* output) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index == 0) {
		float max = res[id1] > res[id2] ? res[id1] : res[id2];
		max = res[id3] > max ? res[id3] : max;
		*output = max;
	}
}

__global__ void updateKappaGPU(float* Kappa1, float* Kappa2, float* Llimit, float* Qtot, int nLine)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int l = index; l < nLine; l += step)
	{
		float kappaNeg1 = Kappa1[l] < 0 ? Kappa1[l] : 0;
		float kappaNeg2 = Kappa2[l] < 0 ? Kappa2[l] : 0;
		float lim = Llimit[l];
		float Q = Qtot[l];
		Kappa1[l] = kappaNeg1 + lim - Q;
		Kappa2[l] = kappaNeg2 + lim + Q;
	}
}
__global__ void diffKappa(float* tempL1, float* Kappa1, float* Kappa2, int nLine)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int l = index; l < nLine; l += step)
	{
		float Kappa1Abs = Kappa1[l] > 0 ? Kappa1[l] : -Kappa1[l]; //2 * (Kappa1[l] > 0) * Kappa1[l] - Kappa1[l]; // Kappa1[l] > 0 ? Kappa1[l] : -Kappa1[l]
		float Kappa2Abs = Kappa2[l] > 0 ? Kappa2[l] : -Kappa2[l]; //2 * (Kappa2[l] > 0) * Kappa2[l] - Kappa2[l]; // Kappa2[l] > 0 ? Kappa2[l] : -Kappa2[l]
		tempL1[l] = Kappa1Abs - Kappa2Abs;
	}
}


__global__ void updateCpOld(float* Cp, float* Cp1, float* Cp2, float* tempN1, float* nVoisin, const float rho1, const int nAgent) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int n = index; n < nAgent; n += step)
	{
		Cp[n] = Cp1[n] + rho1 * (nVoisin[n] * (Cp2[n] + tempN1[n]));
	}
}

__global__ void updateCp(float* Cp, float* Cp1, float* Cp2, const int nAgent) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int n = index; n < nAgent; n += step)
	{
		Cp[n] = Cp1[n] + Cp2[n];
	}
}

__global__ void updateQpart(float* Qpart, float* alpha, const int nAgent) {
	int index = threadIdx.x;
	int step = blockDim.x;
	int l = blockIdx.x;
	extern __shared__ float shAlpha[];

	for (int n = index; n < nAgent; n += step)
	{
		shAlpha[n] = alpha[l * nAgent + n];
	}
	__syncthreads();

	float s_pre = 0;
	int n_pre = nAgent - 1;
	for (int n = (nAgent - index - 1); n >= 0; n -= step)
	{
		float s = 0;
		for (int j = n_pre; j > n; j--) {
			s += shAlpha[j]; // c'est moche cet acc�s de m�moire partag�e
		}
		s = s + s_pre;
		Qpart[l * nAgent + n] = s;
		s_pre = s;
		n_pre = n;
	}
}
__global__ void updateQpartTrans(float* Qpart, float* alpha, const int N, const int nLine) {
	int index = threadIdx.x;
	int step = blockDim.x;
	int l = blockIdx.x;
	extern __shared__ float shAlpha[];

	for (int n = index; n < N; n += step)
	{
		shAlpha[n] = alpha[n * nLine + l]; // moche
	}
	__syncthreads();
	float s_pre = 0;
	int n_pre = N - 1;
	for (int n = (N - index - 1); n >= 0; n -= step)
	{
		float s = 0;
		for (int j = n_pre; j > n; j--) {
			s += shAlpha[j]; // c'est moche cet acc�s de m�moire partag�e
		}
		s = s + s_pre;
		Qpart[n * nLine + l] = s;
		s_pre = s;
		n_pre = n;
	}
}


__global__ void updateQtot(float* Qtot, float* Qpart, float* alpha, const int nLine, const int nAgent) {


	int thIdx = threadIdx.x + blockIdx.x * blockDim.x;
	int step = blockDim.x * gridDim.x;

	for (int l = thIdx; l < nLine; l += step) {

		Qtot[l] = Qpart[l * nAgent] + alpha[l * nAgent];
	}
}
__global__ void updateQtotTrans(float* Qtot, float* Qpart, float* alpha, const int nLine) {


	int thIdx = threadIdx.x + blockIdx.x * blockDim.x;
	int step = blockDim.x * gridDim.x;

	for (int l = thIdx; l < nLine; l += step) {

		Qtot[l] = Qpart[l] + alpha[l];
	}
}

__global__ void updateAlpha(float* alpha, float* G, float* Pn, const int nLine, const int nAgent)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int step = blockDim.x * gridDim.x;
	int N = nAgent * nLine;

	for (int i = index; i < N; i += step)
	{
		int k = i % nAgent;
		alpha[i] = G[i] * Pn[k];
	}
}
__global__ void updateAlphaTrans(float* alpha, float* GTrans, float* Pn, const int nLine, const int nAgent) {

	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int step = blockDim.x * gridDim.x;
	int N = nAgent * nLine;

	for (int i = index; i < N; i += step)
	{
		int k = i / nLine;
		alpha[i] = GTrans[i] * Pn[k];
	}

}


__global__ void updateResX(float* res, float* Kappa1, float* Kappa2, float* KappaPre1, float* KappaPre2, const int nLine) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int l = index; l < nLine; l += step)
	{
		float k1 = (Kappa1[l] < 0) * Kappa1[l];
		float k2 = (Kappa2[l] < 0) * Kappa2[l];
		float k1pre = (KappaPre1[l] < 0) * KappaPre1[l];
		float k2pre = (KappaPre2[l] < 0) * KappaPre2[l];

		k1 -= k1pre;
		k2 -= k2pre;

		res[l] = k1 * k1 + k2 * k2;
	}
}



__global__ void updatePnGPU(float* Pn, float* Tmoy, float* nVoisin, const int nAgent)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int n = index; n < nAgent; n += step)
	{
		Pn[n] = Tmoy[n] * nVoisin[n];
	}

}



/*
__global__ void updateUAiq(float* UAiq, float* u, float* Aiq, int N, int size) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int n = index; n < size; n += step)
	{
		int row = n / N;
		UAiq[n] = u[row] * Aiq[n];
	}
}

// Ru = U*g + epsi
__global__ void updateRu(float* Ru, float* U, float* g, float epsi, int N, int L2) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int i = index; i < L2; i += step)
	{
		
		Ru[i + N] = U[i] * g[i] + epsi;
	}

}

__global__ void updateV(float* v, float* pas, float* alpha, int offset) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index == 0) {
		*v = *v + (*alpha) * pas[offset];
	}
}

__global__ void updateQt(float* qt, float* Pso, float* Pn, float* etaSO, float rho1, int N) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int i = index; i < N; i += step)
	{

		qt[i] = etaSO[i] - rho1 * (Pso[i] + Pn[i]) / 2;
	}


}*/

__global__ void updatePI(float* PI, float* c, float mu, float valMin, int L) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int i = index; i < L; i += step)
	{
		PI[i] = c[i] < valMin ? mu / valMin : mu / c[i];
	}
	if (index == 0) {
		PI[L] = -c[L] / mu;
	}
}


__global__ void updatePso(float* Pso, float* pas, float* alpha, int N) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int i = index; i < N; i += step)
	{

		Pso[i] = Pso[i] + *alpha * pas[i];
	}

}

__global__ void updateU(float* U, float* pas, float* alpha, int N, int L2) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int i = index; i < L2; i += step)
	{

		U[i] = U[i] + (*alpha) * pas[i + N];
	}

}


__global__ void updateEtaPBp3(float* Bp3, float* etaP, float* nVoisin, float* Pso, float* Pn, float rho, const int nAgent) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int step = blockDim.x * gridDim.x;
	for (int l = index; l < nAgent; l += step)
	{
		
		float m = etaP[l] + 0.5 * rho * (Pso[l] - Pn[l]);
		Bp3[l] = (0.5 * (Pso[l] + Pn[l]) + m / rho) / nVoisin[l];
		etaP[l] = m;
	}


}