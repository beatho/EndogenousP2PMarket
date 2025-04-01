#include "../head/MethodP2PGPU.cuh"
 


MethodP2PGPU::MethodP2PGPU(){

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