#pragma once
#include "../head/MethodP2P.cuh"
#define MAX(X, Y) X * (X >= Y) + Y * (Y > X)


MethodP2P::MethodP2P() : Method()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "method constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1abcd, Fb2, Fb3, Fb5, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilisé pendant la simu
}

MethodP2P::~MethodP2P() 
{
}

void MethodP2P::updateLAMBDA(MatrixCPU* LAMBDA, MatrixCPU* trade, float rho)
{
	MatrixCPU times = MatrixCPU(*trade);
	times.addTrans(trade);
	times.multiply(rho);
	times.multiply(0.5);
	LAMBDA->add(LAMBDA, &times);
}

float MethodP2P::updateRes(MatrixCPU* res, MatrixCPU* Tlocal, MatrixCPU* trade, int iter)
{
	MatrixCPU temp(*Tlocal);
	temp.addTrans(Tlocal);

	MatrixCPU temp2(*Tlocal);
	float resR = temp.max2();
	temp2.subtract(trade);
	float resS = temp2.max2();

	res->set(0, iter, resR);
	res->set(1, iter, resS);

	return resR * (resR > resS) + resS * (resR <= resS);

}

float MethodP2P::updateRes(MatrixCPU* res, MatrixCPU* Tlocal, MatrixCPU* trade, int iter, MatrixCPU* Kappa1, MatrixCPU* Kappa2, MatrixCPU* Kappa1_pre, MatrixCPU* Kappa2_pre)
{
	MatrixCPU temp(*Tlocal);
	temp.addTrans(Tlocal);

	MatrixCPU temp2(*Tlocal);
	float resR = temp.max2();
	temp2.subtract(trade);
	float resS = temp2.max2();

	MatrixCPU tempL(*Kappa1);
	MatrixCPU tempL2(*Kappa2);
	Kappa1_pre->projectNeg();
	Kappa2_pre->projectNeg();
	tempL.projectNeg();
	tempL2.projectNeg();
	tempL.subtract(Kappa1_pre);
	tempL2.subtract(Kappa2_pre);
	tempL.multiplyT(&tempL);
	tempL2.multiplyT(&tempL2);
	tempL.add(&tempL2);

	float resX = _ratioEps * sqrt(tempL.max2());
	
	res->set(0, iter, resR);
	res->set(1, iter, resS);
	res->set(2, iter, resX);

	return MAX(MAX(resX,resS),resR);
}




void MethodP2P::updateLAMBDA(MatrixGPU* LAMBDA, MatrixGPU* trade, float rho, MatrixGPU* tempNN)
{
	tempNN->set(trade);
	tempNN->addTrans(trade);
	tempNN->multiply(rho);
	tempNN->multiply(0.5);
	LAMBDA->add(LAMBDA, tempNN);
}

void MethodP2P::updateKappa(MatrixCPU* Kappa1, MatrixCPU* Kappa2, MatrixCPU* L, MatrixCPU* Qtot)
{
	
	
	//
	Kappa1->projectNeg();
	Kappa1->add(L);
	Kappa1->subtract(Qtot);
	
	Kappa2->projectNeg();
	Kappa2->add(L);
	Kappa2->add(Qtot);
	//
	
}

void MethodP2P::updateKappa(MatrixGPU* Kappa1, MatrixGPU* Kappa2, MatrixGPU* L, MatrixGPU* Qtot)
{
	Kappa1->projectNeg();
	Kappa1->add(L);
	Kappa1->subtract(Qtot);
	Kappa2->projectNeg();
	Kappa2->add(L);
	Kappa2->add(Qtot);
}

void MethodP2P::updateCp2(MatrixCPU* Cp2, float rho1, MatrixCPU* Kappa1, MatrixCPU* Kappa2, MatrixCPU* G, MatrixCPU* tempL1, MatrixCPU* Qpart, MatrixCPU* nVoisin, int nLine, int nAgent)
{
	tempL1->subtractAbs(Kappa1, Kappa2);
	//Cp2->multiplyTrans(G, tempL1, 0);

	float r = 0;
	for (int i = 0; i < nAgent; ++i)
	{
		r = 0;
		for (int k = 0; k < nLine; ++k)
		{
			r +=  G->get(k, i) * (tempL1->get(k, 0) + 2 * Qpart->get(k, i));
		}
		Cp2->set(i, 0, r);
	}

	Cp2->multiply(rho1);
	Cp2->multiplyT(nVoisin);
}

float MethodP2P::updateRes(MatrixCPU* res, MatrixGPU* Tlocal, MatrixGPU* trade, int iter, MatrixGPU* tempNN)
{
	tempNN->subtract(Tlocal, trade);
	
	float resS = tempNN->max2();
	tempNN->set(Tlocal);
	tempNN->addTrans(Tlocal);
	float resR = tempNN->max2();
	

	res->set(0, iter, resR);
	res->set(1, iter, resS);
	

	return resR* (resR > resS) + resS * (resR <= resS);
}





void MethodP2P::updatePn(MatrixCPU* Pn, MatrixCPU* Tmoy, MatrixCPU* nVoisin)
{
	Pn->set(Tmoy);
	Pn->multiplyT(nVoisin);
}
void MethodP2P::updatePn(MatrixGPU* Pn, MatrixGPU* Tmoy, MatrixGPU* nVoisin)
{
	Pn->set(Tmoy);
	Pn->multiplyT(nVoisin);
}
void MethodP2P::updatePn(MatrixCPU* Pn, MatrixCPU* trade)
{
	Pn->sum(trade);
}

void MethodP2P::solveWithMinPower(Simparam* result, const Simparam& sim, const StudyCase& cas)
{
	std::cout << "solveWithMinPower : should not be called" << std::endl;
}




void MethodP2P::resetId()
{
	_id = 0;
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
			s += shAlpha[j]; // c'est moche cet accès de mémoire partagée
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
			s += shAlpha[j]; // c'est moche cet accès de mémoire partagée
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

