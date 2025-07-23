#pragma once
#include <device_launch_parameters.h>
#include <cuda_runtime.h>
#include "MatrixGPU.cuh"
#include "Utilities.cuh"
#include "Utilities.h"


#define NSTEPLOCAL 5
#define NMAXPEERPERTRHREAD 5


////IMPLEMENTATION APPEL KERNEL avec template

template <unsigned int blockSize>
__global__ void updateTradePGPU(float* Tlocal, float* Tlocal_pre, float* Tmoy, float* P, float* MU, float* nVoisin, float at1, float at2, float* Bt1, float* Ct,
	float* matLb, float* matUb, float* Ap1, float* Ap12, float* Cp, float* Pmin, float* Pmax, float* CoresAgentLin, float* CoresLinAgent, int const n) {

	__shared__ float shArr[blockSize];
	unsigned int thIdx = threadIdx.x;
	const int step = blockSize;
	float sum = 0;
	int i = blockIdx.x;
	int nVoisinLocal = nVoisin[i];
	int beginLocal = CoresAgentLin[i];
	int endLocal = beginLocal + nVoisinLocal; // = CoresAgentLin[blockIdx.x + 1]

	for (int j = thIdx + beginLocal; j < endLocal; j += step) {
		float m = Tlocal_pre[j] - Tmoy[i] + P[i] - MU[i]; 
		float r = (Bt1[j] * at1 + m * at2 - Ct[j]) / (at1 + at2); 
		float ub = matUb[j]; 
		float lb = matLb[j]; 
		float t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
		Tlocal[j] = t; 
		sum += t;
	}

	shArr[thIdx] = sum;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduce<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		float moy = shArr[0] / nVoisinLocal;
		Tmoy[blockIdx.x] = moy;
		float muLocal = MU[blockIdx.x];
		float b = moy + muLocal;
		float p = (Ap1[blockIdx.x] * b - Cp[blockIdx.x]) / Ap12[blockIdx.x];
		float Pub = Pmax[blockIdx.x];
		float Plb = Pmin[blockIdx.x];
		p = (Pub - p) * (p > Pub) + (Plb - p) * (p < Plb) + p;
		P[blockIdx.x] = p;
		MU[blockIdx.x] = muLocal + moy - p;
	}

}


template <unsigned int blockSize>
__global__ void updateTradePGPULocal(float* Tlocal, float* Tlocal_pre, float* Tmoy, float* P, float* MU, float* nVoisin, float at1, float at2, float* Bt1, float* Ct,
	float* matLb, float* matUb, float* Ap1, float* Ap12, float* Cp, float* Pmin, float* Pmax, float* CoresAgentLin) {

	//Definition de toutes les variables locales
	int i = blockIdx.x; // c'est aussi l'identifiant de l'agent !
	unsigned int thIdx = threadIdx.x;
	const int step = blockSize;
	// ne change pas
	const float Ap1Local = Ap1[i];
	const float CpLocal = Cp[i];
	const float Ap12Local = Ap12[i];
	const float Pub = Pmax[i];
	const float Plb = Pmin[i];
	const int nVoisinLocal = nVoisin[i];
	const int CoresAgentLinLocal = CoresAgentLin[i];
	const int beginLocal = CoresAgentLinLocal + thIdx;
	const int endLocal = CoresAgentLinLocal + nVoisinLocal;

	const float at1local = at1;
	const float at2local = at2;
	const float at12local = at1local + at2local;


	float Bt1local[NMAXPEERPERTRHREAD];
	float Ctlocal[NMAXPEERPERTRHREAD];
	float matUblocal[NMAXPEERPERTRHREAD];
	float matLblocal[NMAXPEERPERTRHREAD];

	float Tlocallocal[NMAXPEERPERTRHREAD]; // change
	float Tlocalprelocal[NMAXPEERPERTRHREAD]; // change
	float MULOCAL;
	float moy;
	float p;
	float sum;
	float bp;
	float m, r, ub, lb, t;
	// le changement doit être partagé par tous les threads du bloc

	__shared__ float MuShared;
	__shared__ float TMoyShared;
	__shared__ float PShared;
	if (thIdx == 0) {
		MuShared = MU[i];
		TMoyShared = Tmoy[i];
		PShared = P[i];
	}
	int k = 0;
	for (int j = beginLocal; j < endLocal; j += step) {
		Bt1local[k] = Bt1[j];
		Ctlocal[k] = Ct[j];
		matUblocal[k] = matUb[j];
		matLblocal[k] = matLb[j];
		//Tlocalprelocal[k] = Tlocal_pre[j];
		Tlocallocal[k] = Tlocal_pre[j];
		k = k + 1;
	}

	__shared__ float shArr[blockSize];


	//Calcul des itérations
	__syncthreads();
	for (int iter = 0; iter < NSTEPLOCAL; iter++) {

		MULOCAL = MuShared; // tous lisent le même : broadcast !
		moy = TMoyShared;
		p = PShared;
		sum = 0;
		k = 0;
		for (int j = beginLocal; j < endLocal; j += step) {
			Tlocalprelocal[k] = Tlocallocal[k];
			m = Tlocallocal[k] - moy + p - MULOCAL;
			r = (Bt1local[k] * at1local + m * at2local - Ctlocal[k]) / (at12local);
			ub = matUblocal[k];
			lb = matLblocal[k];
			t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
			Tlocallocal[k] = t;
			sum += t;
			k = k + 1;
		}

		shArr[thIdx] = sum;
		__syncthreads();
		if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
		if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
		if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
		if (thIdx < 32) {
			warpReduce<blockSize>(shArr, thIdx);
		}
		__syncthreads();

		if (thIdx == 0) {
			moy = shArr[0] / nVoisinLocal;
			TMoyShared = moy;
			bp = moy + MuShared;
			p = (Ap1Local * bp - CpLocal) / Ap12Local;
			p = (Pub - p) * (p > Pub) + (Plb - p) * (p < Plb) + p;
			PShared = p;
			MuShared = MULOCAL + moy - p;
		}
		__syncthreads();
	}
	//Ecriture des itérations
	__syncthreads();
	k = 0;
	for (int j = beginLocal; j < endLocal; j += step) {
		Tlocal[j] = Tlocallocal[k];
		Tlocal_pre[j] = Tlocalprelocal[k];
		k = k + 1;
	}
	if (thIdx == 0) {
		Tmoy[blockIdx.x] = TMoyShared;// TMoyShared;
		P[blockIdx.x] = PShared;// PShared;
		MU[blockIdx.x] = MuShared;// MuShared;
	}

}


template <unsigned int blockSize>
__global__ void updateTradePGPUShared(float* Tlocal, float* Tlocal_pre, float* Tmoy, float* P, float* MU, float* nVoisin, float at1, float at2, float* Bt1, float* Ct,
	float* matLb, float* matUb, float* Ap1, float* Ap12, float* Cp, float* Pmin, float* Pmax, float* CoresAgentLin) {

	//Definition de toutes les variables locales
	int i = blockIdx.x; // c'est aussi l'identifiant de l'agent !
	unsigned int thIdx = threadIdx.x;
	const int step = blockSize;
	// ne change pas


	float Bt1local[NMAXPEERPERTRHREAD];
	float Ctlocal[NMAXPEERPERTRHREAD];
	float matUblocal[NMAXPEERPERTRHREAD];
	float matLblocal[NMAXPEERPERTRHREAD];

	float Tlocallocal[NMAXPEERPERTRHREAD]; // change
	float Tlocalprelocal[NMAXPEERPERTRHREAD]; // change
	float sum;
	float bp, MULOCAL, moy, p;
	float m, r, ub, lb, t;


	// le changement doit être partagé par tous les threads du bloc

	__shared__ float MuShared;
	__shared__ float TMoyShared;
	__shared__ float PShared;


	// constant et commun à tous les thread d'un bloc
	__shared__ float Ap1Shared;
	__shared__ float CpShared;
	__shared__ float Ap12Shared;
	__shared__ float PmaxShared;
	__shared__ float PminShared;
	__shared__ float nVoisinShared;
	__shared__ float at1Shared;
	__shared__ float at2Shared;
	__shared__ float at12Shared;


	if (thIdx == 0) {
		Ap1Shared = Ap1[i];
		CpShared = Cp[i];
		Ap12Shared = Ap12[i];
		PmaxShared = Pmax[i];
		PminShared = Pmin[i];
		nVoisinShared = nVoisin[i];
		at1Shared = at1;
		at2Shared = at2;
		at12Shared = at1 + at2;
		MuShared = MU[i];
		TMoyShared = Tmoy[i];
		PShared = P[i];
	}
	int k = 0;
	__syncthreads();
	const int CoresAgentLinLocal = CoresAgentLin[i];
	const int beginLocal = CoresAgentLinLocal + thIdx;
	const int endLocal = CoresAgentLinLocal + nVoisinShared;
	for (int j = beginLocal; j < endLocal; j += step) {
		Bt1local[k] = Bt1[j];
		Ctlocal[k] = Ct[j];
		matUblocal[k] = matUb[j];
		matLblocal[k] = matLb[j];
		//Tlocalprelocal[k] = Tlocal_pre[j];
		Tlocallocal[k] = Tlocal_pre[j];
		k = k + 1;
	}

	__shared__ float shArr[blockSize];

	//Calcul des itérations

	for (int iter = 0; iter < NSTEPLOCAL; iter++) {

		MULOCAL = MuShared; // tous lisent le même : broadcast !
		moy = TMoyShared;
		p = PShared;
		sum = 0;
		k = 0;
		for (int j = beginLocal; j < endLocal; j += step) {
			Tlocalprelocal[k] = Tlocallocal[k];
			m = Tlocallocal[k] - moy + p - MULOCAL;
			r = (Bt1local[k] * at1Shared + m * at2Shared - Ctlocal[k]) / (at12Shared);
			ub = matUblocal[k];
			lb = matLblocal[k];
			t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
			Tlocallocal[k] = t;
			sum += t;
			k = k + 1;
		}

		shArr[thIdx] = sum;
		__syncthreads();
		if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
		if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
		if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
		if (thIdx < 32) {
			warpReduce<blockSize>(shArr, thIdx);
		}
		__syncthreads();

		if (thIdx == 0) {
			moy = shArr[0] / nVoisinShared;
			TMoyShared = moy;
			bp = moy + MuShared;
			p = (Ap1Shared * bp - CpShared) / Ap12Shared;
			p = (PmaxShared - p) * (p > PmaxShared) + (PminShared - p) * (p < PminShared) + p; 
			PShared = p;
			MuShared = MULOCAL + moy - p;
		}
		__syncthreads();
	}
	//Ecriture des itérations
	__syncthreads();
	k = 0;
	for (int j = beginLocal; j < endLocal; j += step) {
		Tlocal[j] = Tlocallocal[k];
		Tlocal_pre[j] = Tlocalprelocal[k];
		k = k + 1;
	}
	if (thIdx == 0) {
		Tmoy[blockIdx.x] = TMoyShared;// TMoyShared;
		P[blockIdx.x] = PShared;// PShared;
		MU[blockIdx.x] = MuShared;// MuShared;
	}

}


// problem local en une fois mais propre
template <unsigned int blockSize>
__global__ void updateTradePGPUSharedResidual(float* Tlocal, float* Tlocal_pre, float* Tmoy, float* P, float* MU, float* nVoisin, float at1, float at2, float* Bt1, float* Ct,
	float* matLb, float* matUb, float* Ap1, float* Ap12, float* Cp, float* Pmin, float* Pmax, float* CoresAgentLin, float eps, int nStepL) {

	//Definition de toutes les variables locales
	int i = blockIdx.x; // c'est aussi l'identifiant de l'agent !
	unsigned int thIdx = threadIdx.x;
	const int step = blockSize;
	// ne change pas


	float Bt1local[NMAXPEERPERTRHREAD];
	float Ctlocal[NMAXPEERPERTRHREAD];
	float matUblocal[NMAXPEERPERTRHREAD];
	float matLblocal[NMAXPEERPERTRHREAD];

	float Tlocallocal[NMAXPEERPERTRHREAD]; // change
	float Tlocalprelocal[NMAXPEERPERTRHREAD]; // change
	float sum;
	float bp, MULOCAL, moy, p;
	float m, r, ub, lb, t;
	bool mustContinueLocal;
	// le changement doit être partagé par tous les threads du bloc

	__shared__ float MuShared;
	__shared__ float TMoyShared;
	__shared__ float PShared;


	// constant et commun à tous les thread d'un bloc
	__shared__ float Ap1Shared;
	__shared__ float CpShared;
	__shared__ float Ap12Shared;
	__shared__ float PmaxShared;
	__shared__ float PminShared;
	__shared__ float nVoisinShared;
	__shared__ float at1Shared;
	__shared__ float at2Shared;
	__shared__ float at12Shared;
	__shared__ bool mustContinue;


	if (thIdx == 0) {
		Ap1Shared = Ap1[i];
		CpShared = Cp[i];
		Ap12Shared = Ap12[i];
		PmaxShared = Pmax[i];
		PminShared = Pmin[i];
		nVoisinShared = nVoisin[i];
		at1Shared = at1;
		at2Shared = at2;
		at12Shared = at1 + at2;
		MuShared = MU[i];
		TMoyShared = Tmoy[i];
		PShared = P[i];
		mustContinue = false;
	}
	int k = 0;
	__syncthreads();
	const int CoresAgentLinLocal = CoresAgentLin[i];
	const int beginLocal = CoresAgentLinLocal + thIdx;
	const int endLocal = CoresAgentLinLocal + nVoisinShared;
	float res;
	for (int j = beginLocal; j < endLocal; j += step) {
		Bt1local[k] = Bt1[j];
		Ctlocal[k] = Ct[j];
		matUblocal[k] = matUb[j];
		matLblocal[k] = matLb[j];
		//Tlocalprelocal[k] = Tlocal_pre[j];
		Tlocallocal[k] = Tlocal_pre[j];
		k = k + 1;
	}

	__shared__ float shArr[blockSize];
	__shared__ bool mustcontinueVect[blockSize];

	//Calcul des itérations

	for (int iter = 0; iter < nStepL; iter++) {

		MULOCAL = MuShared; // tous lisent le même : broadcast !
		moy = TMoyShared;
		p = PShared;
		sum = 0;
		mustContinueLocal = false;
		k = 0;

		for (int j = beginLocal; j < endLocal; j += step) {
			Tlocalprelocal[k] = Tlocallocal[k];
			m = Tlocallocal[k] - moy + p - MULOCAL;
			r = (Bt1local[k] * at1Shared + m * at2Shared - Ctlocal[k]) / (at12Shared);
			ub = matUblocal[k];
			lb = matLblocal[k];
			t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
			Tlocallocal[k] = t;
			sum += t;
			res = (t - Tlocalprelocal[k]);
			res = (double)res * res;
			if (res > eps) {
				mustContinueLocal = true;
				//mustContinue = true; // pas de race condition, car l'ordre n'importe pas,
				//mais est ce que cela ne va pas physiquement bloquer ?
			}
			k = k + 1;
		}
		shArr[thIdx] = sum;
		mustcontinueVect[thIdx] = mustContinueLocal;
		__syncthreads();
		if (blockSize >= 512) {
			if (thIdx < 256) {
				shArr[thIdx] += shArr[thIdx + 256];
				mustcontinueVect[thIdx] = mustcontinueVect[thIdx] || mustcontinueVect[thIdx + 256];
			}
			__syncthreads();
		}
		if (blockSize >= 256) {
			if (thIdx < 128) {
				shArr[thIdx] += shArr[thIdx + 128];
				mustcontinueVect[thIdx] = mustcontinueVect[thIdx] || mustcontinueVect[thIdx + 128];
			}
			__syncthreads();
		}
		if (blockSize >= 128) {
			if (thIdx < 64) {
				shArr[thIdx] += shArr[thIdx + 64];
				mustcontinueVect[thIdx] = mustcontinueVect[thIdx] || mustcontinueVect[thIdx + 64];
			}
			__syncthreads();
		}
		if (blockSize >= 64) {
			if (thIdx < 32) {
				warpReduce<blockSize>(shArr, thIdx);
				warpReduceOr<blockSize>(mustcontinueVect, thIdx);
			}
		}
		else if (blockSize >= 32) {
			if (thIdx < 16) {
				warpReduce<blockSize>(shArr, thIdx);
				warpReduceOr<blockSize>(mustcontinueVect, thIdx);
			}
		}

		__syncthreads();

		if (thIdx == 0) {
			mustContinue = mustcontinueVect[0];
			moy = shArr[0] / nVoisinShared;
			TMoyShared = moy;
			bp = moy + MuShared;
			p = (Ap1Shared * bp - CpShared) / Ap12Shared;
			p = (PmaxShared - p) * (p > PmaxShared) + (PminShared - p) * (p < PminShared) + p;
			PShared = p;
			res = p - moy;
			res = (double)res * res;
			if (res > eps) {
				mustContinue = true;
			}
			MuShared = MULOCAL + moy - p;
		}
		__syncthreads();
		if (!mustContinue) {
			break;
		}
		/*else {
			__syncthreads();
			if (thIdx == 0) {
				mustContinue = false;
			}
			__syncthreads();
		}*/
	}
	//Ecriture des itérations
	__syncthreads();
	k = 0;
	for (int j = beginLocal; j < endLocal; j += step) {
		Tlocal[j] = Tlocallocal[k];
		Tlocal_pre[j] = Tlocalprelocal[k];
		k = k + 1;
	}
	if (thIdx == 0) {
		Tmoy[blockIdx.x] = TMoyShared;// TMoyShared;
		P[blockIdx.x] = PShared;// PShared;
		MU[blockIdx.x] = MuShared;// MuShared;
	}

}



template <unsigned int blockSize>
__global__ void updateTradePGPUSharedResidual(float* Tlocal, float* Tlocal_pre, float* Tmoy, float* P, float* MU, float* nVoisin, float at1, float at2, float* Bt1, float* Ct,
	float* matLb, float* matUb, float* Ap1, float* Ap2, float* Ap123, float* Bp2, float* Cp, float* Pmin, float* Pmax, float* CoresAgentLin, float eps, int nStepL) {

	//Definition de toutes les variables locales
	int i = blockIdx.x; // c'est aussi l'identifiant de l'agent !
	unsigned int thIdx = threadIdx.x;
	const int step = blockSize;
	// ne change pas


	float Bt1local[NMAXPEERPERTRHREAD];
	float Ctlocal[NMAXPEERPERTRHREAD];
	float matUblocal[NMAXPEERPERTRHREAD];
	float matLblocal[NMAXPEERPERTRHREAD];

	float Tlocallocal[NMAXPEERPERTRHREAD]; // change
	float Tlocalprelocal[NMAXPEERPERTRHREAD]; // change
	float sum;
	float bp, MULOCAL, moy, p;
	float m, r, ub, lb, t;
	bool mustContinueLocal;
	// le changement doit être partagé par tous les threads du bloc

	__shared__ float MuShared;
	__shared__ float TMoyShared;
	__shared__ float PShared;


	// constant et commun à tous les thread d'un bloc
	__shared__ float Ap1Shared;
	__shared__ float Ap2Shared;
	__shared__ float Bp2Shared;
	__shared__ float CpShared;
	__shared__ float Ap123Shared;
	__shared__ float PmaxShared;
	__shared__ float PminShared;
	__shared__ float nVoisinShared;
	__shared__ float at1Shared;
	__shared__ float at2Shared;
	__shared__ float at12Shared;
	__shared__ bool mustContinue;


	if (thIdx == 0) {
		Ap1Shared = Ap1[i];
		Ap2Shared = Ap2[i];
		Bp2Shared = Bp2[i];
		CpShared = Cp[i];
		Ap123Shared = Ap123[i];
		PmaxShared = Pmax[i];
		PminShared = Pmin[i];
		nVoisinShared = nVoisin[i];
		at1Shared = at1;
		at2Shared = at2;
		at12Shared = at1 + at2;
		MuShared = MU[i];
		TMoyShared = Tmoy[i];
		PShared = P[i];
		mustContinue = false;
	}
	int k = 0;
	__syncthreads();
	const int CoresAgentLinLocal = CoresAgentLin[i];
	const int beginLocal = CoresAgentLinLocal + thIdx;
	const int endLocal = CoresAgentLinLocal + nVoisinShared;
	float res;
	for (int j = beginLocal; j < endLocal; j += step) {
		Bt1local[k] = Bt1[j];
		Ctlocal[k] = Ct[j];
		matUblocal[k] = matUb[j];
		matLblocal[k] = matLb[j];
		//Tlocalprelocal[k] = Tlocal_pre[j];
		Tlocallocal[k] = Tlocal_pre[j];
		k = k + 1;
	}

	__shared__ float shArr[blockSize];
	__shared__ bool mustcontinueVect[blockSize];
	//Calcul des itérations

	for (int iter = 0; iter < nStepL; iter++) {

		MULOCAL = MuShared; // tous lisent le même : broadcast !
		moy = TMoyShared;
		p = PShared;
		mustContinueLocal = false;
		sum = 0;
		k = 0;
		for (int j = beginLocal; j < endLocal; j += step) {
			Tlocalprelocal[k] = Tlocallocal[k];
			m = Tlocallocal[k] - moy + p - MULOCAL;
			r = (Bt1local[k] * at1Shared + m * at2Shared - Ctlocal[k]) / (at12Shared);
			ub = matUblocal[k];
			lb = matLblocal[k];
			t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
			Tlocallocal[k] = t;
			sum += t;
			res = (t - Tlocalprelocal[k]);
			res = (double)res * res;
			if (res > eps) {
				mustContinueLocal = true;
			}
			k = k + 1;
		}

		shArr[thIdx] = sum;
		mustcontinueVect[thIdx] = mustContinueLocal;
		__syncthreads();
		if (blockSize >= 512) {
			if (thIdx < 256) {
				shArr[thIdx] += shArr[thIdx + 256];
				mustcontinueVect[thIdx] = mustcontinueVect[thIdx] || mustcontinueVect[thIdx + 256];
			}
			__syncthreads();
		}
		if (blockSize >= 256) {
			if (thIdx < 128) {
				shArr[thIdx] += shArr[thIdx + 128];
				mustcontinueVect[thIdx] = mustcontinueVect[thIdx] || mustcontinueVect[thIdx + 128];
			}
			__syncthreads();
		}
		if (blockSize >= 128) {
			if (thIdx < 64) {
				shArr[thIdx] += shArr[thIdx + 64];
				mustcontinueVect[thIdx] = mustcontinueVect[thIdx] || mustcontinueVect[thIdx + 64];
			}
			__syncthreads();
		}
		if (blockSize >= 64) {
			if (thIdx < 32) {
				warpReduce<blockSize>(shArr, thIdx);
				warpReduceOr<blockSize>(mustcontinueVect, thIdx);
			}
		}
		else if (blockSize >= 32) {
			if (thIdx < 16) {
				warpReduce<blockSize>(shArr, thIdx);
				warpReduceOr<blockSize>(mustcontinueVect, thIdx);
			}
		}

		__syncthreads();

		if (thIdx == 0) {
			moy = shArr[0] / nVoisinShared;
			mustContinue = mustcontinueVect[0];
			TMoyShared = moy;
			bp = moy + MuShared;
			p = (Ap1Shared * bp + Ap2Shared * Bp2Shared - CpShared) / Ap123Shared;
			p = (PmaxShared - p) * (p > PmaxShared) + (PminShared - p) * (p < PminShared) + p;
			PShared = p;
			res = p - moy;
			res = (double)res * res;
			if (res > eps) {
				mustContinue = true;
			}
			MuShared = MULOCAL + moy - p;
		}
		__syncthreads();
		if (!mustContinue) {
			break;
		}
	}
	//Ecriture des itérations
	__syncthreads();
	k = 0;
	for (int j = beginLocal; j < endLocal; j += step) {
		Tlocal[j] = Tlocallocal[k];
		Tlocal_pre[j] = Tlocalprelocal[k];
		k = k + 1;
	}
	if (thIdx == 0) {
		Tmoy[blockIdx.x] = TMoyShared;// TMoyShared;
		P[blockIdx.x] = PShared;// PShared;
		MU[blockIdx.x] = MuShared;// MuShared;
	}

}




template <unsigned int blockSize>
__global__ void updateTradePGPUSharedResidualCons(float* Tlocal, float* Tlocal_pre, float* Tmoy, float* P, float* MU, float* nVoisin, float at1, float at2, float* Bt1, float* Ct,
	float* matLb, float* matUb, float* Ap1, float* Ap3, float* Ap123, float* Bp3, float* Cp, float* Pmin, float* Pmax, float* CoresAgentLin, float eps, int nStepL) {

	//Definition de toutes les variables locales
	int i = blockIdx.x; // c'est aussi l'identifiant de l'agent !
	unsigned int thIdx = threadIdx.x;
	const int step = blockSize;
	// ne change pas


	float Bt1local[NMAXPEERPERTRHREAD];
	float Ctlocal[NMAXPEERPERTRHREAD];
	float matUblocal[NMAXPEERPERTRHREAD];
	float matLblocal[NMAXPEERPERTRHREAD];

	float Tlocallocal[NMAXPEERPERTRHREAD]; // change
	float Tlocalprelocal[NMAXPEERPERTRHREAD]; // change
	float sum;
	float bp, MULOCAL, moy, p;
	float m, r, ub, lb, t;
	// le changement doit être partagé par tous les threads du bloc

	__shared__ float MuShared;
	__shared__ float TMoyShared;
	__shared__ float PShared;


	// constant et commun à tous les thread d'un bloc
	__shared__ float Ap1Shared;
	__shared__ float Ap3Shared;
	__shared__ float Bp3Shared;
	__shared__ float CpShared;
	__shared__ float Ap123Shared;
	__shared__ float PmaxShared;
	__shared__ float PminShared;
	__shared__ float nVoisinShared;
	__shared__ float at1Shared;
	__shared__ float at2Shared;
	__shared__ float at12Shared;
	__shared__ bool mustContinue;


	if (thIdx == 0) {
		Ap1Shared = Ap1[i];
		CpShared = Cp[i];
		Ap3Shared = Ap3[i];
		Bp3Shared = Bp3[i];
		Ap123Shared = Ap123[i];
		PmaxShared = Pmax[i];
		PminShared = Pmin[i];
		nVoisinShared = nVoisin[i];
		at1Shared = at1;
		at2Shared = at2;
		at12Shared = at1 + at2;
		MuShared = MU[i];
		TMoyShared = Tmoy[i];
		PShared = P[i];
		mustContinue = false;
	}
	int k = 0;
	__syncthreads();
	const int CoresAgentLinLocal = CoresAgentLin[i];
	const int beginLocal = CoresAgentLinLocal + thIdx;
	const int endLocal = CoresAgentLinLocal + nVoisinShared;
	float res;
	for (int j = beginLocal; j < endLocal; j += step) {
		Bt1local[k] = Bt1[j];
		Ctlocal[k] = Ct[j];
		matUblocal[k] = matUb[j];
		matLblocal[k] = matLb[j];
		//Tlocalprelocal[k] = Tlocal_pre[j];
		Tlocallocal[k] = Tlocal_pre[j];
		k = k + 1;
	}

	__shared__ float shArr[blockSize];

	//Calcul des itérations

	for (int iter = 0; iter < nStepL; iter++) {

		MULOCAL = MuShared; // tous lisent le même : broadcast !
		moy = TMoyShared;
		p = PShared;
		sum = 0;
		k = 0;
		for (int j = beginLocal; j < endLocal; j += step) {
			Tlocalprelocal[k] = Tlocallocal[k];
			m = Tlocallocal[k] - moy + p - MULOCAL;
			r = (Bt1local[k] * at1Shared + m * at2Shared - Ctlocal[k]) / (at12Shared);
			ub = matUblocal[k];
			lb = matLblocal[k];
			t = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
			Tlocallocal[k] = t;
			sum += t;
			res = (t - Tlocalprelocal[k]);
			res = (double)res * res;
			if (res > eps) {
				mustContinue = true; // pas de race condition, car l'ordre n'importe pas,
				//mais est ce que cela ne va pas physiquement bloquer ?
			}
			k = k + 1;
		}

		shArr[thIdx] = sum;
		__syncthreads();
		if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
		if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
		if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
		if (thIdx < 32) {
			warpReduce<blockSize>(shArr, thIdx);
		}
		__syncthreads();

		if (thIdx == 0) {
			moy = shArr[0] / nVoisinShared;
			TMoyShared = moy;
			bp = moy + MuShared;
			p = (Ap1Shared * bp + Ap3Shared * Bp3Shared - CpShared) / Ap123Shared;
			p = (PmaxShared - p) * (p > PmaxShared) + (PminShared - p) * (p < PminShared) + p;
			PShared = p;
			res = p - moy;
			res = (double)res * res;
			if (res > eps) {
				mustContinue = true;
			}
			MuShared = MULOCAL + moy - p;
		}
		__syncthreads();
		if (!mustContinue) {
			break;
		}
		else {
			__syncthreads();
			if (thIdx == 0) {
				mustContinue = false;
			}
		}
	}
	//Ecriture des itérations
	__syncthreads();
	k = 0;
	for (int j = beginLocal; j < endLocal; j += step) {
		Tlocal[j] = Tlocallocal[k];
		Tlocal_pre[j] = Tlocalprelocal[k];
		k = k + 1;
	}
	if (thIdx == 0) {
		Tmoy[blockIdx.x] = TMoyShared;// TMoyShared;
		P[blockIdx.x] = PShared;// PShared;
		MU[blockIdx.x] = MuShared;// MuShared;
	}

}


template <unsigned int blockSize>
__global__ void updateCp2b(float* tempN1, float* G, float* Qpart, const int nLine, const int nAgent)
{
	// un bloc par agent
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int n = blockIdx.x;
	float sum = 0;
	__shared__ float shArr[blockSize];
	for (int j = thIdx; j < nLine; j += step) {

		float t = G[j * nAgent + n] * Qpart[j * nAgent + n];
		sum += t;
	}

	shArr[thIdx] = sum;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduce<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		tempN1[n] = 2 * shArr[0];
	}

}

template <unsigned int blockSize>
__global__ void updateCp2bTrans(float* tempN1, float* G, float* Qpart, const int nLine, const int nAgent)
{
	// un bloc par agent
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int n = blockIdx.x;
	float sum = 0;
	__shared__ float shArr[blockSize];
	for (int j = thIdx; j < nLine; j += step) {

		float t = G[n * nLine + j] * Qpart[n * nLine + j];
		sum += t;
	}

	shArr[thIdx] = sum;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduce<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		tempN1[n] = 2 * shArr[0];
	}

}


template <unsigned int blockSize>
__global__ void updateCp2a(float* Cp2, float* diffKappa, float* G, const int nLine, const int nAgent) {
	// un bloc par agent
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int n = blockIdx.x;
	float sum = 0;
	__shared__ float shArr[blockSize];
	for (int j = thIdx; j < nLine; j += step) {

		float t = G[j * nAgent + n] * diffKappa[j];
		sum += t;
	}

	shArr[thIdx] = sum;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduce<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		Cp2[n] = shArr[0];
	}
}

template <unsigned int blockSize>
__global__ void updateCp2aTrans(float* Cp2, float* diffKappa, float* G, const int nLine, const int nAgent) {
	// un bloc par agent
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int n = blockIdx.x;
	float sum = 0;
	__shared__ float shArr[blockSize];
	for (int j = thIdx; j < nLine; j += step) {

		float t = G[n * nLine + j] * diffKappa[j];
		sum += t;
	}

	shArr[thIdx] = sum;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduce<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		Cp2[n] = shArr[0];
	}
}


template <unsigned int blockSize>
__global__ void updateCp2GPU(float* Cp2, float* diffKappa, float* G, float* Qpart, float* nVoisin, float rho1, const int nLine, const int nAgent) {
	// un bloc par agent
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int n = blockIdx.x;
	float sum = 0;
	__shared__ float shArr[blockSize];
	for (int j = thIdx; j < nLine; j += step) {

		float t = G[j * nAgent + n] * (diffKappa[j] + 2 * Qpart[j * nAgent + n]);
		sum += t;
	}

	shArr[thIdx] = sum;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduce<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		Cp2[n] = rho1 * nVoisin[n] * shArr[0];
	}
}

template <unsigned int blockSize>
__global__ void updateCp2GPUTrans(float* Cp2, float* diffKappa, float* G, float* Qpart, float* nVoisin, float rho1, const int nLine, const int nAgent) {
	// un bloc par agent
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int n = blockIdx.x;
	float sum = 0;
	__shared__ float shArr[blockSize];
	for (int j = thIdx; j < nLine; j += step) {
		float Gloc = G[n * nLine + j];
		float dKloc = diffKappa[j];
		float Q = Qpart[n * nLine + j];
		float t = Gloc * (dKloc + 2 * Q);
		sum += t;
	}

	shArr[thIdx] = sum;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduce<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		Cp2[n] = rho1 * nVoisin[n] * shArr[0];
	}
}

template <unsigned int blockSize>
__global__ void updateAPas(float* Apas, float* A, float* pas, int N) {

	// un bloc par ligne
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int l = blockIdx.x;
	float sum = 0;
	__shared__ float shArr[blockSize];
	for (int i = thIdx; i < N; i += step) {

		float t = A[l * N + i] * pas[i];
		sum += t;
	}

	shArr[thIdx] = sum;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduce<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		Apas[l] = shArr[0];
	}
}

template <unsigned int blockSize>
__global__ void updateAlpha(float* alpha, float* u, float* pas, float* c, float* Apas, int N, int L2) {

	// mono-block
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	//int m = blockIdx.x;
	float min = 1;
	__shared__ float sdata[blockSize];
	for (int i = thIdx; i < L2; i += step) {

		float t =  pas[i + N] < 0 ? (- u[i] / pas[i+N]) : 1;
		min = min < t ? min : t;
		t = Apas[i + N] < 0 ? -(c[i] / Apas[i]) : 1;
		min = min < t ? min : t;
	}

	sdata[thIdx] = min;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { 
		sdata[thIdx] = sdata[thIdx] < sdata[thIdx + 256] ? sdata[thIdx] : sdata[thIdx + 256];
	} __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) {
		sdata[thIdx] = sdata[thIdx] < sdata[thIdx + 128] ? sdata[thIdx] : sdata[thIdx + 128];
	} __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) {
		sdata[thIdx] = sdata[thIdx] < sdata[thIdx +  64] ? sdata[thIdx] : sdata[thIdx +  64];
	} __syncthreads(); }
	if (thIdx < 32) {
		warpReduceMin<blockSize>(sdata, thIdx);
	}

	if (thIdx == 0) {
		*alpha = 0.9 * sdata[0];
	}
}

/// PAC
template <unsigned int blockSize>
__global__ void updatePACMuAug(float* Mu, float* Muhat, float* Xhat, float* nVoisin, float* CoresAgentLin, float* indiceNu, float rho, float gamma, float phi) {
	int agent = blockIdx.x;
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	__shared__ float shArr[blockSize];
	int begin = CoresAgentLin[agent];
	int Mn = nVoisin[agent];
	int indiceMu = indiceNu[begin] + agent;
	float t = 0;
	 
	for (int i = thIdx + 1; i < Mn + 1; i += step) {
		t += Xhat[i + begin];
		float muOld = Mu[i + indiceMu];
		float mu = Muhat[i + indiceMu] + rho * gamma * (Xhat[i + begin] + Xhat[i + begin + Mn]);
		Mu[i + indiceMu] = mu;
		Muhat[i + indiceMu] = mu + phi * (mu - muOld);
	}
	shArr[thIdx] = t;
	__syncthreads();
	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduce<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		float pn = Xhat[begin];
		float muOld = Mu[indiceMu];
		float mu = Muhat[indiceMu] + rho * gamma * (pn - shArr[0]);
		Mu[indiceMu] = mu;
		Muhat[indiceMu] = mu + phi * (mu - muOld);
	}
}

template <unsigned int blockSize>
__global__ void updatePACMu(float* Mu, float* Muhat, float* Xhat, float* nVoisin, float* CoresAgentLin, float* indiceNu, float rho, float gamma, float gammahat) {
	int agent = blockIdx.x;
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	__shared__ float shArr[blockSize];
	int begin = CoresAgentLin[agent];
	int Mn = nVoisin[agent];
	int indiceMu = indiceNu[begin] + agent;
	float t = 0;

	for (int i = thIdx + 1; i < Mn + 1; i += step) {
		float Xloc = Xhat[i + begin];
		t += Xloc;
		float muOld = Mu[i + indiceMu];
		float dX = Xloc + Xhat[i + begin + Mn];
		float mu = muOld + rho * gamma * dX;
		Mu[i + indiceMu] =  mu ;
		Muhat[i + indiceMu] = mu + rho * gammahat * dX;
	}
	shArr[thIdx] = t;
	__syncthreads();
	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduce<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		float pn = Xhat[begin];
		float muOld = Mu[indiceMu];
		float dX = pn - shArr[0];
		float mu = muOld + rho * gamma * dX;
		Mu[indiceMu] = mu;
		Muhat[indiceMu] = mu + rho * gammahat * dX;
	}
}


template <unsigned int blockSize>
__global__ void calcConstraintPAC(float* tempM, float* tempM1, float* X, float* nVoisin, float* CoresAgentLin, float* CoresLinTrans, float* CoresIndiceNu) {
	int agent = blockIdx.x;
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	__shared__ float shArr[blockSize];
	int begin = CoresAgentLin[agent];
	int Mn = nVoisin[agent];
	int debutNu = CoresIndiceNu[begin];
	int debutMu = debutNu + agent;
	int fin = Mn;
	float t = 0;

	for (int i = thIdx; i < fin; i += step) {
		int lin = i + begin + 1; // 
		float Xloc = X[lin];
		int linPeer = CoresLinTrans[lin];
		t += Xloc;
		tempM[debutNu + i] = X[lin + Mn] - X[linPeer];
		tempM1[debutMu + 1 + i] = Xloc + X[lin + Mn];/**/
	}
	shArr[thIdx] = t;
	__syncthreads();
	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduce<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		float pn = X[begin];
		tempM1[debutMu] = pn - shArr[0];
	}
}


/// ADMM OPF ---------------------------------------------------------

template <unsigned int blockSize>
__global__ void ComputePFromAgentToBusGPU(float* Pb, float* Pbmin, float* Pbmax, float* CoresSoloBusAgent, float* Pn, float* Pmin, float* Pmax,
	float* CoresAgentBus, float* nAgentByBus, float* beginBus, int  Nagent, int _nBus) {

	__shared__ float shArr[blockSize];
	__shared__ float shArr2[blockSize];
	__shared__ float shArr3[blockSize];
	__shared__ float shArr4[blockSize];
	__shared__ float shArr5[blockSize];
	__shared__ float shArr6[blockSize];
	__shared__ bool mustCompute;

	int thIdx = threadIdx.x;
	int i = blockIdx.x;
	int begin = beginBus[i];
	int end = begin + nAgentByBus[i];

	if (thIdx == 0) {
		mustCompute = nAgentByBus[i] > 1;
		if (nAgentByBus[i] == 1) {
			int	agent = CoresAgentBus[begin];
			CoresSoloBusAgent[i] = agent;
			Pb[i] = Pn[agent];
			Pb[i + _nBus] = Pn[agent + Nagent];
			Pbmin[i] = Pmin[agent];
			Pbmin[i + _nBus] = Pmin[agent + Nagent];
			Pbmax[i] = Pmax[agent];
			Pbmax[i + _nBus] = Pmax[agent + Nagent];
		}
	}
	__syncthreads();
	if (mustCompute) {
		float sum = 0;
		float sum2 = 0;
		float sum3 = 0;
		float sum4 = 0;
		float sum5 = 0;
		float sum6 = 0;
		for (int k = thIdx + begin; k < end; k += blockSize) {
			int	agent = CoresAgentBus[k];
			sum += Pn[agent]; // pas coalescent du tout
			sum2 += Pn[agent + Nagent];
			sum3 += Pmin[agent]; // pas coalescent du tout
			sum4 += Pmin[agent + Nagent];
			sum5 += Pmax[agent]; // pas coalescent du tout
			sum6 += Pmax[agent + Nagent];
		}

		shArr[thIdx] = sum;
		shArr2[thIdx] = sum2;
		shArr3[thIdx] = sum3;
		shArr4[thIdx] = sum4;
		shArr5[thIdx] = sum5;
		shArr6[thIdx] = sum6;
		__syncthreads();
		for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
			if (thIdx < size) {
				shArr[thIdx] += shArr[thIdx + size];
				shArr2[thIdx] += shArr2[thIdx + size];
				shArr3[thIdx] += shArr3[thIdx + size];
				shArr4[thIdx] += shArr4[thIdx + size];
				shArr5[thIdx] += shArr5[thIdx + size];
				shArr6[thIdx] += shArr6[thIdx + size];
			}
			__syncthreads();
		}

		if (thIdx == 0) {
			Pb[i] = shArr[0];
			Pb[i + _nBus] = shArr2[0];
			Pbmin[i] = shArr3[0];
			Pbmin[i + _nBus] = shArr4[0];
			Pbmax[i] = shArr5[0];
			Pbmax[i + _nBus] = shArr6[0];
		}
	}


}

template <unsigned int blockSize>
__global__ void updateY(float* Y, float* H, float* Q, float* sizeOPFADMMbig, float* beginBus, int sizeOPFmax) {

	// un bloc par ligne
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int l = blockIdx.x;
	float sum = 0;
	__shared__ float shArr[blockSize];
	int N = sizeOPFADMMbig[l];
	int indiceBegin = beginBus[l];
	float t;
	for (int i = thIdx; i < N; i += step) {

		t = H[l * sizeOPFmax + i] * Q[indiceBegin + i];
		sum += t;
	}

	shArr[thIdx] = sum;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
	if (thIdx < 32) {
		warpReduce<blockSize>(shArr, thIdx);
	}

	if (thIdx == 0) {
		Y[l] = shArr[0];
	}
}

template <unsigned int blockSize>
__global__ void updateY(float* Y, float* H, float* Q, float* sizeOPFADMMbig, float* beginBus, int sizeOPFmax, int sizeTotal) {

	// un bloc par ligne
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	
	__shared__ float shArr[blockSize];

	for (int l = blockIdx.x; l < sizeTotal; l+= gridDim.x)
	{
		int N = sizeOPFADMMbig[l];
		int indiceBegin = beginBus[l];
		float t = 0;
		float sum = 0;
		for (int i = thIdx; i < N; i += step) {

			t = H[l * sizeOPFmax + i] * Q[indiceBegin + i];
			sum += t;
		}

		shArr[thIdx] = sum;
		__syncthreads();

		if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
		if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
		if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
		if (thIdx < 32) {
			warpReduce<blockSize>(shArr, thIdx);
		}

		if (thIdx == 0) {
			Y[l] = shArr[0];
		}
		__syncthreads();
	}

	
}



__global__ void removeLossAgent(float* _nAgentByBus, float* CoresAgentBusBegin);

__global__ void initVoltageBound(float* VlimReal, float* Vlim, float* constraintLo, float* constraintUp, float* nChild, int nBus);

__global__ void updateQ(float* Q, float* X, float* MU, float _rho, int sizeOPF);


__global__ void updateMUGPU(float* Mu, float* Y, float* X, float rho, int sizeOPF);

__global__ void initPosAgent(float* PosAgent, float* nAgentByBus, float* CoresAgentBusBegin, float* CoresAgentBus);

__global__ void initPQAgent(float* X, float* indiceBusBegin, float* CoresAgentBus, float* nAgentByBus, float* beginBusAgent, float* Pn, int nAgent);
__global__ void initPQAgentV(float* X, float* indiceBusBegin, float* CoresAgentBus, float* nAgentByBus, float* beginBusAgent, float* Pn, int nAgent);

__global__ void initPQV(float* X, float* indiceBusBegin, float* nAgentByBus, float* PnTilde, int nBus);
__global__ void initPQ(float* X, float* indiceBusBegin, float* nAgentByBus, float* PnTilde, int nBus);

template <unsigned int blockSize>
__global__ void ComputeLoss(float* X, float* Pn, float* indiceBusBegin, float* ZsRe, float* ZsIm, int lossType, int nAgent, int nBus) {
	// un seul bloc 
	int thIdx = threadIdx.x;
	int step = blockDim.x;
	int beginLoss = indiceBusBegin[nBus];
	int offset = 0;

	float sum = 0;
	float sum2 = 0;

	__shared__ float shArr[blockSize];
	__shared__ float shArr2[blockSize];

	
	if (lossType == 1) { // CURRENT 
		offset = 1;
		for (int bus = thIdx + 1; bus < nBus; bus += step) {
			int begin = indiceBusBegin[bus];
			sum  += X[begin + 2] * ZsRe[bus - 1];
			sum2 += X[begin + 2] * ZsIm[bus - 1];
		}
	}
	else { // POWER
		offset = nAgent;
		for (int agent = thIdx + 1; agent < nAgent; agent += step) {
			sum  += Pn[agent];
			sum2 += Pn[agent + nAgent];
		}
	}
	shArr[thIdx] = sum;
	shArr2[thIdx] = sum2;
	__syncthreads();

	if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; shArr2[thIdx] += shArr2[thIdx + 256];
	} __syncthreads(); }
	if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; shArr2[thIdx] += shArr2[thIdx + 128];
	} __syncthreads(); }
	if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; shArr2[thIdx] += shArr2[thIdx + 64];
	} __syncthreads(); }
	if (blockSize >= 64) { if (thIdx < 32) { shArr[thIdx] += shArr[thIdx + 32]; shArr2[thIdx] += shArr2[thIdx + 32];
	} __syncthreads(); }
	else if (thIdx < 32) {
		warpReduce<blockSize>(shArr, thIdx);
		warpReduce<blockSize>(shArr2, thIdx);
	}
	if (thIdx == 0) {
		X[beginLoss] = -shArr[0];
		X[beginLoss + offset] = -shArr2[0];
	}

}

__global__ void initDFSPQ(float* X, float* Pb, float* nChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, int nBus);
__global__ void initDFSPQ(float* X, float* nChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, int nBus);

__global__ void defineSizeBig(float* sizeOPFADMMbig, float* nChild, float* CoresBusBegin, float* sizeOPFADMM, float* CoresBusBeginBig, float* nAgentByBus);
__global__ void defineSizeBig(float* sizeOPFADMMbig, float* nChild, float* CoresBusBegin, float* sizeOPFADMM, float* CoresBusBeginBig);
__global__ void defineSizeBig(float* sizeOPFADMMbig, float* nChild, float* CoresBusBegin, float* sizeOPFADMM, float* CoresBusBeginBig, float* nAgentByBus, int lossType, int nBus, int nAgent);


__global__ void divideMultiplyByNagentByBus(float* Apt1, float* Apt2, float* PnTilde, float* PnTmin, float* PnTmax, float* nAgentByBus, float rhol, int nBus);


__global__ void updateXOPFADMM(float* X, float* Chat, float* Vbound, float* nAgentByBus, float* nChild, float* indiceBusBegin, float* CoresChatBegin,
	float* CoresAgentBusBegin, float* CoresAgentBus, float* Cost1, float* Cost2, float* Pmin, float* Pmax, float rho, int nBus, int nAgent, bool Lagrange);


__global__ void updateXOPFADMM(float* X, float* Chat, float* Vbound, float* PnTilde, float* nAgentByBus, float* nChild, float* indiceBusBegin, int nBus, bool Lagrange);
__global__ void updateXEndoMarket(float* X, float* Chat, float* Vbound, float* nChild, float* CoresChatBegin, float* indiceBusBegin, int nBus);

__global__ void updateXOPFADMMCons(float* X, float* Pn, float* Chat, float* Vbound, float* nAgentByBus, float* nChild, float* indiceBusBegin, float* CoresChatBegin,
	float* CoresAgentBusBegin, float* CoresAgentBus, float* Cost1, float* Cost2, float* Pmin, float* Pmax, float rho, int losstype, int nBus, int nAgent, bool Lagrange);
__global__ void updateXPnOPFADMMCons(float* X, float* Pn, float* Chat, float* nAgentByBus, float* indiceBusBegin, float* CoresChatBegin,
	float* CoresAgentBusBegin, float* CoresAgentBus, float* Cost1, float* Cost2, float* Pmin, float* Pmax, float rho, int losstype, int nBus, int nAgent);



__global__ void updateXPn(float* X, float* Pn, float* P, float* nVoisin, float* indiceBusBegin, float* nAgentByBus, float* CoresAgentBusBegin, float* CoresAgentBus, int lossType, int nAgent, int nBus);

__global__ void communicateX(float* X, float* nChild, float* Ancestor, float* Childs, float* indiceBusBegin, float* indiceChildBegin, float* nAgentByBus, int nBus);
__global__ void communicateX(float* X, float* nChild, float* Ancestor, float* Childs, float* indiceBusBegin, float* indiceChildBegin, int nBus);
__global__ void communicateX(float* X, float* nChild, float* Ancestor, float* Childs, float* indiceBusBegin, float* indiceChildBegin, float* nAgentByBus, float* CoresBusAgent, float* PosAgent, int Losstype, int nBus, int nAgent);


__global__ void setPnFromX(float* Pn, float* X, float* indiceBusBegin, float* CoresAgentBus, float* nAgentByBus, float* beginBusAgent, int nAgent);







/// PF ------------------------------------------------------------------------

__global__ void calcWinterCar(float* Pinter, float* Qinter, float* VoltageCar, float* Glin, float* Blin, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B);


__global__ void initECar(float* VoltageRealIm, float v0, float w0, int B);
__global__ void initECar(float* VoltageRealIm, float* E, int B);

__global__ void calcEGPU(float* E, float* VoltageRealIm, int B);
__global__ void calculYGPU(float* Y, float* E, float* Voltage, float* Blin2, float* Glin2, float* CoresLineBus, int B, int L);
__global__ void calculYGPU(float* Y, float* E, float* Blin2, float* Glin2, float* CoresLineBus, int B, int L);
__global__ void calculPhiGPU(float* Phi, float* E, float* Blin2, float* Glin2, float* CoresLineBus, int B, int L);
__global__ void setY(float* Y, float* E, float* Phi, int B, int L);


