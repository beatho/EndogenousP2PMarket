#pragma once 
#include <device_launch_parameters.h>
#include "Utilities.cuh"
#include <cuda_runtime.h>
#include "MatrixGPU.cuh"
#include "Method.h"


class MethodP2PGPU : public Method
{
public:
    MethodP2PGPU();
	virtual ~MethodP2PGPU();
	virtual void setParam(float param) = 0;
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas) = 0;
	//virtual void setTau(float tau) = 0;

	void updateLAMBDA(MatrixGPU* LAMBDA, MatrixGPU* trade, float rho, MatrixGPU* tempNN);
	void updateKappa();
    float calcFc();
	

	float updateRes(int iter);
	virtual float updateResEndo(int iter);

    float calcRes(); // local residuals

	virtual void updatePn();

	void initLinForm(const StudyCase& cas);
	void initSize(const StudyCase& cas);
	void initSimParam(const Simparam& sim); // AFTER initSize
	void initDCEndoGrid(const StudyCase& cas); // init with G transposed
	void initDCEndoMarket();
	void initP2PMarket();
	void initCaseParam(const Simparam& sim, const StudyCase& cas);

	virtual void setResult(Simparam* result, bool casAC);

	virtual void solveWithMinPower(Simparam* result, const Simparam& sim, const StudyCase& cas);
	
	virtual void display();
	

protected:
	// ne change pas avec P0
	clock_t tMarket;

	float _mu = 40.0f;
	float _mu1 = 40.0f;
	float _mul = 40.0f;
	float _tau = 2.0f;
	
	float _rho = 0;
	float _rhol = 0;
	float _rhog = 0;

	int _iterGlobal = 0;
	int _iterG = 0;
	int _iterL = 0;
	int _iterIntern = 0;

	int _stepG = 1;
	int _stepL = 1;
	int _stepIntern = 1;

	float _epsL = 0;
	float _epsG = 0;
	float _epsX = 0;
	float _epsIntern = 0;
	
	bool isAC = false;
	int _nAgent = 0;
	int _nAgentTrue = 0;

	int _nTrade = 0;
	int _nTradeP = 0;
	int _nTradeQ = 0;


	int _blockSize = 256;
	int _numBlocksN = 1;
	int _numBlocksM = 1;
	int _numBlocksL = 1;
	int _numBlocksNL = 1;

	
	float _at1 = 0.0f;
	float _at2 = 0.0f;

	MatrixGPU tempNN; // Matrix temporaire pour aider les calculs
	MatrixGPU tempN1; // plut�t que de re-allouer de la m�moire � chaque utilisation
	
	MatrixGPU tempL1; 
	MatrixGPU tempL2; 
	

	MatrixGPU Tlocal;
	MatrixGPU P; // moyenne des trades
	MatrixGPU Pn; // somme des trades

	MatrixGPU a;
	MatrixGPU Ap1;
	MatrixGPU Ap2; // Mn^2 * (Ap2a + Ap2b)
	MatrixGPU Ap2a; // a
	MatrixGPU Ap2b; // 2 * _rho1 * sum(G^2) 
	MatrixGPU Ap12;
	MatrixGPU Ap3; // rho1*Mn^2 when used
	MatrixGPU Ap123;

	MatrixGPU Bt1;
	MatrixGPU Bt2;
	MatrixGPU Bp1;
	MatrixGPU Ct;
	MatrixGPU matUb;
	
	MatrixGPU nVoisin;
	MatrixGPU tradeLin;
	MatrixGPU Tmoy;
	MatrixGPU LAMBDALin;
	MatrixGPU Tlocal_pre;
	MatrixGPU MU;

	MatrixGPU CoresMatLin;
	MatrixGPU CoresLinAgent;
	MatrixGPU CoresAgentLin;
	MatrixGPU CoresLinVoisin;
	MatrixGPU CoresLinTrans;
	
	// Matrices kept on CPU
	MatrixCPU nVoisinCPU;
	MatrixCPU LAMBDA;
	MatrixCPU trade;
	MatrixCPU resF;
	MatrixCPU resX;



	// change avec P0
	MatrixGPU b;
	MatrixGPU matLb;
	MatrixGPU Cp;
	MatrixGPU Pmin;
	MatrixGPU Pmax;
	
	// Pour le r�seau
	int _nLine = 0;
	int _nBus = 0;
	float _rho1 = 0.0f;

	MatrixGPU Kappa1;
	MatrixGPU Kappa2;
	MatrixGPU Kappa1_pre;
	MatrixGPU Kappa2_pre;
	MatrixGPU Cp1;
	MatrixGPU Cp2;
	MatrixGPU Qpart;
	MatrixGPU Qtot;
	MatrixGPU alpha;

	MatrixGPU G;
	MatrixGPU GTrans;
	MatrixGPU G2; // G^t.*G^t
	MatrixGPU lLimit;


};





// Regroupement des appels kernels pour ne pas devoir d�finir plusieurs fois le m�me avec des noms differents
__global__ void updateLAMBDAGPU(float* LAMBDALin, float* tradeLin, float rho, float* CoresLinTrans, int const N);
__global__ void updateBt1GPU(float* Bt1, float* tradeLin, float rho, float* LAMBDA, float* CoresLinTrans, int const N);
__global__ void updateLAMBDABt1GPU(float* Bt1, float* LAMBDALin, float* tradeLin, float rho, float* CoresLinTrans, int const N);



__global__ void updateDiffGPU(float* tempNN, float* Tlocal, float* CoresLinTrans, int const N);
__global__ void updateResKappa(float* result, float* Kappa1, float* Kappa2, float* Kappapre1, float* Kappapre2, float ratio, int const L);
__global__ void selectResidual(float* res, unsigned int id1, unsigned int id2, unsigned int id3, float* output);
__global__ void updateKappaGPU(float* Kappa1, float* Kappa2, float* Llimit, float* Qtot, int nLine);
__global__ void diffKappa(float* tempL1, float* Kappa1, float* Kappa2, const int nLine);
__global__ void updateCpOld(float* Cp, float* Cp1, float* Cp2, float* tempN1, float* nVoisin, const float rho1, const int nAgent);
__global__ void updateCp(float* Cp, float* Cp1, float* Cp2, const int nAgent);

__global__ void updateQpart(float* Qpart, float* alpha, const int nAgent);
__global__ void updateQpartTrans(float* Qpart, float* alpha, const int N, const int nLine);

__global__ void updateQtot(float* Qtot, float* Qpart, float* alpha, const int nLine, const int nAgent);
__global__ void updateQtotTrans(float* Qtot, float* Qpart, float* alpha, const int nLine);


__global__ void updatePnGPU(float* Pn, float* Tmoy, float* nVoisin, const int nAgent);

__global__ void updateAlpha(float* alpha, float* G, float* Pn, const int nLine, const int nAgent);
__global__ void updateAlphaTrans(float* alpha, float* GTrans, float* Pn, const int nLine, const int nAgent);

__global__ void updateResX(float* res, float* Kappa1, float* Kappa2, float* KappaPre1, float* KappaPre2, const int nLine);




// Methode pour le dual ascent


/*
__global__ void updateQt(float* qt, float* Pso, float* Pn, float* etaSO, float rho1, int N);
__global__ void updateUAiq(float* UAiq, float* u, float* Aiq, int N, int size);
__global__ void updateRu(float* Ru, float* U, float* g, float epsi, int N, int L2);
__global__ void updateV(float* v, float* pas, float* alpha, int offset);*/


__global__ void updatePI(float* PI, float* c, float mu, float valMin, int L);
__global__ void updatePso(float* Pso, float* pas, float* alpha, int N);
__global__ void updateU(float* U, float* pas, float* alpha, int N, int L);

 
__global__ void updateEtaPBp3(float* Bp3, float* etaP, float* nVoisin, float* Pso, float* Pn, float rho, const int nAgent);

