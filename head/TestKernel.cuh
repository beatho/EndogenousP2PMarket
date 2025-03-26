#pragma once
#include <device_launch_parameters.h>
#include <iostream>
#include <math.h>
#include <chrono>
#include <fstream>

#include "MatrixGPU.cuh"
#include "MatrixCPU.h"
#include "StudyCase.h"


float testCalculQpart(int method);

__global__ void calculQpartLineBloc(float* Qpart, float* alpha, const int N);

__global__ void calculQpartAgentBloc(float* Qpart, float* alpha, const int L, const int N);

__global__ void calculQpartLineBlocTrans(float* Qpart, float* alpha, const int N, const int L);

__global__ void calculQpartAgentBlocTrans(float* Qpart, float* alpha, const int L, const int N);

__global__ void calculQpartLineBlocReverse(float* Qpart, float* alpha, const int N);

__global__ void calculQpartLineBlocReverseTrans(float* Qpart, float* alpha, const int N, const int nLine);

__global__ void calculQpartLineBlocReverseBis(float* Qpart, float* alpha, const int N);

__global__ void calculQpartLineBlocReverseBisTrans(float* Qpart, float* alpha, const int N, const int nLine);

float testCalculAlpha(int method);

__global__ void updateAlphaSh(float* alpha, float* G, float* Pn, const int nLine, const int nAgent);

__global__ void updateAlpha2D(float* alpha, float* G, float* Pn, const int nLine, const int nAgent);

__global__ void updateAlpha1D(float* alpha, float* G, float* Pn, const int nLine, const int nAgent);

__global__ void updateAlphaShTrans(float* alpha, float* G, float* Pn, const int nLine, const int nAgent);

__global__ void updateAlpha2DTrans(float* alpha, float* G, float* Pn, const int nLine, const int nAgent);

__global__ void updateAlpha1DTrans(float* alpha, float* G, float* Pn, const int nLine, const int nAgent);


float testCalculCpa(int method);

template <unsigned int blockSize>
__global__ void updateCp2aTest(float* Cp2, float* diffKappa, float* G, const int nLine, const int nAgent);

template <unsigned int blockSize>
__global__ void updateCp2aTestTrans(float* Cp2, float* diffKappa, float* G, const int nLine, const int nAgent);


float testCalculCpb(int method);

template <unsigned int blockSize>
__global__ void updateCp2bTest(float* tempN1, float* G, float* Qpart, const int nLine, const int nAgent);

template <unsigned int blockSize>
__global__ void updateCp2bTestTrans(float* tempN1, float* G, float* Qpart, const int nLine, const int nAgent);

float testCalculQtot(int method);

__global__ void updateQtotTest(float* Qtot, float* Qpart, float* alpha, const int nLine, const int nAgent);

__global__ void updateQtotTestTrans(float* Qtot, float* Qpart, float* alpha, const int nLine);

template <unsigned int blockSize>
__device__ void warpReduceTest(volatile float* sdata, unsigned int tid);



float testCalculCp(int method);


template <unsigned int blockSize>
__global__ void updateCp2Test(float* Cp2, float* diffKappa, float* G, float* Qpart, float* nVoisin, float rho1, const int nLine, const int nAgent);

template <unsigned int blockSize>
__global__ void updateCp2TestTrans(float* Cp2, float* diffKappa, float* G, float* Qpart, float* nVoisin, float rho1, const int nLine, const int nAgent);


float testCalculResX(int method);

__global__ void updateResXTest(float* res, float* Kappa1, float* Kappa2, float* KappaPre1, float* KappaPre2, const int nLine);



float testCalculLAMBDABt1(int method);
__global__ void updateLAMBDAGPUTest(float* LAMBDA, float* trade, float rho, float* CoresLinTrans, int const N);

__global__ void updateBt1GPUTest(float* Bt1, float* tradeLin, float rho, float* LAMBDA, float* CoresLinTrans, int const N);


__global__ void updateLAMBDABt1GPUTest(float* Bt1, float* LAMBDA, float* tradeLin, float rho,  float* CoresLinTrans, int const N);






float testCalculChat(int method, int blockSize, int repartition);
template <unsigned int blockSize>
__global__ void updateChatGPUTest(float* Chat, float* Y, float* MU, float* nChild, float* Ancestor, float* posChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, float _rho, int nBus);
__global__ void updateBpt2Test(float* Bpt2, float* Chat, float* nAgentByBus, int nBus);

template <unsigned int blockSize>
__global__ void updateChatBpt2Test(float* Chat, float* Bpt2, float* Y, float* MU, float* nChild, float* Ancestor, float* posChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, float* nAgentByBus, float _rho, int nBus);
__global__ void updateChatBpt2OneDoAll(float* Chat, float* Bpt2, float* Y, float* MU, float* nChild, float* Ancestor, float* posChild, float* Childs, float* indiceBusBegin, float* indiceChildBegin, float* nAgentByBus, float _rho, int nBus);



float testCalculPnShared(int method, int blockSize, int repartition);

template <unsigned int blockSize>
__global__ void updatePnPGPUSharedResidualTest(float* Pn, float* PnPre, float* PnMoy, float* PnTilde, float* MUL, float* nAgentByBus, float _rhol, float* Ap2, float* Cp, float* Pmin,
	float* Pmax, float* Apt1, float* Apt2, float* Bpt2, float* CoresSoloBusAgent, float* CoresBusAgent, float* CoresBusAgentBegin, float eps, int nIterLMax, int nAgent, int nBus);


template <unsigned int blockSize>
__global__ void updatePnPGPUSharedResidualSameThreadTest(float* Pn, float* PnPre, float* PnMoy, float* PnTilde, float* MUL, float* nAgentByBus, float _rhol, float* Ap2, float* Cp, float* Pmin,
	float* Pmax, float* Apt1, float* Apt2, float* Bpt2, float* CoresSoloBusAgent, float* CoresBusAgent, float* CoresBusAgentBegin, float eps, int nIterLMax, int nAgent, int nBus);


// --------------------------------- Voltage GS -------------------------------------------------------


float testCalculVGS(int method);

__global__ void calculVoltStep2Test(float* VoltageRealIm, float* RMGgrid, float* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, float* CoresTrans, int B);

__global__ void calculVoltStep2bisTest(float* VoltageRealIm, float* RMGgrid, float* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, float* CoresTrans, int B, int BL2);

template <unsigned int blockSize>
__global__ void calculVoltStep1Test(float* VoltageRealIm, float* W0, float* Rgrid, float* Xgrid, float* RMGgrid, float* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B);

template <unsigned int blockSize>
__global__  void calculVoltOneStepTest(float* VoltageRealIm, float* W0, float* Rgrid, float* Xgrid, float* RMGgrid, float* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B);