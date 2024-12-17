#pragma once

#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}


#include <stdio.h>
#include <stdlib.h>

#include "MatrixGPU.cuh"
#include "StudyCase.h"
#include "GPUPF.cuh"
#include "Utilities.cuh"


class GPUPFdistPQ : public GPUPF
{

public:
	GPUPFdistPQ();
	~GPUPFdistPQ();

	virtual void init(const StudyCase& cas, MatrixGPU* Pn);
	bool chekcase();
	virtual void updatePQ(MatrixGPU* PQ);
	
	virtual void solve();
	virtual void calcS();
	virtual void calcW(bool end = false);
	virtual int calcVoltage();
	
	virtual void calcE();
	virtual MatrixGPU getY();
	
	virtual void setE(MatrixGPU* Enew);
	virtual void display2(bool all = true);

private:

	float v0 = 1;
	float w0 = 0;
	int LastBus = 0;


	MatrixGPU ZsRe;
	MatrixGPU ZsIm;
	MatrixGPU Yd;
	MatrixGPU VoltageRealIm;
	MatrixGPU VoltageRealImPre;
	MatrixGPU St;
	MatrixGPU Sf;
	MatrixGPU F; //f(k) = i, the sending node of the branch k
	MatrixGPU nChild;
	MatrixGPU Childs;
	MatrixGPU _indiceChildBegin;
	



};


__global__ void calculStGPU(float* St, float* Voltage, float* W0, float* Yd, int B);
__global__ void calculSGPU(float* St, float* Sf, float* Voltage, float* ZsRe, float* ZsIm, float* nChild, float* Childs, float* indiceChildBegin, int B);

__global__ void updateVoltage(float* Voltage, float* ZsRe, float* ZsIm, float* Sf, float* F, int B, int LastBus);