#pragma once

#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}


#include <stdio.h>
#include <stdlib.h>
#include "MatrixCPU.h"
#include "StudyCase.h"
#include "GPUPF.cuh"


class GPUPFGS : public GPUPF
{
	
public:
	GPUPFGS();
	~GPUPFGS();

	virtual void init(const StudyCase& cas, MatrixGPU* Pn, MatrixGPUD* PnD=nullptr, bool useDouble=false);
	virtual int calcVoltage();
	virtual void calcE();
	virtual void calcW(bool end=false);

	virtual void setE(MatrixGPU* Enew);
	virtual void setE(MatrixGPUD* Enew);

private:

	double v0 = 1;
	double w0 = 0;

	MatrixGPU Rgrid;
	MatrixGPU Xgrid;
	MatrixGPU RMGgrid;
	MatrixGPU RPGgrid;
	MatrixGPU VectorResult;
	MatrixGPU VoltageRealIm;
	
	
	MatrixGPU CoresTrans;

	MatrixGPUD RgridD;
	MatrixGPUD XgridD;
	MatrixGPUD RMGgridD;
	MatrixGPUD RPGgridD;
	MatrixGPUD VectorResultD;
	MatrixGPUD VoltageRealImD;
	

};



__global__ void initEDCar(double* VoltageRealImD, double v0, double w0, int B);
__global__ void initEDCar(double* VoltageRealImD, double* ED, int B);

__global__ void initRX(float* Rgrid, float* Xgrid, float* RMGgrid, float* RPGgrid, float* Glin, float* Blin, float* CoresBusLin, float* nLines);
__global__ void initRXD(double* RgridD, double* XgridD, double* RMGgridD, double* RPGgridD, double* GlinD, double* BlinD, float* CoresBusLin, float* nLines);

template <unsigned int blockSize>
__global__ void calculVoltStep1(float* VoltageRealImD, float* W0, float* Rgrid, float* Xgrid, float* RMGgrid, float* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B);

template <unsigned int blockSize>
__global__ void calculVolDtStep1(double* VoltageRealImD, double* W0, double* Rgrid, double* Xgrid, double* RMGgrid, double* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B);

__global__ void calculVoltStep2(float* VoltageRealIm, float* RMGgrid, float* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, float* CoresTrans, int B);

__global__ void calculVoltDStep2(double* VoltageRealIm, double* RMGgrid, double* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, float* CoresTrans, int B);

__global__ void calculVoltStep2bis(float* VoltageRealIm, float* RMGgrid, float* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, float* CoresTrans, int B, int BL2);

__global__ void calculVoltDStep2bis(double* VoltageRealIm, double* RMGgrid, double* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, float* CoresTrans, int B, int BL2);

__global__ void calcWinterCarD(double* Pinter, double* Qinter, double* E, double* Glin, double* Blin, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B);

__global__ void calcEGPUD(double* ED, double* VoltageRealImD, int B);



