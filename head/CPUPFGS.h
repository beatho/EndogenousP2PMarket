#pragma once

#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}


#include <stdio.h>
#include <stdlib.h>
#include "MatrixCPU.h"
#include "StudyCase.h"
#include "CPUPF.h"


class CPUPFGS : public CPUPF
{
	
public:
	CPUPFGS();
	~CPUPFGS();
	virtual void init(const StudyCase& cas, MatrixCPU* Pn);
	virtual void init(const StudyCase& cas, MatrixCPU* Pn, MatrixCPUD* PnD, bool useDouble);
	virtual void updatePQ(MatrixCPU* PQ);
	virtual void calcW(bool end=false);
	virtual int calcVoltage();
	virtual void calcE();

	virtual void setE(MatrixCPU* Enew);
	virtual void setE(MatrixCPUD* Enew);

private:

	double v0 = 1;
	double w0 = 0;

	MatrixCPU Rgrid;
	MatrixCPU Xgrid;
	MatrixCPU RMGgrid;
	MatrixCPU RPGgrid;
	MatrixCPU VectorResult;
	MatrixCPU VoltageRealIm;


	MatrixCPUD RgridD;
	MatrixCPUD XgridD;
	MatrixCPUD RMGgridD;
	MatrixCPUD RPGgridD;
	MatrixCPUD VectorResultD;
	MatrixCPUD VoltageRealImD;

};


