#pragma once

#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}


#include <stdio.h>
#include <stdlib.h>
#include "MatrixCPU.h"
#include "StudyCase.h"
#include "CPUPF.h"


class CPUPFdist : public CPUPF
{

public:
	CPUPFdist();
	~CPUPFdist();

	virtual void init(const StudyCase& cas, MatrixCPU* Pn);
	bool chekcase();
	virtual void updatePQ(MatrixCPU* PQ);
	
	virtual void solve();
	virtual void calcJ();
	virtual void calcW(bool end = false);
	virtual int calcVoltage();
	
	virtual void calcE();
	
	
	virtual void setE(MatrixCPU* Enew);
	virtual void display2(bool all = true);

private:

	float v0 = 1;
	float w0 = 0;


	MatrixCPU ZsRe;
	MatrixCPU ZsIm;
	MatrixCPU Yd;
	MatrixCPU VoltageRealIm;
	MatrixCPU VoltageRealImPre;
	MatrixCPU Jb;
	MatrixCPU F; //f(k) = i, the sending node of the branch k

};


