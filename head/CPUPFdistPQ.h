#pragma once

#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}


#include <stdio.h>
#include <stdlib.h>
#include "MatrixCPU.h"
#include "StudyCase.h"
#include "CPUPF.h"


class CPUPFdistPQ : public CPUPF
{

public:
	CPUPFdistPQ();
	~CPUPFdistPQ();

	virtual void init(const StudyCase& cas, MatrixCPU* Pn , MatrixCPUD* PnD = nullptr, bool useDouble = false);
	bool chekcase();
	virtual void updatePQ(MatrixCPU* PQ);
	
	virtual void solve();
	virtual void calcS();
	virtual void calcW(bool end = false);
	virtual int calcVoltage();
	
	virtual void calcE();
	virtual MatrixCPU getY();
	
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
	MatrixCPU St;
	MatrixCPU Sf;
	MatrixCPU F; //f(k) = i, the sending node of the branch k

};


