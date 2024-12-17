#pragma once

#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}


#include <stdio.h>
#include <stdlib.h>
#include "MatrixCPU.h"
#include "StudyCase.h"
#include <chrono>

class CPUPF
{
	
public:
	CPUPF();
	~CPUPF();
	virtual void init(const StudyCase& cas, MatrixCPU* PQ);
	virtual void init(const StudyCase& cas, MatrixCPU* Pn, MatrixCPUD* PnD, bool useDouble);
	virtual void solve();
	virtual void updatePQ(MatrixCPU* PQ);
	virtual void calcW(bool end=false);
	virtual void calcJac();
	virtual void calcPhi();
	virtual void calcJacPhiE();
	virtual int calcVoltage();
	virtual void calcE();
	MatrixCPU* calcG();
	MatrixCPU* calcY();


	virtual void setE(MatrixCPU* Enew);
	virtual void setE(MatrixCPUD* Enew);
	void setW(MatrixCPU* Wnew);

	float getPloss();
	float getQloss();
	float getRes();
	int getIter();
	float getP0();
	float getQ0();
	float getTime();
	int getConv();
	MatrixCPU getE();
	virtual MatrixCPU getY(); // (angle, amplitude, flow)

	void display();
	void display2(bool all=true);
	void saveTimeBlock(std::string fileName);

protected:

	bool _useDouble = false;
	int Nagent = 0;
	int Nbus = 0;
	int Nline = 0;
	int Nconstraint = 0;
	int status = 0; // 0 not lauch, -1 failure, 1 converge, 2 not converge
	int B2 = 0; // 2*Nbus =0
	int N2 = 0; // 2*NAgent =0
	double V0 = 1; // voltage at the ref bus
	double theta0 = 0; // voltage angle at the ref bus
	int iterM = 10;
	double epsPF = 0.0005;
	double err = 0;
	int iter = 0;
	clock_t time =0;
	std::string _name = "CPU";

	MatrixCPU I_aug; // agent to bus
	MatrixCPU I; // agent to bus

	MatrixCPU W0;
	MatrixCPU W;
	MatrixCPU dW;
	MatrixCPU E;
	
	MatrixCPU dE;
	MatrixCPU Jac;
	MatrixCPU JacInv;
	MatrixCPU Bgrid;
	MatrixCPU Ggrid;

	MatrixCPU A;
	MatrixCPU P;

	MatrixCPU G;
	MatrixCPU Phi;
	MatrixCPU Y; // [E Phi]

	MatrixCPU tempLN2; // JacPhiE * Jacinv * Iaug
	MatrixCPU tempB2N2; // Jacinv * Iaug
	MatrixCPU JacPhiE; // L * B2
	MatrixCPU CoresLineBus;
	MatrixCPU Ggrid2Bgrid2; // Ggrid^2+Bgrid^2

	// si double
	MatrixCPUD ED;
	MatrixCPUD dED;

	MatrixCPUD BgridD;
	MatrixCPUD GgridD;
	

	
	MatrixCPUD JacD;
	MatrixCPUD JacInvD;
	
	MatrixCPUD W0D;
	MatrixCPUD WD;
	MatrixCPUD dWD;

	MatrixCPUD AD;
	MatrixCPUD PD;


	// si calcul en lineaire
	int BL2 = 0; // Nbus+2*Nline

	MatrixCPU BgridLin;
	MatrixCPU GgridLin;


	MatrixCPUD BgridLinD;
	MatrixCPUD GgridLinD;

	MatrixCPU CoresVoiLin;
	MatrixCPU CoresBusLin;
	MatrixCPU nLines;

	// instrumentation
	MatrixCPU timePerBlock = MatrixCPU(1, 7); // Fb0, Fb1abc, Fb2, Fb3abc, Fb4, Fb5, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	MatrixCPU occurencePerBlock = MatrixCPU(1, 7);; //nb de fois utilisé pendant la simu

	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;


};


