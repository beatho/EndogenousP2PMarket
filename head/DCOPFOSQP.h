#pragma once


#include "MatrixCPU.h"
#include "Simparam.h"
#include "StudyCase.h"
#include "Method.h"


//#include "C:\Program Files\OSQP\osqp\include\osqp.h"
#include <osqp.h>
#include <iostream>
#include <string>
#include <chrono>

/// <summary>
/// But : Ici juste faire un DC-OPf, pas de problème de trade on travaille qu'avec les trade 
/// Change : resolution centralisée avec un OSQP
/// </summary>
/// 
class DCOPFOSQP : public Method 
{
public:
	DCOPFOSQP();
	~DCOPFOSQP();


	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	void init(const Simparam& sim, const StudyCase& cas);
	void updateP0(const StudyCase& cas);
	std::string NAME ="DCOPFOSQP";
	
	
	void display();

private:

	
	
								 
	// ne change pas avec P0
	int _nAgent = 0;
	float _rhog = 0;
	MatrixCPU tempN1;

	MatrixCPU Pn; // somme des trades
	MatrixCPU a;



	// change avec P0
	MatrixCPU b;
	MatrixCPU Pmin;
	MatrixCPU Pmax;
	

	// pour osqp

	OSQPSettings* settings;
	OSQPWorkspace* work;
#ifndef OSQPGPU
	OSQPData* data;
#endif
	
	MatrixCPU Aosqp; // l <A*P<u
	MatrixCPU Hosqp;
	MatrixCPU qTosqp;
	MatrixCPU PsoCPU;
	c_float* Q = nullptr;
	c_float* xResult = nullptr;
	c_float* L = nullptr;
	c_float* U = nullptr;

	int nConst =0;

	// Pour le réseau
	int _nLine=0;
	int _nBus=0;

	MatrixCPU G;
	MatrixCPU lLimit;

};





