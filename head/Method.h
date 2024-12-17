#pragma once
#include "StudyCase.h"
#include "MatrixCPU.h"
#include "MatrixGPU.cuh"
#include "Simparam.h"
//#include "kernelFunction.cuh"


#include <iostream>
#include <string>



enum LossType { POWER, CURRENT };



class Method 
{
public:
	Method();
	virtual ~Method();
	static constexpr float POWERLIMIT = 100000.0f;
	float calcFc(MatrixGPU* cost1, MatrixGPU* cost2, MatrixGPU* trade, MatrixGPU* Pn, MatrixGPU* BETA, MatrixGPU* tempN1, MatrixGPU* tempNN);
	float calcFc(MatrixCPU* cost1, MatrixCPU* cost2, MatrixCPU* trade, MatrixCPU* Pn, MatrixCPU* BETA, MatrixCPU* tempN1, MatrixCPU* tempNN);
	float calcFc(MatrixCPU* cost1, MatrixCPU* cost2, MatrixCPU* Pn, MatrixCPU* tempN1);
	float calcFc(MatrixGPU* cost1, MatrixGPU* cost2, MatrixGPU* Pn, MatrixGPU* tempN2);
	double calcFc(MatrixGPUD* cost1, MatrixGPUD* cost2, MatrixGPUD* Pn, MatrixGPUD* tempN2);
	double calcFc(MatrixCPUD* cost1, MatrixCPUD* cost2, MatrixCPUD* Pn, MatrixCPUD* tempN2);
	
	
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas)=0;
	virtual void updateP0(const StudyCase& cas) = 0;
	virtual void init(const Simparam& sim, const StudyCase& cas) = 0;
	virtual void display() = 0;
	
	virtual void setBestParam(const StudyCase& cas);
	std::string _name = "default";
	void resetId();

protected:
	int _id = 0;
	float _ratioEps = 1;
	MatrixCPU timePerBlock; // Fb0, Fb1abc, Fb2, Fb3abc, Fb4, Fb5, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	MatrixCPU occurencePerBlock; //nb de fois utilisé pendant la simu

};

