#pragma once
#include "StudyCase.h"
#include "MatrixCPU.h"

#include "Simparam.h"
#include "Utilities.h"
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
	
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas)=0;
	virtual void updateP0(const StudyCase& cas) = 0;
	virtual void init(const Simparam& sim, const StudyCase& cas) = 0;
	virtual void display();
	
	virtual void setBestParam(const StudyCase& cas);
	std::string _name = "default";
	void resetId();

protected:
	int _id = 0;
	float _ratioEps = 1;

	MatrixCPU timePerBlock; // Fb0, Fb1abc, Fb2, Fb3abc, Fb4, Fb5, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	MatrixCPU occurencePerBlock; //nb de fois utilis� pendant la simu

};

