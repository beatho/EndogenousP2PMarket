#pragma once
#include "MethodP2P.h"
#include <iostream>
#include <string>
#include <chrono>


class ADMMConst : public MethodP2P
{
public:
	ADMMConst();
	ADMMConst(float rho);
	virtual ~ADMMConst();
	void setParam(float rho);
	void setTau(float tau);
	virtual void solve(Simparam* result, const Simparam& sim, const StudyCase& cas);
	virtual void updateP0(const StudyCase& cas);
	virtual void init(const Simparam& sim, const StudyCase& cas);
	std::string NAME ="ADMMConst";
	void updateGlobalProb();
	void updateLocalProb();
	void updateBt1();
	void updateBt2();
	void updateBp1();
	void updateTl();
	void updateP();
	void updateMU();
	void updateQ();
	void display();
private:
	
};


	




