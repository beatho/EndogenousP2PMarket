#pragma once
#include "MethodP2P.cuh"
#include <iostream>
#include <string>



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
	void updateBt1(MatrixCPU* Bt1, MatrixCPU* trade, float rho, MatrixCPU* LAMBDA);
	void updateBt2(MatrixCPU* Bt2, MatrixCPU* Tlocal, MatrixCPU* Tmoy, MatrixCPU* P, MatrixCPU* MU);
	void updateBp1(MatrixCPU* Bp1, MatrixCPU* MU, MatrixCPU* Tmoy);
	void updateTl(MatrixCPU* Tlocal, float at1, float at2, MatrixCPU* Bt1, MatrixCPU* Bt2, MatrixCPU* Ct, MatrixCPU* matLb, MatrixCPU* matUb);
	float calcRes(MatrixCPU* Tlocal, MatrixCPU* Tlocal_pre, MatrixCPU* Tmoy, MatrixCPU* P);
	void updateP(MatrixCPU* P, MatrixCPU* Ap1, MatrixCPU* Ap2, MatrixCPU* Bp1, MatrixCPU* Cp, MatrixCPU* Pmin, MatrixCPU* Pmax);
	void updateMU(MatrixCPU* MU, MatrixCPU* Tmoy, MatrixCPU* P);
	void updateQ(MatrixCPU* Qpart, MatrixCPU* Qtot, MatrixCPU* alpha, int nAgent, int nLine);
	void display();
private:
	float _rho=0;

};


	




