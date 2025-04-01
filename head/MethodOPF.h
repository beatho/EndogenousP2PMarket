#pragma once

#include "Method.h"




class MethodOPF : public Method
{
public:
	virtual void solveConsensus(float eps, MatrixCPU* PSO) = 0;
	virtual void initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO) = 0;
	virtual void updateConsensus(MatrixCPU* Pmarket) = 0;

	virtual float getPLoss() = 0;
	virtual float getQLoss() = 0;
	virtual MatrixCPU getPb() = 0;
	virtual MatrixCPU getPhi() = 0;
	virtual MatrixCPU getE() = 0;
	double calcFc(MatrixCPUD* cost1, MatrixCPUD* cost2, MatrixCPUD* Pn, MatrixCPUD* tempN1)
	{
		tempN1->set(cost1);
		tempN1->multiply(0.5);
		tempN1->multiplyT(Pn);
		tempN1->add(cost2);
		tempN1->multiplyT(Pn);
		return tempN1->sum();
	};
	float calcFc(MatrixCPU* cost1, MatrixCPU* cost2, MatrixCPU* Pn, MatrixCPU* tempN1)
	{
		tempN1->set(cost1);
		tempN1->multiply(0.5);
		tempN1->multiplyT(Pn);
		tempN1->add(cost2);
		tempN1->multiplyT(Pn);

		return tempN1->sum();
	};
	virtual int feasiblePoint() {
		return 0;
	};
	void setLagrange(bool lagrange) {
		Lagrange = lagrange;
	};
	float getMurhoVar() {
		return _mu;
	};
	float getTaurhoVar() {
		return _tau;
	};
protected:
	float _mu = 100;
	bool Lagrange = false; // true  false
	float _tau = 2;

};

