#pragma once

#include <device_launch_parameters.h>
#include "MethodOPF.h"
#include "MatrixGPU.cuh"




class MethodOPFGPU : public MethodOPF
{
public:
	virtual void solveConsensus(float eps, MatrixGPU* PSO) = 0;
	virtual void updateConsensus(MatrixGPU* Pmarket) = 0;
	virtual void solveConsensus(float eps, MatrixCPU* PSO) = 0;
	virtual void initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO) = 0;
	virtual void updateConsensus(MatrixCPU* Pmarket) = 0;
	
	float calcFc(MatrixGPU* cost1, MatrixGPU* cost2, MatrixGPU* Pn, MatrixGPU* tempN1)
	{
		tempN1->set(cost1);
		tempN1->multiply(0.5);
		tempN1->multiplyT(Pn);
		tempN1->add(cost2);
		tempN1->multiplyT(Pn);
		return tempN1->sum();
	};
	double calcFc(MatrixGPUD* cost1, MatrixGPUD* cost2, MatrixGPUD* Pn, MatrixGPUD* tempN2)
	{
		tempN2->set(cost1);
		tempN2->multiply(0.5);
		tempN2->multiplyT(Pn);
		tempN2->add(cost2);
		tempN2->multiplyT(Pn);
		return tempN2->sum();
	};

protected:
	int _blockSize = 256;
	int _blockSizeSmall = 32;
};

