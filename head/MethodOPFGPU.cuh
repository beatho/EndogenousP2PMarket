#pragma once

#include "MethodOPF.h"




class MethodOPFGPU : public MethodOPF
{
public:
	virtual void solveConsensus(float eps, MatrixGPU* PSO) = 0;
	virtual void updateConsensus(MatrixGPU* Pmarket) = 0;
protected:
	int _blockSize = 256;
	int _blockSizeSmall = 32;
};

