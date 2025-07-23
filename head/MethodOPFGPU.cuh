#pragma once

#include <device_launch_parameters.h>
#include "MethodOPF.h"
#include "MatrixGPU.cuh"




class MethodOPFGPU : public Method
{
public:
	MethodOPFGPU() : Method(){
		timePerBlock = MatrixCPU(1, 9, 0); // Fb0, Fb1 , Fb2, Fb3, Fb5, Fb6 Fb0'
		occurencePerBlock = MatrixCPU(1, 9, 0); //nb de fois utilis� pendant la simu
	};
	virtual void solveConsensus(float eps, MatrixCPU* PSO){};
	virtual void initConsensus(const Simparam& sim, const StudyCase& cas, float rhoSO){};
	virtual void updateConsensus(MatrixCPU* Pmarket){};

	virtual void solveConsensus(float eps, MatrixGPU* PSO) = 0;
	virtual void updateConsensus(MatrixGPU* Pmarket) = 0;
		
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

	virtual MatrixCPU getPb()=0;
	virtual MatrixCPU getPhi()=0;
	virtual MatrixCPU getE()=0;
	virtual float getPLoss()
	{
		float _Ploss = 0;
		switch (losstype)
		{
		case LossType::POWER:
			_Ploss = -Pn.sum(1, _nAgent);
			break;
		case LossType::CURRENT:
			for (int i = 1; i < _nBus; i++) 
			{
				int begin = _indiceBusBeginCPU.get(i, 0);
				_Ploss -= Y.get(begin + indli, 0, false) * ZsRe.get(i - 1, 0, false);
			}
			break;
		}
		return _Ploss;
	}

	virtual float getQLoss()
	{
		float _Ploss = 0;
		switch (losstype)
		{
		case LossType::POWER:
			_Ploss = -Pn.sum(_nAgent+1, 2*_nAgent);
			break;
		case LossType::CURRENT:
			for (int i = 1; i < _nBus; i++) 
			{
				int begin = _indiceBusBeginCPU.get(i, 0);
				_Ploss -= Y.get(begin + indli, 0, false) * ZsIm.get(i - 1, 0, false);
			}
			break;
		}
		return _Ploss;
	}

protected:
	int _blockSize = 256;
	int _blockSizeSmall = 32;

	LossType losstype = LossType::POWER;
	bool Lagrange = false; // true  false
	bool isCurrentLimited = false;

	float _mu = 100;
	float _tau = 2;

	int _nAgent = 0;
	int _nBus = 0;

	MatrixCPU _indiceBusBeginCPU; // size : nBus

	MatrixGPU Pn;
	MatrixGPU Y;// (Pi, Qi, li, vi, pn..., qn..., vai, Pci ..., Qci... , lci...) !!!!!
	MatrixGPU ZsRe;
	MatrixGPU ZsIm;
	
	const int indPi  = 0;
	const int indQi  = 1;
	const int indli  = 2;
	const int indvi  = 3;
	const int indvai = 4;
	const int indpi  = 5;
	const int indqi  = 6;
	const int indChatpi  = 4;
	const int indChatqi  = 5;
};

