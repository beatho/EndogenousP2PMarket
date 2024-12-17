#pragma once
#include "../head/Method.h"
#define MAX(X, Y) X * (X >= Y) + Y * (Y > X)


Method::Method()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "method constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	//timePerBlock = MatrixCPU(1, 11, 0); // Fb0, Fb1abc, Fb2, Fb3, Fb4, Fb5, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	//occurencePerBlock = MatrixCPU(1, 11, 0); //nb de fois utilisé pendant la simu
	
}

Method::~Method()
{
}



float Method::calcFc(MatrixCPU* cost1, MatrixCPU* cost2, MatrixCPU* trade, MatrixCPU* Pn, MatrixCPU* BETA, MatrixCPU* tempN1, MatrixCPU* tempNN)
{

	float fc = 0;
	tempN1->set(cost1);
	tempN1->multiply(0.5);
	tempN1->multiplyT(Pn);
	tempN1->add(cost2);
	tempN1->multiplyT(Pn);
	
	fc = fc + tempN1->sum();
	tempNN->set(trade);
	tempNN->multiplyT(BETA);
	

	fc = fc + tempNN->sum();

	return fc;

}

float Method::calcFc(MatrixCPU* cost1, MatrixCPU* cost2, MatrixCPU* Pn, MatrixCPU* tempN1)
{
	float fc = 0;
	tempN1->set(cost1);
	tempN1->multiply(0.5);
	tempN1->multiplyT(Pn);
	tempN1->add(cost2);
	tempN1->multiplyT(Pn);

	return tempN1->sum();
}

float Method::calcFc(MatrixGPU* cost1, MatrixGPU* cost2, MatrixGPU* Pn, MatrixGPU* tempN1)
{
	float fc = 0;
	tempN1->set(cost1);
	tempN1->multiply(0.5);
	tempN1->multiplyT(Pn);
	tempN1->add(cost2);
	tempN1->multiplyT(Pn);

	return tempN1->sum();
}

double Method::calcFc(MatrixGPUD* cost1, MatrixGPUD* cost2, MatrixGPUD* Pn, MatrixGPUD* tempN2)
{
	double fc = 0;
	tempN2->set(cost1);
	tempN2->multiply(0.5);
	tempN2->multiplyT(Pn);
	tempN2->add(cost2);
	tempN2->multiplyT(Pn);

	return tempN2->sum();
}

double Method::calcFc(MatrixCPUD* cost1, MatrixCPUD* cost2, MatrixCPUD* Pn, MatrixCPUD* tempN2)
{
	double fc = 0;
	tempN2->set(cost1);
	tempN2->multiply(0.5);
	tempN2->multiplyT(Pn);
	tempN2->add(cost2);
	tempN2->multiplyT(Pn);

	return tempN2->sum();
}

void Method::setBestParam(const StudyCase& cas)
{
	// not implemented
	throw std::runtime_error("not implemented for all methods");

}

float Method::calcFc(MatrixGPU* cost1, MatrixGPU* cost2, MatrixGPU* trade, MatrixGPU* Pn, MatrixGPU* BETA, MatrixGPU* tempN1, MatrixGPU* tempNN)
{
	
	tempN1->set(cost1);
	
	tempN1->multiply(0.5);
	tempN1->multiplyT(Pn);
	
	tempN1->add(cost2);
	
	tempN1->multiplyT(Pn);
	
	float fc = tempN1->sum();
	

	tempNN->set(trade);
	
	tempNN->multiplyT(BETA);
	
	fc = fc + tempNN->sum();



	//std::cout << "fc " << fc << std::endl;
	return fc;

}



void Method::resetId()
{
	_id = 0;
}


