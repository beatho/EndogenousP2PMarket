#pragma once
#include "../head/Method.h"
 


Method::Method()
{
#if DEBUG_CONSTRUCTOR
	std::cout << "method constructor" << std::endl;
#endif // DEBUG_CONSTRUCTOR
	//timePerBlock = MatrixCPU(1, 11, 0); // Fb0, Fb1abc, Fb2, Fb3, Fb4, Fb5, Fb0'
	// si les sous ensemble ne sont pas accessible, tout est dans le premier.
	//occurencePerBlock = MatrixCPU(1, 11, 0); //nb de fois utilis� pendant la simu
	
}

Method::~Method()
{
}



void Method::setBestParam(const StudyCase& cas)
{
	// not implemented
	throw std::runtime_error("not implemented for all methods");

}


void Method::display(){
	std::cout << "Method's name : " << _name << std::endl; 
}

void Method::resetId()
{
	_id = 0;
}


