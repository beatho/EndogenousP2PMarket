#include "../head/TestMatrixCPU.h"



int testMatrix()
{
	int n = 1;
	if (!testMConstru()) return n;
	n++;
	if (!testMConstru2()) return n;
	n++;
	if (!testMConstru3()) return n;
	n++;
	if (!testMSet1()) return n;
	n++;
	if (!testMSet2()) return n;
	n++;
	if (!testMGet1()) return n;
	n++;
	if (!testMGet2()) return n;
	n++;
	if (!testMEquality1()) return n;
	n++;
	if (!testMEquality2()) return n;
	n++; // 10
	if (!testMEquality3()) return n;
	n++;
	if (!testMAdd1()) return n;
	n++;
	if (!testMAdd2()) return n;
	n++;
	if (!testMAdd3()) return n;
	n++;
	if (!testMAdd4()) return n;
	n++;
	if (!testMAdd5()) return n;
	n++;
	if (!testMAdd6()) return n;
	n++;
	if (!testMAddVect1()) return n;
	n++;
	if (!testMAddVect2()) return n;
	n++;
	if (!testMAddTrans1()) return n;
	n++;//20
	if (!testMAddTrans2()) return n;
	n++;
	if (!testMSubstract1()) return n;
	n++; 
	if (!testMSubstract2()) return n;
	n++;
	if (!testMSubstract3()) return n;
	n++;
	if (!testMSubstract4()) return n;
	n++;
	if (!testMSubstractVect1()) return n;
	n++;
	if (!testMSubstractVect2()) return n;
	n++;
	if (!testMSubstractTrans1()) return n;
	n++;
	if (!testMSubstractTrans2()) return n;
	n++;
	if (!testMMultiply1()) return n;
	n++;// 30
	if (!testMMultiply2()) return n;
	n++;
	if (!testMMultiply3()) return n;
	n++; 
	if (!testMMultiplyT1()) return n;
	n++;
	if (!testMMultiplyT2()) return n;
	n++;
	if (!testMMultiplyT3()) return n;
	n++;
	if (!testMMultiplyT4()) return n;
	n++;
	if (!testMDivide1()) return n;
	n++;
	if (!testMDivide2()) return n;
	n++;
	if (!testMDivide3()) return n;
	n++;
	if (!testMDivide4()) return n;
	n++;// 40
	if (!testMDivide5()) return n;
	n++;
	if (!testMDivideGJ1()) return n;
	n++; 
	if (!testMDivideGJ2()) return n;
	n++;
	if (!testMmoy1()) return n;
	n++; 
	if (!testMmoy2()) return n;
	n++;
	if (!testMmoy3()) return n;
	n++;
	if (!testMmoy4()) return n;
	n++;
	if (!testMSum1()) return n;
	n++;
	if (!testMSum2()) return n;
	n++;
	if (!testMSum3()) return n;
	n++;// 50
	if (!testMSum4()) return n;
	n++; 
	if (!testMDistance1()) return n;
	n++;
	if (!testMDistance2()) return n;
	n++;
	if (!testMProject1()) return n;
	n++; 
	if (!testMProject2()) return n;
	n++; 
	if (!testMProjectPos()) return n;
	n++;
	if (!testMProjectNeg()) return n;
	n++;
	if (!testMMax()) return n;
	n++;
	if (!testMSwap()) return n;
	n++;// 60
	if (!testMSetFromFile()) return n;
	n++; 
	if (!testMMinAbs()) return n;
	n++;
	if (!testMSetRand()) return n;
	n++; 
	if (!testMSort()) return n;
	n++; 
	#ifdef EIGEN
		if (!testMSolveSys()) return n;
		n++;
		if (!testMDivideEigen()) return n;
		n++;
		if (!testMSetFromEigen()) return n;
		n++;
		if (!testMToEigen()) return n;
		n++;//10
	#endif
	#ifdef OSQP
		if (!testMToCSC1()) return n;
		n++;
		if (!testMToCSC2()) return n;
		n++;
		if (!testMCSC()) return n;
		n++;
	#endif
	
	return 0;
}

bool testMConstru()
{
	std::cout << "default constructor" << std::endl;
	MatrixCPU mempty;
	mempty.display();
	
	std::cout << "null constructor" << std::endl;
	MatrixCPU mnull(3,3);
	mnull.display();
	
	std::cout << "1 constructor" << std::endl;
	MatrixCPU mones(3, 2, 1);
	mones.display();
	
	std::cout << "copy constructor" << std::endl;
	MatrixCPU mones2(mones);
	mones2.display();
	
	std::cout << "empty constructor" << std::endl;
	MatrixCPU vide(3, 0);
	MatrixCPU vide2(0, 3);
	MatrixCPU vide3(0, 0);
	

	return true;
}
bool testMConstru2() {
	std::cout << "2 times constructor" << std::endl;
	MatrixCPU mones3;
	mones3 = MatrixCPU(3, 2, 1);
	mones3.display();
	MatrixCPU mones(3, 2, 1);
	return mones3.isEqual(&mones);
}
bool testMConstru3() {
	std::cout << "2 times constructor bis" << std::endl;
	MatrixCPU* mones3;
	MatrixCPU* mones;
	mones = new MatrixCPU(2, 3, 1);
	mones3 = new MatrixCPU(2, 3, 0);
	mones3->set(mones);
	mones3->display();
	
	bool result = mones3->isEqual(mones);

	DELETEB(mones);
	DELETEB(mones3);

	return result;
}

bool testMSet1()
{
	float value = 4.5;
	int i = 1;
	int j = 2;
	int n = 3;
	MatrixCPU mnull(n, n,1);
	mnull.set(i, j, value); 
	return  (mnull.get(i, j) == value);
}
bool testMSet2()
{
	int n = 3;
	float value = 4;
	int i = 1;
	int j = n;
	MatrixCPU mnull(n, n);
	try
	{
		mnull.set(i, j, value); 
	}
	catch (std::out_of_range& )
	{
		return true;
	}
	return false;
}

bool testMSetFromFile()
{
	MatrixCPU m(2, 3, 1);
	m.set(1, 2, 2);
	MatrixCPU m1(2, 3);
	MatrixCPU m11(2, 3);
	MatrixCPU m2(1, 3);
	MatrixCPU m3(3, 3);
	std::string fileName1 = "test/test.txt";
	std::string fileName2 = "test/test.csv";

	m1.setFromFile(fileName1);
	m1.display();
	m11.setFromFile(fileName2);
	m11.display();

	m2.setFromFile(fileName1);
	m2.display();
	m2.setFromFile(fileName2);
	m2.display();
	m3.setFromFile(fileName1);
	m3.display();
	m3.setFromFile(fileName1);
	m3.display();

	return (m1.isEqual(&m11) && m.isEqual(&m1));
}

bool testMSetRand()
{
	int column = 10;
	int line = 10;
	MatrixCPU mat(line, column);
	MatrixCPU mat2(line, column);
	int nValue = 4;
	int divide = nValue / 0.001;
	mat.setRand(nValue, divide);
	mat.display();
	mat2.setRand(nValue, divide);
	mat2.display();
	return true;
}


bool testMGet1()
{
	float value = 4.0;
	int column = 3;
	int line = 2;
	MatrixCPU mones(line, column, value);

	for (int i = 0;i < line;i++) {
		for (int j = 0;j < column;j++) {
			if (value != mones.get(i, j)) return false;
		}
	}
	return true;
}
bool testMGet2()
{
	int n = 3;
	float value = 4;
	int i = 1;
	int j = n;
	MatrixCPU m(n, n, value);

	try
	{
		m.get(i, j); 
	}
	catch (std::out_of_range& )
	{
		return true;
	}
	return false;
}

bool testMEquality1()
{
	int line = 3; 
	int column = line+1;
	float value = 1.5;
	MatrixCPU m1(line, column, value);
	MatrixCPU m2(line, column,value+1);	
	return !m1.isEqual(&m2);
}
bool testMEquality2()
{
	int line = 3;
	int column = line + 1;
	float value = 1.5;
	MatrixCPU m1(line, column, value);
	MatrixCPU m2(m1);

	return ((m2.isEqual(&m1))&& (m1.isEqual(&m2)));

}
bool testMEquality3()
{
	int line = 3;
	int column = line + 1;
	float value = 1.5;
	MatrixCPU m1(line, column, value);
	MatrixCPU m3(column, line,value);
	
	try
	{
		m3.isEqual(&m1);
	}
	catch (const std::invalid_argument& )
	{
		try
		{
			m1.isEqual(&m3);
		}
		catch (const std::invalid_argument& )
		{
			return true;
		}
		return false;
	}
	return false;

}

/*
  void add(float c);
  void add(MatrixCPU* m1, MatrixCPU* m2); // m = m1 +m2;
  void add(MatrixCPU* m1, float c); // m = m1 + c;
  void add(MatrixCPU* m1);  // m = m + m1;
*/
bool testMAdd1()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;
	float value3 = value + value2;
	
	
	MatrixCPU m1(line, column, value);
	MatrixCPU m2(line, column, value2);
	MatrixCPU m3(line, column, value3); 
	MatrixCPU temp(line, column);

	temp.add(&m1, &m2);
	return temp.isEqual(&m3);

}
bool testMAdd2() {
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;
	float value3 = value + value2;
	MatrixCPU m1(line, column, value);
	MatrixCPU m3(line, column, value3); 
	MatrixCPU temp(line, column);
	temp.add(&m1, value2);
	
	return temp.isEqual(&m3);
}
bool testMAdd3()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;
	float value3 = value + value2;
	
	MatrixCPU m1(line, column, value);
	MatrixCPU m2(line, column, value2);
	MatrixCPU m3(line, column, value3); 
	
	m1.add(&m2);
	return m1.isEqual(&m3);

}
bool testMAdd4()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = -1;
	float value3 = value + value2;
	MatrixCPU m1(line, column, value);
	MatrixCPU m3(line, column, value3);
	
	m1.add(value2);
	
	return m1.isEqual(&m3);
}
bool testMAdd5()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;


	MatrixCPU m1(line, column, value);
	MatrixCPU m2(line, column, value2);
	MatrixCPU m3(column, line); 

	try
	{
		m3.add(&m1, &m2);
	}
	catch (std::invalid_argument& )
	{
		return true;
	}
	return false;
}
bool testMAdd6()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;


	MatrixCPU m1(line, column, value);
	MatrixCPU m2(column, line, value2);


	try
	{
		m1.add(&m2);
	}
	catch (std::invalid_argument& )
	{
		return true;
	}
	return false;
}

//void addVector(MatrixCPU* v);
bool testMAddVect1()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;
	MatrixCPU m1(line, column, value);
	MatrixCPU m4(column, line, value2);
	MatrixCPU vect(1, line, value2);
	
	try {
		m1.addVector(&m4); 
	}
	catch (std::invalid_argument& ) {
		try
		{
			m1.addVector(&vect); 
		}
		catch (std::invalid_argument&)
		{
			return true;
		}
		return false;
	}
	return false;
	
}
bool testMAddVect2()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;
	MatrixCPU m1(line, column, value);
	MatrixCPU m11(line, column, value);
	MatrixCPU m2(line, column, value2);
	MatrixCPU vect(1, column, value2 - value);
	MatrixCPU vect1(line, 1, value - value2);

	m1.addVector(&vect);
	if (!m1.isEqual(&m2)) return false;
	m2.addVector(&vect1);
	if (!m2.isEqual(&m11)) return false;

	return true;
}
//void addTrans(MatrixCPU* m1); // m = m + tm1;

bool testMAddTrans1()
{
	int line = 2;
	int column = line + 1;
	int i = column - 2;
	int j = line - 1;
	float value = 1.5;
	float value2 = 1;
	float value3 = value + value2;
	float value4 = 4;
	float value5 = value + value4;

	MatrixCPU m1(line, column, value);
	MatrixCPU m2(column, line, value2);
	m2.set(i, j, value4);
	MatrixCPU m3(line, column, value3); 
	m3.set(j, i, value5);
	

	m1.addTrans(&m2);
	return m1.isEqual(&m3);

}
bool testMAddTrans2() {
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;
	
	MatrixCPU m1(line, column, value);
	MatrixCPU m2(line, column, value2);
	
	try
	{
		m1.addTrans(&m2);
	}
	catch (std::invalid_argument& )
	{
		return true;
	}

	return false;
}

/*
void subtract(MatrixCPU* m);
void subtract(MatrixCPU* m1, MatrixCPU* m2);
*/

bool testMSubstract1()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;
	float value3 = value - value2;

	MatrixCPU m1(line, column, value);
	MatrixCPU m2(line, column, value2);
	MatrixCPU m3(line, column, value3); 
	MatrixCPU temp(line, column);

	temp.subtract(&m1, &m2);
	return temp.isEqual(&m3);

}
bool testMSubstract2()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;
	float value3 = value - value2;

	MatrixCPU m1(line, column, value);
	MatrixCPU m2(line, column, value2);
	MatrixCPU m3(line, column, value3); 

	m1.subtract(&m2);
	return m1.isEqual(&m3);

}
bool testMSubstract3()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;

	MatrixCPU m1(line, column, value);
	MatrixCPU m2(line, column, value2);
	MatrixCPU m3(column, line);

	try
	{
		m3.subtract(&m1, &m2);
	}
	catch (std::invalid_argument& )
	{
		return true;
	}
	return false;
}
bool testMSubstract4()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;


	MatrixCPU m1(line, column, value);
	MatrixCPU m2(column, line, value2);


	try
	{
		m1.subtract(&m2);
	}
	catch (std::invalid_argument& )
	{
		return true;
	}
	return false;
}

//void subtractVector(MatrixCPU* v);
bool testMSubstractVect1()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;
	MatrixCPU m1(line, column, value);
	MatrixCPU m4(column, line, value2);
	MatrixCPU vect(1, line, value2);

	try {
		m1.subtractVector(&m4);
	}
	catch (std::invalid_argument& ) {
		try
		{
			m1.subtractVector(&vect);
		}
		catch (std::invalid_argument& )
		{
			return true;
		}
		return false;
	}
	return false;

}
bool testMSubstractVect2()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;
	MatrixCPU m1(line, column, value);
	MatrixCPU m11(line, column, value);
	MatrixCPU m2(line, column, value2);
	MatrixCPU vect(1, column, value - value2);
	MatrixCPU vect1(line, 1, value2 - value);

	m1.subtractVector(&vect);
	if (!m1.isEqual(&m2)) return false;
	m2.subtractVector(&vect1);
	if (!m2.isEqual(&m11)) return false;

	return true;
}
//void subtractTrans(MatrixCPU* m);
bool testMSubstractTrans1()
{
	int line = 2;
	int column = line + 1;
	int i = column - 2;
	int j = line - 1;
	float value = 1.5;
	float value2 = 1;
	float value3 = value - value2;
	float value4 = 4;
	float value5 = value - value4;

	MatrixCPU m1(line, column, value);
	MatrixCPU m2(column, line, value2);
	m2.set(i, j, value4);
	MatrixCPU m3(line, column, value3); 
	m3.set(j, i, value5);


	m1.subtractTrans(&m2);
	return m1.isEqual(&m3);

}
bool testMSubstractTrans2() {
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;

	MatrixCPU m1(line, column, value);
	MatrixCPU m2(line, column, value2);

	try
	{
		m1.subtractTrans(&m2);
	}
	catch (std::invalid_argument& )
	{
		return true;
	}

	return false;
}

/*
void multiply(MatrixCPU* m1, MatrixCPU* m2);
void multiply(float c);
*/
bool testMMultiply1()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;
	float value3 = column*value * value2;


	MatrixCPU m1(line, column, value);
	MatrixCPU m2(column, line, value2);
	MatrixCPU m3(line, line, value3); 
	MatrixCPU temp(line, line);

	temp.multiply(&m1, &m2);
	return temp.isEqual(&m3);
}
bool testMMultiply2()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;
	float value3 = column * value * value2;


	MatrixCPU m1(line, column, value);
	MatrixCPU m2(line, column, value2);
	MatrixCPU m3(line, line, value3); 
	MatrixCPU temp(column, column);
	MatrixCPU temp2(line, column);

	try {
		temp.multiplyT(&m1, &m2);
	}
	catch (std::invalid_argument& ) {
		try
		{
			temp2.multiplyT(&m1,&m3);
		}
		catch (std::invalid_argument& )
		{
			return true;
		}
		return false;
	}
	return false;
}
bool testMMultiply3()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;
	float value3 = value * value2;
	MatrixCPU m1(line, column, value);
	MatrixCPU m3(line, column, value3);
	m1.multiply(value2);

	return m1.isEqual(&m3);
}
/*
void multiplyT(MatrixCPU* m);
void multiplyT(MatrixCPU* m1, MatrixCPU* m2);
*/
bool testMMultiplyT1()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;
	float value3 = value * value2;

	MatrixCPU m1(line, column, value);
	MatrixCPU m2(line, column, value2);
	MatrixCPU m3(line, column, value3); 
	MatrixCPU temp(line, column);

	temp.multiplyT(&m1, &m2);
	return temp.isEqual(&m3);
}
bool testMMultiplyT2()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;
	float value3 = value * value2;

	MatrixCPU m1(line, column, value);
	MatrixCPU m2(line, column, value2);
	MatrixCPU m3(line, column, value3); 

	m1.multiplyT(&m2);
	return m1.isEqual(&m3);
}
bool testMMultiplyT3()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;
	float value3 = value * value2;

	MatrixCPU m1(line, column, value);
	MatrixCPU m2(column, line, value2);
	MatrixCPU m3(line, line, value3); 
	MatrixCPU temp(column, column);
	MatrixCPU temp2(line, line);

	try {
		temp.multiplyT(&m1, &m2);
	}
	catch (std::invalid_argument& ) {
		try
		{
			temp.multiplyT(&m1, &m3);
		}
		catch (std::invalid_argument& )
		{
			return true;
		}
		return false;
	}
	return false;
}
bool testMMultiplyT4()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;


	MatrixCPU m1(line, column, value);
	MatrixCPU m2(line, line, value2); 
	

	try {
		m1.multiplyT(&m2);
	}
	catch (std::invalid_argument& ) {
		
		return true;
	}
	return false;
}
/*
	void divide(float c);
	void divideT(MatrixCPU* m);
*/
bool testMDivide1()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1.2;
	float value3 = value /value2;


	MatrixCPU m1(line, column, value);
	MatrixCPU m2(line, column, value2);
	MatrixCPU m3(line, column, value3); 
	

	m1.divideT(&m2);
	return m1.isEqual(&m3);
}
bool testMDivide2()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;


	MatrixCPU m1(line, column, value);
	MatrixCPU m2(column, line, value2);
	

	try
	{
		m1.divideT(&m2);
	}
	catch (std::invalid_argument& )
	{
		return true;
	}
	return false;
}

bool testMDivide3()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 1;


	MatrixCPU m1(line, column, value);
	MatrixCPU m2(line, column, value2);
	m2.set(line - 1, column - 1, 0);

	try
	{
		m1.divideT(&m2);
	}
	catch (std::domain_error& )
	{
		return true;
	}
	return false;
}
bool testMDivide4()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = -1.2;
	float value3 = value /value2;
	MatrixCPU m1(line, column, value);
	MatrixCPU m3(line, column, value3);

	m1.divide(value2);

	return m1.isEqual(&m3);
}
bool testMDivide5()
{
	int line = 2;
	int column = line + 1;
	float value = 1.5;
	float value2 = 0;


	MatrixCPU m1(line, column, value);


	try
	{
		m1.divide(value2);
	}
	catch (std::domain_error&)
	{
		return true;
	}
	return false;
}

bool testMDivideGJ1()
{
	int n = 3;
	MatrixCPU ident(n, n);
	ident.setEyes(1);

	MatrixCPU invert(n, n);
	invert.invertGaussJordan(&ident);

	ident.display();
	invert.display();

	if (!ident.isEqual(&invert)) return false;


	MatrixCPU m1(n, n);
	m1.set(0, 0, 2);
	m1.set(0, 1, -1);
	m1.set(1, 1, -1);
	m1.set(1, 2, 2);
	m1.set(2, 0, -1);
	m1.set(2, 1, 2);
	m1.set(2, 2, 1);

	MatrixCPU m2(n, n);
	MatrixCPU m22(n, n);
	m22.set(0, 0, 5.0 / 8);
	m22.set(0, 1, -1.0 / 8);
	m22.set(0, 2, 1.0 / 4);
	m22.set(1, 0, 1.0 / 4);
	m22.set(1, 1, -1.0 / 4);
	m22.set(1, 2, 1.0 / 2);
	m22.set(2, 0, 1.0 / 8);
	m22.set(2, 1, 3.0 / 8);
	m22.set(2, 2, 1.0 / 4);

	m2.invertGaussJordan(&m1);

	m2.display();
	m22.display();

	MatrixCPU temp(n, n);
	temp.multiply(&m2, &m1);
	temp.display();

	return m2.isEqual(&m22) && temp.isEqual(&ident);
}

bool testMDivideGJ2()
{
	int n = 3;
	MatrixCPU temp1(n, n + 1);
	MatrixCPU temp2(n, n);
	temp2.set(0, 0, 1);
	temp2.set(2, 2, 1);


	try
	{
		temp1.invertGaussJordan(&temp1); // not square matrix
	}
	catch (std::invalid_argument&)
	{
		try
		{
			temp1.invertGaussJordan(&temp2); // not same size
		}
		catch (std::invalid_argument&)
		{
			try
			{
				temp2.invertGaussJordan(&temp2); // not invertible
			}
			catch (std::invalid_argument&)
			{
				return true;
			}
			return false;
		}
		return false;
	}
	return false;
}


// void Moy(MatrixCPU* m, MatrixCPU* nb, int sens=0); 
bool testMmoy1()
{
	int row = 3;
	int column = 4;
	float value = 1.5;
	MatrixCPU m1(row, column, 0);
	m1.set(0, 0, value);
	m1.set(0, 3, value);
	m1.set(1, 0, value);
	m1.set(1, 1, value);
	m1.set(1, 3, value);
	m1.set(2, 2, value);
	MatrixCPU nb1col(1, column, 0);
	nb1col.set(0, 0, 2);
	nb1col.set(0, 1, 1);
	nb1col.set(0, 2, 1);
	nb1col.set(0, 3, 2);
	MatrixCPU temp1(1, column);
	MatrixCPU m3(1, column, value);
	temp1.Moy(&m1, &nb1col, 1);
	return temp1.isEqual(&m3);
}
bool testMmoy2()
{
	int row = 3;
	int column = 4;
	float value = 1.5;
	MatrixCPU m1(row, column, 0);
	m1.set(0, 0, value);
	m1.set(0, 3, value);
	m1.set(1, 0, value);
	m1.set(1, 1, value);
	m1.set(1, 3, value);
	m1.set(2, 2, value);
	MatrixCPU nb2li(row, 1, 2);
	MatrixCPU temp2(row, 1);
	MatrixCPU m5(row, 1, value);
	m5.set(1, 0, value * 3 / 2);
	m5.set(2, 0, value / 2);

	temp2.Moy(&m1, &nb2li, 0);
	return temp2.isEqual(&m5);	
}
bool testMmoy3()
{
	int row = 3;
	int column = 4;
	float value = 1.5;
	MatrixCPU m1(row, column, 0);
	m1.set(0, 0, value);
	m1.set(0, 3, value);
	m1.set(1, 0, value);
	m1.set(1, 1, value);
	m1.set(1, 3, value);
	m1.set(2, 2, value);
	MatrixCPU nb2li(row, 1, 2);
	MatrixCPU temp2(row, 1);
	try
	{
		temp2.Moy(&m1, &nb2li, 1); 
	}
	catch (std::invalid_argument& )
	{
		return true;
	}
	return false;
}
bool testMmoy4()
{
	int row = 3;
	int column = 4;
	float value = 1.5;
	MatrixCPU m1(row, column, 0);
	m1.set(0, 0, value);
	m1.set(0, 3, value);
	m1.set(1, 0, value);
	m1.set(1, 1, value);
	m1.set(1, 3, value);
	m1.set(2, 2, value);
	MatrixCPU nb2li(row, 1, 2);
	MatrixCPU temp2(row, 1);
	MatrixCPU temp1(1, column);
	try
	{
		temp1.Moy(&m1, &nb2li, 0);
	}
	catch (std::invalid_argument& )
	{
		try
		{
			temp2.Moy(&m1, &m1, 0); 
		}
		catch (std::invalid_argument& )
		{
			return true;
		}
		return false;
	}
	return false;
}

/*
 float sum() const;
void sum(MatrixCPU* m, int sens = 0);
*/
bool testMSum1()
{
	int row = 3;
	int column = 4;
	float value = 1.5;
	float value2 = row * column * value;
	MatrixCPU m1(row, column, value);
	return value2 == m1.sum();
}
bool testMSum2()
{
	int row = 3;
	int column = 4;
	float value = 1.5;
	float value2 = row * column * value;
	MatrixCPU m1(row, column, value);
	MatrixCPU temp1(1, column);
	MatrixCPU temp2(row, 1);
	MatrixCPU m3(1, column, value * row);
	MatrixCPU m4(row, 1, value * column);

	temp1.sum(&m1, 1);
	if (!temp1.isEqual(&m3)) return false;
	temp2.sum(&m1, 0);
	if (!temp2.isEqual(&m4)) return false;

	return true;
}
bool testMSum3()
{
	int row = 3;
	int column = 4;
	float value = 1.5;
	float value2 = row * column * value;
	MatrixCPU m1(row, column, value);
	MatrixCPU m2(column, column, value);
	MatrixCPU temp1(1, column);
	MatrixCPU temp2(row, 1);
	MatrixCPU m3(1, column, value * row);
	MatrixCPU m4(row, 1, value * column);
	try
	{
		temp2.sum(&m2, 0); 
	}
	catch (std::invalid_argument& )
	{
		try
		{
			m1.sum(&m2, 0); 
		}
		catch (std::invalid_argument& )
		{
			return true;
		}
		return false;
	}
	return false;
}
bool testMSum4()
{
	int row = 3;
	int column = 4;
	float value = 1.5;
	float value2 = row * column * value;
	MatrixCPU m1(row, column, value);
	MatrixCPU m2(column, column, value);
	MatrixCPU temp1(1, column);
	MatrixCPU temp2(row, 1);
	try
	{
		temp1.sum(&m1, 0); 
	}
	catch (std::invalid_argument& )
	{
		try
		{
			temp2.sum(&m1, 1);
		}
		catch (std::invalid_argument& )
		{
			return true;
		}
		return false;
	}
	return false;
}

/*
* float distance2(MatrixCPU* m1) const; 
  float distance2() const; 
*/
bool testMDistance1()
{
	int row = 3;
	int column = 4;
	float value = 1.5;
	float value2 = 4;
	float d1 = sqrtf(row * column * (value * value));
	float d2 = sqrtf(row * column * ((value - value2) * (value - value2)));
	MatrixCPU m1(row, column, value);
	MatrixCPU m2(row, column, value2);

	if (m1.distance2(&m2) != d2) return false;
	if (m1.distance2() != d1) return false;

	return true;
}
bool testMDistance2()
{
	int row = 3;
	int column = 4;
	float value = 1.5;
	float value2 = 4;
	MatrixCPU m1(row, column, value);
	MatrixCPU m3(column, row, value2);

	try
	{
		m1.distance2(&m3); 
	}
	catch (std::invalid_argument& )
	{
		return true;
	}
	return false;



	return true;
}

//void project(MatrixCPU* Lb, MatrixCPU* Ub);
bool testMProject1()
{
	int row = 3;
	int column = 4;
	float value = 1.5;
	float value2 = 0.5;
	float value3 = 3;
	MatrixCPU m1(row, column, value);
	m1.set(0, 2, value2 - 1);
	m1.set(1, 1, value3 + 1);
	MatrixCPU m2(row, column, value);
	m2.set(0, 2, value2);
	m2.set(1, 1, value3);
	MatrixCPU m22(m2);
	MatrixCPU lb(row, column, value2);
	MatrixCPU ub(row, column, value3);





	m1.project(&lb, &ub);
	if (!m1.isEqual(&m22)) return false;
	m2.project(&lb, &ub);
	if (!m2.isEqual(&m22)) return false;
	m2.project(&lb, &lb);
	if (!m2.isEqual(&lb)) return false;

	return true;
}
bool testMProject2()
{
	int row = 3;
	int column = 4;
	float value = 1.5;
	float value2 = 0.5;
	float value3 = 3;
	MatrixCPU m1(row, column, value);
	MatrixCPU m3(column, row, value);
	MatrixCPU lb(row, column, value2);
	MatrixCPU ub(row, column, value3);

	try
	{
		m1.project(&m3, &ub); 
	}
	catch (std::invalid_argument& )
	{
		try
		{
			m1.project(&ub, &lb); // ub>lb
		}
		catch (std::invalid_argument&)
		{
			return true;
		}
		return false;
	}
	return false;

}

bool testMProjectPos()
{
	int row = 2;
	int column = 2;
	float value = 1.5;
	float value2 = -0.5;
	
	MatrixCPU m1(row, column,value2);
	m1.set(0, 1, 0);
	m1.set(1, 1, value);
	MatrixCPU m2(row, column);
	m2.set(1, 1, value);


	m1.projectPos();


	return m1.isEqual(&m2);
}

bool testMProjectNeg()
{
	int row = 2;
	int column = 2;
	float value = -1.5;
	float value2 = 0.5;

	MatrixCPU m1(row, column, value2);
	m1.set(0, 1, 0);
	m1.set(1, 1, value);
	MatrixCPU m2(row, column);
	m2.set(1, 1, value);


	m1.projectNeg();


	return m1.isEqual(&m2);
}

bool testMMax()
{
	int row = 3;
	int column = 4;
	float value = 1.5;
	float value2 = -4;
	float res = 0;
	MatrixCPU m(row, column, value);
	MatrixCPU m2;
	m.set(row - 1, column - 1, value2);

	if (m2.max2() != 0) return false;
	
	res = m.max2();
	return (res == fabs(value2));

}


bool testMSwap()
{
	int line = 4;
	int col = 5;
	float value1 = 2;
	float value2 = -4;
	MatrixCPU m1(line, col, value1);
	MatrixCPU m11(line, col, value1);
	MatrixCPU m2(line, col, value2);
	MatrixCPU m22(line, col, value2);

	
	m1.swap(&m2);
	


	return ((m1.isEqual(&m22)) && (m2.isEqual(&m11)));
}

bool testMMinAbs()
{
	int line = 4;
	int col = 5;
	float value1 = 2;
	float value2 = -4;
	MatrixCPU indice(1, 2);
	MatrixCPU matrix(line, col, 0);

	
	matrix.set(line - 1, col - 1, value1);
	matrix.set(line - 1, col - 2, value2);
	matrix.display();
	float min1 = 0;
	float min2 = value1;

	float min11 = matrix.minAbs(0, line, 0, col, &indice);
	float min22 = matrix.minAbs(0, line, 0, col, &indice, false);

	

	return (min1==min11 && min2==min22);
}

bool testMSort()
{
	int column = 5;
	int line = 4;
	int dim = 1;
	MatrixCPU test(line, column,-1);
	for (int i = 0; i < column; i++) {
		test.set(0, i, i);
	}
	test.set(1, 0, 1);
	test.set(1, 1, 4);
	test.set(1, 2, 2);
	test.set(1, 3, 5);
	test.set(1, 4, 2);

	test.display();
	test.sort(1, 1);
	test.display();

	for (int i = 0; i < column-1; i++) {
		if (test.get(dim, i) > test.get(dim, i + 1)) {
			return false;
		}
	}

	return true;
}





#ifdef EIGEN
	bool testMSetFromEigen()
	{
		Eigen::MatrixXd eMat(2, 2);
		eMat(0, 0) = 1;
		eMat(0, 1) = 2;
		eMat(1, 0) = -1;
		eMat(1, 1) = -4;

		


		MatrixCPU m1(2, 2);
		MatrixCPU m11(2, 2);
		m11.set(0, 0, 1);
		m11.set(0, 1, 2);
		m11.set(1, 0, -1);
		m11.set(1, 1, -4);
		m1.set(&eMat);

		

		return m1.isEqual(&m11);
	}

	bool testMToEigen()
	{
		Eigen::MatrixXd eMat(2, 2);
		MatrixCPU m1(2, 2);
		MatrixCPU m11(2, 2);
		m11.set(0, 0, 1);
		m11.set(0, 1, 2);
		m11.set(1, 0, -1);
		m11.set(1, 1, -4);
		m11.toEigenMatrix(&eMat);
		m1.set(&eMat);


		return m1.isEqual(&m11);
	}

	bool testMDivideEigen()
	{
		int n = 3;
		MatrixCPU ident(n, n);
		ident.setEyes(1);

		MatrixCPU invert(n, n);
		invert.invertEigen(&ident);

		ident.display();
		invert.display();

		if (!ident.isEqual(&invert)) return false;


		MatrixCPU m1(n, n);
		m1.set(0, 0, 2);
		m1.set(0, 1, -1);
		m1.set(1, 1, -1);
		m1.set(1, 2, 2);
		m1.set(2, 0, -1);
		m1.set(2, 1, 2);
		m1.set(2, 2, 1);

		MatrixCPU m2(n, n);
		MatrixCPU m22(n, n);
		m22.set(0, 0, 5.0 / 8);
		m22.set(0, 1, -1.0 / 8);
		m22.set(0, 2, 1.0 / 4);
		m22.set(1, 0, 1.0 / 4);
		m22.set(1, 1, -1.0 / 4);
		m22.set(1, 2, 1.0 / 2);
		m22.set(2, 0, 1.0 / 8);
		m22.set(2, 1, 3.0 / 8);
		m22.set(2, 2, 1.0 / 4);

		m2.invertGaussJordan(&m1);



		MatrixCPU temp(n, n);
		temp.multiply(&m2, &m1);


		return m2.isEqual(&m22) && temp.isEqual(&ident);
	}

	bool testMSolveSys()
	{
		Eigen::Matrix3f M;
		Eigen::Vector3f b;
		M << 1, 2, 3, 4, 5, 6, 7, 8, 10;
		b << 3, 3, 4;
		std::cout << "Here is the matrix A:\n" << M << std::endl;
		std::cout << "Here is the vector b:\n" << b << std::endl;
		Eigen::Vector3f x = M.colPivHouseholderQr().solve(b);
		std::cout << "The solution is:\n" << x << std::endl;

		MatrixCPU Mm(3, 3);
		MatrixCPU bm(3, 1);
		MatrixCPU xm(3, 1);
		MatrixCPU xm2(3, 1);
		MatrixCPU Am(3, 3);
		MatrixCPU P(4, 1);
		for (int i = 0; i < 3; i++) {
			bm.set(i, 0, b(i));
			for (int j = 0; j < 3; j++) {
				Mm.set(i, j, M(i,j));
			}
		}
		bm.display();
		xm.solveSysEigen(&Mm, &bm);
		Mm.LUPFactorization(&Am, &P);
		xm2.solveSys(&Am, &P, &bm);
		

		for (int i = 0; i < 3; i++) {
			if (abs(xm.get(i, 0) - x(i))>0.00001) {
				std::cout << xm.get(i, 0) << " " << x(i) << " "<< abs(xm.get(i, 0) - x(i)) << " "<< (abs(xm.get(i, 0) - x(i)) > 0.00001) << std::endl;
				return false;
			}
		}
		for (int i = 0; i < 3; i++) {
			if (abs(xm2.get(i, 0) - x(i)) > 0.00001) {
				std::cout << xm2.get(i, 0) << " " << x(i) << " " << abs(xm2.get(i, 0) - x(i)) << " " << (abs(xm2.get(i, 0) - x(i)) > 0.00001) << std::endl;
				return false;
			}
		}
		return true;
	}
#endif

#ifdef OSQP
	
bool testMToCSC1()
{
	MatrixCPU P(2, 2);
	P.set(0, 0, 4);
	P.set(0, 1, 1);
	P.set(1, 0, 1);
	P.set(1, 1, 2);
	c_int P_n = P.getNNullHalf();
	c_float* Pdata = new c_float[P_n];
	c_int* Pidx = new c_int[P_n];
	c_int* Pptr = new c_int[P_n];
	P.toCSCHalf(Pdata, Pidx, Pptr);
	c_float P_x[3] = { 4.0, 1.0, 2.0, };
	c_int P_nnz = 3;
	c_int P_i[3] = { 0, 0, 1, };
	c_int P_p[3] = { 0, 1, 3, };

	if(P_nnz != P_n) return false;
	for (int i = 0;i < P_nnz; i++) {
		if ((Pidx[i] != P_i[i]) || (Pptr[i] != P_p[i]) || (Pdata[i] != P_x[i])) {
			return false;
		}
	}
	DELETEA(Pdata);
	DELETEA(Pidx);
	DELETEA(Pptr);


	return true;

	
}

bool testMToCSC2()
{
	MatrixCPU A(3, 2);
	A.set(0, 0, 1);
	A.set(0, 1, 1);
	A.set(1, 0, 1);
	A.set(2, 1, 1);
	c_int An = A.getNNull();
	int col = A.getNCol();
	c_float* Adata = new c_float[An];
	c_int* Aidx = new c_int[An];
	c_int* Aptr = new c_int[col+1];
	
	A.toCSC(Adata, Aidx, Aptr);


	c_float A_x[4] = { 1.0, 1.0, 1.0, 1.0, };
	c_int A_nnz = 4;
	c_int A_i[4] = { 0, 1, 0, 2, };
	c_int A_p[3] = { 0, 2, 4, };
	if (A_nnz != An) return false;
	for (int i = 0;i < A_nnz; i++) {
		if ((Aidx[i] != A_i[i]) || (Adata[i] != A_x[i])) {
			return false;
		}
	}
	for (int i = 0;i < 3; i++) {
		if (A_p[i] != Aptr[i]) {
			return false;
		}
	}
	DELETEA(Adata);
	DELETEA(Aidx);
	DELETEA(Aptr);


	return true;
}
bool testMCSC()
{
	/*
	11   0    0  14   0  16
    0    22   0   0  25  26
    0    0   33  34   0  36
    41   0   43  44   0  46
	
	*/
	MatrixCPU m(4, 6);
	m.set(0, 0, 11);
	m.set(0, 3, 14);
	m.set(0, 5, 16);
	m.set(1, 1, 22);
	m.set(1, 4, 25);
	m.set(1, 5, 26);
	m.set(2, 2, 33);
	m.set(2, 3, 34);
	m.set(2, 5, 36);
	m.set(3, 0, 41);
	m.set(3, 2, 43);
	m.set(3, 3, 44);
	m.set(3, 5, 46);
	m.display();

	float col[7] = { 1, 3, 4, 6, 9, 10, 14 };
	float row[13] = { 1, 4, 2, 3, 4, 1, 3, 4, 2, 1, 2, 3, 4 };
	float entry[13] = { 11.0, 41.0 , 22.0, 33.0, 43.0, 14.0, 34.0, 44.0, 25.0, 16.0, 26.0, 36.0, 46.0 };
		


	c_int A_nnz = m.getNNull();
	int Acol = m.getNCol();
	c_float* Adata = new c_float[A_nnz];
	c_int* Aidx = new c_int[A_nnz];
	c_int* Aptr = new c_int[Acol + 1];
	m.toCSC(Adata, Aidx, Aptr);
	//std::cout << A_nnz << " " << Acol << std::endl;
	


	if (Acol != 6) return false;
	if (A_nnz != 13) return false;

	/*for (int i = 0; i < A_nnz; i++) {
		std::cout << Adata[i] << " " << entry[i] << " " << Aidx[i] + 1 << " " << row[i] << std::endl;
	}
	for (int i = 0; i < Acol + 1; i++) {
		std::cout << Aptr[i] + 1 << " " << col[i] << std::endl;
	}*/


	for (int i = 0; i < A_nnz; i++) {
		if (Adata[i] != entry[i] || (Aidx[i] + 1) != row[i]) return false;
	}
	for (int i = 0; i < Acol+1; i++) {
		if ((Aptr[i]+ 1 ) != col[i]) return false;
	}

	return true;
}

#endif