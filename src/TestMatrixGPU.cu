#include "../head/TestMatrixGPU.cuh"


int testMGPU()
{
	int n = 1;
	if (!testMGPUConstr1()) return n;
	n++;
    if (!testMGPUConstr2()) return n;
    n++;
    if (!testMGPUConstr3()) return n;
    n++;
    if (!testMGPUSet1()) return n;
    n++;
    if (!testMGPUSet2()) return n;
    n++;
    if (!testMGPUSet3()) return n;
    n++;
    if (!testMGPUSet3()) return n;
    n++;
    if (!testMGPUSetForce()) return n;
    n++;
    if (!testMGPUSetGPU()) return n;
    n++;
    if (!testMGPUSetTrans()) return n;
    std::cout << "*************  10   *****************" << std::endl;
    n++;//10
    if (!testMGPUSetBloc()) return n;
    n++; 
    if (!testMGPUTranferG1()) return n;
    n++; 
    if (!testMGPUTranferG2()) return n;
    n++;
    if (!testMGPUTranferC1()) return n;
    n++;
    if (!testMGPUTranferC2()) return n;
    n++;
    if (!testMGPUConv()) return n;
    n++;
    if (!testMGPUAdd1()) return n;
    n++;
    if (!testMGPUAdd2()) return n;
    n++;
    if (!testMGPUAdd3()) return n;
    n++;
    if (!testMGPUAdd4()) return n;
    std::cout << "*************  20   *****************" << std::endl;
    n++;//20
    if (!testMGPUAdd5()) return n;
    n++; 
    if (!testMGPUAdd6()) return n;
    n++;
    if (!testMGPUAddVect1()) return n;
    n++;
    if (!testMGPUAddVect2()) return n;
    n++;
    if (!testMGPUAddTrans1()) return n;
    n++; 
    if (!testMGPUAddTrans2()) return n;
    n++;
    if (!testMGPUSubstract1()) return n;
    n++;
    if (!testMGPUSubstract2()) return n;
    n++;
    if (!testMGPUSubstract3()) return n;
    n++; 
    if (!testMGPUSubstract4()) return n;
    std::cout << "*************  30   *****************" << std::endl;
    n++;//30 
    if (!testMGPUSubstractVect1()) return n;
    n++; 
    if (!testMGPUSubstractVect2()) return n;
    n++;
    if (!testMGPUSubstractTrans1()) return n;
    n++; 
    if (!testMGPUSubstractTrans2()) return n;
    n++;
    if (!testMGPUMultiply()) return n;
    n++; 
    if (!testMGPUMultiply2()) return n;
    n++;
    if (!testMGPUMultiplyT1()) return n;
    n++;
    if (!testMGPUMultiplyT2()) return n;
    n++;
    if (!testMGPUMultiplyT3()) return n;
    n++;
    if (!testMGPUMultiplyT4()) return n;
    std::cout << "*************  40   *****************" << std::endl;
    n++; //40
    if (!testMGPUMultiplyVect()) return n;
    n++;
    if (!testMGPUMultiplyMat()) return n;
    n++;
    if (!testMGPUMultiplyLinearOp()) return n;
    n++;
    if (!testMGPUDivide1()) return n;
    n++; 
    if (!testMGPUDivide2()) return n;
    n++; 
    if (!testMGPUDivide3()) return n;
    n++; 
    if (!testMGPUDivide4()) return n;
    n++;
    if (!testMGPUDivide5()) return n;
    n++;
    if (!testMGPUmoy1()) return n;
    n++; 
    if (!testMGPUmoy2()) return n;
    std::cout << "*************  50   *****************" << std::endl;
    n++;// 50
    if (!testMGPUmoy3()) return n;
    n++; 
    if (!testMGPUmoy4()) return n;
    n++;
    if (!testMGPUProject1()) return n;
    n++;
    if (!testMGPUProject2()) return n;
    n++;
    if (!testMGPUProjectPos()) return n;
    n++; 
    if (!testMGPUProjectNeg()) return n;
    n++; 
    if (!testMGPUSum1()) return n;
    n++; 
    if (!testMGPUSum2()) return n;
    n++;
    if (!testMGPUSum3()) return n;
    n++; 
    if (!testMGPUSum4()) return n;
    std::cout << "*************  60   *****************" << std::endl;
    n++; // 60
    if (!testMGPUSumPartial()) return n;
    n++;
    if (!testMGPUSwap()) return n;
    n++;
    if (!testMGPUDistance()) return n;
    n++;
    if (!testMGPUDistance2()) return n;
    n++;
    if (!testMGPUMax()) return n;
    n++;
    if (!testMGPUMax2()) return n;
    n++; 
    if (!testMGPUMax3()) return n;
    n++;  
    if (!testMGPUDivideGJ1()) return n;
    n++;
    if (!testMGPUDivideGJ2()) return n;
    n++;
    if (!testMGPUSolveSys()) return n;
    n++;

    return 0;

}

bool testMGPUConstr1(){
    std::cout<< "default constructor"<<std::endl;

    MatrixGPU m;

    return true;

}
bool testMGPUConstr2() {
    std::cout << "param constructor" << std::endl;

    int line = 2;
    int column = line + 1;
    float value = 20;
    MatrixGPU m(line,column,value);

    m.display();

    return true;

}
bool testMGPUConstr3() {
    std::cout << "copy constructor" << std::endl;

    int line = 2;
    int column = line + 1;
    float value = 20;
    MatrixCPU m(line, column, value);
    MatrixGPU m2(m);

    m2.display();

    return true;

}


bool testMGPUSet1()
{
    float value2 = 1;
    float value = 4.5;
    int i = 1;
    int j = 2;
    int n = 3;
    MatrixGPU mnull(n, n,value2);
    mnull.set(i, j, value);
    return  (mnull.get(i, j) == value);
}
bool testMGPUSet2()
{
    int n = 3;
    float value = 4;
    int i = 1;
    int j = n;
    MatrixGPU mnull(n, n);
    try
    {
        mnull.set(i, j, value); 
    }
    catch (std::out_of_range&)
    {
        return true;
    }
    return false;
}
bool testMGPUSet3()
{
    float value = 4.5;
    int n = 3;
    MatrixGPU m1(n, n);
    MatrixGPU m11(n, n);
    MatrixGPU m2(n, n, value);
    MatrixGPU m22(n, n, value);
    m1.set(&m2);
    m11.transferGPU();
    m22.transferGPU();
    

    

    m11.set(&m22);
    m11.transferCPU();
    m22.transferCPU();



    return  ((m1.isEqual(&m2)) && (m11.isEqual(&m22)) && (m1.isEqual(&m22)));
}
bool testMGPUSet4()
{
    float value = 4.5;
    int n = 3;
    MatrixGPU m1(n, n);
    MatrixGPU m2(n-1, n, value);
    MatrixGPU m22(n, n, value);
    m22.transferGPU();
    try
    {
        m1.set(&m2); 
    }
    catch (std::invalid_argument&)
    {
        try
        {
            m1.set(&m22); 
        }
        catch (std::invalid_argument&)
        {
            return true;
        }
        return false;
    }
    return false;
}

bool testMGPUSetForce()
{
    float value2 = 1;
    float value = 4.5;
    int i = 1;
    int j = 2;
    int n = 3;
    MatrixGPU mnull(n, n, value2);
    MatrixGPU mnull2(n, n, value2, true);
    mnull2.set(i, j, value, true);
    mnull.set(i, j, value);

    mnull.display();
    mnull2.display(true);

    return  (mnull.get(i, j) == value) && (mnull2.get(i, j, false) == value);
}

bool testMGPUSetGPU()
{
    float value2 = 1;
    float value = 4.5;
    int i = 1;
    int j = 2;
    int n = 3;
    MatrixGPU m(n, n, value2, 1);
    m.set(i, j, value, true);
    float test = m.get(i, j, false);
    m.transferCPU();
    return  (m.get(i, j) == value) && (test ==value);
}

bool testMGPUSetTrans()
{
    int row = 4;
    int col = 5;
    float value1 = 2;
    float value2 = -1;
    MatrixGPU m(row, col, value1);
    MatrixGPU mTrans1(col, row, value1);
    MatrixGPU mTrans2(col, row, 0);
    MatrixGPU mTrans3(col, row, 0, 1);

    m.set(row - 2, col - 2, value2);
    mTrans1.set(col - 2, row - 2, value2);
    mTrans2.setTrans(&m);
    m.transferGPU();
    mTrans3.setTrans(&m);
    mTrans3.transferCPU();

   


    return (mTrans1.isEqual(&mTrans2) && mTrans2.isEqual(&mTrans3));
}

bool testMGPUSetBloc()
{
    int row = 4;
    int col = 5;
    float value1 = 2;
    float value2 = -1;
    MatrixGPU m1(3*row, 4*col, 0, 1);
    MatrixGPU m11(3 * row, 4 * col);
    MatrixGPU m2(row, col, value1);
    
    int iBegin = 0;
    int iEnd = row;
    int jBegin = 0;
    int jEnd = col;

    int iBegin2 = row + 1;
    int iEnd2 = 2*row + 1;
    int jBegin2 = 2 * col + 1;
    int jEnd2 = 3 * col + 1;

    m11.setBloc(iBegin, iEnd, jBegin, jEnd, &m2);
    m11.setBloc(iBegin2, iEnd2, jBegin2, jEnd2, &m2, value2);

    m2.transferGPU();

    m1.setBloc(iBegin, iEnd, jBegin, jEnd, &m2);
    m1.setBloc(iBegin2, iEnd2, jBegin2, jEnd2, &m2, value2);

    m1.transferCPU();

    /*std::cout << "--------------------------" << std::endl;
    m1.display();
    std::cout << "--------------------" << std::endl;
    m11.display();
    std::cout << "--------------------------" << std::endl;*/


    return m1.isEqual(&m11);
}



bool testMGPUTranferG1() {

    int line = 2;
    int column = line + 1;
    float value = 20;
    MatrixGPU m(line, column, value);
    m.transferGPU();
    m.display();

    return true;
}
bool testMGPUTranferG2() {
    int line = 2;
    int column = line + 1;
    float value = 20;
    MatrixGPU m(line, column, value);
    m.transferGPU();

    try
    {
        m.transferGPU();
    }
    catch (const std::domain_error&)
    {
        return true;
    }
    return false;
}
bool testMGPUTranferC1() {
    int line = 2;
    int column = line + 1;
    float value = 20;
    MatrixGPU m(line, column, value);
    m.transferGPU();
    m.transferCPU();
    m.display();
    return true;

}
bool testMGPUTranferC2() {
    int line = 2;
    int column = line + 1;
    float value = 20;
    MatrixGPU m(line, column, value);
    m.transferGPU();
    m.transferCPU();
    try
    {
        m.transferCPU();
    }
    catch (const std::domain_error&)
    {
        return true;
    }
    return false;
}

bool testMGPUConv()
{
    int line = 2;
    int column = line + 1;
    float value = 20;
    MatrixGPU m(line, column, value);
    MatrixCPU m1(line, column, value);
    MatrixCPU m2;
    m.toMatCPU(m2);
    MatrixCPU m3;
    m.transferGPU();
    m.toMatCPU(m3);
    return (m1.isEqual(&m2) && m1.isEqual(&m3));
}

bool testMGPUAdd1()
{
    int line = 100;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;
    float value3 = value + value2;

    MatrixGPU m1(line, column, value);
    MatrixGPU m2(line, column, value2);
    MatrixGPU m3(line, column, value3); 
    MatrixGPU temp(line, column);
    temp.add(&m1, &m2);
    if (!temp.isEqual(&m3)) return false;
    m1.transferGPU();
    m2.transferGPU();
    temp.transferGPU();
    temp.add(&m1, &m2);
    temp.transferCPU();


    return temp.isEqual(&m3);

}
bool testMGPUAdd2() {
    int line = 100;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;
    float value3 = value + value2;
    MatrixGPU m1(line, column, value);
    MatrixGPU m3(line, column, value3);
    MatrixGPU temp(line, column);
    temp.add(&m1, value2);
    if (!temp.isEqual(&m3)) return false;
    m1.transferGPU();
    temp.transferGPU();
    temp.add(&m1, value2);
    temp.transferCPU();
    return temp.isEqual(&m3);
    
}
bool testMGPUAdd3()
{
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;
    float value3 = value + value2;

    MatrixGPU m1(line, column, value);
    MatrixGPU m11(line, column, value);
    MatrixGPU m2(line, column, value2);
    MatrixGPU m3(line, column, value3); 
    m1.add(&m2);
    
    if (!m1.isEqual(&m3)) return false;

    m11.transferGPU();
    m2.transferGPU();
    m11.add(&m2);
    m11.transferCPU();
    
    return m11.isEqual(&m3);



}
bool testMGPUAdd4()
{
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = -1;
    float value3 = value + value2;
    MatrixGPU m1(line, column, value);
    MatrixGPU m11(line, column, value);
    MatrixGPU m3(line, column, value3);
    m1.add(value2);

    if (!m1.isEqual(&m3)) return false;

    m11.transferGPU();
    m11.add(value2);
    m11.transferCPU();

    return m11.isEqual(&m3);
}
bool testMGPUAdd5()
{
    int line = 100;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;


    MatrixGPU m1(line, column, value);
    MatrixGPU m2(line, column, value2);
    MatrixGPU m3(column, line);
    MatrixGPU m4(line, column);
    m1.transferGPU();
    m2.transferGPU();
    m3.transferGPU();

    try
    {
        m3.add(&m1, &m2);
    }
    catch (std::invalid_argument&)
    {
        try
        {
            m4.add(&m1, &m2);
        }
        catch (std::invalid_argument&)
        {
            return true;
        }
        return false;
    }
    return false;
}
bool testMGPUAdd6()
{
    int line = 100;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;


    MatrixGPU m1(line, column, value);
    MatrixGPU m2(column, line, value2);
    MatrixGPU m3(line, column);
    m1.transferGPU();
    m2.transferGPU();
    try
    {
        m1.add(&m2);
    }
    catch (std::invalid_argument&)
    {
        try
        {
            m1.add(&m3);
        }
        catch (std::invalid_argument&)
        {
            return true;
        }
        return false;
    }
    return false;
}

bool testMGPUAddVect1()
{
    int line = 100;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;
    MatrixGPU m1(line, column, value);
    MatrixGPU m4(column, line, value2);
    MatrixGPU vect(line, 1, value2);
    vect.transferGPU();
    try {
        m1.addVector(&m4);
    }
    catch (std::invalid_argument&) {
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
bool testMGPUAddVect2()
{
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;
    MatrixGPU m1(line, column, value);
    MatrixGPU m11(line, column, value);
    MatrixGPU m2(line, column, value2);
    MatrixGPU vect(1, column, value2 - value);
    MatrixGPU vect1(line, 1, value2 - value);

    m1.transferGPU();
    m11.transferGPU();
    vect.transferGPU();
    vect1.transferGPU();

    m1.addVector(&vect); 
    m1.transferCPU();
    if (!m1.isEqual(&m2)) return false;
    m11.addVector(&vect1); 
    m11.transferCPU();
    if (!m2.isEqual(&m11)) return false;

    return true;
}

bool testMGPUAddTrans1()
{
    int line = 4;
    int column = line + 1;
    int i = 1;
    int j = 2;
    float value = 2;
    float value2 = 1;
    float value3 = value + value2;
    float value4 = 4.5;
    float value5 = value + value4;

    MatrixGPU m1(line, column, value);
    MatrixGPU m2(column, line, value2);
    m2.set(i, j, value4);
    MatrixGPU m3(line, column, value3); 
    m3.set(j, i, value5);
    m1.transferGPU();
    m2.transferGPU();
    m1.addTrans(&m2);
    m1.transferCPU();
    return m1.isEqual(&m3);
}
bool testMGPUAddTrans2() {
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;

    MatrixGPU m1(line, column, value);
    MatrixGPU m2(line, column, value2);
    MatrixGPU m3(column, line, value2);
    m3.transferGPU();
    try
    {
        m1.addTrans(&m2);
    }
    catch (std::invalid_argument&)
    {
        try
        {
            m1.addTrans(&m3);
        }
        catch (std::invalid_argument&)
        {
            return true;
        }
        return false;
    }

    return false;
}


bool testMGPUSubstract1()
{
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;
    float value3 = value - value2;

    MatrixGPU m1(line, column, value);
    MatrixGPU m2(line, column, value2);
    MatrixGPU m3(line, column, value3); 
    MatrixGPU temp(line, column);
    m1.transferGPU();
    m2.transferGPU();
    temp.transferGPU();
    temp.subtract(&m1, &m2);
    temp.transferCPU();
    return temp.isEqual(&m3);

}
bool testMGPUSubstract2()
{
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;
    float value3 = value - value2;

    MatrixGPU m1(line, column, value);
    MatrixGPU m2(line, column, value2);
    MatrixGPU m3(line, column, value3); 
    m1.transferGPU();
    m2.transferGPU();
    m1.subtract(&m2);
    m1.transferCPU();
    return m1.isEqual(&m3);

}
bool testMGPUSubstract3()
{
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;

    MatrixGPU m1(line, column, value);
    MatrixGPU m2(line, column, value2);
    MatrixGPU m3(column, line);
    MatrixGPU m4(line, column);
    m4.transferGPU();
    try
    {
        m3.subtract(&m1, &m2);
    }
    catch (std::invalid_argument&)
    {
        try
        {
            m4.subtract(&m1, &m2);
        }
        catch (std::invalid_argument&)
        {
            return true;
        }
        return false;
    }
    return false;
}
bool testMGPUSubstract4()
{
    int line = 100;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;


    MatrixGPU m1(line, column, value);
    MatrixGPU m2(column, line, value2);
    MatrixGPU m3(line, column);
    m1.transferGPU();
    m2.transferGPU();
    try
    {
        m1.subtract(&m2);
    }
    catch (std::invalid_argument&)
    {
        try
        {
            m1.subtract(&m3);
        }
        catch (std::invalid_argument&)
        {
            return true;
        }
        return false;
    }
    return false;
}

bool testMGPUSubstractVect1()
{
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;
    MatrixGPU m1(line, column, value);
    MatrixGPU m4(column, line, value2);
    MatrixGPU vect(1, line, value2);

    try {
        m1.subtractVector(&m4);
    }
    catch (std::invalid_argument&) {
        try
        {
            m1.subtractVector(&vect);
        }
        catch (std::invalid_argument&)
        {
            return true;
        }
        return false;
    }
    return false;

}
bool testMGPUSubstractVect2()
{
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;
    MatrixGPU m1(line, column, value);
    MatrixGPU m11(line, column, value);
    MatrixGPU m2(line, column, value2);
    MatrixGPU vect(1, column, value - value2);
    MatrixGPU vect1(line, 1, value - value2);

    m1.transferGPU();
    vect.transferGPU();
    vect1.transferGPU();
    m11.transferGPU();
    m1.subtractVector(&vect);
    m1.transferCPU();
    if (!m1.isEqual(&m2)) return false;

    m11.subtractVector(&vect1);
    m11.transferCPU();
    if (!m2.isEqual(&m11)) return false;

    return true;
}

bool testMGPUSubstractTrans1()
{
    int line = 4;
    int column = line + 1;
    int i = 1;
    int j = 2;
    float value = 2;
    float value2 = 1;
    float value3 = value - value2;
    float value4 = 4.5;
    float value5 = value - value4;

    MatrixGPU m1(line, column, value);
    MatrixGPU m2(column, line, value2);
    m2.set(i, j, value4);
    MatrixGPU m3(line, column, value3); 
    m3.set(j, i, value5);
    m1.transferGPU();
    m2.transferGPU();
    m1.subtractTrans(&m2);
    m1.transferCPU();
    return m1.isEqual(&m3);

}
bool testMGPUSubstractTrans2() {
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;

    MatrixGPU m1(line, column, value);
    MatrixGPU m2(line, column, value2);
    MatrixGPU m3(column, line, value2);
    m3.transferGPU();
    try
    {
        m1.subtractTrans(&m2);
    }
    catch (std::invalid_argument&)
    {
        try
        {
            m1.subtractTrans(&m3);
        }
        catch (std::invalid_argument&)
        {
            return true;
        }
        return false;
    }

    return false;
}

bool testMGPUMultiply()
{
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;
    float value3 = value * value2;
    MatrixGPU m1(line, column, value);
    MatrixGPU m3(line, column, value3);
    m1.transferGPU();
    m1.multiply(value2);
    m1.transferCPU();
    return m1.isEqual(&m3);
}

bool testMGPUMultiply2()
{
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    int line = 4;
    int column = line + 1;
    float value = 1.5;
    float value2 = 2;
    float value3 = -3;

    MatrixGPU result(line, 1, 0, 1);

    MatrixGPU result2(line, 1, 0);
    MatrixGPU Mat(line, column, value);

    MatrixGPU vect(column, 1, value2);

    Mat.set(line - 2, column - 2, value3);
    vect.set(column - 3, 0, value3);

    try
    {
        result2.multiply(&Mat, &vect);

    }
    catch (const std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return false;
    }

    Mat.transferGPU();
    vect.transferGPU();
   
    try
    {
        t1 = std::chrono::high_resolution_clock::now();
        result.multiply(&Mat, &vect);
        t2 = std::chrono::high_resolution_clock::now();
        std::cout << "temps de calcul de multiply su GPU " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return false;
    }
    
    
    result.transferCPU();
   
    bool testresult = result2.isEqual(&result);

    if (!testresult) {
        result.display();
        result2.display();
    }


    return testresult;
}

bool testMGPUMultiplyT1()
{
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;
    float value3 = value * value2;

    MatrixGPU m1(line, column, value);
    MatrixGPU m2(line, column, value2);
    MatrixGPU m3(line, column, value3); 
    MatrixGPU temp(line, column);

    m1.transferGPU();
    m2.transferGPU();
    temp.transferGPU();

    temp.multiplyT(&m1, &m2);
    temp.transferCPU();

    return temp.isEqual(&m3);
}
bool testMGPUMultiplyT2()
{
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;
    float value3 = value * value2;

    MatrixGPU m1(line, column, value);
    MatrixGPU m2(line, column, value2);
    MatrixGPU m3(line, column, value3); 
    m1.transferGPU();
    m2.transferGPU();
    
    m1.multiplyT(&m2);
    m1.transferCPU();

    return m1.isEqual(&m3);
}
bool testMGPUMultiplyT3()
{
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;
    float value3 = value * value2;

    MatrixGPU m1(line, column, value);
    MatrixGPU m2(column, line, value2);
    MatrixGPU m3(line, column, value3);
    MatrixGPU temp(line, column);
    m1.transferGPU();
    m2.transferGPU();
    temp.transferGPU();

    try {
        temp.multiplyT(&m1, &m2);
    }
    catch (std::invalid_argument&) {
        try
        {
            temp.multiplyT(&m1, &m3);
        }
        catch (std::invalid_argument&)
        {
            return true;
        }
        return false;
    }
    return false;
}
bool testMGPUMultiplyT4()
{
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;


    MatrixGPU m1(line, column, value);
    MatrixGPU m2(line, line, value2);
    MatrixGPU m3(line, column, value);
    m1.transferGPU();
    m2.transferGPU();

    try {
        m1.multiplyT(&m2);
    }
    catch (std::invalid_argument&) {

        try
        {
            m1.multiplyT(&m3);
        }
        catch (std::invalid_argument&)
        {
            return true;
        }
        return false;
    }
    return false;
}

bool testMGPUMultiplyVect()
{
    int nRow = 4;
    int nCol = 6;
    float value1 = 4;
    float value2 = 72.36;
    float value3 = -42.78;
    float value4 = 19.832;

    MatrixGPU y(nRow, 1, 0, 1);
    MatrixGPU yCPU(nRow, 1);

    MatrixGPU A(nRow, nCol, value1);
    MatrixGPU x(nCol, 1, value2);
    A.set(nRow - 1, nCol - 2, value3);
    x.set(nCol - 3, 0, value4);

    yCPU.multiply(&A, &x);
    A.transferGPU();
    x.transferGPU();

    y.multiply(&A, &x);
    y.transferCPU();


    return yCPU.isEqual(&y);
}

bool testMGPUMultiplyMat()
{
    int nRow = 4;
    int nCol = 6;
    int common = 5;
    float value1 = 4;
    float value2 = 72.36;
    float value3 = -42.78;
    float value4 = 19.832;

    MatrixGPU y(nRow, nCol, 0, 1);
    MatrixGPU yCPU(nRow, nCol);

    MatrixGPU A(nRow, common, value1);
    MatrixGPU x(common, nCol, value2);
    A.set(nRow - 1, common - 2, value3);
    x.set(common - 3, nCol - 4, value4);

    yCPU.multiplyMat(&A, &x);
    A.transferGPU();
    x.transferGPU();

    y.multiplyMat(&A, &x);
    y.transferCPU();


    return yCPU.isEqual(&y);
}

bool testMGPUMultiplyLinearOp()
{
    int nRow = 4;
    int nCol = 6;
    float value1 = 4;
    float value2 = 72.36;
    float value3 = -42.78;
    float value4 = 19.832;
    float value5 = 12.654;
    float value6 = 93.47;

    MatrixGPU y(nRow, 1, 0, 1);
    MatrixGPU yCPU(nRow, 1);

    MatrixGPU A(nRow, nCol, value1);
    MatrixGPU x(nCol, 1, value2);
    MatrixGPU b(nCol, 1, value5);
    A.set(nRow - 1, nCol - 2, value3);
    x.set(nCol - 3, 0, value4);
    b.set(nCol - 1, 0, value6);

    yCPU.linearOperation(&A, &x, &b);
    A.transferGPU();
    x.transferGPU();
    b.transferGPU();

    y.linearOperation(&A, &x, &b);
    y.transferCPU();


    return yCPU.isEqual(&y);
}



bool testMGPUDivide1()
{
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1.2;
    float value3 = value / value2;


    MatrixGPU m1(line, column, value);
    MatrixGPU m2(line, column, value2);
    MatrixGPU m3(line, column, value3); 
    m1.transferGPU();
    m2.transferGPU();
    m1.divideT(&m2);
    m1.transferCPU();
    return m1.isEqual(&m3);
}
bool testMGPUDivide2()
{
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;


    MatrixGPU m1(line, column, value);
    MatrixGPU m2(column, line, value2);
    MatrixGPU m3(line, column, value);
    m3.transferGPU();
    try
    {
        m1.divideT(&m2);
    }
    catch (std::invalid_argument&)
    {
        try
        {
            m1.divideT(&m3);
        }
        catch (std::invalid_argument&)
        {
            return true;
        }
        return false;
    }
    return false;
}
bool testMGPUDivide3()
{
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;


    MatrixGPU m1(line, column, value);
    MatrixGPU m2(line, column, value2);
    m2.set(line - 1, column - 1, 0);
    
    m1.transferGPU();
    m2.transferGPU();

    m1.divideT(&m2);
    m1.transferCPU();

    return (m1.get(line - 1, column - 1) == std::numeric_limits<float>::infinity());
    
}
bool testMGPUDivide4()
{
    int line = 2;
    int column = line + 1;
    float value = 1.5;
    float value2 = -1.2;
    float value3 = value / value2;
    MatrixGPU m1(line, column, value);
    MatrixGPU m3(line, column, value3);

    m1.transferGPU();

    m1.divide(value2);

    m1.transferCPU();

    return m1.isEqual(&m3);
}
bool testMGPUDivide5()
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

bool testMGPUmoy1()
{
    int ligne = 3;
    int colonne = 4;
    float value = 1.5;
    MatrixGPU m1(ligne, colonne, 0);
    m1.set(0, 0, value);
    m1.set(0, 3, value);
    m1.set(1, 0, value);
    m1.set(1, 1, value);
    m1.set(1, 3, value);
    m1.set(2, 2, value);
    MatrixGPU nb1col(1, colonne, 0);
    nb1col.set(0, 0, 2);
    nb1col.set(0, 1, 1);
    nb1col.set(0, 2, 1);
    nb1col.set(0, 3, 2);
    MatrixGPU temp1(1, colonne);
    MatrixGPU m3(1, colonne, value);


    temp1.transferGPU();
    m1.transferGPU();
    nb1col.transferGPU();
    temp1.Moy(&m1, &nb1col, 1);
    temp1.transferCPU();

    return temp1.isEqual(&m3);
}
bool testMGPUmoy2()
{
    int ligne = 3;
    int colonne = 4;
    float value = 1.5;
    MatrixGPU m1(ligne, colonne, 0);
    m1.set(0, 0, value);
    m1.set(0, 3, value);
    m1.set(1, 0, value);
    m1.set(1, 1, value);
    m1.set(1, 3, value);
    m1.set(2, 2, value);
    MatrixGPU nb2li(ligne, 1, 2);
    MatrixGPU temp2(ligne, 1);
    MatrixGPU m5(ligne, 1, value);
    m5.set(1, 0, value * 3 / 2);
    m5.set(2, 0, value / 2);

    temp2.transferGPU();
    m1.transferGPU();
    nb2li.transferGPU();
    temp2.Moy(&m1, &nb2li, 0);
    temp2.transferCPU();
    

    return temp2.isEqual(&m5);
}
bool testMGPUmoy3()
{
    int ligne = 3;
    int colonne = 4;
    float value = 1.5;
    MatrixGPU m1(ligne, colonne, 0);
    m1.set(0, 0, value);
    m1.set(0, 3, value);
    m1.set(1, 0, value);
    m1.set(1, 1, value);
    m1.set(1, 3, value);
    m1.set(2, 2, value);
    MatrixGPU nb2li(ligne, 1, 2);
    MatrixGPU temp2(ligne, 1);
    try
    {
        temp2.Moy(&m1, &nb2li, 1); 
    }
    catch (std::invalid_argument&)
    {
        return true;
    }
    return false;
}
bool testMGPUmoy4()
{
    int ligne = 3;
    int colonne = 4;
    float value = 1.5;
    MatrixGPU m1(ligne, colonne, 0);
    m1.set(0, 0, value);
    m1.set(0, 3, value);
    m1.set(1, 0, value);
    m1.set(1, 1, value);
    m1.set(1, 3, value);
    m1.set(2, 2, value);
    MatrixGPU nb2li(ligne, 1, 2);
    MatrixGPU temp2(ligne, 1);
    MatrixGPU temp1(1, colonne);
    try
    {
        temp1.Moy(&m1, &nb2li, 0); 
    }
    catch (std::invalid_argument&)
    {
        try
        {
            temp2.Moy(&m1, &m1, 0); 
        }
        catch (std::invalid_argument&)
        {
            return true;
        }
        return false;
    }
    return false;
}

bool testMGPUProject1()
{
    int ligne = 3;
    int colonne = 4;
    float value = 1.5;
    float value2 = 0.5;
    float value3 = 3;
    MatrixGPU m1(ligne, colonne, value);
    m1.set(0, 2, value2 - 1);
    m1.set(1, 1, value3 + 1);
    MatrixGPU m2(ligne, colonne, value);
    m2.set(0, 2, value2);
    m2.set(1, 1, value3);
    MatrixGPU m22(m2);
    MatrixGPU Lb(ligne, colonne, value2);
    MatrixGPU Ub(ligne, colonne, value3);

    m1.transferGPU();
    m2.transferGPU();
    Lb.transferGPU();
    Ub.transferGPU();

    m1.project(&Lb, &Ub);
    m1.transferCPU();
    if (!m1.isEqual(&m22)) return false;
    m2.project(&Lb, &Ub);
    m2.transferCPU();
    if (!m2.isEqual(&m22)) return false;

    return true;
}
bool testMGPUProject2()
{
    int ligne = 3;
    int colonne = 4;
    float value = 1.5;
    float value2 = 0.5;
    float value3 = 3;
    MatrixGPU m1(ligne, colonne, value);
    MatrixGPU m3(colonne, ligne, value);
    MatrixGPU lb(ligne, colonne, value2);
    MatrixGPU ub(ligne, colonne, value3);

    lb.transferGPU();
    try
    {
        m1.project(&m3, &ub); 
    }
    catch (std::invalid_argument&)
    {
        try
        {
            m1.project(&lb, &ub); 
        }
        catch (std::invalid_argument&)
        {
            return true;
        }
        return false;
    }
    return false;

}

bool testMGPUProjectPos()
{
    int row = 2;
    int column = 2;
    float value = 1.5;
    float value2 = -0.5;

    MatrixGPU m1(row, column, value2);
    m1.set(0, 1, 0);
    m1.set(1, 1, value);
    MatrixGPU m2(row, column);
    m2.set(1, 1, value);
    m1.transferGPU();
    

    m1.projectPos();
    m1.transferCPU();

    return m1.isEqual(&m2);
}

bool testMGPUProjectNeg()
{
    int row = 2;
    int column = 2;
    float value = -1.5;
    float value2 = 0.5;

    MatrixGPU m1(row, column, value2);
    m1.set(0, 1, 0);
    m1.set(1, 1, value);
    MatrixGPU m2(row, column);
    m2.set(1, 1, value);


    m1.transferGPU();


    m1.projectNeg();
    m1.transferCPU();


    return m1.isEqual(&m2);
}



bool testMGPUSum1()
{
    int ligne = 100;
    int colonne = 200;
    float value = 1.5;
    float value3 = 5;
    float value2 = (ligne * colonne -1) * value+value3;
    MatrixGPU m1(ligne, colonne, value);
    MatrixGPU m11(ligne, colonne, value);
    m1.set(ligne - 1, colonne - 1, value3);
    m11.set(ligne - 1, colonne - 1, value3);
    m1.transferGPU();
    float value4 = m11.sum();
    float value5 = m1.sum();

    std::cout << "testMGPUSum1 : " << value2 << " " << value4 << " " << value5 << std::endl;
    return ((value2 == value4) && (value2==value5));
}
bool testMGPUSum2()
{
    int ligne = 3;
    int colonne = 4;
    float value = 1.5;
    float value2 = ligne * colonne * value;
    MatrixGPU m1(ligne, colonne, value);
    MatrixGPU temp2(ligne, 1);
    MatrixGPU m4(ligne, 1, value * colonne);

    m1.transferGPU();
    temp2.transferGPU();

    temp2.sum(&m1);
    temp2.transferCPU();

    return temp2.isEqual(&m4);
}
bool testMGPUSum3()
{
    int ligne = 3;
    int colonne = 4;
    float value = 1.5;
    float value2 = ligne * colonne * value;
    MatrixCPU m1(ligne, colonne, value);
    MatrixCPU m2(colonne, colonne, value);
    MatrixCPU temp2(ligne, 1);

    try
    {
        temp2.sum(&m2, 0); 
    }
    catch (std::invalid_argument&)
    {
        try
        {
            m1.sum(&m2, 0); 
        }
        catch (std::invalid_argument&)
        {
            return true;
        }
        return false;
    }
    return false;
}
bool testMGPUSum4()
{
    int ligne = 3;
    int colonne = 4;
    float value = 1.5;
    float value2 = ligne * colonne * value;
    MatrixGPU m1(ligne, colonne, value);
    MatrixGPU temp2(ligne, 1);

    temp2.transferGPU();
    try
    {
        temp2.sum(&m1); 
    }
    catch (std::invalid_argument&)
    {
        return true;
    }
    return false;
}

bool testMGPUSumPartial()
{
    int ligne = 200;
    int colonne = 100;
    float value = 1.5;
    float value3 = 5;
    float value2 = (ligne * colonne) * value /2 ;
    float value22 = ((ligne * colonne) / 2 - 1)* value + value3;
    MatrixGPU m1(ligne, colonne, value);
    MatrixGPU m11(ligne, colonne, value);
    m1.set(ligne - 1, colonne - 1, value3);
    m11.set(ligne - 1, colonne - 1, value3);
    m1.transferGPU();
    float value4 = m11.sum(0, ligne * colonne / 2);
    float value5 = m11.sum(ligne * colonne / 2, ligne * colonne);

    float value42 = m1.sum(0, ligne * colonne / 2);
    float value52 = m1.sum(ligne * colonne / 2, ligne * colonne);

    std::cout << "testMGPUSumPartial : " << value2  << " " << value4 << " " << value42 << std::endl;
    std::cout << "testMGPUSumPartial : " << value22 << " " << value5 << " " << value52 << std::endl;
    return ((value2 == value4) && (value2 == value42) && (value22 == value5) && (value22 == value52));
}

bool testMGPUSwap() {
    int line = 100;
    int column = line + 1;
    float value = 1.5;
    float value2 = 1;
    float value3 = value + value2;

    MatrixGPU m1(line, column, value);
    MatrixGPU m11(line, column, value);
    MatrixGPU m2(line, column, value2);
    MatrixGPU m22(line, column, value2);

    m1.transferGPU();
    m2.transferGPU();
    m1.swap(&m2);

    m1.transferCPU();
    m2.transferCPU();

    return ((m1.isEqual(&m22)) && (m2.isEqual(&m11)));

}

bool testMGPUDistance()
{
    int colonne = 8;
    int ligne = 8;
    float value1 = 2;
    MatrixGPU m(ligne, colonne, value1);
    float value2 = sqrtf(ligne * colonne * value1 * value1);
    float value4 = m.distance2();

    if (value4 != value2) return false;

    m.transferGPU();
    float value3 = m.distance2();


    return (value2==value3);
}

bool testMGPUDistance2()
{
    int colonne = 2345;
    int ligne = 1234;
    float value1 = 2;
    float value5 = -5;
    MatrixGPU m(ligne, colonne, value1);
    MatrixGPU m2(ligne, colonne, value5);
    float value2 = sqrtf(ligne * colonne * (value1-value5) * (value1 - value5));
    float value4 = m.distance2(&m2);


    std::cout << "distance " << value2 << " " << value4 << std::endl;


    if (value4 != value2) return false;

    m.transferGPU();
    m2.transferGPU();
    float value3 = m.distance2(&m2);

    std::cout << "distance " << value2 << " " << value3 << std::endl;
    return (value2 == value3);
}

bool testMGPUMax()
{
    int colonne = 1;
    int ligne = 1;
    float value1 = 2;
    float value2 = 4;
    MatrixGPU m(ligne, colonne, value1);
    m.set(ligne - 1, colonne - 1, value2);

    float value3 = m.max2();
    


    if (value3 != value2) return false;

    m.transferGPU();

    float value4 = m.max2();

    std::cout << "testMGPUMax : " << value2 << " " << value3 << " " << value4 << std::endl;
    return (value2 == value4);
}

bool testMGPUMax2()
{
    int colonne = 1000;
    int ligne = 1000;
    float value1 = 2;
    float value2 = 4;
    float value3 = -4;
    MatrixGPU m(ligne, colonne, value1);
    m.set(ligne - 2, colonne - 2, value2);
    MatrixGPU m2(ligne, colonne, value3);
    

    float value4 = m.max2(&m2);


    if (value4 != (value2-value3)) return false;

    m.transferGPU();
    m2.transferGPU();
    value4 = m.max2(&m2);

    
    return ((value2 - value3) == value4);
}


bool testMGPUMax3()
{
    int colonne = 1000;
    int ligne = 1000;
    float value1 = 2;
    float value2 = -4;
    MatrixGPU m(ligne, colonne, value1);
    m.set(ligne - 2, colonne - 2, value2);


    int pos1 = (ligne - 2) * colonne + colonne - 2;
    int pos2 = 0;
    float value3 = m.max2(&pos2);
    

    if ((abs(value3) != abs(value2))||(pos1 != pos2)) return false;

    m.transferGPU();

    float value4 = m.max2();

   
    return ((abs(value2) == abs(value4)) && (pos1 == pos2));
}

bool testMGPUDivideGJ1()
{
    int n = 3;
    MatrixGPU ident(n, n, 0, 1);
    ident.setEyes(1);

    MatrixGPU invert(n, n, 0, 1);
    invert.invertGaussJordan(&ident);
    
    ident.transferCPU();
    invert.transferCPU();

    if (!ident.isEqual(&invert)) return false;
   

    MatrixGPU m1(n, n);
    m1.set(0, 0, 2);
    m1.set(0, 1, -1);
    m1.set(1, 1, -1);
    m1.set(1, 2, 2);
    m1.set(2, 0, -1);
    m1.set(2, 1, 2);
    m1.set(2, 2, 1);
    m1.display();
    m1.transferGPU();

    MatrixGPU m2(n, n, 0,1);
    MatrixGPU m22(n, n);
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

    m2.transferCPU();

    m22.display();
    m2.display();
    //
    MatrixCPU m4(5, 5);
    
    m4.set(0, 0, 1.5);
    m4.set(0, 4, 1);
    m4.set(1, 1, 1.5);
    m4.set(1, 2, 1);
    m4.set(1, 3, -1);
    m4.set(1, 4, 1);
    m4.set(2, 1, -0.00000002829688838801303063519299030303955078125);
    m4.set(3, 1, 0.00000002829688838801303063519299030303955078125);
    m4.set(2, 2, -0.800027906894683837890625);
    m4.set(3, 3, -0.799972116947174072265625);
    m4.set(4, 0, 1);
    m4.set(4, 1, 1);
    m4.display();

    MatrixGPU m4GPU(m4,1);

    MatrixCPU m3(5,5);
    MatrixGPU m3GPU(5, 5, 0, 1);
    MatrixCPU m3CPU(5, 5);
    
    m3.invertGaussJordan(&m4);
    m3GPU.invertGaussJordan(&m4GPU);
    
    m3GPU.toMatCPU(m3CPU);

    m3CPU.display();
    m3.display();

    return m2.isEqual(&m22) && m3CPU.isEqual(&m3);
}

bool testMGPUDivideGJ2()
{
    int n = 3;
    MatrixGPU temp1(n, n + 1);
    MatrixGPU temp2(n, n);
    temp2.set(0, 0, 1);
    temp2.set(2, 2, 1);

    //temp2.transferGPU();
    //temp1.transferGPU();

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

bool testMGPUSolveSys()
{
    Eigen::Matrix3f M;
    Eigen::Vector3f b;
    M << 1, 2, 3, 4, 5, 6, 7, 8, 10;
    b << 3, 3, 4;
    std::cout << "Here is the matrix M:\n" << M << std::endl;
    std::cout << "Here is the vector b:\n" << b << std::endl;
    Eigen::Vector3f x = M.colPivHouseholderQr().solve(b);
    std::cout << "The solution is:\n" << x << std::endl;

    MatrixGPU Mm(3, 3);
    MatrixGPU bm(3, 1);
    MatrixGPU xm(3, 1, 0, 1);
    MatrixGPU Am(3, 3, 0, 1);
    MatrixGPU P(4, 1, 0, 1);

    MatrixGPU Mm2(3, 3);
    MatrixGPU bm2(3, 1);
    MatrixGPU xm2(3, 1);
    MatrixGPU Am2(3, 3);
    MatrixGPU P2(4, 1);

    for (int i = 0; i < 3; i++) {
        bm.set(i, 0, b(i));
        bm2.set(i, 0, b(i));
        for (int j = 0; j < 3; j++) {
            Mm.set(i, j, M(i, j));
            Mm2.set(i, j, M(i, j));
        }
    }
    /*std::cout << "M sur CPU" << std::endl;
    Mm2.display();
    std::cout << "M sur GPU" << std::endl;
    Mm.display();
    std::cout << "b sur CPU" << std::endl;
    bm2.display();
    std::cout << "b sur GPU" << std::endl;
    bm.display();*/

    Mm.transferGPU();
    bm.transferGPU();

    Mm.LUPFactorization(&Am, &P);
    Mm2.LUPFactorization(&Am2, &P2);
    std::cout << "A sur CPU" << std::endl;
    Am2.display();
    std::cout << "A sur GPU" << std::endl;
    Am.display(true);
    std::cout << "P sur CPU" << std::endl;
    P2.display();
    std::cout << "P sur GPU" << std::endl;
    P.display(true);/**/

    
    xm.solveSys(&Am, &P, &bm);
    xm2.solveSys(&Am2, &P2, &bm2);
   
    xm.transferCPU();
    

    /*std::cout << "solution sur CPU" << std::endl;
    xm2.display();
    std::cout << "solution sur GPU" << std::endl;
    xm.display();*/

    for (int i = 0; i < 3; i++) {
        if (abs(xm.get(i, 0) - x(i)) > 0.00001) {
            std::cout << xm.get(i, 0) << " " << x(i) << " " << std::endl;
            return false;
        }
    }
    for (int i = 0; i < 3; i++) {
        if (abs(xm2.get(i, 0) - x(i)) > 0.00001) {
            std::cout << xm2.get(i, 0) << " " << x(i) << " " << std::endl;
            return false;
        }
    }


    return true;
}
