#pragma once

#include "MatrixCPU.h"
#include <iostream>
#include <stdexcept>

int testMatrix();


bool testMConstru();
bool testMConstru2();
bool testMConstru3();

bool testMEquality1();
bool testMEquality2();
bool testMEquality3();

bool testMGet1();
bool testMGet2();

bool testMSet1();
bool testMSet2();
bool testMSetFromFile();
bool testMSetRand();
bool testMSetFromEigen();
bool testMToEigen();


bool testMAdd1();
bool testMAdd2();
bool testMAdd3();
bool testMAdd4();
bool testMAdd5();
bool testMAdd6();

bool testMAddVect1();
bool testMAddVect2();

bool testMAddTrans1();
bool testMAddTrans2();

bool testMSubstract1();
bool testMSubstract2();
bool testMSubstract3();
bool testMSubstract4();

bool testMSubstractVect1();
bool testMSubstractVect2();

bool testMSubstractTrans1();
bool testMSubstractTrans2();

bool testMMultiply1();
bool testMMultiply2();
bool testMMultiply3();
bool testMMultiplyT1();
bool testMMultiplyT2();
bool testMMultiplyT3();
bool testMMultiplyT4();

bool testMDivide1();
bool testMDivide2();
bool testMDivide3();
bool testMDivide4();
bool testMDivide5();

bool testMDivideGJ1();
bool testMDivideGJ2();

#ifdef EIGEN
    bool testMDivideEigen();
    bool testMSolveSys();
#endif

bool testMmoy1();
bool testMmoy2();
bool testMmoy3();
bool testMmoy4();

bool testMSum1();
bool testMSum2();
bool testMSum3();
bool testMSum4();

bool testMDistance1();
bool testMDistance2();

bool testMProject1();
bool testMProject2();
bool testMProjectPos();
bool testMProjectNeg();

bool testMMax();
bool testMToCSC1();
bool testMToCSC2();

bool testMSwap();

bool testMMinAbs();

bool testMSort();

bool testMCSC();