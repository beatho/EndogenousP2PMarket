#pragma once

#include "MatrixGPU.cuh"
#include "MatrixCPU.h"
#include <iostream>
#include <stdexcept>


int testMGPU(); // host

bool testMGPUConstr1();
bool testMGPUConstr2();
bool testMGPUConstr3();

bool testMGPUSet1();
bool testMGPUSet2();
bool testMGPUSet3();
bool testMGPUSet4();
bool testMGPUSetForce();
bool testMGPUSetGPU();
bool testMGPUSetTrans();
bool testMGPUSetBloc();


bool testMGPUTranferG1();
bool testMGPUTranferG2();
bool testMGPUTranferC1();
bool testMGPUTranferC2();

bool testMGPUConv();

bool testMGPUAdd1();
bool testMGPUAdd2();
bool testMGPUAdd3();
bool testMGPUAdd4();
bool testMGPUAdd5();
bool testMGPUAdd6();

bool testMGPUAddVect1();
bool testMGPUAddVect2();

bool testMGPUAddTrans1();
bool testMGPUAddTrans2();

bool testMGPUSubstract1();
bool testMGPUSubstract2();
bool testMGPUSubstract3();
bool testMGPUSubstract4();

bool testMGPUSubstractVect1();
bool testMGPUSubstractVect2();

bool testMGPUSubstractTrans1();
bool testMGPUSubstractTrans2();

bool testMGPUMultiply();
bool testMGPUMultiply2();

bool testMGPUMultiplyT1();
bool testMGPUMultiplyT2();
bool testMGPUMultiplyT3();
bool testMGPUMultiplyT4();

bool testMGPUMultiplyVect();
bool testMGPUMultiplyMat();
bool testMGPUMultiplyLinearOp();

bool testMGPUDivide1();
bool testMGPUDivide2();
bool testMGPUDivide3();
bool testMGPUDivide4();
bool testMGPUDivide5();

bool testMGPUmoy1();
bool testMGPUmoy2();
bool testMGPUmoy3();
bool testMGPUmoy4();

bool testMGPUProject1();
bool testMGPUProject2();
bool testMGPUProjectPos();
bool testMGPUProjectNeg();

bool testMGPUSum1();
bool testMGPUSum2();
bool testMGPUSum3();
bool testMGPUSum4();
bool testMGPUSumPartial();

bool testMGPUSwap();

bool testMGPUDistance();
bool testMGPUDistance2();

bool testMGPUMax();
bool testMGPUMax2();
bool testMGPUMax3();

bool testMGPUDivideGJ1();
bool testMGPUDivideGJ2();


bool testMGPUSolveSys();