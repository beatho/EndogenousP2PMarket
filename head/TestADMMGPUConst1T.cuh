#pragma once
#include "ADMMGPUConst1T.cuh"
#include "TestKernel.cuh"
#include <math.h>
#include <chrono>

int testADMMGPUConst1T();


bool testADMMGPUConst1TContruct1();
bool testADMMGPUConst1TContruct2();
bool testADMMGPUConst1TContruct3();

bool testADMMGPUConst1TSolve1();
bool testADMMGPUConst1TSolve2();
bool testADMMGPUConst1TSolve3();


bool testADMMGPUConst1TLAMBDA();
bool testADMMGPUConst1TKappa();
bool testADMMGPUConst1TBt1();
bool testADMMGPUConst1TCP();
bool testADMMGPUConst1TCpb();

bool testADMMGPUConst1TUpdateRes();
bool testADMMGPUConst1TCalcRes();

bool testADMMGPUConst1TTradeP();

bool testADMMGPUConst1Talpha();
bool testADMMGPUConst1TQ();




