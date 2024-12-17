#pragma once
#include "ADMMGPUConst1.cuh"
#include "TestKernel.cuh"
#include <math.h>
#include <chrono>

int testADMMGPUConst1();
void testADMMGPUConst1Time(int test=0);

bool testADMMGPUConst1Contruct1();
bool testADMMGPUConst1Contruct2();
bool testADMMGPUConst1Contruct3();

bool testADMMGPUConst1Solve1();
bool testADMMGPUConst1Solve2();
bool testADMMGPUConst1Solve3();


bool testADMMGPUConst1LAMBDA();
bool testADMMGPUConst1Kappa();
bool testADMMGPUConst1Bt1();
bool testADMMGPUConst1CP();
bool testADMMGPUConstCpb();

bool testADMMGPUConst1UpdateRes();
bool testADMMGPUConst1CalcRes();

bool testADMMGPUConst1TradeP();

bool testADMMGPUConst1alpha();
bool testADMMGPUConst1Q();


void testADMMGPUConst1TimeLAMBDA();
void testADMMGPUConst1TimeBt1();
void testADMMGPUConst1TimeTradeP();
void testADMMGPUConst1TimeUpdateRes();
void testADMMGPUConst1TimeCalcRes();

