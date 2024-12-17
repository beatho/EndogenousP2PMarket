#pragma once
#include "ADMMGPUConst3.cuh"
#include <math.h>
#include <chrono>

int testADMMGPUConst3();
void testADMMGPUConst3Time(int test=0);

bool testADMMGPUConst3Contruct1();
bool testADMMGPUConst3Contruct2();
bool testADMMGPUConst3Contruct3();

bool testADMMGPUConst3Solve1();
bool testADMMGPUConst3Solve2();
bool testADMMGPUConst3LAMBDA();
bool testADMMGPUConst3Bt1();

bool testADMMGPUConst3UpdateRes();
bool testADMMGPUConst3CalcRes();

bool testADMMGPUConst3TradeP();

void testADMMGPUConst3TimeLAMBDA();
void testADMMGPUConst3TimeBt1();
void testADMMGPUConst3TimeTradeP();
void testADMMGPUConst3TimeUpdateRes();
void testADMMGPUConst3TimeCalcRes();

