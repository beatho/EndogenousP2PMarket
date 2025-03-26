#pragma once
#include <device_launch_parameters.h>
#include "ADMMGPUConst2.cuh"
#include <math.h>
#include <chrono>

int testADMMGPUConst2();
void testADMMGPUConst2Time(int test=0);

bool testADMMGPUConst2Contruct1();
bool testADMMGPUConst2Contruct2();
bool testADMMGPUConst2Contruct3();

bool testADMMGPUConst2Solve1();
bool testADMMGPUConst2Solve2();
bool testADMMGPUConst2LAMBDA();
bool testADMMGPUConst2Bt1();

bool testADMMGPUConst2UpdateRes();
bool testADMMGPUConst2CalcRes();

bool testADMMGPUConst2TradeP();

void testADMMGPUConst2TimeLAMBDA();
void testADMMGPUConst2TimeBt1();
void testADMMGPUConst2TimeTradeP();
void testADMMGPUConst2TimeUpdateRes();
void testADMMGPUConst2TimeCalcRes();

