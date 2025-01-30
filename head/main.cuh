#pragma once




// utilile
#include "MatrixCPU.h"
#include "MatrixCPUD.h"
#include "MatrixGPU.cuh"
#include "MatrixGPUD.cuh"
#include "StudyCase.h"
#include "Simparam.h"
#include "System.cuh"
#include "System.h"

//Market
#include "ADMMMarket.h"
#include "ADMMMarketOpenMP.h"
#ifdef OSQP
    #include "OSQP.h" 
    #include "OSQPCentralized2.h"
    #include <osqp.h>
#endif

#include "PAC.h" 
#include "PACOpenMP.h" 
#include "PACGPU.cuh"
#include "ADMMMarketGPU.cuh"


//PF
#include "CPUPF.h"
#include "GPUPF.cuh"
#include "CPUPFGS.h"
#include "GPUPFGS.cuh"
#include "CPUPFdist.h"
#include "CPUPFdistPQ.h"
#include "GPUPFdistPQ.cuh"

//OPF
#include "OPFADMM.h"
#include "OPFADMM2.h"
#include "OPFADMMGPU.cuh"
#include "OPFADMMGPU2.cuh"
#include "OPFADMMCons.h" // fonction cout pour les pertes !!!
#include "OPFADMMConsGPU.cuh" // fonction cout pour les pertes !!!

// Market Endo
#include "MarEndoCons.h"
#include "EndoPF.h"
#include "MarketEndoDirect.h"
#include "EndoPFGPU.cuh"
#include "MarketEndoDirectGPU.cuh"
#include "MarEndoConsGPU.cuh"


#include <stdio.h>
#include <iostream>
#include <time.h>
#include <cuda_runtime.h>
#include <algorithm>
#include <iterator>
#include <random>
#include <vector>
#include <cudaProfiler.h>





// fichier de test
#include "TestMatrixCPU.h"
#include "TestAgent.h"
#include "TestStudyCase.h"
#include "TestADMMConst.h"
#include "TestADMMConst1.h"

#include "TestADMMGPUConst1.cuh"
#include "TestADMMGPUConst1T.cuh"
#include "TestADMMGPUConst2.cuh"
#include "TestADMMGPUConst3.cuh"
#include "TestPAC.h"
#include "TestPACConst.h"
#include "TestSimparam.h"
#include "TestSysteme.h"
#include "TestMatrixGPU.cuh"
#include "TestMatrixCPU.h"
#include "TestKernel.cuh"
#include "TestUtilities.cuh"




// test on the global and local step impact
void SimulationTempStepOpti();

void comparaisonArticle();
void testExemple();


// Performance test, alea case or europeen case.
void SimuTemporal();

void SimuTemporal(std::string name);
void SimuTemporalRho(std::string name = "Europe");
void SimuTemporalTestFeeder();
void SimuTemporalLlimit(std::string name = "Europe");
void SimuTemporalConvergence(std::string name = "Europe");
void SimuTemporalWOConstraint(std::string name);



void SimuCompare();
void SimuCompareParra();
void SimuStatPowerTech();
void SimuStat();
void SimuStatOPF();
void SimuSensiStudyCase();


// PF
void SimuStatPFSGE();
void SimuStatPFTransport();
void SimuStatPFCompare();
// fin PF

// OPF
void SimuStatOPFCompare();
void SimuCompareISGT();
// fin OPF

// Market Endo
void SimuStatMarketEndo();
void SimuStatMarketEndoAC();
void SimuStatADMMGPU();
void SimuTemporalTestFeederEndo();
void SimuTemporalTestFeederEndoAll();
void SimuStatMarketEndoACAgent();
void SimuStatMarketEndoGrid();

void testADMMGPUtemp();


void testCPUPF();
void testOPF();
void testADMMACConst();
void testMarket();
void testMarketEndo();

void getErrorSensiPwerLine();

//utilitaire
float pow10(int n);
int getNFileline(std::string nameFile);
std::string generateDate(int year, int month, int day, int hour);

/*
On ne peut inclure qu'un seul osqp � la fois, (car m�me nom)
C:\Program Files \OSQP\osqp\include\osqp


*/