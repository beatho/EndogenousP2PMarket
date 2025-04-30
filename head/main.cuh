// utilile

#include "MatrixGPU.cuh"
#include "MatrixGPUD.cuh"
#include "System.cuh"


//Market

#include "PACGPU.cuh"
#include "ADMMMarketGPU.cuh"


//PF
#include "GPUPFdistPQ.cuh"
#include "GPUPF.cuh"
#include "GPUPFGS.cuh"
//OPF

#include "OPFADMMGPU.cuh"
#include "OPFADMMGPU2.cuh"
#include "OPFADMMConsGPU.cuh" // fonction cout pour les pertes !!!

// Market Endo
#include "EndoPFGPU.cuh"
#include "MarketEndoDirectGPU.cuh"
#include "MarEndoConsGPU.cuh"

#include <cudaProfiler.h>
#include <cuda_runtime.h>





// fichier de test
#include "TestADMMGPUConst1.cuh"
#include "TestADMMGPUConst1T.cuh"
#include "TestADMMGPUConst2.cuh"
#include "TestADMMGPUConst3.cuh"
#include "TestKernel.cuh"
#include "TestUtilities.cuh"
#include "TestMatrixGPU.cuh"

#include "TestMatrixCPU.h"
#include "TestAgent.h"
#include "TestStudyCase.h"
#include "TestADMMConst.h"


#include "TestPAC.h"
#include "TestPACConst.h"
#include "TestSimparam.h"
#include "TestSysteme.h"

#include "TestMatrixCPU.h"

#include "Utilities.h"





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