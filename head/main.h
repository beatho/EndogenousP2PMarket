#pragma once

// Standard

#include <stdio.h>
#include <iostream>
#include <time.h>

#include <algorithm>
#include <iterator>
#include <random>
#include <vector>



#include "MatrixCPU.h"
#include "MatrixCPUD.h"
#include "StudyCase.h"
#include "Simparam.h"
#include "System.h"
#include "Utilities.h"

// Market

#include "ADMMMarket.h"
#include "ADMMMarketOpenMP.h"
#ifdef OSQP
    #include "OSQP.h" 
    #include "OSQPCentralized2.h"
    #include <osqp.h>
#endif

#include "PAC.h" 
#include "PACOpenMP.h" 

//  PF
#include "CPUPF.h"
#include "CPUPFGS.h"
#include "CPUPFdist.h"
#include "CPUPFdistPQ.h"

// OPF
#include "OPFADMM.h"
#include "OPFADMM2.h"
#include "OPFADMMCons.h" // fonction cout pour les pertes !!!

// Endo
#include "MarEndoCons.h"
#include "EndoPF.h"
#include "MarketEndoDirect.h"


// test
#include "TestMatrixCPU.h"
#include "TestAgent.h"
#include "TestStudyCase.h"
#include "TestADMMConst.h"


#include "TestPAC.h"
#include "TestPACConst.h"
#include "TestSimparam.h"
#include "TestSysteme.h"

#include "TestMatrixCPU.h"

int main2();

void testMarketEndo(int numCase, std::string caseName);