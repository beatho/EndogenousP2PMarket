

#include <device_launch_parameters.h>
#include "Utilities.cuh"
#include "Utilities.h"
#include "MatrixGPU.cuh"




int testUtilities();

void compareCPUGPU();



bool testcoefPolynome3From4to2coef1();
bool testcoefPolynome3From4to2coef2();
bool testcoefPolynome3From4to2coef3();

bool testresolveRealPolynome3without2term1();
bool testresolveRealPolynome3without2term2();
bool testresolveRealPolynome3without2term3();


bool testresolveRealPolynome4without2term();
bool testresolveRealPolynome4without2term2();

bool testresolveRealPolynome4without2termLagrange();
bool testresolveRealPolynome4without2term2Lagrange();

bool testresolveRealPolynome3without2termGPU();
bool testresolveRealPolynome4without2termGPU();
bool testresolveRealPolynome4without2termGPULagrange();

#ifdef EIGEN

bool testPolyEigen3();
bool testPolyEigen4();

bool testresolveRealPolynome3without2termEigen();
bool testresolveRealPolynome4without2termEigen();
#endif

bool testresolveRealPolynome3Newton1();
bool testresolveRealPolynome3Newton2();
bool testresolveRealPolynome3Newton3();



bool testresolveRealPolynome4Newton1();
bool testresolveRealPolynome4Newton2();


bool testresolveRealPolynome3Laguerre1();
bool testresolveRealPolynome3Laguerre2();
bool testresolveRealPolynome3Laguerre3();

bool testresolveRealPolynome3Halley1();
bool testresolveRealPolynome3Halley2();
bool testresolveRealPolynome3Halley3();

bool testresolveRealPolynome4Halley1();
bool testresolveRealPolynome4Halley2();


bool testresolveRealPolynome3GPU();
bool testresolveRealPolynome4GPU();