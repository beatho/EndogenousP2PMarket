#pragma once

#include <iostream>
#include <math.h>

#ifdef Eigen
    #include <unsupported/Eigen/Polynomials>
#endif
#include <chrono>


#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

#define CHECK_CUDA_ERROR(val) check((val), #val, __FILE__, __LINE__)
#define CHECK_LAST_CUDA_ERROR() checkLast(__FILE__, __LINE__)

#define BILLION 1000000000

#define MILLION 1000000

float min(float a, float b);

// am�lioration : passer MatrixCPU en param�tre (au moins c'est safe)

int resolveRealPolynome3without2term(double* root, double* coef);
int resolveRealPolynome3Newton(double* root, double* coef, double init = 0);
int resolveRealPolynome3Laguerre(double* root, double* coef, double init = 0);
int resolveRealPolynome3Halley(double* root, double* coef, double init = 0);

int resvolveRealPolynome4without2term(double* root, double* coef);
int resvolveRealPolynome4without2termLagrange(double* root, double* coef);
int resvolveRealPolynome4without2term(double* root, double* coef, bool Lagrange); // selectionne la m�thode

int resolveRealPolynome4Newton(double* root, double* coef, double init = 0);
int resolveRealPolynome4Halley(double* root, double* coef, double init = 0);


void coefPolynome3From4to2coef(double* coef4, double* coef2);

double findAntpoly3Neg(double p, double q);
double findAntpoly3Pos(double p, double q);


// EIGEN
#ifdef EIGEN
    int resolveRealPolynome3without2termEigen(double* root, double* coef);

    int resvolveRealPolynome4without2termEigen(double* root, double* coef);

#endif

