#pragma once


#include <cuda_runtime.h>

#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

#define CHECK_CUDA_ERROR(val) check((val), #val, __FILE__, __LINE__)
#define CHECK_LAST_CUDA_ERROR() checkLast(__FILE__, __LINE__)

template <typename T>
void check(T err, const char* const func, const char* const file, const int line);
void checkLast(const char* const file, const int line);

__device__ int resolveRealPolynome3without2termGPU( double* root,  double p,  double q);

__device__ int resvolveRealPolynome4without2termGPU( double* root,  double b,  double d,  double e);
__device__ int resvolveRealPolynome4without2termGPULagrange(double* root, double b, double d, double e);
__device__ int resvolveRealPolynome4without2termGPU(double* root, double b, double d, double e, bool Lagrange);

__device__ void coefPolynome3From4to2coefGPU(double* coef4, double* coef2);

__global__ void resolveSeveralRealPolynome3termGPU(double* nRoot, double* roots, double* coefs, int nPoly);

__global__ void resolveSeveralRealPolynome4WO2termGPU(double* nRoot, double* roots, double* coefs, int nPoly);
__global__ void resolveSeveralRealPolynome4WO2termGPULagrange(double* nRoot, double* roots, double* coefs, int nPoly);


__device__ int resolveRealPolynome3GPU(double* root, double b, double c, double d);

__device__ int resvolveRealPolynome4GPU(double* root, double b, double c, double d, double e);


__global__ void resolveSeveralRealPolynome3GPU(double* nRoot, double* roots, double* coefs, int nPoly);

__global__ void resolveSeveralRealPolynome4GPU(double* nRoot, double* roots, double* coefs, int nPoly);

__host__ __device__ double findAntpoly3Neg(double p, double q);
__host__ __device__ double findAntpoly3Pos(double p, double q);

