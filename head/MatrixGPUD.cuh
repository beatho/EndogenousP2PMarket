#pragma once

#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif


#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <curand.h>
#include <curand_kernel.h>

#include <cuda_runtime.h>

#include <device_launch_parameters.h>


#include "MatrixCPU.h"

#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}

class MatrixGPUD  
{
    int _row;
    int _column;
    bool _GPU = false;
    int _N;
    
    double rand1();
    double* _preallocation = nullptr;
    double* _preallocationFloat = nullptr; // host memory pinned just one double 
   

public:
    double* _matrixCPU = nullptr;
    double* _matrixGPU = nullptr;
    const int _blockSize = 256;
    int _numBlocks = 0;
    bool preallocation = false;
    MatrixGPUD();
    MatrixGPUD(int line, int column, double value = 0.0, bool pos=false);
    MatrixGPUD(const MatrixCPUD& m, bool pos=false);
    MatrixGPUD(const MatrixCPU& m, bool pos = false);
    MatrixGPUD(const MatrixGPUD& m);
    MatrixGPUD& operator= (const MatrixGPUD& m);
    MatrixGPUD& operator= (const MatrixCPUD& m);

    void preallocateReduction();

    void transferGPU();
    void transferCPU();

    
    double get(int i, int j, bool verbose=true) const;
    int getNCol() const;
    int getNLin() const;
    void getCol(MatrixGPUD* col, int numCol, int offset=0);
    bool getPos() const;
    bool dim(MatrixGPUD* m) const;
    bool isEqual(MatrixGPUD* m, double pre = 0.000001) const;
    void toMatCPU(MatrixCPUD& m) const;


    void setSize(int row, int column);
    
    void set(int i, int j, double value, bool force=false);
    void setEyes(double value);
    void setEyes(MatrixGPUD* m);
    void setRand(double eps); 
    void set(MatrixGPUD* m, bool synchrone=true, cudaStream_t stream=0);
    void set(MatrixCPUD* m);
    void setTrans(MatrixGPUD* m);
    void setBloc(int iBegin, int iEnd, int jBegin, int jEnd, MatrixGPUD* m);
    void setBloc(int iBegin, int iEnd, int jBegin, int jEnd, MatrixGPUD* m, double factor);
    void setBloc(int iBegin, int iEnd, int jBegin, int jEnd, MatrixCPUD* m);
    void swap(MatrixGPUD* m);
    void replace(double previous, double newValue);

    void add(double c);
    void add(MatrixGPUD* m1, MatrixGPUD* m2); // m = m1 +m2;
    void add(MatrixGPUD* m1, double c); // m = m1 + c;
    void add(MatrixGPUD* m1);  // m = m + m1;
    void addVector(MatrixGPUD* v);
    void addTrans(MatrixGPUD* m1); // m = m + tm1;

    void subtract(MatrixGPUD* m);
    void subtract(MatrixGPUD* m1, MatrixGPUD* m2);
    void subtractVector(MatrixGPUD* v);
    void subtractTrans(MatrixGPUD* m);
    void multiply(double c);
    void multiplyMat(MatrixGPUD* A, MatrixGPUD* B); // real multiplication A*B
    void multiply(MatrixGPUD* Mat, MatrixGPUD* vect, bool trans = 0); // par defaut vecteur colonne 
    void linearOperation(MatrixGPUD* A, MatrixGPUD* x, MatrixGPUD* b, bool trans = 0); // Ax + b , vecteur colonne
    void multiplyT(MatrixGPUD* m);
    void multiplyT(MatrixGPUD* m1, MatrixGPUD* m2);
    void divide(double c);
    void divideT(MatrixGPUD* m);

    void invertGaussJordan(MatrixGPUD* mToInvert);
    void LUPFactorization(MatrixGPUD* A, MatrixGPUD* P);

    void solveSysUpper(MatrixGPUD* U); // Ux = y
    void solveSysLower(MatrixGPUD* L, MatrixGPUD* b, MatrixGPUD* P); // Ly = Pb
    void solveSys(MatrixGPUD* A, MatrixGPUD* P, MatrixGPUD* b); // LUx = Pb , A = (L-I) + U
  
    void Moy(MatrixGPUD* m, MatrixGPUD* nb, int sens = 0); // moy = moy(m) en consid�rant nb terme par ligne (sens =0) ou par colonne (sens =1)
    void project(MatrixGPUD* Lb, MatrixGPUD* Ub);
    void projectNeg(); // min(m,0)
    void projectPos(); // max(m,0)
    double sum() const;
    double sum(int begin, int end);
    void sum(MatrixGPUD* m); // vecteur colonne
    double distance2();
    double distance2(MatrixGPUD* m);
    double max2() const; // renvoie la norme infini 
    double max2(int* indice); // renvoie la norme infini et son emplacement (en 1D)
    double max2(MatrixGPUD* m) const;
    void display(bool force=false);
    void displayBloc(int iBegin, int iEnd, int jBegin, int jEnd, bool force = false);
    
    void swapLine(int line1, int line2);

    void saveCSV(const std::string& filename, std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app, int trans = 0) const;
    void saveCSVForce(const std::string& filename, std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app, int trans = 0);

    ~MatrixGPUD();
};


__global__ void setup_kernelD(curandState* state);


__global__ void generate_kernel(curandState* my_curandstate, double* result, double eps, const unsigned int n );


__global__ void setGPU(double* mat1, double* mat2, int N);
__global__ void setGPU(double* mat1, const double value, int N);
__global__ void setGPUunique(double* mat1, const double value, int pos);
__global__ void setTransGPU(double* mat1, double* matToTrans, const int column, const int row);
__global__ void setColGPU(double* mat1, double* mat2, const int numCol, const int column, const int row, const int offset);
__global__ void setEyesGPU(double* mat2, const double value, const int col, const int row);
__global__ void setEyesGPU(double* mat2, double* mat1, const int col, const int row);

__global__ void replaceGPU(double* mat, const double previous, const double newValue, const int N);


__global__ void SetBlocGPU(double* out, double* in, int ibegin, int iend, int jbegin, int jend, int col);
__global__ void SetBlocGPU(double* out, double* in, int ibegin, int iend, int jbegin, int jend, int col, double factor);

__global__ void addGPU(double* mat, double c, int N);
__global__ void addGPU(double* mat1, double* mat2, int N);
__global__ void addGPU(double* mat1, double* mat2, double* mat3, int N);
__global__ void addGPU(double* mat1, double* mat2, double c, int N);
__global__ void addVectorGPU1(double* mat1, double* vect, const int n, int N);
__global__ void addVectorGPU2(double* mat1, double* vect, const int n, int N);
__global__ void addTransGPU(double* out, double* mat1, double* mat2, const int c, const int l, int N);

__global__ void substractGPU(double* mat1, double* mat2, int N);
__global__ void substractGPU(double* mat1, double* mat2, double* mat3, int N);
__global__ void substractVectorGPU1(double* mat1, double* vect, const int n, int N);
__global__ void substractVectorGPU2(double* mat1, double* vect, const int n, int N);
__global__ void substractTransGPU(double* out, double* mat1, double* mat2, const int c, const int l, int N);
__global__ void multiplyGPU(double* mat, const double c, int N);
__global__ void multiplyTGPU(double* mat1, double* mat2, int N);
__global__ void multiplyTGPU(double* mat1, double* mat2, double* mat3, int N);
__global__ void divideGPU(double* mat, const double c, int N);
__global__ void divideGPU(double* mat1, double* mat2, int N);
__global__ void moyGPU1(double* res, double* mat1, double* nb, const int n, const int N);
__global__ void moyGPU2(double* res, double* mat1, double* nb, const int n, const int N);
__global__ void projectGPU(double* mat1, double* lb, double* ub, int N);
__global__ void projectGPUNeg(double* mat1, int N);
__global__ void projectGPUPos(double* mat1, int N);

__global__ void sumGPU(double* res, double* mat1, const int line, const int column);
__global__ void sumGPU2(double* res, double* mat1, const int column);

template <unsigned int blockSize>
__device__ void warpReduce(volatile double* sdata, unsigned int tid) {
    if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
    if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
    if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
    if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
    if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
    if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
};


__device__ int sumCommSingleWarp(volatile double* shArr);

template <unsigned int blockSize>
__global__ void sumMonoBlock(double* g_idata, double* g_odata, unsigned int n);

template <unsigned int blockSize>
__global__ void SumMultiBlock(double* g_idata, double* g_odata, unsigned int begin, unsigned int end);

template <unsigned int blockSize>
__global__ void distanceMultiBlock(double* g_idata, double* g_odata, unsigned int n);

template <unsigned int blockSize>
__global__ void distanceMultiBlock(double* g_idata, double* g_idata2, double* g_odata, unsigned int n);

template <unsigned int blockSize>
__global__ void SumMultiBlock(double* g_idata, double* g_odata, unsigned int n);

template <unsigned int blockSize>
__global__ void SumEachRow(double* g_idata, double* g_odata, const int nCol);

__device__ double warpReduceMax(volatile double* sdata);

__device__ void warpReduceMax(volatile double* sdata, volatile int* pos);

template <unsigned int blockSize>
__global__ void maxMonoBlock(double* g_idata, double* g_odata, unsigned int n);

template <unsigned int blockSize>
__global__ void maxMonoBlock(double* g_idata, double* g_odata, unsigned int n, int* pos);

template <unsigned int blockSize>
__global__ void maxMultiBlock(double* g_idata, double* g_odata, unsigned int n);

template <unsigned int blockSize>
__global__ void maxMultiBlock(double* g_idata, double* g_odata, unsigned int n, int* pos);

template <unsigned int blockSize>
__global__ void maxMultiBlock(double* g_idata, double* g_idata2, double* g_odata, unsigned int n);


template <unsigned int blockSize>
__device__ void warpReduceMin(volatile double* sdata, unsigned int tid) {
    if (blockSize >= 64) sdata[tid] = sdata[tid] < sdata[tid + 32] ? sdata[tid] : sdata[tid + 32];
    if (blockSize >= 32) sdata[tid] = sdata[tid] < sdata[tid + 16] ? sdata[tid] : sdata[tid + 16];
    if (blockSize >= 16) sdata[tid] = sdata[tid] < sdata[tid +  8] ? sdata[tid] : sdata[tid +  8];
    if (blockSize >=  8) sdata[tid] = sdata[tid] < sdata[tid +  4] ? sdata[tid] : sdata[tid +  4];
    if (blockSize >=  4) sdata[tid] = sdata[tid] < sdata[tid +  2] ? sdata[tid] : sdata[tid +  2];
    if (blockSize >=  2) sdata[tid] = sdata[tid] < sdata[tid +  1] ? sdata[tid] : sdata[tid +  1];
};




template <unsigned int blockSize>
__global__ void multiplyGPU(double* result, double* Mat, double* vect, int N) {

    // un bloc par ligne
    int thIdx = threadIdx.x;
    int step = blockDim.x;
    int l = blockIdx.x;
    double sum = 0;
    __shared__ double shArr[blockSize];
    for (int i = thIdx; i < N; i += step) {

        double t = Mat[l * N + i] * vect[i];
        sum += t;
    }

    shArr[thIdx] = sum;
    __syncthreads();

    if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
    if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
    if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
    if (thIdx < 32) {
        warpReduce<blockSize>(shArr, thIdx);
    }

    if (thIdx == 0) {
        result[l] = shArr[0];
    }
}


template <unsigned int blockSize>
__global__ void linearOpGPU(double* result, double* A, double* x, double* b, int N) {

    // un bloc par element de result
    int thIdx = threadIdx.x;
    int step = blockDim.x;
    int l = blockIdx.x;
    double sum = 0;
    __shared__ double shArr[blockSize];
    for (int i = thIdx; i < N; i += step) {

        double t = A[l * N + i] * x[i];
        sum += t;
    }

    shArr[thIdx] = sum;
    __syncthreads();

    if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
    if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
    if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
    if (thIdx < 32) {
        warpReduce<blockSize>(shArr, thIdx);
    }

    if (thIdx == 0) {
        result[l] = shArr[0] + b[l];
    }
}
/// for GJ inversion

__global__ void normalisationGJ(double* mat, const int row, const int nCol, const double factor);

__global__ void swapLineGJ(double* mat, const int row1, const int row2, const int nCol);

__global__ void eliminationGJ(double* mat,double* matAug, const int r, const int nRow, const int nCol);


__global__ void initPermMatr(double* P, const int N);
__global__ void updatePermMatr(double* P, const int line1, const int line2, const int N);
__global__ void updateLUPFactorization(double* A, const int col, const int N);

__global__ void setPermute(double* y, double* b, double* P, const int N);
__global__ void solveLowSys(double* A, double* y, const int iter, const int N);
__global__ void solveUpSys(double* A, double* y, const int iter, const int N);

__global__ void solveSysGPU(double* A, double* y, const int N);