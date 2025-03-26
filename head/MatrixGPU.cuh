#pragma once

#include <device_launch_parameters.h>
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
#include "MatrixGPUD.cuh"
#include "Utilities.cuh"
#include "Utilities.h"


#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}

#define FULL_MASK 0xffffffff

class MatrixGPU  
{
    int _row;
    int _column;
    bool _GPU = false;
    int _N;
    
    float rand1();
    float* _preallocation = nullptr;
    float* _preallocationFloat = nullptr; // host memory pinned just one float 
   

public:
    float* _matrixCPU = nullptr;
    float* _matrixGPU = nullptr;
    const int _blockSize = 256;
    int _numBlocks = 0;
    bool preallocation = false;
    MatrixGPU();
    MatrixGPU(int line, int column, float value = 0.0, bool pos=false);
    MatrixGPU(const MatrixCPU& m, bool pos=false);
    MatrixGPU(const MatrixGPU& m);
    MatrixGPU(const MatrixGPUD& m);
    
    MatrixGPU& operator= (const MatrixGPU& m);
    MatrixGPU& operator= (const MatrixGPUD& m);
    MatrixGPU& operator= (const MatrixCPU& m);

    void preallocateReduction();

    void transferGPU();
    void transferCPU();

    
    float get(int i, int j, bool verbose=true) const;
    int getNCol() const;
    int getNLin() const;
    void getCol(MatrixGPU* col, int numCol, int offset=0);
    bool getPos() const;
    bool dim(MatrixGPU* m) const;
    bool isEqual(MatrixGPU* m, float pre = 0.000001) const;
    void toMatCPU(MatrixCPU& m) const;
    void toMatGPUD(MatrixGPUD& m) const;


 
    
    void set(int i, int j, float value, bool force=false);
    void setEyes(float value);
    void setEyes(MatrixGPU* m);
    void setRand(float eps); 
    void set(MatrixGPU* m, bool synchrone=true, cudaStream_t stream=0);
    void set(MatrixCPU* m);
    void set(double value);
    void setTrans(MatrixGPU* m);
    void setBloc(int iBegin, int iEnd, int jBegin, int jEnd, MatrixGPU* m);
    void setBloc(int iBegin, int iEnd, int jBegin, int jEnd, MatrixGPU* m, float factor);
    void setBloc(int iBegin, int iEnd, int jBegin, int jEnd, MatrixCPU* m);
    void setBloc(int iBegin, int iEnd, int jBegin, int jEnd, float value);
    void swap(MatrixGPU* m);
    void replace(float previous, float newValue);

    void add(float c);
    void add(MatrixGPU* m1, MatrixGPU* m2); // m = m1 +m2;
    void add(MatrixGPU* m1, float c); // m = m1 + c;
    void add(MatrixGPU* m1);  // m = m + m1;
    void addVector(MatrixGPU* v);
    void addTrans(MatrixGPU* m1); // m = m + tm1;

    void subtract(MatrixGPU* m);
    void subtract(MatrixGPU* m1, MatrixGPU* m2);
    void subtractVector(MatrixGPU* v);
    void subtractTrans(MatrixGPU* m);
    void multiply(float c);
    void multiplyMat(MatrixGPU* A, MatrixGPU* B); // real multiplication A*B
    void multiply(MatrixGPU* Mat, MatrixGPU* vect, bool trans = 0); // par defaut vecteur colonne 
    void MultiplyMatTransVec(MatrixGPU* MatToTrans, MatrixGPU* vect, bool rowVector = 0); // par defaut vecteur colonne
    void linearOperation(MatrixGPU* A, MatrixGPU* x, MatrixGPU* b, bool trans = 0); // Ax + b , vecteur colonne
    void multiplyT(MatrixGPU* m);
    void multiplyT(MatrixGPU* m1, MatrixGPU* m2);
    void divide(float c);
    void divideT(MatrixGPU* m);

    void invertGaussJordan(MatrixGPU* mToInvert);
    void LUPFactorization(MatrixGPU* A, MatrixGPU* P);

    void solveSysUpper(MatrixGPU* U); // Ux = y
    void solveSysLower(MatrixGPU* L, MatrixGPU* b, MatrixGPU* P); // Ly = Pb
    void solveSys(MatrixGPU* A, MatrixGPU* P, MatrixGPU* b); // LUx = Pb , A = (L-I) + U
  
    void Moy(MatrixGPU* m, MatrixGPU* nb, int sens = 0); // moy = moy(m) en consid�rant nb terme par ligne (sens =0) ou par colonne (sens =1)
    void project(MatrixGPU* Lb, MatrixGPU* Ub);
    void projectNeg(); // min(m,0)
    void projectPos(); // max(m,0)
    float sum() const;
    float sum(int begin, int end);
    void sum(MatrixGPU* m); // vecteur colonne
    float distance2();
    float distance2(MatrixGPU* m);
    float max2() const; // renvoie la norme infini 
    float max2(int* indice); // renvoie la norme infini et son emplacement (en 1D)
    float max2(MatrixGPU* m) const;
    void display(bool force=false);
    void displayBloc(int iBegin, int iEnd, int jBegin, int jEnd, bool force = false);
    
    void swapLine(int line1, int line2);

    void saveCSV(const std::string& filename, std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app, int trans = 0) const;
    void saveCSVForce(const std::string& filename, std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app, int trans = 0);

    ~MatrixGPU();
};


__global__ void setup_kernel(curandState* state);


__global__ void generate_kernel(curandState* my_curandstate, float* result, float eps, const unsigned int n );


__global__ void setGPU(float* mat1, float* mat2, int N);
__global__ void setGPUFD(float* mat1, double* mat2, int N);
__global__ void setGPUDF(double* mat1, float* mat2, int N);
__global__ void setGPU(float* mat1, const float value, int N);
__global__ void setGPUunique(float* mat1, const float value, int pos);
__global__ void setTransGPU(float* mat1, float* matToTrans, const int column, const int row);
__global__ void setColGPU(float* mat1, float* mat2, const int numCol, const int column, const int row, const int offset);
__global__ void setEyesGPU(float* mat2, const float value, const int col, const int row);
__global__ void setEyesGPU(float* mat2, float* mat1, const int col, const int row);

__global__ void replaceGPU(float* mat, const float previous, const float newValue, const int N);


__global__ void SetBlocGPU(float* out, float* in, int ibegin, int iend, int jbegin, int jend, int col);
__global__ void SetBlocGPU(float* out, float value, int ibegin, int iend, int jbegin, int jend, int col);
__global__ void SetBlocGPU(float* out, float* in, int ibegin, int iend, int jbegin, int jend, int col, float factor);

__global__ void addGPU(float* mat, float c, int N);
__global__ void addGPU(float* mat1, float* mat2, int N);
__global__ void addGPU(float* mat1, float* mat2, float* mat3, int N);
__global__ void addGPU(float* mat1, float* mat2, float c, int N);
__global__ void addVectorGPU1(float* mat1, float* vect, const int n, int N);
__global__ void addVectorGPU2(float* mat1, float* vect, const int n, int N);
__global__ void addTransGPU(float* out, float* mat1, float* mat2, const int c, const int l, int N);

__global__ void substractGPU(float* mat1, float* mat2, int N);
__global__ void substractGPU(float* mat1, float* mat2, float* mat3, int N);
__global__ void substractVectorGPU1(float* mat1, float* vect, const int n, int N);
__global__ void substractVectorGPU2(float* mat1, float* vect, const int n, int N);
__global__ void substractTransGPU(float* out, float* mat1, float* mat2, const int c, const int l, int N);
__global__ void multiplyGPU(float* mat, const float c, int N);
__global__ void multiplyTGPU(float* mat1, float* mat2, int N);
__global__ void multiplyTGPU(float* mat1, float* mat2, float* mat3, int N);
__global__ void divideGPU(float* mat, const float c, int N);
__global__ void divideGPU(float* mat1, float* mat2, int N);
__global__ void moyGPU1(float* res, float* mat1, float* nb, const int n, const int N);
__global__ void moyGPU2(float* res, float* mat1, float* nb, const int n, const int N);
__global__ void projectGPU(float* mat1, float* lb, float* ub, int N);
__global__ void projectGPUNeg(float* mat1, int N);
__global__ void projectGPUPos(float* mat1, int N);

__global__ void sumGPU(float* res, float* mat1, const int line, const int column);
__global__ void sumGPU2(float* res, float* mat1, const int column);

template <unsigned int blockSize>
__inline__ __device__ void warpReduce(volatile float* sdata, unsigned int tid) {
    if (blockSize >= 64) {
        if (tid < 32) {
            sdata[tid] += sdata[tid + 32];
        }
    }
    __syncwarp();
    if (blockSize >= 32) {
        if (tid < 16) {
            sdata[tid] += sdata[tid + 16];
        }
    }
    __syncwarp();
    if (blockSize >= 16) {
        if (tid < 8) {
            sdata[tid] += sdata[tid +  8];
        }
    }
    __syncwarp();
    if (blockSize >= 8) {
        if (tid < 4) {
            sdata[tid] += sdata[tid +  4];
        }
    }
    __syncwarp();
    if (blockSize >= 4) {
        if (tid < 2) {
            sdata[tid] += sdata[tid +  2];
        }
    }
    __syncwarp();
    if (blockSize >= 2) {
        if (tid < 1) {
            sdata[tid] += sdata[tid +  1];
        }
    }
};/**/


template <unsigned int blockSize>
__inline__ __device__ void warpReduceOr(volatile bool* mustcontinueVect, unsigned int tid) {
    if (blockSize >= 64) {
        if (tid < 32) {
            mustcontinueVect[tid] = mustcontinueVect[tid] || mustcontinueVect[tid + 32];
        }
    }
    __syncwarp();
    if (blockSize >= 32) {
        if (tid < 16) {
            mustcontinueVect[tid] = mustcontinueVect[tid] || mustcontinueVect[tid + 16];
        }
    }
    __syncwarp();
    if (blockSize >= 16) {
        if (tid < 8) {
            mustcontinueVect[tid] = mustcontinueVect[tid] || mustcontinueVect[tid +  8];
        }
    }
    __syncwarp();
    if (blockSize >= 8) {
        if (tid < 4) {
            mustcontinueVect[tid] = mustcontinueVect[tid] || mustcontinueVect[tid +  4];
        }
    }
    __syncwarp();
    if (blockSize >= 4) {
        if (tid < 2) {
            mustcontinueVect[tid] = mustcontinueVect[tid] || mustcontinueVect[tid +  2];
        }
    }
    __syncwarp();
    if (blockSize >= 2) {
        if (tid < 1) {
            mustcontinueVect[tid] = mustcontinueVect[tid] || mustcontinueVect[tid +  1];
        }
    }
};/**/

/*template <unsigned int blockSize>
__device__ void warpReduce(volatile float* sdata, unsigned int tid) {
    float val = sdata[tid];
    if (blockSize >= 64) {
        val += sdata[tid + 32];
        __syncwarp();
    }
    unsigned mask = __ballot_sync(FULL_MASK, tid < 16);
    if (blockSize >= 32) val += __shfl_down_sync(mask, val, 16);
    if (blockSize >= 16) val += __shfl_down_sync(mask, val,  8);
    if (blockSize >= 8)  val += __shfl_down_sync(mask, val,  4);
    if (blockSize >= 4)  val += __shfl_down_sync(mask, val,  2);
    if (blockSize >= 2)  val += __shfl_down_sync(mask, val,  1);
    if (tid == 0) sdata[0] = val;
};
*/
__device__ int sumCommSingleWarp(volatile float* shArr);

template <unsigned int blockSize>
__global__ void sumMonoBlock(float* g_idata, float* g_odata, unsigned int n);

template <unsigned int blockSize>
__global__ void distanceMultiBlock(float* g_idata, float* g_odata, unsigned int n);

template <unsigned int blockSize>
__global__ void distanceMultiBlock(float* g_idata, float* g_idata2, float* g_odata, unsigned int n);

template <unsigned int blockSize>
__global__ void SumMultiBlock(float* g_idata, float* g_odata, unsigned int n);

template <unsigned int blockSize>
__global__ void SumMultiBlock(float* g_idata, float* g_odata, unsigned int begin, unsigned int end);

template <unsigned int blockSize>
__global__ void SumEachRow(float* g_idata, float* g_odata, const int nCol);

/*
unsigned mask = __ballot_sync(FULL_MASK, threadIdx.x < NUM_ELEMENTS);
if (threadIdx.x < NUM_ELEMENTS) {
    val = input[threadIdx.x];
    for (int offset = 16; offset > 0; offset /= 2)
        val += __shfl_down_sync(mask, val, offset);
    …
}

*/

/*template <unsigned int blockSize>
__device__ float warpReduceMax(volatile float* sdata, unsigned int tid) {
    unsigned mask = __ballot_sync(FULL_MASK, tid < 16);
    float val = sdata[tid];
    if (blockSize >= 64) val = val > sdata[tid + 32] ? val : sdata[tid + 32];
    __syncwarp();
    if (blockSize >= 32) val = val > __shfl_down_sync(mask, val, 16) ? val : __shfl_down_sync(mask, val, 16);
    if (blockSize >= 16) val = val > __shfl_down_sync(mask, val,  8) ? val : __shfl_down_sync(mask, val,  8);
    if (blockSize >= 8)  val = val > __shfl_down_sync(mask, val,  4) ? val : __shfl_down_sync(mask, val,  4);
    if (blockSize >= 4)  val = val > __shfl_down_sync(mask, val,  2) ? val : __shfl_down_sync(mask, val,  2);
    if (blockSize >= 2)  val = val > __shfl_down_sync(mask, val,  1) ? val : __shfl_down_sync(mask, val,  1);
    if (tid == 0) sdata[0] = val;
};*/

template <unsigned int blockSize>
__inline__ __device__ float warpReduceMax(volatile float* sdata, unsigned int tid) {
    if (blockSize >= 64) {
        if (tid < 32) {
            sdata[tid] = sdata[tid] > sdata[tid + 32] ? sdata[tid] : sdata[tid + 32];
        }
    }
    __syncwarp();
    if (blockSize >= 32) {
        if (tid < 16) {
            sdata[tid] = sdata[tid] > sdata[tid + 16] ? sdata[tid] : sdata[tid + 16];
        }
    }
    __syncwarp();
    if (blockSize >= 16) {
        if (tid < 8) {
            sdata[tid] = sdata[tid] > sdata[tid + 8] ? sdata[tid] : sdata[tid + 8];
        }
    }
    __syncwarp();
    if (blockSize >= 8) {
        if (tid < 4) {
            sdata[tid] = sdata[tid] > sdata[tid + 4] ? sdata[tid] : sdata[tid + 4];
        }
    }
    __syncwarp();
    if (blockSize >= 4) {
        if (tid < 2) {
            sdata[tid] = sdata[tid] > sdata[tid + 2] ? sdata[tid] : sdata[tid + 2];
        }
    }
    __syncwarp();
    if (blockSize >= 2) {
        if (tid < 1) {
            sdata[tid] = sdata[tid] > sdata[tid + 1] ? sdata[tid] : sdata[tid + 1];
        }
    }
    /*if (blockSize >= 64) sdata[tid] = sdata[tid] > sdata[tid + 32] ? sdata[tid] : sdata[tid + 32];
    __syncwarp();
    if (blockSize >= 32) sdata[tid] = sdata[tid] > sdata[tid + 16] ? sdata[tid] : sdata[tid + 16];
    __syncwarp();
    if (blockSize >= 16) sdata[tid] = sdata[tid] > sdata[tid + 8] ? sdata[tid] : sdata[tid + 8];
    __syncwarp();
    if (blockSize >= 8) sdata[tid] = sdata[tid] > sdata[tid + 4] ? sdata[tid] : sdata[tid + 4];
    __syncwarp();
    if (blockSize >= 4) sdata[tid] = sdata[tid] > sdata[tid + 2] ? sdata[tid] : sdata[tid + 2];
    __syncwarp();
    if (blockSize >= 2) sdata[tid] = sdata[tid] > sdata[tid + 1] ? sdata[tid] : sdata[tid + 1];
    __syncwarp();*/
};

template <unsigned int blockSize>
__device__ void warpReduceMaxPos(volatile float* sdata, volatile int* pos);

template <unsigned int blockSize>
__global__ void maxMonoBlock(float* g_idata, float* g_odata, unsigned int n);

template <unsigned int blockSize>
__global__ void maxMonoBlock(float* g_idata, float* g_odata, unsigned int n, int* pos);

template <unsigned int blockSize>
__global__ void maxMultiBlock(float* g_idata, float* g_odata, unsigned int n);

template <unsigned int blockSize>
__global__ void maxMultiBlock(float* g_idata, float* g_odata, unsigned int n, int* pos);

template <unsigned int blockSize>
__global__ void maxMultiBlock(float* g_idata, float* g_idata2, float* g_odata, unsigned int n);


template <unsigned int blockSize>
__device__ void warpReduceMin(volatile float* sdata, unsigned int tid) {
    if (blockSize >= 64) sdata[tid] = sdata[tid] < sdata[tid + 32] ? sdata[tid] : sdata[tid + 32];
    __syncwarp();
    if (blockSize >= 32) sdata[tid] = sdata[tid] < sdata[tid + 16] ? sdata[tid] : sdata[tid + 16];
    __syncwarp();
    if (blockSize >= 16) sdata[tid] = sdata[tid] < sdata[tid +  8] ? sdata[tid] : sdata[tid +  8];
    __syncwarp();
    if (blockSize >=  8) sdata[tid] = sdata[tid] < sdata[tid +  4] ? sdata[tid] : sdata[tid +  4];
    __syncwarp();
    if (blockSize >=  4) sdata[tid] = sdata[tid] < sdata[tid +  2] ? sdata[tid] : sdata[tid +  2];
    __syncwarp();
    if (blockSize >=  2) sdata[tid] = sdata[tid] < sdata[tid +  1] ? sdata[tid] : sdata[tid +  1];
    __syncwarp();
};




template <unsigned int blockSize>
__global__ void multiplyGPU(float* result, float* Mat, float* vect, int N) {

    // un bloc par ligne
    int thIdx = threadIdx.x;
    int step = blockDim.x;
    int l = blockIdx.x;
    float sum = 0;
    __shared__ float shArr[blockSize];
    for (int i = thIdx; i < N; i += step) {

        float t = Mat[l * N + i] * vect[i];
        sum += t;
    }

    shArr[thIdx] = sum;
    __syncthreads();

    if (blockSize >= 512) { if (thIdx < 256) { shArr[thIdx] += shArr[thIdx + 256]; } __syncthreads(); }
    if (blockSize >= 256) { if (thIdx < 128) { shArr[thIdx] += shArr[thIdx + 128]; } __syncthreads(); }
    if (blockSize >= 128) { if (thIdx < 64) { shArr[thIdx] += shArr[thIdx + 64]; } __syncthreads(); }
    if (blockSize >= 64) {
        if (thIdx < 32) {
            warpReduce<blockSize>(shArr, thIdx);
        }
    }
    if (blockSize >= 32 && blockSize < 64) {
        if (thIdx < 16) {
            warpReduce<blockSize>(shArr, thIdx);
        }
    }
        
    __syncthreads();
    if (thIdx == 0) {
        result[l] = shArr[0];
    }
}

template <unsigned int blockSize>
__global__ void multiplyGPUMatVectTrans(float* result, float* Mat, float* vect, int N) {

    // un bloc par ligne
    int thIdx = threadIdx.x;
    int step = blockDim.x;
    int l = blockIdx.x;
    float sum = 0;
    __shared__ float shArr[blockSize];
    for (int i = thIdx; i < N; i += step) {

        float t = Mat[i * N + l] * vect[i]; // accès matrice pas coalescent !!!
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
__global__ void linearOpGPU(float* result, float* A, float* x, float* b, int N) {

    // un bloc par element de result
    int thIdx = threadIdx.x;
    int step = blockDim.x;
    int l = blockIdx.x;
    float sum = 0;
    __shared__ float shArr[blockSize];
    for (int i = thIdx; i < N; i += step) {

        float t = A[l * N + i] * x[i];
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

__global__ void normalisationGJ(float* mat, const int row, const int nCol, const float factor);

__global__ void swapLineGJ(float* mat, const int row1, const int row2, const int nCol);

__global__ void eliminationGJ(float* mat,float* matAug, const int r, const int nRow, const int nCol);


__global__ void initPermMatr(float* P, const int N);
__global__ void updatePermMatr(float* P, const int line1, const int line2,  const int N);
__global__ void updateLUPFactorization(float* A, const int col, const int N);

__global__ void setPermute(float* y, float* b, float* P, const int N);
__global__ void solveLowSys(float* A, float* y, const int iter, const int N);
__global__ void solveUpSys(float* A, float* y, const int iter, const int N);

__global__ void solveSysGPU(float* A, float* y, const int N);