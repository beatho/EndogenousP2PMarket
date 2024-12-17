#pragma once


#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <Eigen\Dense>
#include <glob_opts.h>
#include <osqp.h>
#include "MatrixCPUD.h"


#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}

class MatrixCPU  
{
    int _row = 0;
    int _column = 0;
    float rand1();

public:
    float* _matrixCPU = nullptr;
   
    MatrixCPU();
    MatrixCPU(int line, int column, float value = 0.0);
    MatrixCPU(const MatrixCPU & m);
    MatrixCPU(const MatrixCPUD & m);
    MatrixCPU& operator= (const MatrixCPU& m);
    MatrixCPU& operator= (const MatrixCPUD& m);

    void toMatCPUD(MatrixCPUD& m) const;

   
    inline float get(int i, int j) const;
    double getD(int i, int j) const;
    int getNCol() const;
    int getNLin() const;
    void getLin(MatrixCPU* line, int i) const;
    bool dim(MatrixCPU* m) const;
    bool isEqual(MatrixCPU* m, float pre = 0.000001) const;
    int getNNull() const;
    int getNNullHalf() const;
    void swap(MatrixCPU* m);
    void getBloc(MatrixCPU* dest, int iBegin, int iEnd, int jBegin, int jEnd);
    
   
    void setSize(int l, int c);
    inline void set(int i, int j, float value);
    void set(float value);
    void set(MatrixCPU* m);
    void setTrans(MatrixCPU* m);
    void set(Eigen::MatrixXd* eigenMatrix);
    void setBloc(int iBegin, int iEnd, int jBegin, int jEnd, MatrixCPU* m);
    void setEyes(float v);
    void setEyes(MatrixCPU* vect);
    void setRand(int eps, int divide); // rand()%eps/divide
    void setRand1(float eps);
    void setFromFile(std::string file,int entete=0);
    void add(float c);
    void add(MatrixCPU* m1, MatrixCPU* m2); // m = m1 +m2;
    void add(MatrixCPU* m1, float c); // m = m1 + c;
    void add(MatrixCPU* m1);  // m = m + m1;
    void increment(int i, int j, float add); //m[i,j] = m[i,j] + add
    void addVector(MatrixCPU* v);
    void addTrans(MatrixCPU* m1); // m = m + tm1;
    
    void subtract(MatrixCPU* m); // this - m
    void subtract(MatrixCPU* m1, MatrixCPU* m2); // m1 - m2
    void subtractRow(int row1, int row2, float factor); // row1 - factor*row2
    void subtractAbs(MatrixCPU* m1, MatrixCPU* m2);
    void subtractVector(MatrixCPU* v);
    void subtractTrans(MatrixCPU* m);

    void multiply(MatrixCPU* m1, MatrixCPU* m2);
    void multiplyTrans(MatrixCPU* m1, MatrixCPU* m2, int secondToTrans=1);
    void multiplyDiag(MatrixCPU* m1, MatrixCPU* m2); // multiply but only keep the diagonal
    void multiply(float c);
    void multiplyT(MatrixCPU* m);
    void multiplyT(MatrixCPU* m1, MatrixCPU* m2);
    void multiplyTVector(MatrixCPU* m, MatrixCPU* v, int sens=0);
    void divide(float c);
    void divideT(MatrixCPU* m);
    
    void invertGaussJordan(MatrixCPU* mToInvert);
    void invertEigen(MatrixCPU* mToInvert);
    void LDLFactorization(MatrixCPU* L, MatrixCPU* D);
    void LUPFactorization(MatrixCPU* A, MatrixCPU* P);

    void solveSysUpper(MatrixCPU* U, MatrixCPU* y); // Ux = y
    void solveSysLower(MatrixCPU* L, MatrixCPU* b , MatrixCPU* P); // Ly = Pb
    void solveSys(MatrixCPU* A, MatrixCPU* P, MatrixCPU* b); // LUx = Pb , A = (L-I) + U
    void solveSysEigen(MatrixCPU* M, MatrixCPU* b); // Mx = b


    void MultiplyMatVec(MatrixCPU* m, MatrixCPU* vect, int sens = 0);
    void MultiplyMatTransVec(MatrixCPU* mToTrans, MatrixCPU* vect, int sensRow = 0);
    void MultiplyMatMat(MatrixCPU* m1, MatrixCPU* m2);


    
    void project(MatrixCPU* Lb, MatrixCPU* Ub);
    void projectNeg(); // min(m,0)
    void projectPos(); // max(m,0)


    float distance2(MatrixCPU* m1) const; // ||m-m1||_2
    float distance2() const; // ||m||_2 
    float min2() const;
    float min2Nnull(float eps = 0) const; // non nul
    float max2() const; // ||m||_inf 
    float max2(MatrixCPU* m2); // ||m-m2||_inf
    float maxAbs(int iBegin, int iEnd, int jBegin, int jEnd, MatrixCPU* indices = nullptr);
    float minAbs(int iBegin, int iEnd, int jBegin, int jEnd, MatrixCPU* indices = nullptr, bool Null=true);
    void Moy(MatrixCPU* m, MatrixCPU* nb, int sens=0); 
    float sum() const;
    float sumD() const;
    void sum(MatrixCPU* m, int sens = 0);
    void sumT(MatrixCPU* m, int sens = 0); // if matrix and vector not in the same sens
    void sort(int dim = 0, int sens = 0);
    void sortLine(int line1, int line2, int dim);
    void sortColumn(int col1, int col2, int dim);
    void fusionLine(int col1, int col2, int dim);
    void fusionColumn(int col1, int col2, int dim);
    void swapLine(int line1, int line2);
    void swapColumn(int col1, int col2);

    void RelativeEror(MatrixCPU* MatRef, MatrixCPU* Mat); 


    c_float* toCFloat();
    c_int toCSC(c_float* data, c_int* idx, c_int* ptr);
    c_int toCSCHalf(c_float* data, c_int* idx, c_int* ptr);
    void toEigenMatrix(Eigen::MatrixXd* eigenMatrix);

    void display() const;
    void displayBloc(int iBegin, int iEnd, int jBegin, int jEnd) const;
    void saveCSV(const std::string& filename, std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app, int trans=0, std::string separator = ";") const;

    static void saveTabMatCSV(MatrixCPU* tab, int sizeTab, const std::string& filename, std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app, int trans = 0, std::string separator = ";");
    ~MatrixCPU();
};

