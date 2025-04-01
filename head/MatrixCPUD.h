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
#ifdef EIGEN
    #include <Eigen\Dense>
#endif

#ifdef OSQP
    #include <glob_opts.h>
    #include <osqp.h>
#endif


#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}

class MatrixCPUD  
{
    int _row = 0;
    int _column = 0;
    double rand1();

public:
    double* _matrixCPU = nullptr;
   
    MatrixCPUD();
    MatrixCPUD(int line, int column, double value = 0.0);
    MatrixCPUD(const MatrixCPUD & m);
    MatrixCPUD& operator= (const MatrixCPUD& m);

   
    double get(int i, int j) const;
    
    int getNCol() const;
    int getNLin() const;
    void getLin(MatrixCPUD* line, int i) const;
    bool dim(MatrixCPUD* m) const;
    bool isEqual(MatrixCPUD* m, double pre = 0.000001) const;
    int getNNull() const;
    int getNNullHalf() const;
    void swap(MatrixCPUD* m);
    void getBloc(MatrixCPUD* dest, int iBegin, int iEnd, int jBegin, int jEnd);
    
   
    void setSize(int l, int c);
    void set(int i, int j, double value);
    void set(MatrixCPUD* m);
    void set(double value);
    void setTrans(MatrixCPUD* m);
    void setBloc(int iBegin, int iEnd, int jBegin, int jEnd, MatrixCPUD* m);
    void setEyes(double v);
    void setEyes(MatrixCPUD* vect);
    void setRand(int eps, int divide); // rand()%eps/divide
    void setRand1(double eps);
    void setFromFile(std::string file,int entete=0);
    void add(double c);
    void add(MatrixCPUD* m1, MatrixCPUD* m2); // m = m1 +m2;
    void add(MatrixCPUD* m1, double c); // m = m1 + c;
    void add(MatrixCPUD* m1);  // m = m + m1;
    void increment(int i, int j, double add); //m[i,j] = m[i,j] + add
    void addVector(MatrixCPUD* v);
    void addTrans(MatrixCPUD* m1); // m = m + tm1;
    
    void subtract(MatrixCPUD* m);
    void subtract(MatrixCPUD* m1, MatrixCPUD* m2);
    void subtractRow(int row1, int row2, double factor); // row1 - factor*row2
    void subtractAbs(MatrixCPUD* m1, MatrixCPUD* m2);
    void subtractVector(MatrixCPUD* v);
    void subtractTrans(MatrixCPUD* m);
    void multiply(MatrixCPUD* m1, MatrixCPUD* m2);
    void multiplyTrans(MatrixCPUD* m1, MatrixCPUD* m2, int numToTRans=1);
    void multiplyDiag(MatrixCPUD* m1, MatrixCPUD* m2); // multiply but only keep the diagonal
    void multiply(double c);
    void multiplyT(MatrixCPUD* m);
    void multiplyT(MatrixCPUD* m1, MatrixCPUD* m2);
    void multiplyTVector(MatrixCPUD* m, MatrixCPUD* v, int sens=0);
    void divide(double c);
    void divideT(MatrixCPUD* m);
    
    void invertGaussJordan(MatrixCPUD* mToInvert);
    void LDLFactorization(MatrixCPUD* L, MatrixCPUD* D);
    void LUPFactorization(MatrixCPUD* A, MatrixCPUD* P);

    void solveSysUpper(MatrixCPUD* U, MatrixCPUD* y); // Ux = y
    void solveSysLower(MatrixCPUD* L, MatrixCPUD* b , MatrixCPUD* P); // Ly = Pb
    void solveSys(MatrixCPUD* A, MatrixCPUD* P, MatrixCPUD* b); // LUx = Pb , A = (L-I) + U
    void solveSysGaussSeidel(MatrixCPUD* M, MatrixCPUD* b); // Mx = b


    void MultiplyMatVec(MatrixCPUD* m, MatrixCPUD* vect, int sens = 0);
    void MultiplyMatMat(MatrixCPUD* m1, MatrixCPUD* m2);


    
    void project(MatrixCPUD* Lb, MatrixCPUD* Ub);
    void projectNeg(); // min(m,0)
    void projectPos(); // max(m,0)


    double distance2(MatrixCPUD* m1) const; // ||m-m1||_2
    double distance2() const; // ||m||_2 
    double max2() const; // ||m||_inf 
    double maxAbs(int iBegin, int iEnd, int jBegin, int jEnd, MatrixCPUD* indices = nullptr);
    double minAbs(int iBegin, int iEnd, int jBegin, int jEnd, MatrixCPUD* indices = nullptr, bool Null=true);
    void Moy(MatrixCPUD* m, MatrixCPUD* nb, int sens=0); 
    double sum() const;
    double sumD() const;
    void sum(MatrixCPUD* m, int sens = 0);
    void sumT(MatrixCPUD* m, int sens = 0); // if matrix and vector not in the same sens
    void sort(int dim = 0, int sens = 0);
    void sortLine(int line1, int line2, int dim);
    void sortColumn(int col1, int col2, int dim);
    void fusionLine(int col1, int col2, int dim);
    void fusionColumn(int col1, int col2, int dim);
    void swapLine(int line1, int line2);
    void swapColumn(int col1, int col2);

    void RelativeEror(MatrixCPUD* MatRef, MatrixCPUD* Mat); 


    void display() const;
    void displayBloc(int iBegin, int iEnd, int jBegin, int jEnd) const;
    void saveCSV(const std::string& filename, std::ios_base::openmode mode= std::ios_base::in | std::ios_base::out, int trans=0) const;
    ~MatrixCPUD();


    // def Eigen
#ifdef EIGEN
    void set(Eigen::MatrixXd* eigenMatrix);
    void toEigenMatrix(Eigen::MatrixXd* eigenMatrix);
    void solveSysEigen(MatrixCPU* M, MatrixCPU* b); // Mx = b
    void invertEigen(MatrixCPU* mToInvert);
#endif
// def Eigen
#ifdef OSQP
    c_float* toCFloat();
    c_int toCSC(c_float* data, c_int* idx, c_int* ptr);
    c_int toCSCHalf(c_float* data, c_int* idx, c_int* ptr);
#endif



};

