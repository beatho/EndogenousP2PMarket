#pragma once
#include "MatrixCPU.h"

/*
Pour créer les paramètres :
  - taille 1 : - IterG, iterL, iterIntern, ,
               - stepG, stepL, stepIntern,
               - epsG, epsL, epsX, epsInter, 
               - rho, rhoL, rho1
               - Sbase, Vbase
               - N, B, L, 
               - iterFinal, temps, fc
  - taille N*N : lambda, trade,
  - taille N : Pn, MU
  - taille L : _delta1, _delta2
  - taille iter : _resF
*/

enum indParam {iterG_ind, iterL_ind, iterItern_ind, stepG_ind, stepL_ind, stepIntern_ind, epsG_ind,
     epsL_ind, epsX_ind, epsIntern_ind, rho_ind, rhoL_ind, rho1_ind, SbaseP_ind, VbaseP_ind, finParam_ind};
enum indSizes {nAgentP_ind, nBusP_ind, nLineP_ind, finSizes_ind};
enum indRes   {iterF_ind, temps_ind, fc_ind, iterMax_ind, stepG2_ind, resR_ind, resS_ind, resX_ind, finRes_ind};


class ParamInterface
{
private:
    MatrixCPU _sizes;
    MatrixCPU _param;
    MatrixCPU _trade;
    MatrixCPU _Pn;
    MatrixCPU _MU;
    MatrixCPU _lambda;
    MatrixCPU _delta;
public:
    ParamInterface(int N, int B, int L);
    void setIter(int iterG, int iterL, int iterIntern=1);
    void setStep(int stepG, int stepL, int stepInter=1);
    void setEps(float epsG, float epsL, float epsX=0, float epsInter=0);
    void setRho(float rho, float rhoL=0, float rho1=0);
    void setVbase(float Vbase);
    void setSbase(float Sbase);
    void initProbleme(MatrixCPU trade, MatrixCPU Pn);
    void initDual(MatrixCPU lambda, MatrixCPU MU);
    void initDelta(MatrixCPU delta1, MatrixCPU delta2);
    MatrixCPU getSize();
    MatrixCPU getParam();
    MatrixCPU getTrade();
    MatrixCPU getLambda();
    MatrixCPU getDelta();
    MatrixCPU getPn();
    MatrixCPU getMU();
    ~ParamInterface();
};
class ResultInterface
{
private:
    MatrixCPU _sizes;
    MatrixCPU _results;
    MatrixCPU _lambda;
    MatrixCPU _trade;
    MatrixCPU _Pn;
    MatrixCPU _MU;
    MatrixCPU _delta;
    MatrixCPU _resF;
public:
    ResultInterface(int N, int B, int L);
    ~ResultInterface();
    void setProbleme(MatrixCPU trade, MatrixCPU Pn);
    void setDual(MatrixCPU lambda, MatrixCPU MU);
    void setDelta(MatrixCPU delta1, MatrixCPU delta2);
    void changeIterStep(int iter, int step);
    //iterF_ind, temps_ind, fc_ind, iterMax_ind, stepG2_ind, resR_ind, resS_ind, resX_ind, finRes_ind};
    void setResult(int iterF, float temps, float fc, MatrixCPU resF);
    MatrixCPU getTrade();
    MatrixCPU getLambda();
    MatrixCPU getDelta();
    MatrixCPU getPn();
    MatrixCPU getMU();
    MatrixCPU getResF();
    MatrixCPU getResults();
};


