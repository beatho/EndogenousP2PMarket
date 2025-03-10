#include "../head/ParamInterface.h"

/*
Pour créer les paramètres :
  - taille 1 : - IterG, iterL, iterIntern, ,
               - stepG, stepL, stepIntern,
               - epsG, epsL, epsX, epsInter, 
               - rho, rhoL, rho1
               - N, B, L, 
               - iterFinal, temps, fc
  - taille N*N : lambda, trade,
  - taille N : Pn, MU
  - taille L : _delta1, _delta2
  - taille iter : _resF

*/
ParamInterface::ParamInterface(int N, int B, int L)
{
    _sizes = MatrixCPU(1,3);
    _sizes.set(0, nAgentP_ind, N);
    _sizes.set(0, nBusP_ind, B);
    _sizes.set(0, nLineP_ind, L);

    _param  = MatrixCPU(1, finParam_ind);
    _lambda = MatrixCPU(2*N, N);
    _trade  = MatrixCPU(2*N, N);
    _Pn     = MatrixCPU(2*N, 1);
    _MU     = MatrixCPU(2*N, 1);
    _delta  = MatrixCPU(L  , 2);


    _param.set(0, SbaseP_ind, 1);
    _param.set(0, VbaseP_ind, 1);

}
ParamInterface::~ParamInterface(){}

void ParamInterface::setIter(int iterG, int iterL, int iterIntern){
    _param.set(0, iterG_ind, iterG);
    _param.set(0, iterL_ind, iterL);
    _param.set(0, iterItern_ind, iterIntern);
}
void ParamInterface::setStep(int stepG, int stepL, int stepInter){
    _param.set(0, stepG_ind, stepG);
    _param.set(0, stepL_ind, stepL);
    _param.set(0, stepIntern_ind, stepInter);
}
void ParamInterface::setEps(float epsG, float epsL, float epsX, float epsInter){
    _param.set(0, epsG_ind, epsG);
    _param.set(0, epsL_ind, epsL);
    _param.set(0, epsIntern_ind, epsInter);
    _param.set(0, epsX_ind, epsX);
}
void ParamInterface::setRho(float rho, float rhoL, float rho1){
    _param.set(0, rho_ind, rho);
    _param.set(0, rhoL_ind, rhoL);
    _param.set(0, rho1_ind, rho1);
}
void ParamInterface::setSbase(float Sbase){
    _param.set(0, SbaseP_ind, Sbase);
}
void ParamInterface::setVbase(float Vbase){
    _param.set(0, VbaseP_ind, Vbase);
}
void ParamInterface::initProbleme(MatrixCPU trade, MatrixCPU Pn){
    _trade = trade;
    _Pn = Pn;
}
void ParamInterface::initDual(MatrixCPU lambda, MatrixCPU MU){
    _lambda = lambda;
    _MU = MU;
}
void ParamInterface::initDelta(MatrixCPU delta1, MatrixCPU delta2){
    for(int i=0; i<_sizes.get(0,2); i++){
        _delta.set(i, 0, delta1.get(i,0));
        _delta.set(i, 0, delta2.get(i,0));
    }
}


MatrixCPU ParamInterface::getSize(){return _sizes;}
MatrixCPU ParamInterface::getParam(){return _param;}
MatrixCPU ParamInterface::getTrade(){return _trade;}
MatrixCPU ParamInterface::getLambda(){return _lambda;}
MatrixCPU ParamInterface::getDelta(){return _delta;}
MatrixCPU ParamInterface::getPn(){return _Pn;}
MatrixCPU ParamInterface::getMU(){return _MU;}

/*
MatrixCPU _sizes;
    MatrixCPU _results;
    MatrixCPU _lambda;
    MatrixCPU _trade;
    MatrixCPU _Pn;
    MatrixCPU _MU;
    MatrixCPU _delta;

*/


ResultInterface::ResultInterface(int N, int B, int L)
{
    _sizes = MatrixCPU(1,3);
    _sizes.set(0,0,N);
    _sizes.set(0,1,B);
    _sizes.set(0,2,L);

    _results = MatrixCPU(1, finRes_ind);

    _lambda = MatrixCPU(2*N, N);
    _trade  = MatrixCPU(2*N, N);
    _Pn     = MatrixCPU(2*N, 1);
    _MU     = MatrixCPU(2*N, 1);
    _delta  = MatrixCPU(L  , 2);
}

ResultInterface::~ResultInterface(){}

void ResultInterface::setProbleme(MatrixCPU trade, MatrixCPU Pn){
    _trade = trade;
    _Pn = Pn;
}
void ResultInterface::setDual(MatrixCPU lambda, MatrixCPU MU){
    _lambda = lambda;
    _MU = MU;
}
void ResultInterface::setDelta(MatrixCPU delta1, MatrixCPU delta2){
    for(int i=0; i<_sizes.get(0,2); i++){
        _delta.set(i, 0, delta1.get(i,0));
        _delta.set(i, 0, delta2.get(i,0));
    }
}
void ResultInterface::changeIterStep(int iter, int step){
    if(iter){
        _sizes.set(0, iterMax_ind, iter);
    }
    if(step){
        _sizes.set(0, stepG2_ind, step);
    }
    int _iter = _sizes.get(0, iterMax_ind);
    int _step = _sizes.get(0, stepG2_ind);
    if(_step){
        _resF = MatrixCPU(3, _iter/_step + 1);
    }
    
}
void ResultInterface::setResult(int iterF, float temps, float fc, MatrixCPU resF){
    _results.set(0, iterF_ind, iterF);
    _results.set(0, temps_ind, temps);
    _results.set(0, fc_ind, fc);

    _resF = resF;
    int indiceFinal =  (_sizes.get(0, iterMax_ind) - 1) / _sizes.get(0, stepG2_ind);
    _results.set(0, resR_ind, _resF.get(0, indiceFinal));
    _results.set(0, resS_ind, _resF.get(1, indiceFinal));
    _results.set(0, resX_ind, _resF.get(2, indiceFinal));
}


MatrixCPU ResultInterface::getTrade(){return _trade;}
MatrixCPU ResultInterface::getLambda(){return _lambda;}
MatrixCPU ResultInterface::getDelta(){return _delta;}
MatrixCPU ResultInterface::getPn(){return _Pn;}
MatrixCPU ResultInterface::getMU(){return _MU;}
MatrixCPU ResultInterface::getResF(){return _resF;}
MatrixCPU ResultInterface::getResults(){return _results;}