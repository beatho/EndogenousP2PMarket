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
ParamInterface::ParamInterface(int N, int B, int L, int Lconst)
{
    _sizes = MatrixCPU(1,finSizes_ind);
    _sizes.set(0, nAgentP_ind, N);
    _sizes.set(0, nBusP_ind, B);
    _sizes.set(0, nLineP_ind, L);
    _sizes.set(0, nLineCons_ind, Lconst);

    _param  = MatrixCPU(1, finParam_ind);
    _param.set(0, SbaseP_ind, 1);
    _param.set(0, VbaseP_ind, 1);

    _lambda = MatrixCPU(2*(N + 1), N + 1);
    _trade  = MatrixCPU(2*(N + 1), N + 1);
    _Pn     = MatrixCPU(2*(N + 1), 1);
    _MU     = MatrixCPU(2*(N + 1), 1);
    _delta  = MatrixCPU(Lconst, 2);


   

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
    for(int i=0; i<_sizes.get(0, nLineCons_ind); i++){
        _delta.set(i, 0, delta1.get(i,0));
        _delta.set(i, 1, delta2.get(i,0));
    }
}


MatrixCPU ParamInterface::getSize(){return _sizes;}
MatrixCPU ParamInterface::getParam(){return _param;}
MatrixCPU ParamInterface::getTrade(){return _trade;}
MatrixCPU ParamInterface::getLambda(){return _lambda;}
MatrixCPU ParamInterface::getDelta(){return _delta;}
MatrixCPU ParamInterface::getPn(){return _Pn;}
MatrixCPU ParamInterface::getMU(){return _MU;}


void ParamInterface::display(int type){

    std::cout << "Simulation's parameters : " << std::endl;
    std::cout << "means default value if = 0" << std::endl;
    if(type == 1){
        
        std::cout << "Agents' count : " << _sizes.get(0, nAgentP_ind) << std::endl;
        std::cout << "rho / rho1 / rhol : " << _param.get(0, rho_ind) << " / " << _param.get(0, rho1_ind) << " / " 
        << _param.get(0, rhoL_ind) << std::endl;
        std::cout << "k_max / j_max / j_intern : " << _param.get(0, iterMax_ind) << " / " << _param.get(0,iterL_ind) 
        << " / " << _param.get(0, iterItern_ind) << std::endl;
        std::cout << "eps_g / eps_x / eps_l / eps_intern: " << _param.get(0, epsG_ind) << " / " << _param.get(0, epsX_ind) 
        << " / " << _param.get(0,epsL_ind) << " / " << _param.get(0, epsIntern_ind) << std::endl;
        
        std::cout << "StepG / StepL / StepIntern : " << _param.get(0, stepG_ind) << " / " << _param.get(0, stepL_ind) 
        << " / " << _param.get(0, stepIntern_ind) << std::endl;
    } else if (type==2){
        //enum indParam {iterG_ind, iterL_ind, iterItern_ind, stepG_ind, stepL_ind, stepIntern_ind, epsG_ind,
            //epsL_ind, epsX_ind, epsIntern_ind, rho_ind, rhoL_ind, rho1_ind, SbaseP_ind, VbaseP_ind, finParam_ind};
        //enum indSizes {nAgentP_ind, nBusP_ind, nLineP_ind, nLineCons_ind, finSizes_ind};
        std::cout << "sizes : nAgent, nBus, nLines, nLinesCosntraints " <<std::endl;
        _sizes.display();
        std::cout << "param : iterG, iterL, iterInt, stepG, stepL, stepInt, epsG, espL, epsX, epsInt, rho, rhoL, rhoX, Sbase, Vbase";
        _param.display();
    } else{
        std::cout << "Agents' count : " << _sizes.get(0, nAgentP_ind) << std::endl;
        std::cout << "rho / rho1 / rhol : " << _param.get(0, rho_ind) << " / " << _param.get(0, rho1_ind) << " / " 
        << _param.get(0, rhoL_ind) << std::endl;
        std::cout << "k_max / j_max / j_intern : " << _param.get(0, iterMax_ind) << " / " << _param.get(0,iterL_ind) 
        << " / " << _param.get(0, iterItern_ind) << std::endl;
        std::cout << "eps_g / eps_x / eps_l / eps_intern: " << _param.get(0, epsG_ind) << " / " << _param.get(0, epsX_ind) 
        << " / " << _param.get(0,epsL_ind) << " / " << _param.get(0, epsIntern_ind) << std::endl;
        std::cout << "StepG / StepL / StepIntern : " << _param.get(0, stepG_ind) << " / " << _param.get(0, stepL_ind) 
        << " / " << _param.get(0, stepIntern_ind) << std::endl;

        std::cout << "initialisation, trade, lambda, Pn, Mu, delta" << std::endl;
        _trade.display();
        _lambda.display();
        _Pn.display();
        _MU.display();
        _delta.display();

    }
   


}

/*
MatrixCPU _sizes;
    MatrixCPU _results;
    MatrixCPU _lambda;
    MatrixCPU _trade;
    MatrixCPU _Pn;
    MatrixCPU _MU;
    MatrixCPU _delta;

*/


ResultInterface::ResultInterface(int N, int B, int L, int Lconst)
{
    _sizes = MatrixCPU(1, finSizes_ind);
    _sizes.set(0, nAgentP_ind, N);
    _sizes.set(0, nBusP_ind, B);
    _sizes.set(0, nLineP_ind, L);
    _sizes.set(0, nLineCons_ind, Lconst);

    _results = MatrixCPU(1, finRes_ind);

    _lambda = MatrixCPU(2*(N+1), (N+1));
    _trade  = MatrixCPU(2*(N+1), (N+1));
    _Pn     = MatrixCPU(2*(N+1), 1);
    _MU     = MatrixCPU(2*(N+1), 1);
    _delta  = MatrixCPU(Lconst, 2);
}

ResultInterface::~ResultInterface(){}

void ResultInterface::setProbleme(MatrixCPU trade, MatrixCPU Pn){
    int M = Pn.getNLin();
    int N = trade.getNCol(); //_sizes.get(0, nAgentP_ind);
    int offset = 0;
    if(M == N){ // cas DC
        offset = 1;
    }
    /*else if( M == 2*N){ // cas AC agent des pertes compris
        // nothing
    } else if(M == 2*(N + 1)){ // cas AC avec ajout agent des pertes
        N = N + 1;
    } else{
        std::cout << "[WARNING] : Pn has not an expected size, unnexpected error can occur" <<std::endl;
    }*/
    for(int i = 0; i <M ; i++){
        _Pn.set(i + offset, 0, Pn.get(i,0));
        for(int j=0; j<N; j++){
            _trade.set(i + offset, j + offset, trade.get(i,j));
        }
    }
}

void ResultInterface::setDual(MatrixCPU lambda, MatrixCPU MU){
    int M = MU.getNLin(); // N ou 2*(N + 1)
    int N = lambda.getNCol(); // N ou N + 1
    int offset = 0;
    if(M == N){ // cas DC
        offset = 1;
    }
    for(int i = 0; i <M ; i++){
        _MU.set(i + offset, 0, MU.get(i,0));
        for(int j=0; j<N; j++){
            _lambda.set(i + offset, j + offset, lambda.get(i,j));
        }
    }
    //std::cout<< "affichage dans result offset " << offset<< " M =  " << M << " N " << N <<std::endl;
    //_lambda.display();
}
void ResultInterface::setDelta(MatrixCPU delta1, MatrixCPU delta2){
    for(int i=0; i<_sizes.get(0,nLineCons_ind); i++){
        _delta.set(i, 0, delta1.get(i,0));
        _delta.set(i, 1, delta2.get(i,0));
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
void ResultInterface::setResult(int iterF, int stepG, float temps, float fc, MatrixCPU resF){
    
    _resF = resF;

    _results.set(0, iterF_ind, iterF);
    _results.set(0, temps_ind, temps);
    _results.set(0, stepG2_ind, stepG);
    
    if(_resF.getNCol()>0){
        _results.set(0, fc_ind, fc);
        int indiceFinal =  (iterF - 1) / stepG;
        _results.set(0, resR_ind, _resF.get(0, indiceFinal));
        _results.set(0, resS_ind, _resF.get(1, indiceFinal));
        _results.set(0, resX_ind, _resF.get(2, indiceFinal));
    } else {
        _results.set(0, resR_ind, fc);
        _results.set(0, fc_ind, 0);
    }
    
}

void ResultInterface::setvarPhysic(MatrixCPU Pb, MatrixCPU Phi, MatrixCPU E){
    _Pb = Pb;
    _Phi = _Phi;
    _E = E;
}

MatrixCPU ResultInterface::getTrade(){return _trade;}
MatrixCPU ResultInterface::getLambda(){return _lambda;}
MatrixCPU ResultInterface::getDelta(){return _delta;}
MatrixCPU ResultInterface::getPn(){return _Pn;}
MatrixCPU ResultInterface::getMU(){return _MU;}
MatrixCPU ResultInterface::getResF(){return _resF;}
MatrixCPU ResultInterface::getResults(){return _results;}
MatrixCPU ResultInterface::getPb(){return _Pb;}
MatrixCPU ResultInterface::getPhi(){return _Phi;}
MatrixCPU ResultInterface::getE(){return _E;}

void ResultInterface::display(StudyCaseInterface* _case, int type){
   
   MatrixCPU Agents = _case->getAgentCase();
   int nAgent = _case->getN();
   MatrixCPU a(nAgent + 1, 1);
   MatrixCPU b(nAgent + 1, 1);
   MatrixCPU Pmax(2*(nAgent + 1), 1);
   MatrixCPU Pmin(2*(nAgent + 1), 1);
   for(int n=0; n<nAgent; n++){
       a.set(n + 1, 0, Agents.get(n, a_ind));
       b.set(n + 1, 0, Agents.get(n, b_ind));
       Pmax.set(n + 1, 0, Agents.get(n, Pmax_ind));
       Pmin.set(n + 1, 0, Agents.get(n, Pmin_ind));
       Pmax.set(nAgent + n + 2, 0, Agents.get(n, Qmax_ind));
       Pmin.set(nAgent + n + 2, 0, Agents.get(n, Qmin_ind));
   }

    if(type== 0){ // default
        std::cout << " Simulation result : " << std::endl;
        std::cout << "Study Case of agent, buses, lines and constraint line" << std::endl;
        _sizes.display();
        //enum indRes   {iterF_ind, temps_ind, fc_ind, iterMax_ind, stepG2_ind, resR_ind, resS_ind, resX_ind, finRes_ind};
        std::cout << "Results : iter, time, fc, iterMax, stepG, resR, resS, resX " <<std::endl;
        _results.display();

    } else if (type == 1){ // Market
        std::cout << "===============================================================|" << std::endl;
        std::cout << "        Market Simulation result :  System Summary             |" << std::endl;
        std::cout << "===============================================================|" << std::endl;
        std::cout << "Agents' count : " << _sizes.get(0, nAgentP_ind) << std::endl;
        std::cout << "f_c : " << _results.get(0, fc_ind) << std::endl;
        std::cout << "iter : " << _results.get(0, iterF_ind) << std::endl;
        std::cout << "Residuals : symetrie " << _results.get(0, resR_ind) << " convergence " << _results.get(0, resS_ind) << " grid " << _results.get(0, resX_ind) << std::endl;
        std::cout << "Computation time : " << _results.get(0, temps_ind) << std::endl;
        std::cout << std::endl << std::endl;
        std::cout << std::endl << std::endl;
        std::cout << "========================================================================================================|" << std::endl;
        std::cout << "      Agent Data                                                                                        |" << std::endl;
        std::cout << "========================================================================================================|" << std::endl;
        std::cout << " Agent  |  Cost    |  Cost    |          Power Injection           |           Power Injection          |" << std::endl;
        std::cout << "  #     |   a (pu) |   b (pu) |  P (pu)  | Pmin (pu)  | Pmax (pu)  |  Q (pu)   | Qmin (pu)  | Qmax (pu) |" << std::endl;
        std::cout << "--------|----------|----------|----------|------------|------------|-----------|------------|-----------|" << std::endl;

        std::cout << std::setw(8) << 0 << "|" << std::setw(9) <<  0 << " |" << std::setw(10)
                << 0 << "|" << std::setw(10) << _Pn.get(0, 0) << "|" << std::setw(12)
                << "-inf" << "|" << std::setw(12) << Pmax.get(0, 0)
                << "|" << std::setw(11) << _Pn.get(nAgent + 1, 0) << "|" << std::setw(12) << "-inf"
                << "|" << std::setw(11) << "inf" << "|" << std::endl;
        for (int n = 1; n <  _sizes.get(0, nAgentP_ind) + 1; n++) {
            std::cout << std::setw(8) << n << "|" << std::setw(9) << a.get(n, 0) << " |" << std::setw(10)
                << b.get(n, 0) << "|" << std::setw(10) << _Pn.get(n, 0) << "|" << std::setw(12)
                << Pmin.get(n, 0) << "|" << std::setw(12) << Pmax.get(n, 0)
                << "|" << std::setw(11) << _Pn.get(n + nAgent, 0) << "|" << std::setw(12) << Pmin.get(n + nAgent, 0)
                << "|" << std::setw(11) << Pmax.get(n + nAgent, 0) << "|" << std::endl;
        }


	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "                      END PRINT                                                                         |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;

        if(_sizes.get(0, nAgentP_ind) < 20){
            std::cout << " Trade : " << std::endl;
            _trade.display();
            std::cout << " Lambda : " << std::endl;
            _lambda.display();
        }

    } else if (type == 2){ // PF
        int Nbus =  _sizes.get(0, nBusP_ind);
        std::cout << "The error of this state is " << _results.get(0, resR_ind) << std::endl;
        std::cout << "===============================================================|" << std::endl;
        std::cout << "      System Summary                                           |" << std::endl;
        std::cout << "===============================================================|" << std::endl;
        std::cout << "Buses            " << Nbus << std::endl;
        std::cout << "Branches         " << _sizes.get(0, nLineP_ind) << std::endl;

        std::cout << "===============================================================|" << std::endl;
        std::cout << "      Bus Data                                                 |" << std::endl;
        std::cout << "===============================================================|" << std::endl;
        std::cout << " Bus |          Voltage        |  Power = Generation  + Load   |" << std::endl;
        std::cout << "  #  |    Mag(pu) |  Ang(deg)  |    P (pu)     |     Q (pu)    |" << std::endl;
        std::cout << "-----|------------|------------|---------------|---------------|" << std::endl;

        //std::cout << 0 << "      " << E.get(Nbus, 0) << "             " << E.get(0, 0) * (abs(E.get(0, 0)) > 0.0001) * 180 / 3.1415 << "              " << (abs(W.get(0, 0)) > 0.0001) * W.get(0, 0) << "         " << (abs(W.get(Nbus, 0)) > 0.0001) * W.get(Nbus, 0) << std::endl;

        float seuil = 0.0001;
        
        std::cout << std::setw(5) << 0 << "|" << std::setw(11) << _E.get(Nbus, 0) << "*|" << std::setw(11) << _E.get(0, 0) * (abs(_E.get(0, 0)) > seuil) * 180 / 3.1415
            << "*|" << std::setw(15) << (abs(_Pb.get(0, 0)) > seuil) * _Pb.get(0, 0) << "|" << std::setw(15) << (abs(_Pb.get(Nbus, 0)) > seuil) * _Pb.get(Nbus, 0)
            << "|" << std::endl;
        for (int b = 1; b < Nbus; b++) {
            //std::cout.width(10);
            //std::cout << b << "      " << E.get(b + Nbus, 0) << "        " << E.get(b, 0) * (abs(E.get(b, 0)) > 0.0001) * 180 / 3.1415 << "          " << (abs(W.get(b, 0)) > 0.0001) * W.get(b, 0) << "         " << (abs(W.get(b + Nbus, 0)) > 0.0001) * W.get(b + Nbus, 0) << std::endl;
            std::cout << std::setw(5) << b << "|" << std::setw(11) << _E.get(b + Nbus, 0) << " |" << std::setw(11)
                << _E.get(b, 0) * (abs(_E.get(b, 0)) > seuil) * 180 / 3.1415 << " |" << std::setw(15)
                << (abs(_Pb.get(b, 0)) > seuil) * _Pb.get(b, 0) << "|" << std::setw(15) << (abs(_Pb.get(b + Nbus, 0)) > seuil) * _Pb.get(b + Nbus, 0)
                << "|" << std::endl;

        }
        std::cout << std::endl << std::endl;
        std::cout << " Phi " << std::endl;
        _Phi.display();


        std::cout << "===============================================================================================|" << std::endl;
        std::cout << "                      END PRINT                                                                |" << std::endl;
        std::cout << "===============================================================================================|" << std::endl;
    } else if (type == 3){ // EndoMarket
        int Nbus =  _sizes.get(0, nBusP_ind);
        std::cout << "===============================================================|" << std::endl;
        std::cout << "        Market Simulation result :  System Summary             |" << std::endl;
        std::cout << "===============================================================|" << std::endl;
        std::cout << "Agents' count :  " << _sizes.get(0, nAgentP_ind) << std::endl;
        std::cout << "Buses            " <<  Nbus << std::endl;
        std::cout << "Branches         " << _sizes.get(0, nLineP_ind) << std::endl;
        std::cout << "f_c : " << _results.get(0, fc_ind) << std::endl;
        std::cout << "iter : " << _results.get(0, iterF_ind) << std::endl;
        std::cout << "Residuals : symetrie " << _results.get(0, resR_ind) << " convergence " << _results.get(0, resS_ind) << " grid " << _results.get(0, resX_ind) << std::endl;
        std::cout << "Computation time : " << _results.get(0, temps_ind) << std::endl;
        std::cout << std::endl << std::endl;
        std::cout << std::endl << std::endl;
	
        std::cout << std::endl << std::endl;
        std::cout << std::endl << std::endl;
        std::cout << "========================================================================================================|" << std::endl;
        std::cout << "      Agent Data                                                                                        |" << std::endl;
        std::cout << "========================================================================================================|" << std::endl;
        std::cout << " Agent  |  Cost    |  Cost    |          Power Injection           |           Power Injection          |" << std::endl;
        std::cout << "  #     |   a (pu) |   b (pu) |  P (pu)  | Pmin (pu)  | Pmax (pu)  |  Q (pu)   | Qmin (pu)  | Qmax (pu) |" << std::endl;
        std::cout << "--------|----------|----------|----------|------------|------------|-----------|------------|-----------|" << std::endl;

        std::cout << std::setw(8) << 0 << "|" << std::setw(9) <<  0 << " |" << std::setw(10)
                << 0 << "|" << std::setw(10) << _Pn.get(0, 0) << "|" << std::setw(12)
                << Pmin.get(0, 0) << "|" << std::setw(12) << Pmax.get(0, 0)
                << "|" << std::setw(11) << _Pn.get(nAgent, 0) << "|" << std::setw(12) << "-inf"
                << "|" << std::setw(11) << 0 << "|" << std::endl;
        for (int n = 1; n <  _sizes.get(0, nAgentP_ind) + 1; n++) {
            
            std::cout << std::setw(8) << n << "|" << std::setw(9) << a.get(n, 0) << " |" << std::setw(10)
                << b.get(n, 0) << "|" << std::setw(10) << _Pn.get(n, 0) << "|" << std::setw(12)
                << Pmin.get(n, 0) << "|" << std::setw(12) << Pmax.get(n, 0)
                << "|" << std::setw(11) << _Pn.get(n + nAgent, 0) << "|" << std::setw(12) << Pmin.get(n + nAgent, 0)
                << "|" << std::setw(11) << Pmax.get(n + nAgent, 0) << "|" << std::endl;
        }
        std::cout << "===============================================================|" << std::endl;
        std::cout << "      Bus Data                                                 |" << std::endl;
        std::cout << "===============================================================|" << std::endl;
        std::cout << " Bus |          Voltage        |  Power = Generation  + Load   |" << std::endl;
        std::cout << "  #  |    Mag(pu) |  Ang(deg)  |    P (pu)     |     Q (pu)    |" << std::endl;
        std::cout << "-----|------------|------------|---------------|---------------|" << std::endl;

        //std::cout << 0 << "      " << E.get(Nbus, 0) << "             " << E.get(0, 0) * (abs(E.get(0, 0)) > 0.0001) * 180 / 3.1415 << "              " << (abs(W.get(0, 0)) > 0.0001) * W.get(0, 0) << "         " << (abs(W.get(Nbus, 0)) > 0.0001) * W.get(Nbus, 0) << std::endl;

        float seuil = 0.0001;
        
        std::cout << std::setw(5) << 0 << "|" << std::setw(11) << _E.get(Nbus, 0) << "*|" << std::setw(11) << _E.get(0, 0) * (abs(_E.get(0, 0)) > seuil) * 180 / 3.1415
            << "*|" << std::setw(15) << (abs(_Pb.get(0, 0)) > seuil) * _Pb.get(0, 0) << "|" << std::setw(15) << (abs(_Pb.get(Nbus, 0)) > seuil) * _Pb.get(Nbus, 0)
            << "|" << std::endl;
        for (int b = 1; b < Nbus; b++) {
            //std::cout.width(10);
            //std::cout << b << "      " << E.get(b + Nbus, 0) << "        " << E.get(b, 0) * (abs(E.get(b, 0)) > 0.0001) * 180 / 3.1415 << "          " << (abs(W.get(b, 0)) > 0.0001) * W.get(b, 0) << "         " << (abs(W.get(b + Nbus, 0)) > 0.0001) * W.get(b + Nbus, 0) << std::endl;
            std::cout << std::setw(5) << b << "|" << std::setw(11) << _E.get(b + Nbus, 0) << " |" << std::setw(11)
                << _E.get(b, 0) * (abs(_E.get(b, 0)) > seuil) * 180 / 3.1415 << " |" << std::setw(15)
                << (abs(_Pb.get(b, 0)) > seuil) * _Pb.get(b, 0) << "|" << std::setw(15) << (abs(_Pb.get(b + Nbus, 0)) > seuil) * _Pb.get(b + Nbus, 0)
                << "|" << std::endl;

        }
        std::cout << std::endl << std::endl;
        std::cout << " Phi " << std::endl;
        _Phi.display();

	std::cout << "========================================================================================================|" << std::endl;
	std::cout << "                      END PRINT                                                                         |" << std::endl;
	std::cout << "========================================================================================================|" << std::endl;

   
    }

}