#include "../head/GPUPFdistPQ.cuh"


GPUPFdistPQ::GPUPFdistPQ() {}
GPUPFdistPQ::~GPUPFdistPQ() {}

void GPUPFdistPQ::init(const StudyCase& cas, MatrixGPU* PQ)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
    timePerBlock = MatrixCPU(1, 9); // Fb0 : init, Fb1ab : Flu, Fb2abc: Tension , FB3 : puissance, Fb4 erreur, Fb0 mise à jour

    occurencePerBlock = MatrixCPU(1, 9);; //nb de fois utilisé pendant la simu

#endif // INSTRUMENTATION
    

    //std::cout << "init PF GPUPFdistPQ" <<std::endl;
    //PQ->display(true);
    Nagent = cas.getNagent();
    Nbus = cas.getNBus();
    B2 = 2 * Nbus;
    N2 = 2 * Nagent;
   
    Nline = cas.getNLine(true); // ne doit pas être réduit ici !!!
    BL2 = Nbus + 2 * Nline;
    Nconstraint = B2 + Nline;
    //std::cout << Nline << " " << Nbus << std::endl;
    iterM = 30;
    iter = 0;
    epsPF = 0.00005;
    status = 0;
    numBlock = ceil((Nbus + _blockSize - 1) / _blockSize);
    //std::cout << numBlock << " " << _blockSize << std::endl;
    if (Nbus > 1024) {
        throw std::invalid_argument("too much bus, must change the computation of the voltage and S");
    }
    //CHECK_LAST_CUDA_ERROR();
    V0 = cas.getV0();
    theta0 = cas.gettheta0();

    //std::cout << "V0 :" << V0 << " theta0 " << theta0 << std::endl;
    v0 = V0 * cos(theta0);
    w0 = V0 * sin(theta0);
    _name = "Power summation method GPU"; // meilleure covergence quand c'est beaucoup chargé

    W0 = MatrixGPU(B2, 1, 0, 1);
    CoresAgentBus = MatrixGPU(cas.getCoresAgentBusLin(), 1);
    CoresAgentBusBegin = MatrixGPU(cas.getCoresAgentBusLinBegin(), 1);
    NagentByBus = MatrixGPU(cas.getNagentByBus(), 1);

    removeLossAgent << <1, 1 >> > (NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU);

    //CHECK_LAST_CUDA_ERROR();

    calculW0Bis(PQ);
    //CHECK_LAST_CUDA_ERROR();
    /*std::cout << " W0 : " << std::endl;
    W0.display(true);*/
    Y = MatrixGPU(2 * Nbus + Nline, 1, 0, 1);
    W = MatrixGPU(B2, 1, 0, 1);
   
    W.preallocateReduction(); // calcul des pertes
   
    dW = MatrixGPU(B2, 1, 0, 1);
    dW.preallocateReduction();
   
    _Pintermediate = MatrixGPU(BL2, 1, 0, 1);
    
    _Qintermediate = MatrixGPU(BL2, 1, 0, 1);
    E = MatrixGPU(B2, 1, 0, 1);
   
    VoltageRealIm = MatrixGPU(B2, 1, 0, 1);
    
    VoltageRealImPre = MatrixGPU(B2, 1, 0, 1);
    
    VoltageRealImPre.preallocateReduction();
    

    initE << <numBlock, _blockSize >> > (E._matrixGPU, theta0, V0, Nbus);
    initECar << <numBlock, _blockSize >> > (VoltageRealIm._matrixGPU, v0, w0, Nbus);
    
    //VoltageRealIm.display(true);

    // W0[2 * N] : puissance active et réactive au noeud (I*[P Q])
    // W[2 * N] : puissance obtenue par calcul à partir de E
    // dW[2 * N] : derive de puissance
    // E[2 * N] : angle puis tension [O et 1] pour l'init ?
    // dE[2 * N] : derive de angle puis tension
   
  
    CoresLineBus = MatrixGPU(cas.getCoresLineBus(true));
    _CoresVoiLin = MatrixGPU(cas.getCoresVoiLin(), 1);
    _CoresBusLin = MatrixGPU(cas.getCoresBusLin(), 1);
    _nLines = MatrixGPU(cas.getNLines(), 1);
    
    //CHECK_LAST_CUDA_ERROR();


    //Bgrid = cas.getLineSuceptance();
    //Ggrid = cas.getLineReactance();
    _Blin = MatrixGPU(cas.getBlin(), 1);
    _Glin = MatrixGPU(cas.getGlin(), 1);
    _Blin2 = MatrixGPU(cas.getBlin2(), 1);
    _Glin2 = MatrixGPU(cas.getGlin2(), 1);


    // specificite algo
    // CoresLineBus.display();
    //std::cout << Nbus << " " << Nline << std::endl;
    ZsRe = MatrixGPU(cas.getZsRe(), 1);
    ZsIm = MatrixGPU(cas.getZsImag(), 1);
    Yd = MatrixGPU(cas.getYd(), 1);
    chekcase();

    F = MatrixGPU(Nbus, 1, -1); // F_i = bus antécédent de i
    nChild = MatrixGPU(Nbus, 1, 0);
    CoresLineBusGPU = MatrixGPU(2, Nline);
    F.set(1, 0, 0);
    if (Nbus != (Nline + 1)) {
        std::cout << "Warning this is not a distribution network, F not set" << std::endl;
    
    }
    else {
        for (int lold = 0; lold < Nline; lold++) {
            int busTo = CoresLineBus.get(lold, 1);
            int busFrom = CoresLineBus.get(lold, 0);
            F.set(busTo, 0, busFrom);
            nChild.set(busFrom, 0, nChild.get(busFrom, 0) + 1);
            CoresLineBusGPU.set(0, lold, busFrom);
            CoresLineBusGPU.set(1, lold, busTo);
        }
    } 
    LastBus = cas.getLastBus(); 
    //CoresLineBus.display();
    //std::cout << "LastBus " << LastBus << std::endl;
    //CHECK_LAST_CUDA_ERROR();
    int debutChild = 0;
    MatrixCPU nChildTemp(Nbus, 1, 0);
    _indiceChildBegin = MatrixGPU(Nline, 1);
    Childs = MatrixGPU(Nbus, 1);
    for (int i = 0; i < Nbus; i++) {
        if (i > 0) {
            _indiceChildBegin.set(i - 1, 0, debutChild);

            int Ai = F.get(i, 0);
            Childs.set(_indiceChildBegin.get(Ai, 0) + nChildTemp.get(Ai, 0), 0, i);
            nChildTemp.increment(Ai, 0, 1);
            debutChild += nChild.get(i - 1, 0);
        }

    }
  
    //CoresLineBus.transferGPU();
    CoresLineBusGPU.transferGPU();
    F.transferGPU();
    nChild.transferGPU();
    Childs.transferGPU();
    _indiceChildBegin.transferGPU();
    
   
    //CHECK_LAST_CUDA_ERROR();
   
    St = MatrixGPU(2 * Nline, 1, -1, 1);
    Sf = MatrixGPU(2 * Nline, 1, -1, 1);
  
    /*St.display(true);
    Sf.display(true);
    
    ZsRe.display();
    ZsIm.display();
    Yd.display();*/
    
    //CHECK_LAST_CUDA_ERROR();

   
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 0, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION
   //std::cout << " fin init" << std::endl;

}

bool GPUPFdistPQ::chekcase()
{
    bool transfertToDo = false;
    if (CoresLineBus.getPos()) {
        transfertToDo = true;
        CoresLineBus.transferCPU();
    }
    if (Nbus != (Nline + 1)) {
        std::cout << "wrong number of line "<< Nline << "against "<< Nbus << std::endl;
        return false;
    }
    for (int i = 0; i < Nline; i++) {
        if (CoresLineBus.get(i, 1) != (i + 1)) {
            std::cout << "wrong numerotation of line " << CoresLineBus.get(i, 1) << "against " << (i + 1) << std::endl;
            return false;
        }
        if (CoresLineBus.get(i, 0) > CoresLineBus.get(i, 1)) {
            std::cout << "wrong numeoration of bus " << CoresLineBus.get(i, 0) << "against " << CoresLineBus.get(i, 1) << std::endl;
            return false;
        }
    }
    if (ZsRe.getNLin() == 0 || Yd.getNLin() == 0 || ZsIm.getNLin() == 0) {
        std::cout << "matrice non defined, ZsRe, Zs Im, Yd" << std::endl;
        ZsRe.display();
        ZsIm.display();
        Yd.display();
        return false;
    }

    if (transfertToDo) {
        CoresLineBus.transferGPU();
    }
    //std::cout << "checkcase OK " << std::endl;

    return true;
}

void GPUPFdistPQ::updatePQ(MatrixGPU* PQ)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    calculW0Bis(PQ);
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 8, 1);
#endif
}

void GPUPFdistPQ::solve()
{
    //std::cout << "solve distPQ" << std::endl;
    time = clock();
    err = 2 * epsPF;
    iter = 0;
    //std::cout << epsPF << " " << iterM << std::endl;
    status = 1;
    while (err > epsPF && iter < iterM) {
       
        calcS();
        
        //Jb.display();
        //VoltageRealImPre.display();

#ifdef INSTRUMENTATION
        t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION      
        calcVoltage();
        
#ifdef INSTRUMENTATION
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 3, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 3, 1);
        t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
        

        err = VoltageRealIm.distance2(&VoltageRealImPre);
        
#ifdef INSTRUMENTATION
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 7, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 7, 1);
#endif // INSTRUMENTATION
        VoltageRealImPre.set(&VoltageRealIm);
       
        /*if (err < epsPF) {
            calcW();
            dW.subtract(&W0, &W); // dW = W0 - W
            err = dW.max2(); //err = ||dW||
        }*/

        iter++;
        //std::cout << err << " ";
       
    }
    //
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    calcW(true);
   // W.display(true);
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 6, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 6, 1);
#endif // INSTRUMENTATION
    /*std::cout << "tension bus entree puis sortie" << std::endl;
    W0.display(true);
    W.display(true);*/

    if (iter >= iterM) {
        status = 2;
        if (err > 100 * epsPF) {
            status = -1;
        }
        //std::cout << "fin solve " << iter << " " << err << std::endl;
    }

    time = clock() - time;


}

void GPUPFdistPQ::calcS()
{
    // step 2
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    //std::cout << "step 2" << std::endl;
    // Set receiving end branch ?ow equal to the sum of the demand at receiving end (s^k_d) and the power drawn in the admittance(y^k_d) connected to bus k
   
    calculStGPU << <numBlock, _blockSize >> > (St._matrixGPU, VoltageRealIm._matrixGPU, W0._matrixGPU, Yd._matrixGPU, Nbus);
    
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 1, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 1, 1);
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    //std::cout << "step 3" << std::endl;
    // step 3
    //Backward sweep: Perform current summation, starting from the branch with the biggest index and heading towards the branch
    //whose index is equal to 1. The current of branch k is added to the current of the branch whose index is equal to i = f(k)
    
    
    // : 117
    
    calculSGPU << <1, Nline, (Nline * (sizeof(bool) + sizeof(int)) + 2 * Nline * sizeof(float)) >> > (St._matrixGPU, Sf._matrixGPU, VoltageRealIm._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, nChild._matrixGPU, Childs._matrixGPU, _indiceChildBegin._matrixGPU, Nbus);
   
    
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 2, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 2, 1);
#endif // INSTRUMENTATION
    //Jb.display();

}

void GPUPFdistPQ::calcW(bool end)
{
    
    calcWinterCar << <Nbus, _blockSize, B2 * sizeof(float) >> > (_Pintermediate._matrixGPU, _Qintermediate._matrixGPU, VoltageRealIm._matrixGPU, _Glin._matrixGPU, _Blin._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);

    if (!end) { // pendant simu, la puissance à ce noeud est libre
        switch (_blockSize) {
        case 512:
            calcWGPU<512> << <Nbus, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 256:
            calcWGPU<256> << <Nbus, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 128:
            calcWGPU<128> << <Nbus, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 64:
            calcWGPU< 64> << <Nbus, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 32:
            calcWGPU< 32> << <Nbus, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 16:
            calcWGPU< 16> << <Nbus, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case  8:
            calcWGPU<  8> << <Nbus, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case  4:
            calcWGPU<  4> << <Nbus, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case  2:
            calcWGPU<  2> << <Nbus, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case  1:
            calcWGPU<  1> << <Nbus, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        }
    }
    else {
        switch (_blockSize) {
        case 512:
            calcWGPU<512> << <Nbus, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 256:
            calcWGPU<256> << <Nbus, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 128:
            calcWGPU<128> << <Nbus, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 64:
            calcWGPU< 64> << <Nbus, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 32:
            calcWGPU< 32> << <Nbus, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 16:
            calcWGPU< 16> << <Nbus, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case  8:
            calcWGPU<  8> << <Nbus, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case  4:
            calcWGPU<  4> << <Nbus, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case  2:
            calcWGPU<  2> << <Nbus, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case  1:
            calcWGPU<  1> << <Nbus, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        }

    }

}



int GPUPFdistPQ::calcVoltage()
{
   // std::cout << "step 4" << std::endl;
    //Forward sweep: The receiving end bus voltages are calculated with known branch currents and sending bus voltages. 
    //Sf.display(true);
   updateVoltage << <1, Nline, Nbus* (8*sizeof(bool) + 2 * sizeof(float)) >> > (VoltageRealIm._matrixGPU, ZsRe._matrixGPU, ZsIm._matrixGPU, Sf._matrixGPU, F._matrixGPU, Nbus, LastBus);
   
    return 0;
}


void GPUPFdistPQ::calcE()
{
    calcEGPU << <numBlock, _blockSize >> > (E._matrixGPU, VoltageRealIm._matrixGPU, Nbus);

}

MatrixGPU GPUPFdistPQ::getY()
{
    
    calcE();
    calculYGPU <<<numBlock, _blockSize >> > (Y._matrixGPU, E._matrixGPU, VoltageRealIm._matrixGPU, _Blin2._matrixGPU, _Glin2._matrixGPU, CoresLineBusGPU._matrixGPU, Nbus, Nline);
    
    return Y;
}

void GPUPFdistPQ::setE(MatrixGPU* Enew)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    E = *Enew;
    if (!E.getPos()) {
        E.transferGPU();
    }
    initECar << <numBlock, _blockSize >> > (VoltageRealIm._matrixGPU, E._matrixGPU, Nbus);
    //CHECK_LAST_CUDA_ERROR();
    VoltageRealImPre.set(&VoltageRealIm);
    
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION
}

void GPUPFdistPQ::display2(bool all)
{
    std::cout.precision(3);
    float errV = err;
    W.transferCPU();
    W0.transferCPU();
    E.transferCPU();
    dW.transferCPU();

    if (iter == 0) {
        std::cout << "algorithm not launch" << std::endl;
        calcW(true);
       
        float temp = W.get(0, 0);
        float temp2 = W.get(Nbus, 0);
        W.set(0, 0, W0.get(0, 0));
        W.set(Nbus, 0, W0.get(Nbus, 0));
        dW.subtract(&W0, &W); // dW = W0 - W
        err = dW.max2(); //err = ||dW||
        W.set(0, 0, temp);
        W.set(Nbus, 0, temp2);
        
    }
    else if (iter < iterM) {
        std::cout << "method " << _name << " converged in " << iter << " iterations." << std::endl;
        std::cout << "Converged in " << (float)time / CLOCKS_PER_SEC << " seconds" << std::endl;
        
        std::cout << " Computation with float simple precision" << std::endl;
        float temp = W.get(0, 0);
       
        float temp2 = W.get(Nbus, 0);
        W.set(0, 0, W0.get(0, 0));
        W.set(Nbus, 0, W0.get(Nbus, 0));
        dW.subtract(&W0, &W); // dW = W0 - W
        err = dW.max2(); //err = ||dW||
        W.set(0, 0, temp);
        W.set(Nbus, 0, temp2);
        

    }
    else {
        std::cout << "method " << _name << " not converged in " << iter << " iterations." << std::endl;
        std::cout << "time taken " << (float)time / CLOCKS_PER_SEC << " seconds" << std::endl;
    }
    std::cout << "The power error of this state is " << err << std::endl;
    std::cout << "The tension error of this state is " << errV << std::endl;
    std::cout << "===============================================================|" << std::endl;
    std::cout << "      System Summary                                           |" << std::endl;
    std::cout << "===============================================================|" << std::endl;
    std::cout << "Buses            " << Nbus << std::endl;
    std::cout << "Branches         " << Nline << std::endl;
    std::cout << "Ploss            " << getPloss() << std::endl;
    std::cout << "Qloss            " << getQloss() << std::endl;


    std::cout << std::endl << std::endl;
    if (all) {

        std::cout << "===============================================================================================|" << std::endl;
        std::cout << "      Bus Data                                                                                 |" << std::endl;
        std::cout << "===============================================================================================|" << std::endl;
        std::cout << " Bus |          Voltage        |  Power = Generation  + Load   |  Init = Generation  + Load    |" << std::endl;
        std::cout << "  #  |    Mag(pu) |  Ang(deg)  |    P (pu)     |     Q (pu)    |    P (pu)     |     Q (pu)    |" << std::endl;
        std::cout << "-----|------------|------------|---------------|---------------|---------------|---------------|" << std::endl;

        //std::cout << 0 << "      " << E.get(Nbus, 0) << "             " << E.get(0, 0) * (abs(E.get(0, 0)) > 0.0001) * 180 / 3.1415 << "              " << (abs(W.get(0, 0)) > 0.0001) * W.get(0, 0) << "         " << (abs(W.get(Nbus, 0)) > 0.0001) * W.get(Nbus, 0) << std::endl;

        float seuil = 0.0001;
       
        
        std::cout << std::setw(5) << 0 << "|" << std::setw(11) << E.get(Nbus, 0) << "*|" << std::setw(11) << E.get(0, 0) * (abs(E.get(0, 0)) > seuil) * 180 / 3.1415
            << "*|" << std::setw(15) << (abs(W.get(0, 0)) > seuil) * W.get(0, 0) << "|" << std::setw(15) << (abs(W.get(Nbus, 0)) > seuil) * W.get(Nbus, 0)
            << "|" << std::setw(15) << W0.get(0, 0) << "|" << std::setw(15)
            << W0.get(Nbus, 0) << "|" << std::endl;
        for (int b = 1; b < Nbus; b++) {
            //std::cout.width(10);
            //std::cout << b << "      " << E.get(b + Nbus, 0) << "        " << E.get(b, 0) * (abs(E.get(b, 0)) > 0.0001) * 180 / 3.1415 << "          " << (abs(W.get(b, 0)) > 0.0001) * W.get(b, 0) << "         " << (abs(W.get(b + Nbus, 0)) > 0.0001) * W.get(b + Nbus, 0) << std::endl;
            std::cout << std::setw(5) << b << "|" << std::setw(11) << E.get(b + Nbus, 0) << " |" << std::setw(11)
                << E.get(b, 0) * (abs(E.get(b, 0)) > seuil) * 180 / 3.1415 << " |" << std::setw(15)
                << (abs(W.get(b, 0)) > seuil) * W.get(b, 0) << "|" << std::setw(15) << (abs(W.get(b + Nbus, 0)) > seuil) * W.get(b + Nbus, 0)
                << "|" << std::setw(15) << W0.get(b, 0) << "|" << std::setw(15)
                << W0.get(b + Nbus, 0) << "|" << std::endl;

        }
        
    }
    else {
        float seuil = 0.0001;
        std::cout << " Bus |          Voltage        |  Power = Generation  + Load   |  Init = Generation  + Load    |" << std::endl;
        std::cout << "  #  |    Mag(pu) |  Ang(deg)  |    P (pu)     |     Q (pu)    |    P (pu)     |     Q (pu)    |" << std::endl;
        std::cout << "-----|------------|------------|---------------|---------------|---------------|---------------|" << std::endl;
        
        
         std::cout << std::setw(5) << 0 << "|" << std::setw(11) << E.get(Nbus, 0) << "*|" << std::setw(11) << E.get(0, 0) * (abs(E.get(0, 0)) > seuil) * 180 / 3.1415
                << "*|" << std::setw(15) << (abs(W.get(0, 0)) > seuil) * W.get(0, 0) << "|" << std::setw(15) << (abs(W.get(Nbus, 0)) > seuil) * W.get(Nbus, 0)
                << "|" << std::setw(15) << W0.get(0, 0) << "|" << std::setw(15)
                << W0.get(Nbus, 0) << "|" << std::endl;
        
    }

    std::cout << "===============================================================================================|" << std::endl;
    std::cout << "                      END PRINT                                                                |" << std::endl;
    std::cout << "===============================================================================================|" << std::endl;
}



/*
for (int l = 0; l < Nline; l++) {
        int k = l + 1; // bus to
        float vRe = VoltageRealImPre.get(k, 0);
        float vIm = VoltageRealImPre.get(k + Nbus, 0);
        float p = -W0.get(k, 0);
        float q = -W0.get(k + Nbus, 0);
        float y = l > 0 ? Yd.get(l, 0) : 0;

        float SRe =  p;
        float SIm =  q - y * (vRe * vRe + vIm * vIm);
        St.set(l, 0, SRe);
        St.set(l + Nline, 0, SIm);
    }

*/
__global__ void calculStGPU(float* St, float* Voltage, float* W0, float* Yd, int B) {
    int thIdx = threadIdx.x + blockDim.x * blockIdx.x;
    int size = gridDim.x * blockDim.x;
    for (int k = thIdx + 1; k < B; k += size) { // bus to
        int l = k - 1; // line
        float vRe = Voltage[k];
        float vIm = Voltage[k + B];
        float p = -W0[k];
        float q = -W0[k + B];
        float y = Yd[l];

        float SRe = p;
        float SIm = q - y * (vRe * vRe + vIm * vIm);
        St[l] = SRe;
        St[l + B - 1] = SIm;

    }
}



/*
for (int k = 1; k < Nbus; k++) {
        // branch l entre le bus i=F(k) et k = l + 1;

        int i = F.get(k, 0); // busFrom
        int l = k - 1; // line


        float vRe = VoltageRealIm.get(i, 0);
        float vIm = VoltageRealIm.get(i  + Nbus, 0);
        float V2 = vRe * vRe + vIm * vIm;

        float zRe = ZsRe.get(l, 0);
        float zIm = ZsIm.get(l, 0);
        float SRe = Sf.get(l, 0);
        float SIm = Sf.get(l + Nline, 0);

        float vRe2 = vRe - (zRe * (SRe * vRe + SIm * vIm) + zIm * (SIm * vRe - SRe * vIm)) / V2;
        float vIm2 = vIm - (zIm * (SRe * vRe + SIm * vIm) - zRe * (SIm * vRe - SRe * vIm)) / V2;

        VoltageRealIm.set(k, 0, vRe2);
        VoltageRealIm.set(k + Nbus, 0, vIm2);
}
*/

__global__ void updateVoltage(float* Voltage, float* ZsRe, float* ZsIm, float* Sf, float* F, int B, int LastBus) {
    
    int line = threadIdx.x;
    int size = blockDim.x;
    // un seul block
    extern __shared__ float globalMemory[];
    float* VoltageSh = globalMemory;
    bool* hasFinished = (bool*) &globalMemory[2 * B];

    __shared__ bool notfinished;
    bool mustCompute = false;
    int L = B - 1;
    int bus = line + 1;
    if (line == 0) {
        notfinished = true;
    }
    for (int k = line; k < 2 * B; k += size) {
        VoltageSh[k] = Voltage[k];
    }
    __syncthreads();
    if (line < L) {
        hasFinished[line] = false;
        int busFrom = F[bus];
        float zRe = ZsRe[line];
        float zIm = ZsIm[line];
        float SRe = Sf[line];
        float SIm = Sf[line + L];
        mustCompute = (busFrom == 0);
    
        while (notfinished) {
            if (mustCompute) { // divergent mais on n'y peut rien
                float vRe = VoltageSh[busFrom];
                float vIm = VoltageSh[busFrom + B];
                float V2 = vRe * vRe + vIm * vIm;


                float vRe2 = vRe - (zRe * (SRe * vRe + SIm * vIm) + zIm * (SIm * vRe - SRe * vIm)) / V2;
                float vIm2 = vIm - (zIm * (SRe * vRe + SIm * vIm) - zRe * (SIm * vRe - SRe * vIm)) / V2;

                VoltageSh[bus] = vRe2;
                VoltageSh[bus + B] = vIm2;


                hasFinished[line] = true;
                if (line == LastBus - 1) {
                    notfinished = false;
                }
            }
            __syncthreads();
            // trouver qui doit tourner � la prochaine boucle
            mustCompute = !(hasFinished[line]) && hasFinished[busFrom - 1];
           /* if (line == 0) {
                notfinished = false;
            }
            __syncthreads();
                // tous ecrive la même chose
            if (!hasFinished[line]) {
                notfinished = true;
            }*/
            
                __syncthreads();
            }
        }
    __syncthreads();
    for (int k = line; k < 2 * B; k += size) {
        Voltage[k] = VoltageSh[k];
    }
}


/*

for (int l = Nline-1; l >= 0; l--) {
        int k = l + 1; // busTo
        int i = F.get(k, 0); // busFrom
        int lprev = i - 1;

        float SRe = St.get(l,0);
        float SIm = St.get(l + Nline, 0);
        float vRe = VoltageRealImPre.get(k, 0);
        float vIm = VoltageRealImPre.get(k + Nbus, 0);

        float SfRe = SRe + ZsRe.get(l, 0) * (SRe * SRe + SIm * SIm) / (vRe * vRe + vIm * vIm);
        float SfIm = SIm + ZsIm.get(l, 0) * (SRe * SRe + SIm * SIm) / (vRe * vRe + vIm * vIm);

        Sf.set(l, 0, SfRe);
        Sf.set(l + Nline, 0, SfIm);
        if (lprev > -1) {
            St.increment(lprev, 0, Sf.get(l, 0));
            St.increment(lprev + Nline, 0, Sf.get(l + Nline, 0));
        }

    }
  */



__global__ void calculSGPU(float* St, float* Sf, float* Voltage, float* ZsRe, float* ZsIm, float* nChild, float* Childs, float* indiceChildBegin, int B) {

    int L = B - 1;
    extern __shared__ float  globalMemory2[];
    float* SfSh        = (float*) globalMemory2;
    int*   ChildsSh    = (int*)  &globalMemory2[2 * L];
    bool*  hasfinished = (bool*) &ChildsSh[L];
   

  /**/


    __shared__ bool notfinished;

    int line = threadIdx.x;
    int step = blockDim.x;
    
    if (line == 0) {
        notfinished = true;
    }/**/
    for (int l = line; l < L; l += step) {
        hasfinished[l] = false;
        ChildsSh[l] = Childs[l];
        SfSh[l] = 0.0f;
        SfSh[l + L] = 0;
    } 

    __syncthreads();
   

     
   
   
    if (line < L) {
        //hasfinished[line] = false;
        int bus = line + 1;
        int indiceChild = (bus < (B - 1)) ? indiceChildBegin[bus] : 0;
        int nb = nChild[bus];
        bool mustCompute = (nb == 0);
        float vRe = Voltage[bus];
        float vIm = Voltage[bus + B];
        float vNorm = (vRe * vRe + vIm * vIm);
        float ZRe = ZsRe[line];
        float ZIm = ZsIm[line];
        float StRe = St[line];
        float StIm = St[line + L];
   
        __syncthreads();/**/
        while (notfinished) {
            if (mustCompute) { // divergent mais on n'y peut rien
                for (int i = 0; i < nb; i++) { // calcul St 
                    int c = ChildsSh[indiceChild + i];
                    int lineChild = c - 1;
                    StRe += SfSh[lineChild];
                    StIm += SfSh[lineChild + L];
                }
                St[line] = StRe;
                St[line + L] = StIm;

                // calcul Sf
           
                float StNorm = (StRe * StRe + StIm * StIm);

                float SfRe = StRe + ZRe * StNorm / vNorm;
                float SfIm = StIm + ZIm * StNorm / vNorm;

                SfSh[line] = SfRe;
                SfSh[line + L] = SfIm;
                
          
                hasfinished[line] = true;
                
                
                if (line == 0) {
                    notfinished = false;
                }
            }
            __syncthreads();
            // trouver qui doit tourner � la prochaine boucle
            mustCompute = !(hasfinished[line]);

            for (int i = 0; i < nb; i++) {
                int c = ChildsSh[indiceChild + i];
                int lineChild = c - 1;
                mustCompute = (mustCompute && hasfinished[lineChild]); // il suffit qu'un enfant n'a pas fini pour que cela soit false
            }
            __syncthreads();
        }
    }
   
    for (int l = line; l < L; l += step) {
        Sf[l] = SfSh[l];
        Sf[l + L] = SfSh[l + L];
    }

   

}




/*for (int b = 0; b < B2; b++) {
        Y.set(b, 0, E.get(b, 0));
    }
    int line = 0;

    for (int i = 0; i < Nbus; i++) {
        int k = CoresBusLin.get(i, 0);

        float ei = VoltageRealIm.get(i, 0);
        float fi = VoltageRealIm.get(Nbus + i, 0);

        for (int voisin = k + 1; voisin < (k + nLines.get(i, 0)); voisin++) {
            int j = CoresVoiLin.get(voisin, 0);
            float ej = VoltageRealIm.get(j, 0);
            float fj = VoltageRealIm.get(j + Nbus, 0);
            if (j > i) {
                float B = BgridLin.get(voisin, 0);
                float G = GgridLin.get(voisin, 0);
                float Pij = (ei * ei + fi * fi - ei * ej - fi * fj) * G + (ei * fj - ej * fi) * B;

                Y.set(2 * Nbus + line, 0, Pij);

                line++;
            }
        }
    }*/



