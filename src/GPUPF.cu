#include "../head/GPUPF.cuh"


GPUPF::GPUPF(){}
GPUPF::~GPUPF(){
   

}
void GPUPF::init(const StudyCase& cas, MatrixGPU* PQ)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
    timePerBlock = MatrixCPU(1, 9); // Fb0 : init, Fb1ab : Flu, Fb2abc: Tension , FB3 : puissance, Fb4 erreur, Fb0 mise � jour

    occurencePerBlock = MatrixCPU(1, 9);; //nb de fois utilis� pendant la simu
#endif // INSTRUMENTATION
    
   // std::cout << "init PF NR GPU simple" << std::endl;
    std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
    Nagent = cas.getNagent();
    Nbus = cas.getNBus();
    B2 = 2 * Nbus;
    N2 = 2 * Nagent;
    Nline = cas.getNLine(true); // ne doit pas �tre r�duit ici !!!
    BL2 = Nbus + 2 * Nline;
    Nconstraint = B2 + Nline;
    iterM = 20;
    iter = 0;
    V0 = cas.getV0();
    theta0 = cas.gettheta0();
    I = MatrixGPU(cas.getCoresBusAgentLin(), 1);
    status = 0;

    CoresAgentBus = MatrixGPU(cas.getCoresAgentBusLin(), 1);
    CoresAgentBusBegin = MatrixGPU(cas.getCoresAgentBusLinBegin(), 1);
    NagentByBus = MatrixGPU(cas.getNagentByBus(), 1);
    removeLossAgent << <1, 1 >> > (NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU);
    //I.display(true);
    _name = "Newton";
    numBlock = Nbus;
    _useDouble = false;


    Bgrid = cas.getLineSuceptance();
    Ggrid = cas.getLineReactance();



    Y = MatrixGPU(2 * Nbus + Nline, 1, 0, 1);
    CoresLineBus = MatrixGPU(cas.getCoresLineBus(true));
    _CoresVoiLin = MatrixGPU(cas.getCoresVoiLin(), 1);
    _CoresBusLin = MatrixGPU(cas.getCoresBusLin(), 1);
    _nLines = MatrixGPU(cas.getNLines(), 1);
    CoresLineBusGPU = MatrixGPU(2, Nline);

    for (int lold = 0; lold < Nline; lold++) {
        int busTo = CoresLineBus.get(lold, 1);
        int busFrom = CoresLineBus.get(lold, 0);
        CoresLineBusGPU.set(0, lold, busFrom);
        CoresLineBusGPU.set(1, lold, busTo);
    }
    CoresLineBusGPU.transferGPU();
    _Blin2 = MatrixGPU(cas.getBlin2(), 1);
    _Glin2 = MatrixGPU(cas.getGlin2(), 1);
    Phi = MatrixGPU(Nline, 1, 0, 1);
    //_nLines.display(true);


    
    _Blin = MatrixGPU(cas.getBlin(), 1);
    _Glin = MatrixGPU(cas.getGlin(), 1);

    W = MatrixGPU(B2, 1, 0, 1);
    _Pintermediate = MatrixGPU(BL2, 1, 0, 1);
    _Qintermediate = MatrixGPU(BL2, 1, 0, 1);
    dW = MatrixGPU(B2, 1, 0, 1);
    dW.preallocateReduction();
    E = MatrixGPU(B2, 1, 0, 1);
    dE = MatrixGPU(B2, 1, 0, 1);
    Jac = MatrixGPU(B2, B2);
    Jac.set(0, 0, 1);
    Jac.set(Nbus, Nbus, 1);
    Jac.transferGPU();
    JacInv = MatrixGPU(B2, B2, 0, 1);

    initE << <numBlock, _blockSize >> > (E._matrixGPU, theta0, V0, Nbus);
    //E.display(true);
    /*std::cout << " Bgrid : " << std::endl;
    Bgrid.display();
    std::cout << " Ggrid : " << std::endl;
    Ggrid.display();*/
    W0 = MatrixGPU(B2, 1, 0, 1);

    calculW0Bis(PQ);

    /*std::cout << " PQ : " << std::endl;
    PQ->display(true);
    std::cout << " W0 : " << std::endl;
    W0.display(true);

    CoresAgentBus.display(true);
    CoresAgentBusBegin.display(true);
    NagentByBus.display(true);
    std::cout << "N: " << Nagent << " B= " << Nbus << std::endl;*/


    //std::cout << " W0 : " << std::endl;
    //W0.display(true);
    calcW();

    dW.subtract(&W0, &W);
    A = MatrixGPU(B2, B2, 0, 1);
    P = MatrixGPU(B2 + 1, 1, 0, 1);

    
#ifdef INSTRUMENTATION
    cudaDeviceSynchronize();
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 0, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION

    //std::cout << " fin init" << std::endl;

}
void GPUPF::init(const StudyCase& cas, MatrixGPU* PQ, MatrixGPUD* PQD, bool useDouble)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
    timePerBlock = MatrixCPU(1, 9); // Fb0 : init, Fb1ab : Flu, Fb2abc: Tension , FB3 : puissance, Fb4 erreur, Fb0 mise � jour

    occurencePerBlock = MatrixCPU(1, 9);; //nb de fois utilis� pendant la simu
#endif // INSTRUMENTATION
   // std::cout << "init PF NR GPU" << std::endl;
    std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
    Nagent = cas.getNagent();
    Nbus = cas.getNBus();
    B2 = 2 * Nbus;
    N2 = 2 * Nagent;
    Nline = cas.getNLine(true); // ne doit pas �tre r�duit ici !!!
    BL2 = Nbus + 2 * Nline;
    Nconstraint = B2 + Nline;
    iterM = 20;
    iter = 0;
    V0 = cas.getV0();
    theta0 = cas.gettheta0();
    I = MatrixGPU(cas.getCoresBusAgentLin(),1);
    status = 0;

    CoresAgentBus = MatrixGPU(cas.getCoresAgentBusLin(), 1);
    CoresAgentBusBegin = MatrixGPU(cas.getCoresAgentBusLinBegin(), 1);
    NagentByBus = MatrixGPU(cas.getNagentByBus(), 1);
    //I.display(true);
    _name = "Newton";
    numBlock = Nbus;
    _useDouble = useDouble;
   

    Bgrid = cas.getLineSuceptance();
    Ggrid = cas.getLineReactance();

    
     
    Y = MatrixGPU(2 * Nbus + Nline, 1, 0, 1);
    CoresLineBus = MatrixGPU(cas.getCoresLineBus(true));
    _CoresVoiLin = MatrixGPU(cas.getCoresVoiLin(), 1);
    _CoresBusLin = MatrixGPU(cas.getCoresBusLin(), 1);
    _nLines = MatrixGPU(cas.getNLines(), 1);
    CoresLineBusGPU = MatrixGPU(2, Nline);
    
    for (int lold = 0; lold < Nline; lold++) {
        int busTo = CoresLineBus.get(lold, 1);
        int busFrom = CoresLineBus.get(lold, 0);
        CoresLineBusGPU.set(0, lold, busFrom);
        CoresLineBusGPU.set(1, lold, busTo);
    }
    CoresLineBusGPU.transferGPU();
    _Blin2 = MatrixGPU(cas.getBlin2(), 1);
    _Glin2 = MatrixGPU(cas.getGlin2(), 1);
    //_nLines.display(true);


    if (_useDouble) {
        
        _BlinD = MatrixGPUD(cas.getBlinD(), 1);
        _GlinD = MatrixGPUD(cas.getGlinD(), 1);

        WD = MatrixGPUD(B2, 1, 0, 1);
        _PintermediateD = MatrixGPUD(BL2, 1, 0, 1);
        _QintermediateD = MatrixGPUD(BL2, 1, 0, 1);
      

        dWD = MatrixGPUD(B2, 1, 0, 1);
        dWD.preallocateReduction();
        ED = MatrixGPUD(B2, 1, 0, 1);
        initED <<<numBlock,_blockSize>>>(ED._matrixGPU, theta0, V0, Nbus);
        
        //ED.display(true);

        dED = MatrixGPUD(B2, 1, 0, 1);
        JacD = MatrixGPUD(B2, B2);
        if (JacD.getPos()) {
            JacD.transferCPU();
        }
        JacD.set(0, 0, 1);
        JacD.set(Nbus, Nbus, 1);
        JacD.transferGPU();

        JacInvD = MatrixGPUD(B2, B2, 0, 1);
        W0D = MatrixGPUD(B2, 1, 0, 1);

        /*std::cout << " PQ : " << std::endl;
        PQD->display(true);
        std::cout << " W0 : " << std::endl;
        W0D.display(true);
       
        CoresAgentBus.display(true);
        CoresAgentBusBegin.display(true);
        NagentByBus.display(true);

        std::cout << "N: " << Nagent<< " B= "<< Nbus << std::endl;*/
        calculW0DBis(PQD);
        /*ED.display(true);
        _GlinD.display(true);
        _BlinD.display(true);
        std::cout << "------------" << std::endl;
        _CoresVoiLin.display(true);
        _CoresBusLin.display(true);
        _nLines.display(true);
        std::cout << "------------" << std::endl;

        _PintermediateD.display(true);
        _QintermediateD.display(true);*/
        //std::cout << "*******" << std::endl;



        //std::cout << " W0 : " << std::endl;
        //W0D.display(true);

        calcW();

        dWD.subtract(&W0D, &WD);
        AD = MatrixGPUD(B2, B2, 0, 1);
        PD = MatrixGPUD(B2 + 1, 1, 0, 1);
    
    }
    else {
        
        _Blin = MatrixGPU(cas.getBlin(), 1);
        _Glin = MatrixGPU(cas.getGlin(), 1);
           
        W = MatrixGPU(B2, 1, 0, 1);
        _Pintermediate = MatrixGPU(BL2, 1, 0, 1);
        _Qintermediate = MatrixGPU(BL2, 1, 0, 1);
        dW = MatrixGPU(B2, 1, 0, 1);
        dW.preallocateReduction();
        E = MatrixGPU(B2, 1, 0, 1);
        dE = MatrixGPU(B2, 1, 0, 1);
        Jac = MatrixGPU(B2, B2);
        Jac.set(0, 0, 1);
        Jac.set(Nbus, Nbus, 1);
        Jac.transferGPU();
        JacInv = MatrixGPU(B2, B2, 0, 1);

        initE << <numBlock, _blockSize >> > (E._matrixGPU, theta0, V0, Nbus);
        //E.display(true);
        /*std::cout << " Bgrid : " << std::endl;
        Bgrid.display();
        std::cout << " Ggrid : " << std::endl;
        Ggrid.display();*/
        W0 = MatrixGPU(B2, 1, 0, 1);

        calculW0Bis(PQ); 
       
        /*std::cout << " PQ : " << std::endl;
        PQ->display(true);
        std::cout << " W0 : " << std::endl;
        W0.display(true);

        CoresAgentBus.display(true);
        CoresAgentBusBegin.display(true);
        NagentByBus.display(true);
        std::cout << "N: " << Nagent << " B= " << Nbus << std::endl;*/
      
        
        //std::cout << " W0 : " << std::endl;
        //W0.display(true);
        calcW();
    
        dW.subtract(&W0, &W);
        A = MatrixGPU(B2, B2, 0, 1);
        P = MatrixGPU(B2 + 1, 1, 0, 1);
      
    }

   
        
    /*Ggrid2Bgrid2 = MatrixGPU(Nbus, Nbus);
    for (int i = 0; i < Nbus; i++) {
        for (int j = 0; j < Nbus; j++) {
            Ggrid2Bgrid2.set(i, j, sqrt(Ggrid.get(i, j) * Ggrid.get(i, j) + Bgrid.get(i, j) * Bgrid.get(i, j)));
        }
    }

    
    G = MatrixGPU(Nconstraint, N2);
    Phi = MatrixGPU(Nline, 1);
    Y = MatrixGPU(Nconstraint, 1);
    tempLN2 = MatrixGPU(Nline, N2);
    JacPhiE = MatrixGPU(Nline, B2);
    tempB2N2 = MatrixGPU(B2, N2);*/
    
    
#ifdef INSTRUMENTATION
    cudaDeviceSynchronize();
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 0, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION

    //std::cout << " fin init" << std::endl;

}


void GPUPF::solve() {
  
    //std::cout << "solve Newton" << std::endl;
    time = clock();
    err = 2 * epsPF;
    iter = 0;
    int failure = 0;
    status = 1;
    //std::cout << epsPF << " " << iterM << std::endl;
    while (err > epsPF && iter < iterM) {
        
        failure = calcVoltage();
        if (failure) {
            status = -1;
            time = clock() - time;
            //std::cout << "failure ! " << iter << " " << err << std::endl;
            return;
        }
#ifdef INSTRUMENTATION
        t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
        calcW();
#ifdef INSTRUMENTATION
        cudaDeviceSynchronize();
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 6, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 6, 1);
        t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
        if (_useDouble) {
            dWD.subtract(&W0D, &WD); // dW = W0 - W
            err = dWD.max2(); //err = ||dW||
        }
        else {
            dW.subtract(&W0, &W); // dW = W0 - W
            err = dW.max2(); //err = ||dW||
        }
#ifdef INSTRUMENTATION
        cudaDeviceSynchronize();
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 7, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 7, 1);
#endif // INSTRUMENTATION      
        iter++;
        //std::cout << err << " * ";
    }
    //std::cout << std::endl;
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    calcW(true);
#ifdef INSTRUMENTATION
    cudaDeviceSynchronize();
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 6, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 6, 1);
#endif // INSTRUMENTATION

    if (iter >= iterM) {
        status = 2;
        if (err > 100 * epsPF) {
            status = -1;
        }
        //std::cout << "fin solve " << iter<<" " << err << std::endl;
    }
    
    time = clock() - time;
        
}

void GPUPF::updatePQ(MatrixGPU* PQ)
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

void GPUPF::calculW0(MatrixGPU* PQ)
{
    switch (_blockSize) {
    case 512:
        calcW0GPU<512> << <numBlock, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case 256:
        calcW0GPU<256> << <numBlock, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case 128:
        calcW0GPU<128> << <numBlock, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case 64:
        calcW0GPU< 64> << <numBlock, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case 32:
        calcW0GPU< 32> << <numBlock, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case 16:
        calcW0GPU< 16> << <numBlock, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case  8:
        calcW0GPU<  8> << <numBlock, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case  4:
        calcW0GPU<  4> << <numBlock, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case  2:
        calcW0GPU<  2> << <numBlock, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case  1:
        calcW0GPU<  1> << <numBlock, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    }
}

void GPUPF::calculW0D(MatrixGPUD* PQD)
{
    switch (_blockSize) {
    case 512:
        calcW0GPUD<512> << <Nbus, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case 256:
        calcW0GPUD<256> << <Nbus, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case 128:
        calcW0GPUD<128> << <Nbus, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case 64:
        calcW0GPUD< 64> << <Nbus, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case 32:
        calcW0GPUD< 32> << <Nbus, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case 16:
        calcW0GPUD< 16> << <Nbus, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case  8:
        calcW0GPUD<  8> << <Nbus, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case  4:
        calcW0GPUD<  4> << <Nbus, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case  2:
        calcW0GPUD<  2> << <Nbus, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    case  1:
        calcW0GPUD<  1> << <Nbus, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, I._matrixGPU, Nagent, Nbus);
        break;
    }
}

void GPUPF::calculW0Bis(MatrixGPU* PQ)
{
    // prend en compte le premier agent !!! 
    
    switch (_blockSize) {
    case 512:
        calcW0GPUBis<512> << <Nbus, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case 256:
        calcW0GPUBis<256> << <Nbus, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case 128:
        calcW0GPUBis<128> << <Nbus, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case 64:
        calcW0GPUBis< 64> << <Nbus, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case 32:
        calcW0GPUBis< 32> << <Nbus, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case 16:
        calcW0GPUBis< 16> << <Nbus, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case  8:
        calcW0GPUBis<  8> << <Nbus, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case  4:
        calcW0GPUBis<  4> << <Nbus, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case  2:
        calcW0GPUBis<  2> << <Nbus, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case  1:
        calcW0GPUBis<  1> << <Nbus, _blockSize >> > (W0._matrixGPU, PQ->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    }
}

void GPUPF::calculW0DBis(MatrixGPUD* PQD)
{
    switch (_blockSize) {
    case 512:
        calcW0GPUDBis<512> << <numBlock, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case 256:
        calcW0GPUDBis<256> << <numBlock, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case 128:
        calcW0GPUDBis<128> << <numBlock, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case 64:
        calcW0GPUDBis< 64> << <numBlock, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case 32:
        calcW0GPUDBis< 32> << <numBlock, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case 16:
        calcW0GPUDBis< 16> << <numBlock, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case  8:
        calcW0GPUDBis<  8> << <numBlock, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case  4:
        calcW0GPUDBis<  4> << <numBlock, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case  2:
        calcW0GPUDBis<  2> << <numBlock, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    case  1:
        calcW0GPUDBis<  1> << <numBlock, _blockSize >> > (W0D._matrixGPU, PQD->_matrixGPU, CoresAgentBus._matrixGPU, NagentByBus._matrixGPU, CoresAgentBusBegin._matrixGPU, Nagent, Nbus);
        break;
    }
}


void GPUPF::calcW(bool end)
{
    /*
    ED.display(true);
        _GlinD.display(true);
        _BlinD.display(true); 
        std::cout << "------------" << std::endl;
        _CoresVoiLin.display(true);
        _CoresBusLin.display(true);
        _nLines.display(true);
        std::cout << "------------" << std::endl;
    _PintermediateD.display(true);
        _QintermediateD.display(true);
    _Qintermediate.display(true);
    std::cout << "------------" << std::endl;*/ 
    if (_useDouble) {
        
        calcWinterD << <numBlock, _blockSize, B2 * sizeof(double) >> > (_PintermediateD._matrixGPU, _QintermediateD._matrixGPU, ED._matrixGPU, _GlinD._matrixGPU, _BlinD._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
        
        if (!end) { // pendant simu, la puissance � ce noeud est libre
            switch (_blockSize) {
            case 512:
                calcWGPUD<512> << <numBlock, _blockSize >> > (WD._matrixGPU, W0D._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 256:
                calcWGPUD<256> << <numBlock, _blockSize >> > (WD._matrixGPU, W0D._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 128:
                calcWGPUD<128> << <numBlock, _blockSize >> > (WD._matrixGPU, W0D._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 64:
                calcWGPUD< 64> << <numBlock, _blockSize >> > (WD._matrixGPU, W0D._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 32:
                calcWGPUD< 32> << <numBlock, _blockSize >> > (WD._matrixGPU, W0D._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 16:
                calcWGPUD< 16> << <numBlock, _blockSize >> > (WD._matrixGPU, W0D._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case  8:
                calcWGPUD<  8> << <numBlock, _blockSize >> > (WD._matrixGPU, W0D._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case  4:
                calcWGPUD<  4> << <numBlock, _blockSize >> > (WD._matrixGPU, W0D._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case  2:
                calcWGPUD<  2> << <numBlock, _blockSize >> > (WD._matrixGPU, W0D._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case  1:
                calcWGPUD<  1> << <numBlock, _blockSize >> > (WD._matrixGPU, W0D._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            }
        }
        else {
            switch (_blockSize) {
            case 512:
                calcWGPUD<512> << <numBlock, _blockSize >> > (WD._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 256:
                calcWGPUD<256> << <numBlock, _blockSize >> > (WD._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 128:
                calcWGPUD<128> << <numBlock, _blockSize >> > (WD._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 64:
                calcWGPUD< 64> << <numBlock, _blockSize >> > (WD._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 32:
                calcWGPUD< 32> << <numBlock, _blockSize >> > (WD._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 16:
                calcWGPUD< 16> << <numBlock, _blockSize >> > (WD._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case  8:
                calcWGPUD<  8> << <numBlock, _blockSize >> > (WD._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case  4:
                calcWGPUD<  4> << <numBlock, _blockSize >> > (WD._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case  2:
                calcWGPUD<  2> << <numBlock, _blockSize >> > (WD._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case  1:
                calcWGPUD<  1> << <numBlock, _blockSize >> > (WD._matrixGPU, _PintermediateD._matrixGPU, _QintermediateD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            }

        }
    }
    else {
       /* E.display(true);
        _Glin.display(true);
        _Blin.display(true);
        std::cout << "------------" << std::endl;
        _CoresVoiLin.display(true);
        _CoresBusLin.display(true);
        _nLines.display(true);
        std::cout << "------------" << std::endl;*/

        calcWinter << <numBlock, _blockSize, B2 * sizeof(float) >> > (_Pintermediate._matrixGPU, _Qintermediate._matrixGPU, E._matrixGPU, _Glin._matrixGPU, _Blin._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);

        /*_Pintermediate.display(true);
        _Qintermediate.display(true);*/
        if (!end) { // pendant simu, la puissance � ce noeud est libre
            switch (_blockSize) {
            case 512:
                calcWGPU<512> << <numBlock, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 256:
                calcWGPU<256> << <numBlock, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 128:
                calcWGPU<128> << <numBlock, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 64:
                calcWGPU< 64> << <numBlock, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 32:
                calcWGPU< 32> << <numBlock, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 16:
                calcWGPU< 16> << <numBlock, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case  8:
                calcWGPU<  8> << <numBlock, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case  4:
                calcWGPU<  4> << <numBlock, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case  2:
                calcWGPU<  2> << <numBlock, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case  1:
                calcWGPU<  1> << <numBlock, _blockSize >> > (W._matrixGPU, W0._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            }
        }
        else {
            switch (_blockSize) {
            case 512:
                calcWGPU<512> << <numBlock, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 256:
                calcWGPU<256> << <numBlock, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 128:
                calcWGPU<128> << <numBlock, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 64:
                calcWGPU< 64> << <numBlock, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 32:
                calcWGPU< 32> << <numBlock, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case 16:
                calcWGPU< 16> << <numBlock, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case  8:
                calcWGPU<  8> << <numBlock, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case  4:
                calcWGPU<  4> << <numBlock, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case  2:
                calcWGPU<  2> << <numBlock, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            case  1:
                calcWGPU<  1> << <numBlock, _blockSize >> > (W._matrixGPU, _Pintermediate._matrixGPU, _Qintermediate._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
                break;
            }

        }

    }
   
    //W.display(true);
    
}

void GPUPF::calcJac()
{
    if (_useDouble) {
        calcJacGPUD <<<numBlock, _blockSize, B2 * sizeof(double) >>> (JacD._matrixGPU, WD._matrixGPU, ED._matrixGPU, _GlinD._matrixGPU, _BlinD._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
    }
    else {
        calcJacGPU <<<numBlock, _blockSize,  B2 * sizeof(float) >>> (Jac._matrixGPU, W._matrixGPU, E._matrixGPU, _Glin._matrixGPU, _Blin._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
    }
}

void GPUPF::calcPhi()
{
    calcE();
    
    calculPhiGPU << <numBlock, _blockSize >> > (Phi._matrixGPU, E._matrixGPU, _Blin2._matrixGPU, _Glin2._matrixGPU, CoresLineBusGPU._matrixGPU, Nbus, Nline);
}

void GPUPF::calcJacPhiE()
{
    for (int l = 0; l < Nline; l++) { //angle
        int i = CoresLineBus.get(l, 0); //from 
        int i2 = i + Nbus;
        int j = CoresLineBus.get(l, 1); // to
        int j2 = j + Nbus;
        
        JacPhiE.set(l, i2, E.get(j2, 0) * Ggrid2Bgrid2.get(i, j));
        JacPhiE.set(l, j2, E.get(i2, 0) * Ggrid2Bgrid2.get(i, j));

        /*float dTheta_ij = E.get(i, 0) - E.get(j, 0);

        JacPhiE.set(l, i, -E.get(i2, 0) * E.get(j2, 0) * (Ggrid.get(i, j) * sin(dTheta_ij) - Bgrid.get(i, j) * cos(dTheta_ij)));
        JacPhiE.set(l, j,  E.get(i2, 0) * E.get(j2, 0) * (Ggrid.get(i, j) * sin(dTheta_ij) - Bgrid.get(i, j) * cos(dTheta_ij)));
        JacPhiE.set(l, i2, E.get(j2, 0) *                (Ggrid.get(i, j) * cos(dTheta_ij) + Bgrid.get(i, j) * sin(dTheta_ij)));
        JacPhiE.set(l, j2, E.get(i2, 0) *                (Ggrid.get(i, j) * cos(dTheta_ij) + Bgrid.get(i, j) * sin(dTheta_ij)));
        */
    }

}

int GPUPF::calcVoltage()
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    calcJac();
#ifdef INSTRUMENTATION
    cudaDeviceSynchronize();
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 3, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 3, 1);
#endif // INSTRUMENTATION
    


    if (_useDouble) {
#ifdef INSTRUMENTATION
        t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
        try
        {
            JacD.LUPFactorization(&AD, &PD);
        }
        catch (const std::exception&)
        {
            return 1; // failure
        }
        
        /*AD.display(true);
        PD.display(true);
         dWD.display(true);*/
#ifdef INSTRUMENTATION
        cudaDeviceSynchronize();
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 4, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 4, 1);
#endif // INSTRUMENTATION

        dED.solveSys(&AD, &PD, &dWD);
#ifdef INSTRUMENTATION
        cudaDeviceSynchronize();
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 5, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 5, 1);
#endif // INSTRUMENTATION
        //dED.display(true);
        //std::cout << "**********" << std::endl;
       
       /*try
        {
            JacInvD.invertGaussJordan(&JacD);
        }
        catch (const std::exception&)
        {
            JacD.display(true);
            exit(0);
        }
        dED.multiply(&JacInvD, &dWD);*//**/// dE = Jac_inv * dW;
        ED.add(&ED, &dED);// E = E + dE;/*
    }
    else {
        //std::cout << " Jac : " << std::endl;
        /*try
        {
            JacInv.invertGaussJordan(&Jac);
        }
        catch (const std::exception&)
        {
            exit(0);
        }
        dE.multiply(&JacInv, &dW);*/// dE = Jac_inv * dW;
#ifdef INSTRUMENTATION
        t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
        try
        {
            Jac.LUPFactorization(&A, &P);
        }
        catch (const std::exception&)
        {
            return 1;
        }
#ifdef INSTRUMENTATION
        cudaDeviceSynchronize();
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 4, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 4, 1);
#endif // INSTRUMENTATION
        dE.solveSys(&A, &P, &dW);/**/
#ifdef INSTRUMENTATION
        cudaDeviceSynchronize();
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 5, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 5, 1);
#endif // INSTRUMENTATION       
        E.add(&E, &dE);// E = E + dE;
    }
    return 0;
}

void GPUPF::calcE()
{
    // nothing to do
}

MatrixGPU* GPUPF::calcG()
{
    calcJacPhiE();
   
    tempB2N2.multiplyMat(&JacInv, &I_aug);
    
    tempLN2.multiplyMat(&JacPhiE, &tempB2N2);
    
    G.setBloc(0, B2, 0, N2, &tempB2N2);
    
    G.setBloc(B2, Nconstraint, 0, N2, &tempLN2);

    return &G;
}

MatrixGPU GPUPF::getY()
{
    //std::cout << "getY GPUPF" << std::endl;
    //CHECK_LAST_CUDA_ERROR();
    // E.display(true);
    calcPhi();
    //Phi.display();
    //CHECK_LAST_CUDA_ERROR();
    setY << <numBlock, _blockSize >> > (Y._matrixGPU, E._matrixGPU, Phi._matrixGPU, Nbus, Nline);
    //CHECK_LAST_CUDA_ERROR();
    //Y.display(true);
    return Y;
}

void GPUPF::setE(MatrixGPU* Enew)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    E = *Enew;
    if (!E.getPos()) {
        E.transferGPU();
    }
    if (_useDouble) {
        E.toMatGPUD(ED);
    }
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 8, 1);
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    calcW();
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 6, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 6, 1);
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    if (_useDouble) {
        dWD.subtract(&W0D, &WD);
        err = dWD.max2(); //err = ||dW||
    }
    else {
        dW.subtract(&W0, &W);
        err = dW.max2(); //err = ||dW||
    }
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 7, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 7, 1);
#endif // INSTRUMENTATION


    
    
}

void GPUPF::setE(MatrixGPUD* Enew)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    ED = *Enew;    
    if (!ED.getPos()) {
        ED.transferGPU();
    }
    if (!_useDouble) {
        E = ED;
    }
#ifdef INSTRUMENTATION
    cudaDeviceSynchronize();
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 8, 1);
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    calcW();
#ifdef INSTRUMENTATION
    cudaDeviceSynchronize();
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 6, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 6, 1);
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    if (_useDouble) {
        dWD.subtract(&W0D, &WD);
        err = dWD.max2(); //err = ||dW||
    }
    else {
        dW.subtract(&W0, &W);
        err = dW.max2(); //err = ||dW||
    }
#ifdef INSTRUMENTATION
    cudaDeviceSynchronize();
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 7, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 7, 1);
#endif // INSTRUMENTATION
}

void GPUPF::setW(MatrixGPU* Wnew)
{
    W = *Wnew;
}

float GPUPF::getPloss()
{
    float s = 0;
    
    if (_useDouble) {
        s = WD.sum(0, Nbus);
    }
    else {
        s = W.sum(0, Nbus);
    }
   
    return s; // desequilibre !
}

float GPUPF::getQloss()
{
    float s = 0;
    
    if (_useDouble) {
        s = WD.sum(Nbus, B2);

    }
    else {
        s = W.sum(Nbus, B2);
    }
    
    return s; // desequilibre !
}

float GPUPF::getRes()
{
    return err;
}

int GPUPF::getIter()
{
    return iter;
}

float GPUPF::getP0()
{
    if (_useDouble) {
        return WD.get(0, 0, false);
    }
    else {
        return W.get(0, 0, false);
    }
}

float GPUPF::getQ0()
{
    if (_useDouble) {
        return WD.get(Nbus, 0, false);
    }
    else {
        return W.get(Nbus, 0, false);
    }
}

float GPUPF::getTime()
{
    return (float)time / CLOCKS_PER_SEC;
}

int GPUPF::getConv()
{
    return status;
}

MatrixCPU GPUPF::getE()
{
    if (_useDouble) {
        if(ED.max2()==0){
            calcE();
        }
        E = ED;
    }
    if(E.max2()==0){
        calcE();
    }
    MatrixCPU ECPU;
    E.toMatCPU(ECPU);

    return ECPU;
}
MatrixCPU GPUPF::getW()
{
    if (_useDouble) {
        if (WD.get(0,0)== 0)
        {
           calcW(true);
        }
        W = WD;
    }
    if(W.get(0,0)== 0){
        calcW(true);
    }

    MatrixCPU WCPU;
    W.toMatCPU(WCPU);

    return WCPU;
}



void GPUPF::display()
{
    std::cout << "-----------Resultat du PF --------" << std::endl;

    std::cout << "Nombre d'iter " << iter << " precision atteinte " << err << " temps de resolution " << (float)time / CLOCKS_PER_SEC <<  std::endl;
    //std::cout << " Puissance d'entree" << std::endl;
    //W0.display();
    std::cout << " Puissance active-reactive" << std::endl;
    W.display();
    std::cout << " Tension angle-tension" << std::endl;
    E.display();

    
    std::cout << "Pertes actives  " << getPloss() << std::endl;
    std::cout << "Pertes reactive " << getQloss() << std::endl;
    std::cout << "---------------------------" << std::endl;
}

void GPUPF::display2(bool all)
{
    std::cout.precision(3);

    if (_useDouble) {
        WD.transferCPU();
        W0D.transferCPU();
        dWD.transferCPU();
        ED.transferCPU();
    }
    else {
        W.transferCPU();
        W0.transferCPU();
        dW.transferCPU();
        E.transferCPU();
    }


    if (iter == 0) {
        std::cout << "algorithm not launch" << std::endl;
        if (_useDouble) {
            double temp = WD.get(0, 0);
            double temp2 = WD.get(Nbus, 0);
            WD.set(0, 0, W0.get(0, 0));
            WD.set(Nbus, 0, W0.get(Nbus, 0));
            dWD.subtract(&W0D, &WD); // dW = W0 - W
            err = dWD.max2(); //err = ||dW||
            /*for (int b = 0; b < Nbus; b++) {
                std::cout << b << " " << dW.get(b, 0) << " " << dW.get(Nbus + b, 0) << std::endl;
            }*/
            WD.set(0, 0, temp);
            WD.set(Nbus, 0, temp2);
        }
        else {
            float temp = W.get(0, 0);
            float temp2 = W.get(Nbus, 0);
            W.set(0, 0, W0.get(0, 0));
            W.set(Nbus, 0, W0.get(Nbus, 0));
            dW.subtract(&W0, &W); // dW = W0 - W
            err = dW.max2(); //err = ||dW||
            /*for (int b = 0; b < Nbus; b++) {
                std::cout << b << " " << dW.get(b, 0) << " " << dW.get(Nbus + b, 0) << std::endl;
            }*/
            W.set(0, 0, temp);
            W.set(Nbus, 0, temp2);
        }
    }
    else if (iter < iterM) {
        std::cout << "method " << _name << " on GPU converged in " << iter << " iterations." << std::endl;
        std::cout << "Converged in " << (float)time / CLOCKS_PER_SEC << " seconds" << std::endl;
        if (_useDouble) {
            std::cout << " Computation with double precision" << std::endl;
            double temp = WD.get(0, 0);
            double temp2 = WD.get(Nbus, 0);
            WD.set(0, 0, W0D.get(0, 0));
            WD.set(Nbus, 0, W0D.get(Nbus, 0));
            dWD.subtract(&W0D, &WD); // dW = W0 - W
            err = dWD.max2(); //err = ||dW||
            /*for (int b = 0; b < Nbus; b++) {
                std::cout << b << " " << dW.get(b, 0) << " " << dW.get(Nbus + b, 0) << std::endl;
            }*/
            WD.set(0, 0, temp);
            WD.set(Nbus, 0, temp2);
        }
        else {
            std::cout << " Computation with float simple precision" << std::endl;
            float temp = W.get(0, 0);
            float temp2 = W.get(Nbus, 0);
            W.set(0, 0, W0.get(0, 0));
            W.set(Nbus, 0, W0.get(Nbus, 0));
            dW.subtract(&W0, &W); // dW = W0 - W
            err = dW.max2(); //err = ||dW||
            /*for (int b = 0; b < Nbus; b++) {
                std::cout << b << " " << dW.get(b, 0) << " " << dW.get(Nbus + b, 0) << std::endl;
            }*/
            W.set(0, 0, temp);
            W.set(Nbus, 0, temp2);
        }
    }
    else {
        std::cout << "method " << _name << " on GPU not converged in " << iter << " iterations." << std::endl;
        std::cout << "time taken " << (float)time / CLOCKS_PER_SEC << " seconds" << std::endl;
        if (_useDouble) {
            std::cout << " Computation with double precision" << std::endl;
        }
        else {
            std::cout << " Computation with float simple precision, maibe try with double to converge" << std::endl;
        }
    }
    std::cout << "The power error of this state is " << err << std::endl;
    std::cout << "===============================================================|" << std::endl;
    std::cout << "      System Summary                                           |" << std::endl;
    std::cout << "===============================================================|" << std::endl;
    std::cout << "Buses            " << Nbus << std::endl;
    std::cout << "Branches         " << Nline << std::endl;
    std::cout << "Ploss            " << getPloss() << std::endl;
    std::cout << "Qloss            " << getQloss() << std::endl;


    std::cout << std::endl << std::endl;
    std::cout << "===============================================================================================|" << std::endl;
    std::cout << "      Bus Data                                                                                 |" << std::endl;
    std::cout << "===============================================================================================|" << std::endl;
    std::cout << " Bus |          Voltage        |  Power = Generation  + Load   |  Init = Generation  + Load    |" << std::endl;
    std::cout << "  #  |    Mag(pu) |  Ang(deg)  |    P (pu)     |     Q (pu)    |    P (pu)     |     Q (pu)    |" << std::endl;
    std::cout << "-----|------------|------------|---------------|---------------|---------------|---------------|" << std::endl;
    float seuil = 0.0001;
    //std::cout << 0 << "      " << E.get(Nbus, 0) << "             " << E.get(0, 0) * (abs(E.get(0, 0)) > 0.0001) * 180 / 3.1415 << "              " << (abs(W.get(0, 0)) > 0.0001) * W.get(0, 0) << "         " << (abs(W.get(Nbus, 0)) > 0.0001) * W.get(Nbus, 0) << std::endl;
    //ED.display(true);
    if (all) {
        if (_useDouble) {
            std::cout << std::setw(5) << 0 << "|" << std::setw(11) << ED.get(Nbus, 0) << "*|" << std::setw(11) << ED.get(0, 0) * (abs(ED.get(0, 0)) > seuil) * 180 / 3.1415
                << "*|" << std::setw(15) << (abs(WD.get(0, 0)) > seuil) * WD.get(0, 0) << "|" << std::setw(15) << (abs(WD.get(Nbus, 0)) > seuil) * WD.get(Nbus, 0)
                << "|" << std::setw(15) << W0D.get(0, 0) << "|" << std::setw(15)
                << W0D.get(Nbus, 0) << "|" << std::endl;
            for (int b = 1; b < Nbus; b++) {
                //std::cout.width(10);
                //std::cout << b << "      " << E.get(b + Nbus, 0) << "        " << E.get(b, 0) * (abs(E.get(b, 0)) > 0.0001) * 180 / 3.1415 << "          " << (abs(W.get(b, 0)) > 0.0001) * W.get(b, 0) << "         " << (abs(W.get(b + Nbus, 0)) > 0.0001) * W.get(b + Nbus, 0) << std::endl;
                std::cout << std::setw(5) << b << "|" << std::setw(11) << ED.get(b + Nbus, 0) << " |" << std::setw(11)
                    << ED.get(b, 0) * (abs(ED.get(b, 0)) > seuil) * 180 / 3.1415 << " |" << std::setw(15)
                    << (abs(WD.get(b, 0)) > seuil) * WD.get(b, 0) << "|" << std::setw(15) << (abs(WD.get(b + Nbus, 0)) > seuil) * WD.get(b + Nbus, 0)
                    << "|" << std::setw(15) << (abs(W0D.get(b, 0)) > seuil) * W0D.get(b, 0) << "|" << std::setw(15)
                    << W0D.get(b + Nbus, 0) << "|" << std::endl;

            }
        }
        else {
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
    }
    else {
        if (_useDouble) {
            std::cout << std::setw(5) << 0 << "|" << std::setw(11) << ED.get(Nbus, 0) << "*|" << std::setw(11) << ED.get(0, 0) * (abs(ED.get(0, 0)) > seuil) * 180 / 3.1415
                << "*|" << std::setw(15) << (abs(WD.get(0, 0)) > seuil) * WD.get(0, 0) << "|" << std::setw(15) << (abs(WD.get(Nbus, 0)) > seuil) * WD.get(Nbus, 0)
                << "|" << std::setw(15) << W0D.get(0, 0) << "|" << std::setw(15)
                << W0D.get(Nbus, 0) << "|" << std::endl;
        }
        else {
            std::cout << std::setw(5) << 0 << "|" << std::setw(11) << E.get(Nbus, 0) << "*|" << std::setw(11) << E.get(0, 0) * (abs(E.get(0, 0)) > seuil) * 180 / 3.1415
                << "*|" << std::setw(15) << (abs(W.get(0, 0)) > seuil) * W.get(0, 0) << "|" << std::setw(15) << (abs(W.get(Nbus, 0)) > seuil) * W.get(Nbus, 0)
                << "|" << std::setw(15) << W0.get(0, 0) << "|" << std::setw(15)
                << W0.get(Nbus, 0) << "|" << std::endl;
        }
    }




    std::cout << "===============================================================================================|" << std::endl;
    std::cout << "                      END PRINT                                                                |" << std::endl;
    std::cout << "===============================================================================================|" << std::endl;

}

void GPUPF::saveTimeBlock(std::string fileName)
{
    std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
    float factor = 1000000; // go from ns to ms fot the printed time


    if (occurencePerBlock.get(0, 0) != 0) {
        std::cout << "total resolution time :" << timePerBlock.sum() / (1000 * factor) << "s" << std::endl;
        std::cout << " Fb0 : " << timePerBlock.get(0, 0) / factor << "ms and occurence :" << occurencePerBlock.get(0, 0) << std::endl;
        if (occurencePerBlock.get(0, 3) != 0) {
            std::cout << " Fb1a : " << timePerBlock.get(0, 1) / factor << "ms and occurence :" << occurencePerBlock.get(0, 1) << std::endl;
            std::cout << " Fb1b : " << timePerBlock.get(0, 2) / factor << "ms and occurence :" << occurencePerBlock.get(0, 2) << std::endl;
            std::cout << " Fb1c : " << timePerBlock.get(0, 3) / factor << "ms and occurence :" << occurencePerBlock.get(0, 3) << std::endl;
        }
        else {
            std::cout << " Fb1a : " << timePerBlock.get(0, 1) / factor << "ms and occurence :" << occurencePerBlock.get(0, 1) << std::endl;
            std::cout << " Fb1b : " << timePerBlock.get(0, 2) / factor << "ms and occurence :" << occurencePerBlock.get(0, 2) << std::endl;
        }


        std::cout << " Fb2 : " << timePerBlock.get(0, 4) / factor << "ms and occurence :" << occurencePerBlock.get(0, 4) << std::endl;

        std::cout << " Fb3 : " << timePerBlock.get(0, 5) / factor << "ms and occurence :" << occurencePerBlock.get(0, 5) << std::endl;

        if (occurencePerBlock.get(0, 6) > 0) {
            std::cout << " Fb4 : " << timePerBlock.get(0, 6) / factor << "ms and occurence :" << occurencePerBlock.get(0, 6) << std::endl;

        }
    }
    else {
        std::cout << "pas de temps � afficher, ou alors il n'y a pas eut d'initialisation" << std::endl;
    }

    occurencePerBlock.saveCSV(fileName, mode);
    timePerBlock.saveCSV(fileName, mode);
}


template <unsigned int blockSize>
__global__ void calcW0GPU(float* W0, float* PQ, float* Cores, int N, int B) {
    __shared__ float shArr[blockSize];
    __shared__ float shArr2[blockSize];
    __shared__ bool mustCompute;

    int thIdx = threadIdx.x;
    int i = blockIdx.x;


    if (thIdx == 0) {
        mustCompute = false;
    }
    __syncthreads();
    float sum = 0;
    float sum2 = 0;
   
    for (int k = thIdx; k < N; k += blockSize) {
        if (Cores[k] == i) // c'est tr�s divergent, c'est nul
        {
            sum += PQ[k + 1];
            sum2 += PQ[k + 1 + N];
            mustCompute = true;
        } 
        /* ce n'est plus divergent, mais beaucoup plus d'acc�s m�moire...
         sum += Pinter[k] * (Cores[k] == i);
         sum2 += Qinter[k]* (Cores[k] == i);
        
        */
    }
    if (mustCompute) {
        shArr[thIdx] = sum;
        shArr2[thIdx] = sum2;
        __syncthreads();
        for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
            if (thIdx < size) {
                shArr[thIdx] += shArr[thIdx + size];
                shArr2[thIdx] += shArr2[thIdx + size];
            }
            __syncthreads();
        }

        if (thIdx == 0) {
            W0[i] = shArr[0];
            W0[i + B] = shArr2[0];
        }
    }
    else {
        W0[i] = 0;
        W0[i + B] = 0;
    }
    
}



template <unsigned int blockSize>
__global__ void calcW0GPUD(double* W0D, double* PQD, float* Cores, int N, int B) {
    __shared__ double shArr[blockSize];
    __shared__ double shArr2[blockSize];
    __shared__ bool mustCompute;

    int thIdx = threadIdx.x;
    int i = blockIdx.x;

    if (thIdx == 0) {
        mustCompute = false;
    }
    __syncthreads();

    double sum = 0;
    double sum2 = 0;

    for (int k = thIdx; k < N-1; k += blockSize) {
        if ((int) Cores[k] == i) // c'est tr�s divergent, c'est nul
        {
            sum += PQD[k + 1];
            sum2 += PQD[k + 1 + N];
            mustCompute = true;
        }
        /* ce n'est plus divergent, mais beaucoup plus d'acc�s m�moire...
         sum += PQD[k] * (Cores[k] == i);
         sum2 += PQD[k + N]* (Cores[k] == i);

        */
    }
    if (mustCompute) {
        shArr[thIdx] = sum;
        shArr2[thIdx] = sum2;
        __syncthreads();
        for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
            if (thIdx < size) {
                shArr[thIdx] += shArr[thIdx + size];
                shArr2[thIdx] += shArr2[thIdx + size];
            }
            __syncthreads();
        }

        if (thIdx == 0) {
            W0D[i] = shArr[0];
            W0D[i + B] = shArr2[0];
        }
    }
    else {
        W0D[i] = 0;
        W0D[i + B] = 0;
    }
}


template <unsigned int blockSize>
__global__ void calcW0GPUBis(float* W0, float* PQ, float* Cores, float* nAgentByBus, float* beginBus, int N, int B) {
    __shared__ float shArr[blockSize];
    __shared__ float shArr2[blockSize];
    __shared__ bool mustCompute;

    int thIdx = threadIdx.x;
    int i = blockIdx.x;
    int begin = beginBus[i];
    int end = begin + nAgentByBus[i];

    if (thIdx == 0) {
        mustCompute = nAgentByBus[i] > 0;
    }
    __syncthreads();
    if (mustCompute) { 
        float sum = 0;
        float sum2 = 0;
         for (int k = thIdx + begin; k < end; k += blockSize) {
             int indice = Cores[k];
             sum += PQ[indice];
             sum2 += PQ[indice + N];
         }
            
        shArr[thIdx] = sum;
        shArr2[thIdx] = sum2;
        __syncthreads();
        for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
            if (thIdx < size) {
                shArr[thIdx] += shArr[thIdx + size];
                shArr2[thIdx] += shArr2[thIdx + size];
            }
            __syncthreads();
        }

        if (thIdx == 0) {
            W0[i] = shArr[0];
            W0[i + B] = shArr2[0];
        }
    }
}


template <unsigned int blockSize>
__global__ void calcW0GPUDBis(double* W0D, double* PQD, float* Cores, float* nAgentByBus, float* beginBus, int N, int B) {
    __shared__ double shArr[blockSize];
    __shared__ double shArr2[blockSize];
    __shared__ bool mustCompute;

    int thIdx = threadIdx.x;
    int i = blockIdx.x;
    int begin = beginBus[i];
    int end = begin + nAgentByBus[i];

    if (thIdx == 0) {
        mustCompute = nAgentByBus[i] > 0;
    }
    __syncthreads();
    if (mustCompute) {
        double sum = 0;
        double sum2 = 0;
        for (int k = thIdx + begin; k < end; k += blockSize) {
            int indice = Cores[k];
            sum += PQD[indice];
            sum2 += PQD[indice + N];
        }

        shArr[thIdx] = sum;
        shArr2[thIdx] = sum2;
        __syncthreads();
        for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
            if (thIdx < size) {
                shArr[thIdx] += shArr[thIdx + size];
                shArr2[thIdx] += shArr2[thIdx + size];
            }
            __syncthreads();
        }

        if (thIdx == 0) {
            W0D[i] = shArr[0];
            W0D[i + B] = shArr2[0];
        }
    }
}


__global__ void initE(float* E, float theta0, float V0, int B) {


    int thIdx = threadIdx.x + blockDim.x * blockIdx.x;
    int size = gridDim.x * blockDim.x;
    for (int i = thIdx; i < B; i += size) {
        E[i] = theta0;
        E[i + B] = V0;
    }
}

__global__ void initED(double* ED, double theta0, double V0, int B) {


    int thIdx = threadIdx.x + blockDim.x * blockIdx.x;
    int size = gridDim.x * blockDim.x;
    for (int i = thIdx; i < B; i += size) {
        ED[i] = theta0;
        ED[i + B] = V0;
    }
}



template <unsigned int blockSize>
__global__ void calcWGPU(float* W, float* W0, float* Pinter, float* Qinter, float* CoresBusLin, float* nLines, int B) {
    __shared__ float shArr[blockSize];
    __shared__ float shArr2[blockSize];
    int thIdx = threadIdx.x;
    int i = blockIdx.x; // bus !!!

    int begin = CoresBusLin[i];
    int end = begin + nLines[i];

 
    float sum = 0;
    float sum2 = 0;
    if (i == 0) {
        if (thIdx == 0) {
            W[0] = W0[0];
            W[B] = W0[B];
        }
    }
    else {
        for (int i = begin + thIdx; i < end; i += blockSize) {
            sum += Pinter[i];
            sum2 += Qinter[i];
        }
        

        shArr[thIdx] = sum;
        shArr2[thIdx] = sum2;
        __syncthreads();
        for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
            if (thIdx < size) {
                shArr[thIdx] += shArr[thIdx + size];
                shArr2[thIdx] += shArr2[thIdx + size];
            }
            
            __syncthreads();
        }
    
        if (thIdx == 0) {
            W[i] = shArr[0];
            W[i + B] = shArr2[0];
        }
    }
   
    
}

template <unsigned int blockSize>
__global__ void calcWGPU(float* W, float* Pinter, float* Qinter, float* CoresBusLin, float* nLines, int B) {
    __shared__ float shArr[blockSize];
    __shared__ float shArr2[blockSize];
    int thIdx = threadIdx.x;
    int i = blockIdx.x;
    int begin = CoresBusLin[i];
    int end = begin + nLines[i];


    float sum = 0;
    float sum2 = 0;
    for (int i = begin + thIdx; i < end; i += blockSize) {
        sum += Pinter[i];
        sum2 += Qinter[i];
    }


    shArr[thIdx] = sum;
    shArr2[thIdx] = sum2;
    __syncthreads();
    for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
        if (thIdx < size) {
            shArr[thIdx] += shArr[thIdx + size];
            shArr2[thIdx] += shArr2[thIdx + size];
        }

        __syncthreads();
    }
    if (thIdx == 0) {
        W[i] = shArr[0];
        W[i + B] = shArr2[0];
    }
}

__global__ void calcWinter(float* Pinter,float*Qinter, float* E, float* Glin, float* Blin, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B) {


    int index = threadIdx.x;
    int step = blockDim.x;
    int i = blockIdx.x;
    extern __shared__ float shE[];
    int begin = CoresBusLin[i];
    int end = begin + nLines[blockIdx.x];
    int B2 = 2 * B;

    for (int n = index; n < B2; n += step)
    {
        shE[n] = E[n];
    }
    __syncthreads();

    for (int l = begin + index; l < end; l += step) {
        int k = CoresVoiLin[l];
        float g = Glin[l];
        float b = Blin[l];
        float dt = shE[i] - shE[k];
        float cdt = cos(dt);
        float sdt = sin(dt);
        float v = shE[k + B] * shE[i + B];
       

        Pinter[l] = v * (g * cdt + b * sdt);
        Qinter[l] = v * (g * sdt - b * cdt);

    }

}

__global__ void calcJacGPU(float* Jac, float* W, float* E, float* Glin, float* Blin, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B) {


    int index = threadIdx.x;
    int step = blockDim.x;
    int i = blockIdx.x;
    extern __shared__ float shE[];
    int begin = CoresBusLin[i];
    int end = begin + nLines[i];
    int B2 = 2 * B;

    for (int n = index; n < B2; n += step)
    {
        shE[n] = E[n];
    }
    __syncthreads();
    if (i != 0) {
        for (int l = begin + index; l < end; l += step) {
                int k = CoresVoiLin[l];
                float g = Glin[l];
                float b = Blin[l];
                float dt = shE[i] - shE[k];
                float cdt = cos(dt);
                float sdt = sin(dt);
                float vi = shE[i + B];
                float vk = shE[k + B];
                float p = W[i];
                float q = W[i + B];

                int i2 = i + B;
                int k2 = k + B;
        
      
                //
                Jac[i * B2 + k]   = (-q - b * vi * vi) * (i == k) + ( vi * vk * (g * sdt - b * cdt)) * (i != k);
                Jac[i * B2 + k2]  = (p / vi + g * vi)  * (i == k) + ( vi * (g * cdt + b * sdt)) * (i != k);
                Jac[i2 * B2 + k]  = (p - g * vi * vi)  * (i == k) + (-vi * vk * (g * cdt + b * sdt)) * (i != k);
                Jac[i2 * B2 + k2] = (q / vi - b * vi)  * (i == k) + ( vi * (g * sdt - b * cdt)) * (i != k);

        }
    }
    

}





template <unsigned int blockSize>
__global__ void calcWGPUD(double* W, double* W0, double* Pinter, double* Qinter, float* CoresBusLin, float* nLines, int B) {
    __shared__ double shArr[blockSize];
    __shared__ double shArr2[blockSize];
    int thIdx = threadIdx.x;
    int i = blockIdx.x;

    int begin = CoresBusLin[i];
    int end = begin + nLines[i];


    double sum = 0;
    double sum2 = 0;
    if (i == 0) {
        if (thIdx == 0) {
            W[0] = W0[0];
            W[B] = W0[B];
        }
    }
    else {
        for (int j = begin + thIdx; j < end; j += blockSize) {
            sum += Pinter[j];
            sum2 += Qinter[j];
        }


        shArr[thIdx] = sum;
        shArr2[thIdx] = sum2;
        __syncthreads();
        for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
            if (thIdx < size) {
                shArr[thIdx] += shArr[thIdx + size];
                shArr2[thIdx] += shArr2[thIdx + size];
            }

            __syncthreads();
        }

        if (thIdx == 0) {
            W[i] = shArr[0];
            W[i + B] = shArr2[0];
        }
    }


}

template <unsigned int blockSize>
__global__ void calcWGPUD(double* W, double* Pinter, double* Qinter, float* CoresBusLin, float* nLines, int B) {
    __shared__ double shArr[blockSize];
    __shared__ double shArr2[blockSize];
    int thIdx = threadIdx.x;
    int i = blockIdx.x;
    int begin = CoresBusLin[i];
    int end = begin + nLines[i];


    double sum = 0;
    double sum2 = 0;
    for (int i = begin + thIdx; i < end; i += blockSize) {
        sum += Pinter[i];
        sum2 += Qinter[i];
    }


    shArr[thIdx] = sum;
    shArr2[thIdx] = sum2;
    __syncthreads();
    for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
        if (thIdx < size) {
            shArr[thIdx] += shArr[thIdx + size];
            shArr2[thIdx] += shArr2[thIdx + size];
        }

        __syncthreads();
    }
    if (thIdx == 0) {
        W[i] = shArr[0];
        W[i + B] = shArr2[0];
    }
}

__global__ void calcWinterD(double* Pinter, double* Qinter, double* E, double* Glin, double* Blin, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B) {


    int index = threadIdx.x;
    int step = blockDim.x;
    int i = blockIdx.x;
    extern __shared__ double shED[];
    int begin = CoresBusLin[i];
    int end = begin + nLines[blockIdx.x];
    int B2 = 2 * B;

    for (int n = index; n < B2; n += step)
    {
        shED[n] = E[n];
    }
    __syncthreads();

    for (int l = begin + index; l < end; l += step) {
        int k = CoresVoiLin[l];
        double g = Glin[l];
        double b = Blin[l];
        double dt = shED[i] - shED[k];
        double cdt = cos(dt);
        double sdt = sin(dt);
        double v = shED[k + B] * shED[i + B];


        Pinter[l] = v * (g * cdt + b * sdt);
        Qinter[l] = v * (g * sdt - b * cdt);

    }

}

__global__ void calcJacGPUD(double* Jac, double* W, double* E, double* Glin, double* Blin, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B) {


    int index = threadIdx.x;
    int step = blockDim.x;
    int i = blockIdx.x;
    extern __shared__ double shED[];
    int begin = CoresBusLin[i];
    int end = begin + nLines[i];
    int B2 = 2 * B;
    __shared__ double p;
    __shared__ double q;
    if (index == 0) {
        p = W[i];
        q = W[i + B];
    }

    for (int n = index; n < B2; n += step)
    {
        shED[n] = E[n];
    }
    __syncthreads();
    if (i != 0) {
        for (int l = begin + index; l < end; l += step) {
            int k = CoresVoiLin[l];
            double g = Glin[l];
            double b = Blin[l];
            double dt = shED[i] - shED[k];
            double cdt = cos(dt);
            double sdt = sin(dt);
            double vi = shED[i + B];
            double vk = shED[k + B];


            int i2 = i + B;
            int k2 = k + B;


            //
            Jac[i * B2 + k] = (-q - b * vi * vi) * (i == k) + (vi * vk * (g * sdt - b * cdt)) * (i != k);
            Jac[i * B2 + k2] = (p / vi + g * vi) * (i == k) + (vi * (g * cdt + b * sdt)) * (i != k);
            Jac[i2 * B2 + k] = (p - g * vi * vi) * (i == k) + (-vi * vk * (g * cdt + b * sdt)) * (i != k);
            Jac[i2 * B2 + k2] = (q / vi - b * vi) * (i == k) + (vi * (g * sdt - b * cdt)) * (i != k);

        }
    }


}