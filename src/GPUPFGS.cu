#include "../head/GPUPFGS.CUh"


GPUPFGS::GPUPFGS(){}
GPUPFGS::~GPUPFGS(){}

void GPUPFGS::init(const StudyCase& cas, MatrixGPU* PQ, MatrixGPUD* PQD, bool useDouble)
{

#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
    timePerBlock = MatrixCPU(1, 9); // Fb0 : init, Fb1ab : Flu, Fb2abc: Tension , FB3 : puissance, Fb4 erreur, Fb0 mise � jour

    occurencePerBlock = MatrixCPU(1, 9);; //nb de fois utilis� pendant la simu
#endif // INSTRUMENTATION
    Nagent = cas.getNagent();
    Nbus = cas.getNBus();
    B2 = 2 * Nbus;
    N2 = 2 * Nagent;
    Nline = cas.getNLine(true); // ne doit pas �tre r�duit ici !!!
    BL2 = Nbus + 2 * Nline;
    Nconstraint = B2 + Nline;
    iterM = 5000;
    iter = 0;
    V0 = cas.getV0();
    theta0 = cas.gettheta0();
    v0 = V0 * cos(theta0);
    w0 = V0 * sin(theta0);
    _name = "Gauss-Seidel";
    _useDouble = useDouble;
    status = 0;
    
    //I = MatrixGPU(cas.getCoresBusAgentLin(), 1);
    CoresAgentBus = MatrixGPU(cas.getCoresAgentBusLin(), 1);
    CoresAgentBusBegin = MatrixGPU(cas.getCoresAgentBusLinBegin(), 1);
    NagentByBus = MatrixGPU(cas.getNagentByBus(), 1);
    //I.display(true);
    numBlock = Nbus;
    _useDouble = useDouble;


    

    CoresLineBus = MatrixGPU(cas.getCoresLineBus(true));
    _CoresVoiLin = MatrixGPU(cas.getCoresVoiLin());
    _CoresBusLin = MatrixGPU(cas.getCoresBusLin());
    _nLines = MatrixGPU(cas.getNLines());
    CoresTrans = MatrixGPU(BL2, 1, 0);

    int* decompte = new int[Nbus];
    for (int i = 0; i < Nbus; i++) {
        decompte[i] = 0;
    }

    
    for (int i = 0; i < Nbus; i++) {
        int begin = _CoresBusLin.get(i,0);
        for (int l = begin + 1; l < (begin + _nLines.get(i, 0)); l++) { // l = (j, i)
            int j = _CoresVoiLin.get(l, 0);
            CoresTrans.set(l, 0, _CoresBusLin.get(j, 0) + 1 + decompte[j]);
            decompte[j]++;
        }
    }
    CoresLineBusGPU = MatrixGPU(2, Nline);

    for (int lold = 0; lold < Nline; lold++) {
        int busTo = CoresLineBus.get(lold, 1);
        int busFrom = CoresLineBus.get(lold, 0);
        CoresLineBusGPU.set(0, lold, busFrom);
        CoresLineBusGPU.set(1, lold, busTo);
    }
    CoresLineBusGPU.transferGPU();
  
    DELETEA(decompte);

    _CoresVoiLin.transferGPU();
    _CoresBusLin.transferGPU();
    _nLines.transferGPU();
    CoresTrans.transferGPU();

    if (useDouble) {
        //BgridD = cas.getLineSuceptanceD();
        //GgridD = cas.getLineReactanceD();
        _BlinD = MatrixGPUD(cas.getBlinD(), 1);
        _GlinD = MatrixGPUD(cas.getGlinD(), 1);

        WD = MatrixGPUD(B2, 1, 0, 1);
        _PintermediateD = MatrixGPUD(BL2, 1, 0, 1);
        _QintermediateD = MatrixGPUD(BL2, 1, 0, 1);

        dWD = MatrixGPUD(B2, 1, 0, 1);
        dWD.preallocateReduction();
        ED = MatrixGPUD(B2, 1, 0, 1);

        initED << <numBlock, _blockSize >> > (ED._matrixGPU, theta0, V0, Nbus);
        // ED.display(true);

        dED = MatrixGPUD(B2, 1, 0, 1);
       
       
        W0D = MatrixGPUD(B2, 1, 0, 1);
        calculW0DBis(PQD);
        //W0D.display(true);

        /*Ggrid2Bgrid2 = MatrixCPU(Nbus, Nbus);
        for (int i = 0; i < Nbus; i++) {
            for (int j = 0; j < Nbus; j++) {
                Ggrid2Bgrid2.set(i, j, sqrt(GgridD.get(i, j) * GgridD.get(i, j) + BgridD.get(i, j) * BgridD.get(i, j)));
            }
        }*/
        
        RgridD = MatrixGPUD(Nbus, 1, 0, 1);
        XgridD = MatrixGPUD(Nbus, 1, 0, 1);
        RMGgridD = MatrixGPUD(BL2, 1, 0, 1);
        RPGgridD = MatrixGPUD(BL2, 1, 0, 1);
        //VectorResultD = MatrixGPUD(B2, 1, 0, 1);
        VoltageRealImD = MatrixGPUD(B2, 1, 0, 1);

        initEDCar <<<numBlock, _blockSize >> > (VoltageRealImD._matrixGPU, theta0, V0, Nbus);
        //VoltageRealImD.display(true);
        

        initRXD << <numBlock, _blockSize >> > (RgridD._matrixGPU, XgridD._matrixGPU, RMGgridD._matrixGPU, RPGgridD._matrixGPU, _GlinD._matrixGPU, _BlinD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU);
        //initRXD2 << <numBlock, _blockSize >> > (RgridD._matrixGPU, XgridD._matrixGPU, RMGgridD._matrixGPU, RPGgridD._matrixGPU, _GlinD._matrixGPU, _BlinD._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, CoresTrans._matrixGPU);


        calcW();
        //WD.display(true);
        dWD.subtract(&W0D, &WD);

        //VoltageRealImD.display(true);

        
       /*_GlinD.display(true);
       _BlinD.display(true);

       RgridD.display(true);
       XgridD.display(true);
       _CoresBusLin.display(true);
       _CoresVoiLin.display(true);
       RMGgridD.display(true);
       RPGgridD.display(true);*/

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

        initE << <numBlock, _blockSize >> > (E._matrixGPU, theta0, V0, Nbus);
        //E.display(true);

        dE = MatrixGPUD(B2, 1, 0, 1);


        W0 = MatrixGPUD(B2, 1, 0, 1);
        calculW0Bis(PQ);
        //W0.display(true);

        /*Ggrid2Bgrid2 = MatrixCPU(Nbus, Nbus);
        for (int i = 0; i < Nbus; i++) {
            for (int j = 0; j < Nbus; j++) {
                Ggrid2Bgrid2.set(i, j, sqrt(GgridD.get(i, j) * GgridD.get(i, j) + BgridD.get(i, j) * BgridD.get(i, j)));
            }
        }*/

        Rgrid = MatrixGPU(Nbus, 1, 0, 1);
        Xgrid = MatrixGPU(Nbus, 1, 0, 1);
        RMGgrid = MatrixGPU(BL2, 1, 0, 1);
        RPGgrid = MatrixGPU(BL2, 1, 0, 1);
       // VectorResult = MatrixGPU(B2, 1, 0, 1);
        VoltageRealIm = MatrixGPU(B2, 1, 0, 1);

        initECar << <numBlock, _blockSize >> > (VoltageRealIm._matrixGPU, theta0, V0, Nbus);
        initRX << <numBlock, _blockSize >> > (Rgrid._matrixGPU, Xgrid._matrixGPU, RMGgrid._matrixGPU, RPGgrid._matrixGPU, _Glin._matrixGPU, _Blin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU);

        
        /*_Glin.display(true);
        _Blin.display(true);

        Rgrid.display(true);
        Xgrid.display(true);
        _CoresBusLin.display(true);
        _CoresVoiLin.display(true);
        RMGgrid.display(true);
        RPGgrid.display(true);*/

        calcW();
        //W.display(true);
        dW.subtract(&W0, &W);


        //VoltageRealIm.display(true);

    }


    _Blin2 = MatrixGPU(cas.getBlin2(), 1);
    _Glin2 = MatrixGPU(cas.getGlin2(), 1);
    /*std::cout << " Bgrid : " << std::endl;
    Bgrid.display();
    std::cout << " Ggrid : " << std::endl;
    Ggrid.display(); */

    //std::cout << " E : " << std::endl;
    //E.display();
    
	// W0[2 * N] : puissance active et r�active au noeud (I*[P Q])
	// W[2 * N] : puissance obtenue par calcul � partir de E
	// dW[2 * N] : derive de puissance
	// E[2 * N] : angle puis tension [O et 1] pour l'init ?
	// dE[2 * N] : derive de angle puis tension
	// Jac[2 * N][2 * N] : jacobienne
	// Jac_inv[2 * N][2 * N]: inverse de la jacobienne
	
	// B[N][N], G[N][N] : caract�ristique des lignes entre les noeuds i et j

    /*G = MatrixCPU(Nconstraint, N2);
    Phi = MatrixCPU(Nline, 1);
    Y = MatrixCPU(Nconstraint, 1);
    tempLN2 = MatrixCPU(Nline, N2);
    JacPhiE = MatrixCPU(Nline, B2);
    tempB2N2 = MatrixCPU(B2, N2);*/
    



#ifdef INSTRUMENTATION
    cudaDeviceSynchronize();
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 0, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION
    //std::cout << numBlock << " " << _blockSize << std::endl;
    //std::cout << " fin init" << std::endl;

}



int GPUPFGS::calcVoltage()
{
    
    //calcE();
    //E.display();

  
   
    if (_useDouble) {  
        //VoltageRealImD.display(true);
#ifdef INSTRUMENTATION
        t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
      
        switch (_blockSize) {
        case 512:
            calculVolDtStep1<512> << <numBlock, _blockSize , B2 * sizeof(double) >> > (VoltageRealImD._matrixGPU, W0D._matrixGPU, RgridD._matrixGPU, XgridD._matrixGPU, RMGgridD._matrixGPU, RPGgridD._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 256:
            calculVolDtStep1<256> << <numBlock, _blockSize, B2 * sizeof(double) >> > (VoltageRealImD._matrixGPU, W0D._matrixGPU, RgridD._matrixGPU, XgridD._matrixGPU, RMGgridD._matrixGPU, RPGgridD._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 128:
            calculVolDtStep1<128> << <numBlock, _blockSize, B2 * sizeof(double) >> > (VoltageRealImD._matrixGPU, W0D._matrixGPU, RgridD._matrixGPU, XgridD._matrixGPU, RMGgridD._matrixGPU, RPGgridD._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 64:
            calculVolDtStep1< 64> << <numBlock, _blockSize, B2 * sizeof(double) >> > (VoltageRealImD._matrixGPU, W0D._matrixGPU, RgridD._matrixGPU, XgridD._matrixGPU, RMGgridD._matrixGPU, RPGgridD._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 32:
            calculVolDtStep1< 32> << <numBlock, _blockSize, B2 * sizeof(double) >> > (VoltageRealImD._matrixGPU, W0D._matrixGPU, RgridD._matrixGPU, XgridD._matrixGPU, RMGgridD._matrixGPU, RPGgridD._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 16:
            calculVolDtStep1< 16> << <numBlock, _blockSize, B2 * sizeof(double) >> > (VoltageRealImD._matrixGPU, W0D._matrixGPU, RgridD._matrixGPU, XgridD._matrixGPU, RMGgridD._matrixGPU, RPGgridD._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case  8:
            calculVolDtStep1<  8> << <numBlock, _blockSize, B2 * sizeof(double) >> > (VoltageRealImD._matrixGPU, W0D._matrixGPU, RgridD._matrixGPU, XgridD._matrixGPU, RMGgridD._matrixGPU, RPGgridD._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case  4:
            calculVolDtStep1<  4> << <numBlock, _blockSize, B2 * sizeof(double) >> > (VoltageRealImD._matrixGPU, W0D._matrixGPU, RgridD._matrixGPU, XgridD._matrixGPU, RMGgridD._matrixGPU, RPGgridD._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case  2:
            calculVolDtStep1<  2> << <numBlock, _blockSize, B2 * sizeof(double) >> > (VoltageRealImD._matrixGPU, W0D._matrixGPU, RgridD._matrixGPU, XgridD._matrixGPU, RMGgridD._matrixGPU, RPGgridD._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case  1:
            calculVolDtStep1<  1> << <numBlock, _blockSize, B2 * sizeof(double) >> > (VoltageRealImD._matrixGPU, W0D._matrixGPU, RgridD._matrixGPU, XgridD._matrixGPU, RMGgridD._matrixGPU, RPGgridD._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        }
        //cudaError_t c = cudaPeekAtLastError();
        //std::cout << cudaGetErrorString(c) << std::endl;
#ifdef INSTRUMENTATION
        cudaDeviceSynchronize();
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 3, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 3, 1);
        t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
        //VoltageRealImD.display(true);
    
        //VoltageRealImD.display(true);
        
        //calculVoltDStep2 <<< 1, _blockSize, B2 * sizeof(double) >> > (VoltageRealImD._matrixGPU, RMGgridD._matrixGPU, RPGgridD._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, CoresTrans._matrixGPU, Nbus);

        calculVoltDStep2bis <<<1, _blockSize, 2 * (BL2 + Nbus) * sizeof(double) >>> (VoltageRealImD._matrixGPU, RMGgridD._matrixGPU, RPGgridD._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, CoresTrans._matrixGPU, Nbus, BL2);
        
        //VoltageRealImD.display(true);
        //std::cout << "----" << std::endl;
#ifdef INSTRUMENTATION
        cudaDeviceSynchronize();
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 4, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 4, 1);
#endif // INSTRUMENTATION

    }
    else {
        //VoltageRealIm.display(true);
#ifdef INSTRUMENTATION
        t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
        switch (_blockSize) {
        case 512:
            calculVoltStep1<512> << <numBlock, _blockSize, B2 * sizeof(float) >> > (VoltageRealIm._matrixGPU, W0._matrixGPU, Rgrid._matrixGPU, Xgrid._matrixGPU, RMGgrid._matrixGPU, RPGgrid._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 256:
            calculVoltStep1<256> << <numBlock, _blockSize, B2 * sizeof(float) >> > (VoltageRealIm._matrixGPU, W0._matrixGPU, Rgrid._matrixGPU, Xgrid._matrixGPU, RMGgrid._matrixGPU, RPGgrid._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 128:
            calculVoltStep1<128> << <numBlock, _blockSize, B2 * sizeof(float) >> > (VoltageRealIm._matrixGPU, W0._matrixGPU, Rgrid._matrixGPU, Xgrid._matrixGPU, RMGgrid._matrixGPU, RPGgrid._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 64:
            calculVoltStep1< 64> << <numBlock, _blockSize, B2 * sizeof(float) >> > (VoltageRealIm._matrixGPU, W0._matrixGPU, Rgrid._matrixGPU, Xgrid._matrixGPU, RMGgrid._matrixGPU, RPGgrid._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 32:
            calculVoltStep1< 32> << <numBlock, _blockSize, B2 * sizeof(float) >> > (VoltageRealIm._matrixGPU, W0._matrixGPU, Rgrid._matrixGPU, Xgrid._matrixGPU, RMGgrid._matrixGPU, RPGgrid._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case 16:
            calculVoltStep1< 16> << <numBlock, _blockSize, B2 * sizeof(float) >> > (VoltageRealIm._matrixGPU, W0._matrixGPU, Rgrid._matrixGPU, Xgrid._matrixGPU, RMGgrid._matrixGPU, RPGgrid._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case  8:
            calculVoltStep1<  8> << <numBlock, _blockSize, B2 * sizeof(float) >> > (VoltageRealIm._matrixGPU, W0._matrixGPU, Rgrid._matrixGPU, Xgrid._matrixGPU, RMGgrid._matrixGPU, RPGgrid._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case  4:
            calculVoltStep1<  4> << <numBlock, _blockSize, B2 * sizeof(float) >> > (VoltageRealIm._matrixGPU, W0._matrixGPU, Rgrid._matrixGPU, Xgrid._matrixGPU, RMGgrid._matrixGPU, RPGgrid._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case  2:
            calculVoltStep1<  2> << <numBlock, _blockSize, B2 * sizeof(float) >> > (VoltageRealIm._matrixGPU, W0._matrixGPU, Rgrid._matrixGPU, Xgrid._matrixGPU, RMGgrid._matrixGPU, RPGgrid._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
        case  1:
            calculVoltStep1<  1> << <numBlock, _blockSize, B2 * sizeof(float) >> > (VoltageRealIm._matrixGPU, W0._matrixGPU, Rgrid._matrixGPU, Xgrid._matrixGPU, RMGgrid._matrixGPU, RPGgrid._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);
            break;
    }
#ifdef INSTRUMENTATION
        cudaDeviceSynchronize();
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 3, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 3, 1);
        t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

       // VoltageRealIm.display(true);

        
        //calculVoltStep2 << <1, _blockSize, B2 * sizeof(float)>>> (VoltageRealIm._matrixGPU, RMGgrid._matrixGPU, RPGgrid._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, CoresTrans._matrixGPU, Nbus);

        calculVoltStep2bis <<<1, _blockSize, 2 * (BL2 + Nbus) * sizeof(float) >> > (VoltageRealIm._matrixGPU, RMGgrid._matrixGPU, RPGgrid._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, CoresTrans._matrixGPU, Nbus, BL2);
        
        //VoltageRealIm.display(true);
        //std::cout << "----" << std::endl;
#ifdef INSTRUMENTATION
        cudaDeviceSynchronize();
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 4, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 4, 1);
#endif // INSTRUMENTATION

    }
    return 0;

}


void GPUPFGS::calcE()
{
    if (_useDouble) {
        calcEGPUD << <numBlock, _blockSize >> > (ED._matrixGPU, VoltageRealImD._matrixGPU, Nbus);
    }
    else {
        calcEGPU << <numBlock, _blockSize >> > (E._matrixGPU, VoltageRealIm._matrixGPU, Nbus);
    }
   

}

void GPUPFGS::calcW(bool end)
{
    if (_useDouble) {

        calcWinterCarD << <numBlock, _blockSize, B2 * sizeof(double) >> > (_PintermediateD._matrixGPU, _QintermediateD._matrixGPU, VoltageRealImD._matrixGPU, _GlinD._matrixGPU, _BlinD._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);

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

        calcWinterCar << <numBlock, _blockSize, B2 * sizeof(float) >> > (_Pintermediate._matrixGPU, _Qintermediate._matrixGPU, VoltageRealIm._matrixGPU, _Glin._matrixGPU, _Blin._matrixGPU, _CoresVoiLin._matrixGPU, _CoresBusLin._matrixGPU, _nLines._matrixGPU, Nbus);

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
           // _Pintermediate.display(true);
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
            //W.display(true);
        }

    }



}

void GPUPFGS::setE(MatrixGPU* Enew)
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
        initEDCar << <numBlock, _blockSize >> > (VoltageRealImD._matrixGPU, ED._matrixGPU, Nbus);
    }
    else {
        initECar << <numBlock, _blockSize >> > (VoltageRealIm._matrixGPU, E._matrixGPU, Nbus);
    }
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 8, 1);
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    calcW();
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
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
    timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 7, 1);
#endif // INSTRUMENTATION

}

void GPUPFGS::setE(MatrixGPUD* Enew)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    ED = *Enew;
    if (!ED.getPos()) {
        ED.transferGPU();
    }
    if (_useDouble) {
        initEDCar << <numBlock, _blockSize >> > (VoltageRealImD._matrixGPU, ED._matrixGPU, Nbus);
    }
    else {
        E = ED;
        initECar << <numBlock, _blockSize >> > (VoltageRealIm._matrixGPU, E._matrixGPU, Nbus);
    }
    
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 8, 1);
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    calcW();
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
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
    timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 7, 1);
#endif // INSTRUMENTATION�
}



__global__ void initEDCar(double* VoltageRealImD, double* ED, int B) {


    int thIdx = threadIdx.x + blockDim.x * blockIdx.x;
    int size = gridDim.x * blockDim.x;
    for (int i = thIdx; i < B; i += size) {
        double V0 = ED[i + B];
        double theta0 = ED[i];

        VoltageRealImD[i] = V0 * cos(theta0);
        VoltageRealImD[i + B] = V0 * sin(theta0);
    }
}
__global__ void initEDCar(double* VoltageRealImD, double v0, double w0, int B) {


    int thIdx = threadIdx.x + blockDim.x * blockIdx.x;
    int size = gridDim.x * blockDim.x;
    for (int i = thIdx; i < B; i += size) {

        VoltageRealImD[i] = v0;
        VoltageRealImD[i + B] = w0;
    }
}



__global__ void initRX(float* Rgrid, float* Xgrid, float* RMGgrid, float* RPGgrid, float* Glin, float* Blin, float* CoresBusLin, float* nLines) {
    int index = threadIdx.x;
    int step = blockDim.x;
    int i = blockIdx.x;
    
    int begin = CoresBusLin[i];
    int end = begin + nLines[i];
   


    float norm = Glin[begin] * Glin[begin] + Blin[begin] * Blin[begin];
    float r =    Glin[begin] / norm;
    float x =   -Blin[begin] / norm;

    for (int l = begin + index + 1; l < end; l += step) {
       
        float m = Glin[l] * r - Blin[l] * x;
        float n = Blin[l] * r + Glin[l] * x;

        RMGgrid[l] = m;
        RPGgrid[l] = n;

    } 
    if (index == 0) {
        Rgrid[i] = r;
        Xgrid[i] = x;
    }
}
__global__ void initRXD(double* RgridD, double* XgridD, double* RMGgridD, double* RPGgridD, double* GlinD, double* BlinD, float* CoresBusLin, float* nLines) {
    int index = threadIdx.x;
    int step = blockDim.x;
    int i = blockIdx.x;

    int begin = CoresBusLin[i];
    int end = begin + nLines[i];


    double norm = (GlinD[begin] * GlinD[begin] + BlinD[begin] * BlinD[begin]);
    double r = GlinD[begin] / norm; // Re(1/Y)
    double x = -BlinD[begin] / norm; // Im(1/Y)

  
    
    for (int l = begin + index + 1; l < end; l += step) {

        float m = GlinD[l] * r - BlinD[l] * x;
        float n = BlinD[l] * r + GlinD[l] * x;

        RMGgridD[l] = m; //Re(Yij/Yii)
        RPGgridD[l] = n; //Im(Yij/Yii)

    }  
    if (index == 0) {
        RgridD[i] = r;
        XgridD[i] = x;
    }
}





/*
 for (int i = 1; i < Nbus; i++) {
            double vi = VoltageRealImD.get(i, 0);
            double wi = VoltageRealImD.get(i + Nbus, 0);
            double norm = vi * vi + wi * wi;
            double c = (W0D.get(i, 0) * vi + W0D.get(i + Nbus, 0) * wi) / norm;
            double d = (W0D.get(i, 0) * wi - W0D.get(i + Nbus, 0) * vi) / norm;
            double sum1 = c * RgridD.get(i, 0) - d * XgridD.get(i, 0);
            double sum2 = d * RgridD.get(i, 0) + c * XgridD.get(i, 0);
            
        }

*/
template <unsigned int blockSize>
__global__ void calculVoltStep1(float* VoltageRealIm, float* W0, float* Rgrid, float* Xgrid, float* RMGgrid, float* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B) {
    
    __shared__ float shArr[blockSize];
    __shared__ float shArr2[blockSize];
    extern __shared__ float shE[];
    int thIdx = threadIdx.x;
    int i = blockIdx.x;
    int step = blockSize;
    int begin = CoresBusLin[i];
    int end = begin + nLines[i];
    int B2 = 2 * B;

    if (i != 0) {
        for (int n = thIdx; n < B2; n += step)
        {
            shE[n] = VoltageRealIm[n];
        }
        __syncthreads();
        float sum = 0;
        float sum2 = 0;
        for (int l = begin + thIdx + 1; l < end; l += step) {
            int k = CoresVoiLin[l];
            if (k > i) {
                sum  -= (RMGgrid[l] * shE[k] - RPGgrid[l] * shE[k + B]);
                sum2 -= (RPGgrid[l] * shE[k] + RMGgrid[l] * shE[k + B]);
            }
        }


        shArr[thIdx] = sum;
        shArr2[thIdx] = sum2;
        __syncthreads();

        if (blockSize >= 512) { 
            if (thIdx < 256) {
                shArr[thIdx] += shArr[thIdx + 256]; 
                shArr2[thIdx] += shArr2[thIdx + 256];
            } 
        __syncthreads();
        }
        if (blockSize >= 256) {
            if (thIdx < 128) {
                shArr[thIdx] += shArr[thIdx + 128]; 
                shArr2[thIdx] += shArr2[thIdx + 128];
            } 
            __syncthreads();
        }
        if (blockSize >= 128) {
            if (thIdx < 64) {
                shArr[thIdx] += shArr[thIdx + 64]; 
                shArr[thIdx] += shArr2[thIdx + 64];
            } __syncthreads(); 
        }
        if (thIdx < 32) {
            warpReduce<blockSize>(shArr, thIdx);
            warpReduce<blockSize>(shArr2, thIdx);
        }
        if (thIdx == 0) {
            float vi = shE[i];
            float wi = shE[i + B];
            float r = Rgrid[i];
            float x = Xgrid[i];
            float W0_local = W0[i];
            float W0B_local = W0[i + B];

            float norm = vi * vi + wi * wi;
            float c = (W0_local * vi + W0B_local * wi) / norm;
            float d = (W0_local * wi - W0B_local * vi) / norm;


       
            VoltageRealIm[i] = shArr[0] + c * r - d * x;
            VoltageRealIm[i + B] = shArr2[0] + d * r + c * x;
        }
    }
}


template <unsigned int blockSize>
__global__ void calculVolDtStep1(double* VoltageRealImD, double* W0, double* Rgrid, double* Xgrid, double* RMGgrid, double* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B) {

    __shared__ double shArr[blockSize];
    __shared__ double shArr2[blockSize];
    extern __shared__ double shED[];
    int thIdx = threadIdx.x;
    int i = blockIdx.x;
    int step = blockSize;
    int begin = CoresBusLin[i];
    int end = begin + nLines[i];
    int B2 = 2 * B;

    if (i != 0) {
        for (int n = thIdx; n < B2; n += step)
        {
            shED[n] = VoltageRealImD[n];
        }
        __syncthreads();
        double sum = 0;
        double sum2 = 0;
        for (int l = begin + thIdx; l < end; l += step) {
            int k = CoresVoiLin[l];
            if (k > i) {
                sum += -(RMGgrid[l] * shED[k] - RPGgrid[l] * shED[k + B]);
                sum2 += -(RPGgrid[l] * shED[k] + RMGgrid[l] * shED[k + B]);
            }
           
        }


        shArr[thIdx] = sum;
        shArr2[thIdx] = sum2;
        __syncthreads();

        if (blockSize >= 512) {
            if (thIdx < 256) {
                shArr[thIdx] += shArr[thIdx + 256];
                shArr2[thIdx] += shArr2[thIdx + 256];
            }
            __syncthreads();
        }
        if (blockSize >= 256) {
            if (thIdx < 128) {
                shArr[thIdx] += shArr[thIdx + 128];
                shArr2[thIdx] += shArr2[thIdx + 128];
            }
            __syncthreads();
        }
        if (blockSize >= 128) {
            if (thIdx < 64) {
                shArr[thIdx] += shArr[thIdx + 64];
                shArr[thIdx] += shArr2[thIdx + 64];
            } __syncthreads();
        }
        if (thIdx < 32) {
            warpReduce<blockSize>(shArr, thIdx);
            warpReduce<blockSize>(shArr2, thIdx);
        }
        if (thIdx == 0) {
            double vi = shED[i];
            double wi = shED[i + B];
            double r = Rgrid[i];
            double x = Xgrid[i];
            double W0_local = W0[i];
            double W0B_local = W0[i + B];

            double norm = vi * vi + wi * wi;
            double c = (W0_local * vi + W0B_local * wi) / norm;
            double d = (W0_local * wi - W0B_local * vi) / norm;



            VoltageRealImD[i] = shArr[0] + c * r - d * x;
            VoltageRealImD[i + B] = shArr2[0] + d * r + c * x;
        }
    }
    
}


/*
for (int iter = 0; iter < Nbus-1; iter++) {
            for (int i = iter + 1; i < Nbus; i++) {
                double db1 = RMGgridD.get(i, iter) * VoltageRealImD.get(iter, 0) - RPGgridD.get(i, iter) * VoltageRealImD.get(iter + Nbus, 0);
                double db2 = RPGgridD.get(i, iter) * VoltageRealImD.get(iter, 0) + RMGgridD.get(i, iter) * VoltageRealImD.get(iter + Nbus, 0);

                VoltageRealImD.increment(i, 0, -db1);
                VoltageRealImD.increment(i + Nbus, 0, -db2);
            }
        }
*
*/


__global__ void calculVoltStep2(float* VoltageRealIm, float* RMGgrid, float* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, float* CoresTrans, int B) {

   /* int thIdx = threadIdx.x;
    int i = blockIdx.x;
    int size = blockDim.x;
    int begin = CoresBusLin[i]; // k = CoresBusLin[iter]; !!!
    int end = begin + nLines[i]; // k + nLines[iter]
    // il ne faut pas Ypj/Ypp mais bien Yjp/Yjj, donc il faut savoir quel voisin est p... 


   if(i>iter){
        for (int voisin = thIdx + begin + 1; voisin < end; voisin += size) {
            int p = CoresVoiLin[voisin];
           
            if (p == iter ) { // pour trouver quel indice est p, c'est plut�t nul
                
                float db1 = RMGgrid[voisin] * VoltageRealIm[iter] - RPGgrid[voisin] * VoltageRealIm[iter + B];
                float db2 = RPGgrid[voisin] * VoltageRealIm[iter] + RMGgrid[voisin] * VoltageRealIm[iter + B];

                VoltageRealIm[i] = VoltageRealIm[i] - db1;
                VoltageRealIm[i + B] = VoltageRealIm[i + B] - db2;
            }
        }
   }*/
    
    int thIdx = threadIdx.x;
    //int i = blockIdx.x; un seul bloc
    int size = blockDim.x;

    extern __shared__ float Voltage[];

    for (int k = thIdx; k < 2 * B; k += size) {
        Voltage[k] = VoltageRealIm[k];
    }
    __syncthreads();


    for (int iter = 0; iter < B - 1; iter++) {
        int begin = CoresBusLin[iter]; // k = CoresBusLin[iter]; !!!
        int end = begin + nLines[iter]; // k + nLines[iter]
        float ei = Voltage[iter];
        float fi = Voltage[iter + B];

        for (int l = thIdx + begin + 1; l < end; l += size) { // voisin
            int j = CoresVoiLin[l];
            
            if (j > iter) {
                int lTrans = CoresTrans[l]; // acc�s pas du tout coalescent !!!

                float ri = RMGgrid[lTrans]; // acc�s pas du tout coalescent mais c'est sur la m�moire partag�
                float li = RPGgrid[lTrans];

                float db1 = ri * ei - li * fi;
                float db2 = li * ei + ri * fi;

                Voltage[j] = Voltage[j] - db1;
                Voltage[j + B] = Voltage[j + B] - db2;
            }
        }
        __syncthreads();
    }

    for (int k = thIdx; k < 2 * B; k += size) {
        VoltageRealIm[k] = Voltage[k];
    }
    __syncthreads();


}


__global__ void calculVoltDStep2(double* VoltageRealIm, double* RMGgrid, double* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, float* CoresTrans, int B) {

    int thIdx = threadIdx.x;
    //int i = blockIdx.x; un seul bloc
    int size = blockDim.x;

    extern __shared__ double VoltageD[];

    for (int k = thIdx; k < 2 * B; k += size) {
        VoltageD[k] = VoltageRealIm[k];
    }
    __syncthreads();


    for (int iter = 0; iter < B - 1; iter++) {
        int begin = CoresBusLin[iter]; // k = CoresBusLin[iter]; !!!
        int end = begin + nLines[iter]; // k + nLines[iter]
        double ei = VoltageD[iter];
        double fi = VoltageD[iter + B];
        __syncthreads();
        for (int l = thIdx + begin + 1; l < end; l += size) { // voisin
            int j = CoresVoiLin[l];

            if (j > iter) {
                int lTrans = CoresTrans[l]; // acc�s coalescent !!!

                double ri = RMGgrid[lTrans]; // acc�s pas du tout coalescent mais c'est sur la m�moire partag�
                double li = RPGgrid[lTrans];

                double db1 = ri * ei - li * fi;
                double db2 = li * ei + ri * fi;

                VoltageD[j] = VoltageD[j] - db1;
                VoltageD[j + B] = VoltageD[j + B] - db2;
            }
        }
        __syncthreads();
    }
    for (int k = thIdx; k < 2 * B; k += size) {
        VoltageRealIm[k] = VoltageD[k];
    }
    __syncthreads();
    
}

__global__ void calculVoltStep2bis(float* VoltageRealIm, float* RMGgrid, float* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, float* CoresTrans, int B, int BL2) {

    
    int thIdx = threadIdx.x;
    //int i = blockIdx.x; un seul bloc
    int size = blockDim.x;
   

    extern __shared__ float RI[];
    
    float* Voltage = &RI[2 * BL2];
    for (int l = thIdx; l < BL2; l += size) { // coalecent
        RI[l] = RMGgrid[l];
        RI[l + BL2] = RPGgrid[l];
    }
    for (int k = thIdx; k < 2 * B; k += size) {
        Voltage[k] = VoltageRealIm[k];
    }
    __syncthreads();


    for (int iter = 0; iter < B - 1; iter++) {
        int begin = CoresBusLin[iter]; // k = CoresBusLin[iter]; !!!
        int end = begin + nLines[iter]; // k + nLines[iter]

        float ei = Voltage[iter];
        float fi = Voltage[iter + B];

        for (int l = thIdx + begin + 1; l < end; l += size) { // voisin
            int j = CoresVoiLin[l];

            if (j > iter) {
                int lTrans = CoresTrans[l];

                float ri = RI[lTrans]; // acc�s pas du tout coalescent mais c'est sur la m�moire partag�
                float li = RI[lTrans + BL2];


                float db1 = ri * ei - li * fi;
                float db2 = li * ei + ri * fi;

                Voltage[j] = Voltage[j] - db1;
                Voltage[j + B] = Voltage[j + B] - db2;
            }
        }
        __syncthreads();
    }
    for (int k = thIdx; k < 2 * B; k += size) {
        VoltageRealIm[k] = Voltage[k];
    }
    __syncthreads();

}


__global__ void calculVoltDStep2bis(double* VoltageRealIm, double* RMGgrid, double* RPGgrid, float* CoresVoiLin, float* CoresBusLin, float* nLines, float* CoresTrans, int B, int BL2) {

    int thIdx = threadIdx.x;
    //int i = blockIdx.x; un seul bloc
    int size = blockDim.x;
    

    extern __shared__ double RID[];
    double* Voltage = &RID[2 * BL2];
    for (int l = thIdx; l < BL2; l += size) { // coalecent
        RID[l] = RMGgrid[l];
        RID[l + BL2] = RPGgrid[l];
    }
    for (int k = thIdx; k < 2 * B; k+=size) {
        Voltage[k] = VoltageRealIm[k];
    }
    __syncthreads();

    for (int iter = 0; iter < B - 1; iter++) {
        int begin = CoresBusLin[iter]; // k = CoresBusLin[iter]; !!!
        int end = begin + nLines[iter]; // k + nLines[iter]

        double ei = Voltage[iter];
        double fi = Voltage[iter + B];

        for (int l = thIdx + begin + 1; l < end; l += size) { // voisin
            int j = CoresVoiLin[l];

            if (j > iter) {
                int lTrans = CoresTrans[l];

                double ri = RID[lTrans]; // acc�s pas du tout coalescent mais c'est sur la m�moire partag�
                double li = RID[lTrans + BL2];


                double db1 = ri * ei - li * fi;
                double db2 = li * ei + ri * fi;

                Voltage[j] = Voltage[j] - db1;
                Voltage[j + B] = Voltage[j + B] - db2;
            }
        }
        __syncthreads();
    }

    for (int k = thIdx; k < 2 * B; k += size) {
        VoltageRealIm[k] = Voltage[k];
    }
    __syncthreads();
}



__global__ void calcWinterCarD(double* Pinter, double* Qinter, double* E, double* Glin, double* Blin, float* CoresVoiLin, float* CoresBusLin, float* nLines, int B) {


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
        double a = g * shED[k] - b * shED[k + B];
        double c = b * shED[k] + g * shED[k + B];

        Pinter[l] = shED[i] * a + shED[i + B] * c;
        Qinter[l] = shED[i + B] * a - shED[i] * c;

    }

}




__global__ void calcEGPUD(double* ED, double* VoltageRealImD, int B) {

    int thIdx = threadIdx.x + blockDim.x * blockIdx.x;
    int size = gridDim.x * blockDim.x;

    for (int i = thIdx; i < B; i += size) {
        double Rev = VoltageRealImD[i];
        double Imv = VoltageRealImD[i + B];


        ED[i + B] = sqrt(Rev * Rev + Imv * Imv);
        ED[i] = atan2(Imv, Rev);
    }

}
