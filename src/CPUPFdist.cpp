#include "../head/CPUPFdist.h"


CPUPFdist::CPUPFdist() {}
CPUPFdist::~CPUPFdist() {}

//Fb0: init, Fb1ab 1 2 : Flu, Fb2abc : 3 4 5 Tension, FB3 6 : puissance, Fb4 7 erreur, 8 Fb0 mise � jour

void CPUPFdist::init(const StudyCase& cas, MatrixCPU* PQ, MatrixCPUD * PnD, bool useDouble)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
    timePerBlock = MatrixCPU(1, 9); // Fb0 : init, Fb1ab : Flu, Fb2abc: Tension , FB3 : puissance, Fb4 erreur, Fb0 mise � jour

    occurencePerBlock = MatrixCPU(1, 9);; //nb de fois utilis� pendant la simu
#endif // INSTRUMENTATION
    if(useDouble){
        std::cout << "WARNING : double precision computation is not yet implemented for this method" << std::endl;
    }    


    Nagent = cas.getNagent();
    Nbus = cas.getNBus();
    B2 = 2 * Nbus;
    N2 = 2 * Nagent;
    Nline = cas.getNLine(true); // ne doit pas �tre r�duit ici !!!
    Nconstraint = B2 + Nline;
    iterM = 30;
    iter = 0;
    epsPF = 0.00005;
    status = 0;

    //std::cout << B2 << std::endl;
    
    V0 = cas.getV0();
    theta0 = cas.gettheta0();
    v0 = (float) (V0 * cos(theta0));
    w0 = (float) (V0 * sin(theta0));
    _name = "Current summation method";

    //std::cout << "init " << _name << std::endl;
    W0 = MatrixCPU(B2, 1);
    I = cas.getCoresBusAgentLin(); // I_n = bus de l'agent n, l'agent 0 est pour g�rer les pertes du march�, il n'existe pas vraiment

    for (int n = 1; n < Nagent; n++) {
        int bus = (int) I.get(n, 0);
        W0.increment(bus, 0, PQ->get(n, 0));
        W0.increment(bus + Nbus, 0, PQ->get(n + Nagent, 0));
    }
   
    //std::cout << " W0 : " << std::endl;
    //W0.display();

    W = MatrixCPU(B2, 1);
    dW = MatrixCPU(B2, 1);
    E = MatrixCPU(B2, 1);
    for (int i = 0; i < Nbus; i++) {
        E.set(i, 0, theta0);
    }

    for (int i = Nbus; i < B2; i++) {
        E.set(i, 0, V0);
    }

    dE = MatrixCPU(B2, 1);
  
   
    //std::cout << " E : " << std::endl;
    //E.display();

    // W0[2 * N] : puissance active et r�active au noeud (I*[P Q])
    // W[2 * N] : puissance obtenue par calcul � partir de E
    // dW[2 * N] : derive de puissance
    // E[2 * N] : angle puis tension [O et 1] pour l'init ?
    // dE[2 * N] : derive de angle puis tension
   

   
    CoresVoiLin = cas.getCoresVoiLin();
    CoresBusLin = cas.getCoresBusLin();
    CoresLineBus = cas.getCoresLineBus(true);
    nLines = cas.getNLines();


    
    BgridLin = cas.getBlin();
    GgridLin = cas.getGlin();


    // specificite algo
    //CoresLineBus.display();
    //std::cout << Nbus << " " << Nline << std::endl;
    ZsRe = cas.getZsRe();
    ZsIm = cas.getZsImag();
    Yd = cas.getYd();
    chekcase();

    F = MatrixCPU(Nbus, 1, -1); // F_i = bus ant�c�dent de i
    F.set(1, 0, 0);
    if (Nbus != (Nline + 1)) {
        std::cout << "Warning this is not a distribution network, F not set" << std::endl;
    
    }
    else {
        for (int lold = 0; lold < Nline; lold++) {
            int busTo = (int) CoresLineBus.get(lold, 1);
            F.set(busTo, 0, CoresLineBus.get(lold, 0));
        }
    }
   
    //CoresLineBus.display();
    //F.display();
   
    VoltageRealIm = MatrixCPU(B2, 1);
    VoltageRealImPre = MatrixCPU(B2, 1);
    Jb = MatrixCPU(2 * Nline, 1, -1);
  
    
    /*ZsRe.display();
    ZsIm.display();
    Yd.display();*/
  

    for (int i = 0; i < Nbus; i++) {
        VoltageRealImPre.set(i, 0, v0);
        VoltageRealImPre.set(i + Nbus, 0, w0);
        VoltageRealIm.set(i, 0, v0);
        VoltageRealIm.set(i + Nbus, 0, w0);
    }
   // std::cout << "*****" << std::endl;
    //VoltageRealIm.display();

#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 0, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION
    //std::cout << " fin init" << std::endl;

}

bool CPUPFdist::chekcase()
{
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

    return true;
}

void CPUPFdist::updatePQ(MatrixCPU* PQ)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    W0 = MatrixCPU(B2, 1);
    //PQ->display();
    for (int n = 1; n < Nagent; n++) {
        int bus = (int) I.get(n, 0);
        W0.increment(bus, 0, PQ->get(n, 0));
        W0.increment(bus + Nbus, 0, PQ->get(n + Nagent, 0));
    }
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 8, 1);
#endif
}

void CPUPFdist::solve()
{
    time = clock();
    err = 2 * epsPF;
    iter = 0;
    //std::cout << epsPF << " " << iterM << std::endl;
    status = 1;
    while (err > epsPF && iter < iterM) {

        calcJ();
        //Jb.display();
        //VoltageRealImPre.display();

#ifdef INSTRUMENTATION
        t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION      
        calcVoltage();
#ifdef INSTRUMENTATION
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 3, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 3, 1);
        t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

       
        err = VoltageRealIm.distance2(&VoltageRealImPre);
#ifdef INSTRUMENTATION
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 7, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 7, 1);
#endif // INSTRUMENTATION
        VoltageRealImPre = VoltageRealIm;

        /*if (err < epsPF) {
            calcW();
            dW.subtract(&W0, &W); // dW = W0 - W
            err = dW.max2(); //err = ||dW||
        }*/

        iter++;
        //std::cout << err << " ";
    }
    //std::cout << std::endl;
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    calcW(true);
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 6, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 6, 1);
#endif // INSTRUMENTATION
    /*std::cout << "tension bus entree puis sortie" << std::endl;
    W0.display();
    W.display();*/

    if (iter >= iterM) {
        status = 2;
        if (err > 100 * epsPF) {
            status = -1;
        }
        //std::cout << "fin solve " << iter << " " << err << std::endl;
    }

    time = clock() - time;


}

void CPUPFdist::calcJ()
{
    // step 2
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    //std::cout << "step 2" << std::endl;
    // Set branch current Flow equal to the sum of current of the demand at receiving end s^k_d
    // and the current drawn in the admittance (y^k_d) connected to bus k 
    for (int l = 0; l < Nline; l++) {
        int k = l + 1; // bus to 

        
        //std::cout << "line " << l << " between bus " << F.get(k, 0) << " and " << k << std::endl;
        

        float vRe = VoltageRealImPre.get(k, 0);
        float vIm = VoltageRealImPre.get(k + Nbus, 0);
        
        float p = -W0.get(k, 0);
        
        float q = -W0.get(k + Nbus, 0);
        
        float y = Yd.get(l, 0);
        
        float jRe =  ((p * vRe + q * vIm) / (vRe * vRe + vIm * vIm)) - y * vIm;
        float jIm = -((q * vRe - p * vIm) / (vRe * vRe + vIm * vIm)) + y * vRe;
         
        
        //std::cout << "vRe " << vRe << " vIm " << vIm << " jRe " << jRe << " jIm " << jIm << std::endl; // << " p " << p << " q " << q << " y " << y << std::endl;

        Jb.set(l, 0, jRe);
        
        Jb.set(l + Nline, 0, jIm);
        
    }
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 1, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 1, 1);
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    //std::cout << "step 3" << std::endl;
    // step 3
    //Backward sweep: Perform current summation, starting from the branch with the biggest index and heading towards the branch
    //whose index is equal to 1. The current of branch k is added to the current of the branch whose index is equal to i = f(k)
    for (int l = Nline - 1; l > 0; l--) {
        int k = l + 1; // busTo
        int i = (int) F.get(k, 0); // busFrom

        int lprev = i - 1;
        if (lprev > -1) {
            Jb.increment(lprev, 0, Jb.get(l, 0));
            Jb.increment(lprev + Nline, 0, Jb.get(l + Nline, 0));
        }
        
    }
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 2, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 2, 1);
#endif // INSTRUMENTATION
    //Jb.display();

}

void CPUPFdist::calcW(bool end)
{
    W.set(0.0);
    double p, q;
    
   /* std::cout << "*****" << std::endl;
    GgridLin.display();
    BgridLin.display();
    CoresBusLin.display();
    nLines.display();
    std::cout << "*****" << std::endl;
    ZsRe.display();
    ZsIm.display();

    std::cout << "*****" << std::endl;
    VoltageRealIm.display();*/

    for (int i = 0; i < Nbus; i++) {
        int k = (int) CoresBusLin.get(i, 0);

        for (int voisin = k; voisin < (k + nLines.get(i, 0)); voisin++) {
            int j = (int) CoresVoiLin.get(voisin, 0);

            double a = GgridLin.getD(voisin, 0) * VoltageRealIm.getD(j, 0) - BgridLin.getD(voisin, 0) * VoltageRealIm.getD(j + Nbus, 0);
            double b = BgridLin.getD(voisin, 0) * VoltageRealIm.getD(j, 0) + GgridLin.getD(voisin, 0) * VoltageRealIm.getD(j + Nbus, 0);


            p = VoltageRealIm.getD(i, 0) * a + VoltageRealIm.getD(i + Nbus, 0) * b;
            q = VoltageRealIm.getD(i + Nbus, 0) * a - VoltageRealIm.getD(i, 0) * b;

            W.increment(i, 0, p);
            W.increment(i + Nbus, 0, q);
        }
    }

    if (!end) { // pendant simu, la puissance � ce noeud est libre
        W.set(0, 0, W0.get(0, 0));
        W.set(Nbus, 0, W0.get(Nbus, 0));

    }

}



int CPUPFdist::calcVoltage()
{
    //std::cout << "step 4" << std::endl;
    //Forward sweep: The receiving end bus voltages are calculated with known branch currents and sending bus voltages. 
    for (int k = 1; k < Nbus; k++) {
        // branch l entre le bus i=F(k) et k = l + 1;
        
        int i = (int) F.get(k, 0); // busFrom
        int l = k - 1; // line

       
        float vRe = VoltageRealIm.get(i, 0); 
        float vIm = VoltageRealIm.get(i  + Nbus, 0);
        
        float zRe = ZsRe.get(l, 0);
        float zIm = ZsIm.get(l, 0);
        float jRe = Jb.get(l, 0);
        float jIm = Jb.get(l + Nline, 0);

        float vRe2 = vRe - (zRe * jRe - zIm * jIm);
        float vIm2 = vIm - (zRe * jIm + zIm * jRe);

        VoltageRealIm.set(k, 0, vRe2);
        VoltageRealIm.set(k + Nbus, 0, vIm2);

    }
    //VoltageRealIm.display();
    return 0;
}


void CPUPFdist::calcE()
{
    for (int i = 0; i < Nbus; i++) {
        float Rev = VoltageRealIm.get(i, 0);
        float Imv = VoltageRealIm.get(i + Nbus, 0);
        float Enorm = sqrt(Rev * Rev + Imv * Imv);
        float Eangle = atan2(Imv, Rev);
        E.set(i + Nbus, 0, Enorm);
        E.set(i, 0, Eangle);
    }

}

void CPUPFdist::setE(MatrixCPU* Enew)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    E = *Enew;
    for (int i = 0; i < Nbus; i++) {
        float v = E.get(i + Nbus, 0);
        float theta = E.get(i, 0);
        
        VoltageRealIm.set(i, 0, v * cos(theta));
        VoltageRealIm.set(i + Nbus, 0, v * sin(theta));
        
        VoltageRealImPre.set(i, 0, v * cos(theta));
        VoltageRealImPre.set(i + Nbus, 0, v * sin(theta));

    }
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION
}


void CPUPFdist::setE(MatrixCPUD* Enew)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    E = *Enew;
    for (int i = 0; i < Nbus; i++) {
        float v = E.get(i + Nbus, 0);
        float theta = E.get(i, 0);
        
        VoltageRealIm.set(i, 0, v * cos(theta));
        VoltageRealIm.set(i + Nbus, 0, v * sin(theta));
        
        VoltageRealImPre.set(i, 0, v * cos(theta));
        VoltageRealImPre.set(i + Nbus, 0, v * sin(theta));

    }
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 8, 1);
#endif // INSTRUMENTATION
}


void CPUPFdist::display(bool all)
{
    std::cout.precision(3);
    float errV = (float) err;
    if (iter == 0) {
        std::cout << "algorithm not launch" << std::endl;
        calcW(true);
        if (_useDouble) {
            double temp = WD.get(0, 0);
            double temp2 = WD.get(Nbus, 0);
            WD.set(0, 0, W0.get(0, 0));
            WD.set(Nbus, 0, W0.get(Nbus, 0));
            dWD.subtract(&W0D, &WD); // dW = W0 - W
            err = dWD.max2(); //err = ||dW||
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
            W.set(0, 0, temp);
            W.set(Nbus, 0, temp2);
        }
    }
    else if (iter < iterM) {
        std::cout << "method " << _name << " converged in " << iter << " iterations." << std::endl;
        std::cout << "Converged in " << (float)time / CLOCKS_PER_SEC << " seconds" << std::endl;
        
        std::cout << " Computation with float simple precision" << std::endl;
        float temp = W.get(0, 0);
        if (temp == 0) {
            calcW(true);
            temp = W.get(0, 0);
        }
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

        float seuil = 0.0001f;
       
        
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
        float seuil = 0.0001f;
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