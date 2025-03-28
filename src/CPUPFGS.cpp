#include "../head/CPUPFGS.h"


CPUPFGS::CPUPFGS(){}
CPUPFGS::~CPUPFGS(){}

void CPUPFGS::init(const StudyCase& cas, MatrixCPU* PQ)
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
    Nconstraint = B2 + Nline;
    iterM = 5000; 
    iter = 0;
    V0 = cas.getV0();
    theta0 = cas.gettheta0();
    v0 = V0 * cos(theta0);
    w0 = V0 * sin(theta0);
    _name = "Gauss-Seidel";
    _useDouble = false;
    status = 0;
    //std::cout << "init " << _name << std::endl;
    I = cas.getCoresBusAgentLin();
   
   
   
    //std::cout << " W0 : " << std::endl;
    //W0.display();
    CoresLineBus = cas.getCoresLineBus(true);
    CoresVoiLin = cas.getCoresVoiLin();
    CoresBusLin = cas.getCoresBusLin();
    nLines = cas.getNLines();

 
    W = MatrixCPU(B2, 1);
    dW = MatrixCPU(B2, 1);
    E = MatrixCPU(B2, 1);
    dE = MatrixCPU(B2, 1);
    Jac = MatrixCPU(B2, B2);
    JacInv = MatrixCPU(B2, B2);
    Bgrid = cas.getLineSuceptance();
    Ggrid = cas.getLineReactance();
    BgridLin = cas.getBlin();
    GgridLin = cas.getGlin();

    for (int i = 0; i < Nbus; i++) {
        E.set(i, 0, theta0);
    }

    for (int i = Nbus; i < B2; i++) {
        E.set(i, 0, V0);
    }
    Ggrid2Bgrid2 = MatrixCPU(Nbus, Nbus);
    for (int i = 0; i < Nbus; i++) {
        for (int j = 0; j < Nbus; j++) {
            Ggrid2Bgrid2.set(i, j, sqrt(Ggrid.get(i, j) * Ggrid.get(i, j) + Bgrid.get(i, j) * Bgrid.get(i, j)));
        }
    }


    Rgrid = MatrixCPU(Nbus, 1);
    Xgrid = MatrixCPU(Nbus, 1);
    RMGgrid = MatrixCPU(Nbus, Nbus);
    RPGgrid = MatrixCPU(Nbus, Nbus);
    VectorResult = MatrixCPU(B2, 1);
    VoltageRealIm = MatrixCPU(B2, 1);

    for (int i = 0; i < Nbus; i++) {
        float r =  Ggrid.get(i, i) / (Ggrid.get(i, i) * Ggrid.get(i, i) + Bgrid.get(i, i) * Bgrid.get(i, i));
        float x = -Bgrid.get(i, i) / (Ggrid.get(i, i) * Ggrid.get(i, i) + Bgrid.get(i, i) * Bgrid.get(i, i));
        Rgrid.set(i, 0, r);
        Xgrid.set(i, 0, x);
        for (int j = 0; j < Nbus; j++) {
            float m = Ggrid.get(i, j) * r - Bgrid.get(i, j) * x;
            float n = Bgrid.get(i, j) * r + Ggrid.get(i, j) * x;
            RMGgrid.set(i, j, m);
            RPGgrid.set(i, j, n);
        }

        VoltageRealIm.set(i, 0, v0);
        VoltageRealIm.set(i + Nbus, 0, w0);
    }
    W0 = MatrixCPU(B2, 1);
    for (int n = 1; n < Nagent; n++) {
        int bus = I.get(n, 0);
        W0.increment(bus, 0, PQ->get(n, 0));
        W0.increment(bus + Nbus, 0, PQ->get(n + Nagent, 0));
    }
    calcW();
    //std::cout << " W : " << std::endl;
    //W.display();
    dW.subtract(&W0, &W);

    /*Ggrid.display();
    Bgrid.display();

        
    GgridLin.display();
    BgridLin.display();

    Rgrid.display();
    Xgrid.display();
    CoresBusLin.display();
    CoresVoiLin.display();
    RMGgrid.display();
    RPGgrid.display();*/
    



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

    G = MatrixCPU(Nconstraint, N2);
    Phi = MatrixCPU(Nline, 1);
    Y = MatrixCPU(Nconstraint, 1);
    tempLN2 = MatrixCPU(Nline, N2);
    JacPhiE = MatrixCPU(Nline, B2);
    tempB2N2 = MatrixCPU(B2, N2);
    



#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 0, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION

    //std::cout << " fin init" << std::endl;

}

void CPUPFGS::init(const StudyCase& cas, MatrixCPU* PQ, MatrixCPUD* PQD, bool useDouble)
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
    std::cout << "init " << _name << std::endl;
    I = cas.getCoresBusAgentLin();



    //std::cout << " W0 : " << std::endl;
    //W0.display();
    CoresLineBus = cas.getCoresLineBus(true);
    CoresVoiLin = cas.getCoresVoiLin();
    CoresBusLin = cas.getCoresBusLin();
  
    nLines = cas.getNLines();

    if (useDouble) {
        BgridD = cas.getLineSuceptanceD();
        GgridD = cas.getLineReactanceD();
        BgridLinD = cas.getBlinD();
        GgridLinD = cas.getGlinD();

        if (BgridD.getNCol() == 0) {
            Bgrid = cas.getLineSuceptance();
            Bgrid.toMatCPUD(BgridD);
        }
        if (GgridD.getNCol() == 0) {
            Ggrid = cas.getLineReactance();
            Ggrid.toMatCPUD(GgridD);
        }

        ED = MatrixCPUD(B2, 1);
        dED = MatrixCPUD(B2, 1);
        JacD = MatrixCPUD(B2, B2);
        JacInvD = MatrixCPUD(B2, B2);

        for (int i = 0; i < Nbus; i++) {
            ED.set(i, 0, theta0);
        }

        for (int i = Nbus; i < B2; i++) {
            ED.set(i, 0, V0);
        }

        W0D = MatrixCPUD(B2, 1);
        W0.toMatCPUD(W0D);
        WD = MatrixCPUD(B2, 1);
        dWD = MatrixCPUD(B2, 1);

        Ggrid2Bgrid2 = MatrixCPU(Nbus, Nbus);
        for (int i = 0; i < Nbus; i++) {
            for (int j = 0; j < Nbus; j++) {
                Ggrid2Bgrid2.set(i, j, sqrt(GgridD.get(i, j) * GgridD.get(i, j) + BgridD.get(i, j) * BgridD.get(i, j)));
            }
        }


        RgridD = MatrixCPUD(Nbus, 1);
        XgridD = MatrixCPUD(Nbus, 1);
        RMGgridD = MatrixCPUD(Nbus, Nbus);
        RPGgridD = MatrixCPUD(Nbus, Nbus);
        VectorResultD = MatrixCPUD(B2, 1);
        VoltageRealImD = MatrixCPUD(B2, 1);

        for (int i = 0; i < Nbus; i++) {
            float norm = (GgridD.get(i, i) * GgridD.get(i, i) + BgridD.get(i, i) * BgridD.get(i, i));
            float r = GgridD.get(i, i) / norm;
            float x = -BgridD.get(i, i) / norm;
            RgridD.set(i, 0, r);
            XgridD.set(i, 0, x);
            for (int j = 0; j < Nbus; j++) {
                float m = GgridD.get(i, j) * r - BgridD.get(i, j) * x;
                float n = BgridD.get(i, j) * r + GgridD.get(i, j) * x;
                RMGgridD.set(i, j, m);
                RPGgridD.set(i, j, n);
            }

            VoltageRealImD.set(i, 0, v0);
            VoltageRealImD.set(i + Nbus, 0, w0);
        }
 
        W0D = MatrixCPUD(B2, 1);
        for (int n = 1; n < Nagent; n++) {
            int bus = I.get(n, 0);
            W0D.increment(bus, 0, PQD->get(n, 0));
            W0D.increment(bus + Nbus, 0, PQD->get(n + Nagent, 0));
        }

        calcW();

        dWD.subtract(&W0D, &WD);

    }
    else {
        W = MatrixCPU(B2, 1);
        dW = MatrixCPU(B2, 1);
        E = MatrixCPU(B2, 1);
        dE = MatrixCPU(B2, 1);
        Jac = MatrixCPU(B2, B2);
        JacInv = MatrixCPU(B2, B2);
        Bgrid = cas.getLineSuceptance();
        Ggrid = cas.getLineReactance();
        BgridLin = cas.getBlin();
        GgridLin = cas.getGlin();

        for (int i = 0; i < Nbus; i++) {
            E.set(i, 0, theta0);
        }

        for (int i = Nbus; i < B2; i++) {
            E.set(i, 0, V0);
        }
        Ggrid2Bgrid2 = MatrixCPU(Nbus, Nbus);
        for (int i = 0; i < Nbus; i++) {
            for (int j = 0; j < Nbus; j++) {
                Ggrid2Bgrid2.set(i, j, sqrt(Ggrid.get(i, j) * Ggrid.get(i, j) + Bgrid.get(i, j) * Bgrid.get(i, j)));
            }
        }


        Rgrid = MatrixCPU(Nbus, 1);
        Xgrid = MatrixCPU(Nbus, 1);
        RMGgrid = MatrixCPU(Nbus, Nbus);
        RPGgrid = MatrixCPU(Nbus, Nbus);
        VectorResult = MatrixCPU(B2, 1);
        VoltageRealIm = MatrixCPU(B2, 1);

        for (int i = 0; i < Nbus; i++) {
            float r = Ggrid.get(i, i) / (Ggrid.get(i, i) * Ggrid.get(i, i) + Bgrid.get(i, i) * Bgrid.get(i, i));
            float x = -Bgrid.get(i, i) / (Ggrid.get(i, i) * Ggrid.get(i, i) + Bgrid.get(i, i) * Bgrid.get(i, i));
            Rgrid.set(i, 0, r);
            Xgrid.set(i, 0, x);
            for (int j = 0; j < Nbus; j++) {
                float m = Ggrid.get(i, j) * r - Bgrid.get(i, j) * x;
                float n = Bgrid.get(i, j) * r + Ggrid.get(i, j) * x;
                RMGgrid.set(i, j, m);
                RPGgrid.set(i, j, n);
            }

            VoltageRealIm.set(i, 0, v0);
            VoltageRealIm.set(i + Nbus, 0, w0);
        }
        W0 = MatrixCPU(B2, 1);
        for (int n = 1; n < Nagent; n++) {
            int bus = I.get(n, 0);
            W0.increment(bus, 0, PQ->get(n, 0));
            W0.increment(bus + Nbus, 0, PQ->get(n + Nagent, 0));
        }
        calcW();
        //std::cout << " W : " << std::endl;
        //W.display();
        dW.subtract(&W0, &W);

        /*Ggrid.display();
        Bgrid.display();


        GgridLin.display();
        BgridLin.display();

        Rgrid.display();
        Xgrid.display();
        CoresBusLin.display();
        CoresVoiLin.display();
        RMGgrid.display();
        RPGgrid.display();*/
    }



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

    G = MatrixCPU(Nconstraint, N2);
    Phi = MatrixCPU(Nline, 1);
    Y = MatrixCPU(Nconstraint, 1);
    tempLN2 = MatrixCPU(Nline, N2);
    JacPhiE = MatrixCPU(Nline, B2);
    tempB2N2 = MatrixCPU(B2, N2);




#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 0, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION

    //std::cout << " fin init" << std::endl;

}

void CPUPFGS::updatePQ(MatrixCPU* PQ)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    W0 = MatrixCPU(B2, 1);
    //PQ->display();
    for (int n = 1; n < Nagent; n++) {
        int bus = I.get(n, 0);
        W0.increment(bus, 0, PQ->get(n, 0));
        W0.increment(bus + Nbus, 0, PQ->get(n + Nagent, 0));
    }
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 8, 1);
#endif
}

void CPUPFGS::calcW(bool end)
{
    double p = 0;
    double q = 0;

    if (_useDouble) {

        WD.set(0.0);
        for (int i = 0; i < Nbus; i++) {
            int k = CoresBusLin.get(i, 0);
            for (int voisin = k; voisin < (k + nLines.get(i, 0)); voisin++) {
                int j = CoresVoiLin.get(voisin, 0);
                double a = GgridLinD.get(voisin, 0) * VoltageRealImD.get(j, 0) - BgridLinD.get(voisin, 0) * VoltageRealImD.get(j + Nbus, 0);
                double b = BgridLinD.get(voisin, 0) * VoltageRealImD.get(j, 0) + GgridLinD.get(voisin, 0) * VoltageRealImD.get(j + Nbus, 0);
               
                p = VoltageRealImD.get(i, 0) * a + VoltageRealImD.get(i + Nbus, 0) * b;
                q = VoltageRealImD.get(i + Nbus, 0) * a - VoltageRealImD.get(i, 0) * b;

                WD.increment(i, 0, p);
                WD.increment(i + Nbus, 0, q);
            }
        }
        /**/
       /*double sp = 0;
       double sq = 0;

       for (int i = 0; i < Nbus; i++) {
           for (int j = 0; j < Nbus; j++) {
               double dTheta_ij = ED.get(i, 0) - ED.get(j, 0);
               p = ED.get(Nbus + j, 0) * (GgridD.get(i, j) * cos(dTheta_ij) + BgridD.get(i, j) * sin(dTheta_ij));
               q = ED.get(Nbus + j, 0) * (GgridD.get(i, j) * sin(dTheta_ij) - BgridD.get(i, j) * cos(dTheta_ij));
               sp = sp + p;
               sq = sq + q;
               if (i == 0 && p!=0) {
                   std::cout << p << std::endl;
               }
           }
           WD.set(i, 0, ED.get(Nbus + i, 0) * sp);
           WD.set(i + Nbus, 0, ED.get(Nbus + i, 0) * sq);
           sp = 0;
           sq = 0;


       }*/
       /*double test = 0;
       for (int i = 0; i < Nbus; i++) {

           Gcos[i] = ED.get(Nbus + i, 0) * GgridD.get(i, i);
           Bcos[i] = ED.get(Nbus + i, 0) * BgridD.get(i, i);
           Gsin[i] = 0;
           Bsin[i] = 0;

       }


       for (int l = 0; l < Nline; l++) {
           int i = CoresLineBus.get(l, 0);
           int j = CoresLineBus.get(l, 1);
           // (i,j)
           double dTheta_ij = ED.get(i, 0) - ED.get(j, 0);

           double cdt = cos(dTheta_ij);
           double sdt = sin(dTheta_ij);

           double Gc = GgridD.get(i, j) * cdt;
           double Bc = BgridD.get(i, j) * cdt;
           double Gs = GgridD.get(i, j) * sdt;
           double Bs = BgridD.get(i, j) * sdt;

           Gcos[i] += ED.get(Nbus + j, 0) * Gc;
           Bcos[i] += ED.get(Nbus + j, 0) * Bc;
           Gsin[i] += ED.get(Nbus + j, 0) * Gs;
           Bsin[i] += ED.get(Nbus + j, 0) * Bs;



           //(j,i)
           Gc = GgridD.get(j, i) * cdt;
           Bc = BgridD.get(j, i) * cdt;
           Gs = GgridD.get(j, i) * sdt;
           Bs = BgridD.get(j, i) * sdt;


           p =  ED.get(Nbus + i, 0) * (Gc - Bs);
           q = -ED.get(Nbus + i, 0) * (Gs + Bc);

           Gcos[j] +=  ED.get(Nbus + i, 0) * Gc;
           Bcos[j] +=  ED.get(Nbus + i, 0) * Bc;
           Gsin[j] += -ED.get(Nbus + i, 0) * Gs;
           Bsin[j] += -ED.get(Nbus + i, 0) * Bs;

       }

       for (int i = 0; i < Nbus; i++) {
           WD.set(i, 0, ED.get(Nbus + i, 0) * (Gcos[i] + Bsin[i]));
           WD.set(i + Nbus, 0, ED.get(Nbus + i, 0) * (Gsin[i] - Bcos[i]));
       }*/
        if (!end) { // pendant simu, la puissance � ce noeud est libre
            WD.set(0, 0, W0D.get(0, 0));
            WD.set(Nbus, 0, W0D.get(Nbus, 0));

        }

    }
    else {
       

        W.set(0.0);
        for (int i = 0; i < Nbus; i++) {
            int k = CoresBusLin.get(i, 0);

            for (int voisin = k; voisin < (k + nLines.get(i, 0)); voisin++) {
                int j = CoresVoiLin.get(voisin, 0);
               
                double a = GgridLin.getD(voisin, 0) * VoltageRealIm.get(j, 0) - BgridLin.getD(voisin, 0) * VoltageRealIm.get(j + Nbus, 0);
                double b = BgridLin.getD(voisin, 0) * VoltageRealIm.get(j, 0) + GgridLin.getD(voisin, 0) * VoltageRealIm.get(j + Nbus, 0);


                p = VoltageRealIm.get(i, 0) * a + VoltageRealIm.get(i + Nbus, 0) * b;
                q = VoltageRealIm.get(i + Nbus, 0) * a - VoltageRealIm.get(i, 0) * b;
               
                W.increment(i, 0, p);
                W.increment(i + Nbus, 0, q);
            }
        }
        /*double Gsin = 0;
        double Bsin = 0;
        double Gcos = 0;
        double Bcos = 0;
        double EE = 0;

        for (int i = 0; i < Nbus; i++) {

            p = E.get(Nbus + i, 0) * Ggrid.get(i, i);
            q = - E.get(Nbus + i, 0) * Bgrid.get(i, i);

            W.set(i, 0, E.getD(Nbus + i, 0) * p);
            W.set(i + Nbus, 0, E.getD(Nbus + i, 0) * q);

        }


        for (int l = 0; l < Nline; l++) {
            int i = CoresLineBus.get(l, 0);
            int j = CoresLineBus.get(l, 1);
            double dTheta_ij = E.get(i, 0) - E.get(j, 0); // pour j c'est -dTheta_ij

            Gsin = Ggrid.getD(i, j) * sin(dTheta_ij);
            Bsin = Bgrid.getD(i, j) * sin(dTheta_ij);
            Gcos = Ggrid.getD(i, j) * cos(dTheta_ij);
            Bcos = Bgrid.getD(i, j) * cos(dTheta_ij);
            EE = E.getD(Nbus + i, 0) * E.getD(Nbus + j, 0);
            p = EE * (Gcos + Bsin );
            q = EE * (Gsin - Bcos);

            W.increment(i, 0, p);
            W.increment(i + Nbus, 0, q);

            p = EE * (Gcos - Bsin);
            q = -EE * (Gsin + Bcos);

            W.increment(j, 0, p);
            W.increment(j + Nbus, 0, q);

        }*/
        /*
        double sp = 0;
        double sq = 0;
        for (int i = 0; i < Nbus; i++) {
            for (int j = 0; j < Nbus; j++) {
                double dTheta_ij = E.get(i,0) - E.get(j,0);
                p =  E.getD(Nbus + j, 0) * (Ggrid.getD(i, j) * cos(dTheta_ij) + Bgrid.getD(i, j) * sin(dTheta_ij));
                q =  E.getD(Nbus + j, 0) * (Ggrid.getD(i, j) * sin(dTheta_ij) - Bgrid.getD(i, j) * cos(dTheta_ij));
                sp = sp + p;
                sq = sq + q;
            }
            W.set(i, 0, E.getD(Nbus + i, 0) * sp);
            W.set(i + Nbus, 0, E.getD(Nbus + i, 0) *sq);
            sp = 0;
            sq = 0;

        }*/
        if (!end) { // pendant simu, la puissance � ce noeud est libre
            W.set(0, 0, W0.get(0, 0));
            W.set(Nbus, 0, W0.get(Nbus, 0));

        }
    }
}



int CPUPFGS::calcVoltage()
{
    //VoltageRealIm.display();
    //calcE();
    //E.display();
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
  

    if (_useDouble) {
        for (int i = 1; i < Nbus; i++) {
            double vi = VoltageRealImD.get(i, 0);
            double wi = VoltageRealImD.get(i + Nbus, 0);
            double norm = vi * vi + wi * wi;
            double c = (W0D.get(i, 0) * vi + W0D.get(i + Nbus, 0) * wi) / norm;
            double d = (W0D.get(i, 0) * wi - W0D.get(i + Nbus, 0) * vi) / norm;
            double sum1 = c * RgridD.get(i, 0) - d * XgridD.get(i, 0);
            double sum2 = d * RgridD.get(i, 0) + c * XgridD.get(i, 0);
            int k = CoresBusLin.get(i, 0);
            for (int voisin = k + 1; voisin < (k + nLines.get(i, 0)); voisin++) { //for (int j = i + 1; j < Nbus; j++) {
                int j = CoresVoiLin.get(voisin, 0);
                if (i < j) {
                    sum1 = sum1 - (RMGgridD.get(i, j) * VoltageRealImD.get(j, 0) - RPGgridD.get(i, j) * VoltageRealImD.get(j + Nbus, 0));
                    sum2 = sum2 - (RPGgridD.get(i, j) * VoltageRealImD.get(j, 0) + RMGgridD.get(i, j) * VoltageRealImD.get(j + Nbus, 0));
                }
           
            }
            VoltageRealImD.set(i, 0, sum1);
            VoltageRealImD.set(i + Nbus, 0, sum2);
        }
        
#ifdef INSTRUMENTATION
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 3, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 3, 1);
        t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

        //VoltageRealImD.display();

        for (int iter = 0; iter < Nbus-1; iter++) {
            int k = CoresBusLin.get(iter, 0);
            for (int voisin = k + 1; voisin < (k + nLines.get(iter, 0)); voisin++) {
                int j = CoresVoiLin.get(voisin, 0);
                if (j > iter) {
                    double db1 = RMGgridD.get(j, iter) * VoltageRealImD.get(iter, 0) - RPGgridD.get(j, iter) * VoltageRealImD.get(iter + Nbus, 0);
                    double db2 = RPGgridD.get(j, iter) * VoltageRealImD.get(iter, 0) + RMGgridD.get(j, iter) * VoltageRealImD.get(iter + Nbus, 0);

                    VoltageRealImD.increment(j, 0, -db1);
                    VoltageRealImD.increment(j + Nbus, 0, -db2);
                }

            }
        }
        //VoltageRealImD.display();
#ifdef INSTRUMENTATION
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 4, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 4, 1);
#endif // INSTRUMENTATION

    }
    else {
        //VoltageRealIm.set(0, 0, v0);
        //VoltageRealIm.set(Nbus, 0, w0);
        for (int i = 1; i < Nbus; i++) {
            double vi = VoltageRealIm.get(i, 0);
            double wi = VoltageRealIm.get(i + Nbus, 0);
            double norm = vi * vi + wi * wi;
            double c = (W0.get(i, 0) * vi + W0.get(i + Nbus, 0) * wi) / norm;
            double d = (W0.get(i, 0) * wi - W0.get(i + Nbus, 0) * vi) / norm;
            double sum1 = c * Rgrid.get(i, 0) - d * Xgrid.get(i, 0);
            double sum2 = d * Rgrid.get(i, 0) + c * Xgrid.get(i, 0);
            int k = CoresBusLin.get(i, 0);
            for (int voisin = k + 1; voisin < (k + nLines.get(i, 0)); voisin++) { //for (int j = i + 1; j < Nbus; j++) {
                int j = CoresVoiLin.get(voisin, 0);
                if (i < j) {
                    sum1 = sum1 - (RMGgrid.get(i, j) * VoltageRealIm.get(j, 0) - RPGgrid.get(i, j) * VoltageRealIm.get(j + Nbus, 0));
                    sum2 = sum2 - (RPGgrid.get(i, j) * VoltageRealIm.get(j, 0) + RMGgrid.get(i, j) * VoltageRealIm.get(j + Nbus, 0));
                }
            }
            VoltageRealIm.set(i, 0, sum1);
            VoltageRealIm.set(i + Nbus, 0, sum2);
        }
#ifdef INSTRUMENTATION
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 3, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 3, 1);
        t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

        //VoltageRealIm.display();

        for (int iter = 0; iter < Nbus - 1; iter++) {
            int k = CoresBusLin.get(iter, 0);
            for (int voisin = k + 1; voisin < (k + nLines.get(iter, 0)); voisin++) {
                int j = CoresVoiLin.get(voisin, 0);
                if (j > iter) {
                    double db1 = RMGgrid.get(j, iter) * VoltageRealIm.get(iter, 0) - RPGgrid.get(j, iter) * VoltageRealIm.get(iter + Nbus, 0);
                    double db2 = RPGgrid.get(j, iter) * VoltageRealIm.get(iter, 0) + RMGgrid.get(j, iter) * VoltageRealIm.get(iter + Nbus, 0);

                    VoltageRealIm.increment(j, 0, -db1);
                    VoltageRealIm.increment(j + Nbus, 0, -db2);
                }
            }
        }
        //VoltageRealIm.display();
#ifdef INSTRUMENTATION
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 4, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 4, 1);
#endif // INSTRUMENTATION

    }
    return 0;


    /**/
       /*for (int i = 1; i < Nbus; i++) {
        float vi = VoltageRealIm.get(i, 0);
        float wi = VoltageRealIm.get(i + Nbus, 0);
        float norm = vi * vi + wi * wi;
        float c = (W0.get(i, 0) * vi + W0.get(i + Nbus, 0) * wi) / norm;
        float d = (W0.get(i, 0) * wi - W0.get(i + Nbus, 0) * vi) / norm;
        float sum1 = c * Rgrid.get(i, 0) - d * Xgrid.get(i, 0);
        float sum2 = d * Rgrid.get(i, 0) + c * Xgrid.get(i, 0);
        for (int j = 0; j < Nbus; j++) {
            if (j != i) {
                sum1 -= RMGgrid.get(i, j) * VoltageRealIm.get(j, 0) - RPGgrid.get(i, j) * VoltageRealIm.get(j + Nbus, 0);
                sum2 -= RPGgrid.get(i, j) * VoltageRealIm.get(j, 0) + RMGgrid.get(i, j) * VoltageRealIm.get(j + Nbus, 0);
            }
        }
        for (int j = 0; j < i; j++) {
            sum1 -= RMGgrid.get(i, j) * VoltageRealIm.get(j, 0) - RPGgrid.get(i, j) * VoltageRealIm.get(j + Nbus, 0);
            sum2 -= RPGgrid.get(i, j) * VoltageRealIm.get(j, 0) + RMGgrid.get(i, j) * VoltageRealIm.get(j + Nbus, 0);

        }
        for (int j = i + 1; j < Nbus; j++) {
            sum1 -= (RMGgrid.get(i, j) * VoltageRealIm.get(j, 0) - RPGgrid.get(i, j) * VoltageRealIm.get(j + Nbus, 0));
            sum2 -= (RPGgrid.get(i, j) * VoltageRealIm.get(j, 0) + RMGgrid.get(i, j) * VoltageRealIm.get(j + Nbus, 0));
        }
        VoltageRealIm.set(i, 0, sum1);
        VoltageRealIm.set(i + Nbus, 0, sum2);
    }*/
    //VoltageRealIm.display();
}


void CPUPFGS::calcE()
{
    std::cout << "calculE de CPUPFGS" << std::endl;
    if (_useDouble) {
        for (int i = 0; i < Nbus; i++) {
            float Rev = VoltageRealImD.get(i, 0);
            float Imv = VoltageRealImD.get(i + Nbus, 0);
            float Enorm = sqrt(Rev * Rev + Imv * Imv);
            float Eangle = atan2(Imv, Rev);
            ED.set(i + Nbus, 0, Enorm);
            ED.set(i, 0, Eangle);
        }
    }
    else {
        for (int i = 0; i < Nbus; i++) {
            float Rev = VoltageRealIm.get(i, 0);
            float Imv = VoltageRealIm.get(i + Nbus, 0);
            float Enorm = sqrt(Rev * Rev + Imv * Imv);
            float Eangle = atan2(Imv, Rev);
            E.set(i + Nbus, 0, Enorm);
            E.set(i, 0, Eangle);
        }
    }
   

}

void CPUPFGS::setE(MatrixCPU* Enew)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    E = *Enew;
   
    if (_useDouble) {
        E.toMatCPUD(ED);
        for (int i = 0; i < Nbus; i++) {
            float v = ED.get(i + Nbus, 0);
            float theta = ED.get(i, 0);
            VoltageRealImD.set(i, 0, v * cos(theta));
            VoltageRealImD.set(i + Nbus, 0, v * sin(theta));

        }
    }
    else {
        for (int i = 0; i < Nbus; i++) {
            float v = E.get(i + Nbus, 0);
            float theta = E.get(i, 0);
            VoltageRealIm.set(i, 0, v * cos(theta));
            VoltageRealIm.set(i + Nbus, 0, v * sin(theta));

        }
    }
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 8, (float) std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 8, 1);
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    calcW();
    //W.display();
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

void CPUPFGS::setE(MatrixCPUD* Enew)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    ED = *Enew;
    if (_useDouble) {
        for (int i = 0; i < Nbus; i++) {
            float v = ED.get(i + Nbus, 0);
            float theta = ED.get(i, 0);
            VoltageRealImD.set(i, 0, v * cos(theta));
            VoltageRealImD.set(i + Nbus, 0, v * sin(theta));

        }
    }
    else {
        E = ED;
        for (int i = 0; i < Nbus; i++) {
            float v = E.get(i + Nbus, 0);
            float theta = E.get(i, 0);
            VoltageRealIm.set(i, 0, v * cos(theta));
            VoltageRealIm.set(i + Nbus, 0, v * sin(theta));

        }
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
