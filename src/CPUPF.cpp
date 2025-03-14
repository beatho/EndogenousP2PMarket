#include "../head/CPUPF.h"
#include "..\head\CPUPF.h"

CPUPF::CPUPF(){}
CPUPF::~CPUPF(){
    
 }
//Fb0: init, Fb1ab 1 2 : Flu, Fb2abc : 3 4 5 Tension, FB3 6 : puissance, Fb4 7 erreur, 8 Fb0 mise � jour

/*void CPUPF::init(const StudyCase& cas, MatrixCPU* PQ)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
    timePerBlock = MatrixCPU(1, 9); // Fb0 : init, Fb1ab : Flu, Fb2abc: Tension , FB3 : puissance, Fb4 erreur, Fb0 mise � jour
    
    occurencePerBlock = MatrixCPU(1, 9);; //nb de fois utilis� pendant la simu
#endif // INSTRUMENTATION
    std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
    Nagent = cas.getNagent();
    Nbus = cas.getNBus();
    B2 = 2 * Nbus;
    N2 = 2 * Nagent;
    Nline = cas.getNLine(true); // ne doit pas �tre r�duit ici !!!
    Nconstraint = B2 + Nline;
    BL2 = Nbus + 2 * Nline;
    iterM = 20;
    iter = 0;
    V0 = cas.getV0();
    theta0 = cas.gettheta0();
    I = cas.getCoresBusAgentLin();
    _name = "Newton";
    _useDouble = false;
    status = 0;
   // std::cout << "init " << _name << std::endl;
   
    //I_aug.display();
   

    CoresVoiLin = cas.getCoresVoiLin();
    CoresBusLin = cas.getCoresBusLin();
    CoresLineBus = cas.getCoresLineBus(true);
    nLines = cas.getNLines();

    W = MatrixCPU(B2, 1);
    dW = MatrixCPU(B2, 1);
    E = MatrixCPU(B2, 1);
    dE = MatrixCPU(B2, 1);
    Jac = MatrixCPU(B2, B2);
    JacInv = MatrixCPU(B2, B2);
    A = MatrixCPU(B2, B2);
    P = MatrixCPU(B2 + 1, 1);

    //Bgrid = cas.getLineSuceptance();
    //Ggrid = cas.getLineReactance();
    BgridLin = cas.getBlin();
    GgridLin = cas.getGlin();
        
    for (int i = 0; i < Nbus; i++) {
        E.set(i, 0, theta0);
    }

    for (int i = Nbus; i < B2; i++) {
        E.set(i, 0, V0);
    }
    W0 = MatrixCPU(B2, 1);
    //PQ->display();
    for (int n = 1; n < Nagent; n++) {
        int bus = I.get(n, 0);
        W0.increment(bus, 0, PQ->get(n, 0));
        W0.increment(bus + Nbus, 0, PQ->get(n + Nagent, 0));
    }
    //std::cout << " W0 : " << std::endl;
    //W0.display();
    //std::cout << " Calcul W init : " << std::endl;
    calcW();
    //W.display();
    dW.subtract(&W0, &W);
        
    G = MatrixCPU(Nconstraint, N2);
    Phi = MatrixCPU(Nline, 1);
    Y = MatrixCPU(Nconstraint, 1);
    tempLN2 = MatrixCPU(Nline, N2);
    JacPhiE = MatrixCPU(Nline, B2);
    tempB2N2 = MatrixCPU(B2, N2);
    
   
    // W0[2 * N] : puissance active et r�active au noeud (I*[P Q])
    // W[2 * N] : puissance obtenue par calcul � partir de E
    // dW[2 * N] : derive de puissance
    // E[2 * N] : angle puis tension [O et 1] pour l'init ?
    // dE[2 * N] : derive de angle puis tension
    // Jac[2 * N][2 * N] : jacobienne
    // Jac_inv[2 * N][2 * N]: inverse de la jacobienne

    // B[N][N], G[N][N] : caract�ristique des lignes entre les noeuds i et j
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 0, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION
    //std::cout << " fin init" << std::endl;

}*/


void CPUPF::init(const StudyCase& cas, MatrixCPU* PQ, MatrixCPUD* PQD, bool useDouble)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
    timePerBlock = MatrixCPU(1, 9); // Fb0 : init, Fb1ab : Flu, Fb2abc: Tension , FB3 : puissance, Fb4 erreur, Fb0 mise � jour

    occurencePerBlock = MatrixCPU(1, 9);; //nb de fois utilis� pendant la simu
#endif // INSTRUMENTATION
    std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
    Nagent = cas.getNagent();
    Nbus = cas.getNBus();
    B2 = 2 * Nbus;
    N2 = 2 * Nagent;
    Nline = cas.getNLine(true); // ne doit pas �tre r�duit ici !!!
    Nconstraint = B2 + Nline;
    BL2 = Nbus + 2 * Nline;
    iterM = 10;
    iter = 0;
    V0 = cas.getV0();
    theta0 = cas.gettheta0();
    I = cas.getCoresBusAgentLin();
    _name = "Newton";
    _useDouble = useDouble;
    status = 0;
   
    //I_aug.display();


    CoresVoiLin = cas.getCoresVoiLin();
    CoresBusLin = cas.getCoresBusLin();
    CoresLineBus = cas.getCoresLineBus(true);
    nLines = cas.getNLines();




    if (_useDouble) {

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
        AD = MatrixCPUD(B2, B2);
        PD = MatrixCPUD(B2 + 1, 1);

        for (int i = 0; i < Nbus; i++) {
            ED.set(i, 0, theta0);
        }

        for (int i = Nbus; i < B2; i++) {
            ED.set(i, 0, V0);
        }


        W0D = MatrixCPUD(B2, 1);

        for (int n = 1; n < Nagent; n++) {
            int bus = I.get(n, 0);
            W0D.increment(bus, 0, PQD->get(n, 0));
            W0D.increment(bus + Nbus, 0, PQD->get(n + Nagent, 0));
        }


        //std::cout << " W0D : " << std::endl;
        //W0D.display();
        WD = MatrixCPUD(B2, 1);
        dWD = MatrixCPUD(B2, 1);

        /*Ggrid2Bgrid2 = MatrixCPU(Nbus, Nbus);
        for (int i = 0; i < Nbus; i++) {
            for (int j = 0; j < Nbus; j++) {
                Ggrid2Bgrid2.set(i, j, sqrt(GgridD.get(i, j) * GgridD.get(i, j) + BgridD.get(i, j) * BgridD.get(i, j)));
            }
        }*/

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
        A = MatrixCPU(B2, B2);
        P = MatrixCPU(B2 + 1, 1);

        //Bgrid = cas.getLineSuceptance();
        //Ggrid = cas.getLineReactance();
        BgridLin = cas.getBlin();
        GgridLin = cas.getGlin();

        for (int i = 0; i < Nbus; i++) {
            E.set(i, 0, theta0);
        }

        for (int i = Nbus; i < B2; i++) {
            E.set(i, 0, V0);
        }
        /*std::cout << " Bgrid : " << std::endl;
        BgridLin.display();
        std::cout << " Ggrid : " << std::endl;
        GgridLin.display();*/
        W0 = MatrixCPU(B2, 1);
        //PQ->display();
        for (int n = 1; n < Nagent; n++) {
            int bus = I.get(n, 0);
            W0.increment(bus, 0, PQ->get(n, 0));
            W0.increment(bus + Nbus, 0, PQ->get(n + Nagent, 0));
        }
        //std::cout << " W0 : " << std::endl;
        //W0.display();
        //std::cout << " Calcul W init : " << std::endl;
        calcW();
        //W.display();
        dW.subtract(&W0, &W);

        /*Ggrid2Bgrid2 = MatrixCPU(Nbus, Nbus);
        for (int i = 0; i < Nbus; i++) {
            for (int j = 0; j < Nbus; j++) {
                Ggrid2Bgrid2.set(i, j, sqrt(Ggrid.get(i, j) * Ggrid.get(i, j) + Bgrid.get(i, j) * Bgrid.get(i, j)));
            }
        }*/

    }






    G = MatrixCPU(Nconstraint, N2);
    Phi = MatrixCPU(Nline, 1);
    Y = MatrixCPU(Nconstraint, 1);
    tempLN2 = MatrixCPU(Nline, N2);
    JacPhiE = MatrixCPU(Nline, B2);
    tempB2N2 = MatrixCPU(B2, N2);





    /*std::cout << " E : " << std::endl;
    E.display();
    std::cout << " W : " << std::endl;
    W.display();*/

    // W0[2 * N] : puissance active et r�active au noeud (I*[P Q])
    // W[2 * N] : puissance obtenue par calcul � partir de E
    // dW[2 * N] : derive de puissance
    // E[2 * N] : angle puis tension [O et 1] pour l'init ?
    // dE[2 * N] : derive de angle puis tension
    // Jac[2 * N][2 * N] : jacobienne
    // Jac_inv[2 * N][2 * N]: inverse de la jacobienne

    // B[N][N], G[N][N] : caract�ristique des lignes entre les noeuds i et j
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 0, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 0, 1);
#endif // INSTRUMENTATION
}


void CPUPF::solve() {
  
    time = clock();
    err = 2 * epsPF;
    iter = 0;
    //std::cout << epsPF << " " << iterM << std::endl;
    int failure = 0;
    status = 1;
    while (err > epsPF && iter < iterM) {
        failure = calcVoltage();
        if (failure) {
            time = clock() - time;
            status = -1;
            //std::cout << "failure ! " << iter << " " << err << std::endl;
            return;
        }
#ifdef INSTRUMENTATION
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
            dWD.subtract(&W0D, &WD); // dW = W0 - W
            err = dWD.max2(); //err = ||dW||
        }
        else {
            dW.subtract(&W0, &W); // dW = W0 - W
            err = dW.max2(); //err = ||dW||
        }
#ifdef INSTRUMENTATION
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 7, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 7, 1);
#endif // INSTRUMENTATION

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
    timePerBlock.increment(0, 6, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 6, 1);
#endif // INSTRUMENTATION

    if (iter >= iterM) {
        status = 2;
        //std::cout << "fin solve " << iter<<" " << err << std::endl;
    }
    
    time = clock() - time;
        
}

void CPUPF::updatePQ(MatrixCPU* PQ)
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
    timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 8, 1);
#endif
}

void CPUPF::calcW(bool end)
{
    std::cout.precision(12);
    double p = 0;
    double q = 0;
        
    if (_useDouble) {

        double cdt = 0;
        double sdt = 0;
       
        double EE = 0;
        WD.set(0.0);
        for (int i = 0; i < Nbus; i++) {
            int k = CoresBusLin.get(i, 0);

            for (int voisin = k; voisin < (k + nLines.get(i, 0)); voisin++) {
                int j = CoresVoiLin.get(voisin, 0);
                double dTheta_ij = ED.get(i, 0) - ED.get(j, 0);
                sdt = sin(dTheta_ij);
                cdt = cos(dTheta_ij);
                EE = ED.get(Nbus + i, 0) * ED.get(Nbus + j, 0);
                p = EE * (GgridLinD.get(voisin, 0) * cdt + BgridLinD.get(voisin, 0) * sdt);
                q = EE * (GgridLinD.get(voisin, 0) * sdt - BgridLinD.get(voisin, 0) * cdt);

                WD.increment(i, 0,  p);
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
        double cdt = 0;
        double sdt = 0;

        double EE = 0;
       
        W.set(0.0);
        for (int i = 0; i < Nbus; i++) {
           // std::cout << i << std::endl;
            int k = CoresBusLin.get(i, 0);

            for (int voisin = k; voisin < (k + nLines.get(i, 0)); voisin++) {
                int j = CoresVoiLin.get(voisin, 0);
                double dTheta_ij = E.getD(i, 0) - E.getD(j, 0);
                sdt = sin(dTheta_ij);
                cdt = cos(dTheta_ij);
                EE = E.getD(Nbus + i, 0) * E.getD(Nbus + j, 0);
                p = EE * (GgridLin.getD(voisin, 0) * cdt + BgridLin.getD(voisin, 0) * sdt);
                q = EE * (GgridLin.getD(voisin, 0) * sdt - BgridLin.getD(voisin, 0) * cdt);

                W.increment(i, 0, p);
                W.increment(i + Nbus, 0, q);
            }
        }
        
        if (!end) { // pendant simu, la puissance � ce noeud est libre
            W.set(0, 0, W0.get(0, 0));
            W.set(Nbus, 0, W0.get(Nbus, 0));
        }
    }
    
    
    /*std::cout << "test = " << test << std::endl;
        

    */
    
}

void CPUPF::calcJac()
{
    if(_useDouble){
        for (int i = 1; i < Nbus; i++) {
            int k = CoresBusLin.get(i, 0);

            for (int voisin = k; voisin < (k + nLines.get(i, 0)); voisin++) {
                int j = CoresVoiLin.get(voisin, 0);
                int i2 = i + Nbus;
                int j2 = j + Nbus;
                if (i == j) {
                    JacD.set(i, i, -WD.get(i2, 0) - BgridLinD.get(voisin, 0) * ED.get(i2, 0) * ED.get(i2, 0));
                    JacD.set(i, i2, WD.get(i, 0) / ED.get(i2, 0) + GgridLinD.get(voisin, 0) * ED.get(i2, 0));
                    JacD.set(i2, i, WD.get(i, 0) - GgridLinD.get(voisin, 0) * ED.get(i2, 0) * ED.get(i2, 0));
                    JacD.set(i2, i2, WD.get(i2, 0) / ED.get(i2, 0) - BgridLinD.get(voisin, 0) * ED.get(i2, 0));

                }
                else {
                    double dTheta_ij = ED.get(i, 0) - ED.get(j, 0);
                    double inter = ED.get(i2, 0) * ED.get(j2, 0) * (GgridLinD.get(voisin, 0) * sin(dTheta_ij) - BgridLinD.get(voisin, 0) * cos(dTheta_ij));
                    JacD.set(i, j, inter);

                    inter = ED.get(i2, 0) * (GgridLinD.get(voisin, 0) * cos(dTheta_ij) + BgridLinD.get(voisin, 0) * sin(dTheta_ij));
                    JacD.set(i, j2, inter);
                    inter = -ED.get(i2, 0) * ED.get(j2, 0) * (GgridLinD.get(voisin, 0) * cos(dTheta_ij) + BgridLinD.get(voisin, 0) * sin(dTheta_ij));
                    JacD.set(i2, j, inter);
                    inter = ED.get(i2, 0) * (GgridLinD.get(voisin, 0) * sin(dTheta_ij) - BgridLinD.get(voisin, 0) * cos(dTheta_ij));
                    JacD.set(i2, j2, inter);

                }
            
            
            }
            
            /*for (int j = 0; j < Nbus; j++) {
                int i2 = i + Nbus;
                int j2 = j + Nbus;
                if (i == j) {
                    JacD.set(i, i, -WD.get(i2, 0) - BgridD.get(i, i) * ED.get(i2, 0) * ED.get(i2, 0));
                    JacD.set(i, i2, WD.get(i, 0) / ED.get(i2, 0) + GgridD.get(i, i) * ED.get(i2, 0));
                    JacD.set(i2, i, WD.get(i, 0) - GgridD.get(i, i) * ED.get(i2, 0) * ED.get(i2, 0));
                    JacD.set(i2, i2, WD.get(i2, 0) / ED.get(i2, 0) - BgridD.get(i, i) * ED.get(i2, 0));

                }
                else {
                    double dTheta_ij = ED.get(i, 0) - ED.get(j, 0);
                    double inter = ED.get(i2, 0) * ED.get(j2, 0) * (GgridD.get(i, j) * sin(dTheta_ij) - BgridD.get(i, j) * cos(dTheta_ij));
                    JacD.set(i, j, inter);

                    inter = ED.get(i2, 0) * (GgridD.get(i, j) * cos(dTheta_ij) + BgridD.get(i, j) * sin(dTheta_ij));
                    JacD.set(i, j2, inter);
                    inter = -ED.get(i2, 0) * ED.get(j2, 0) * (GgridD.get(i, j) * cos(dTheta_ij) + BgridD.get(i, j) * sin(dTheta_ij));
                    JacD.set(i2, j, inter);
                    inter = ED.get(i2, 0) * (GgridD.get(i, j) * sin(dTheta_ij) - BgridD.get(i, j) * cos(dTheta_ij));
                    JacD.set(i2, j2, inter);

                }
            }*/
        }
        // pas de noeud PV

        // noeud de ref
     
        JacD.set(0, 0, 1);
        JacD.set(Nbus, Nbus, 1);
    
    }
    else {
        for (int i = 1; i < Nbus; i++) {
            int k = CoresBusLin.get(i, 0);

            for (int voisin = k; voisin < (k + nLines.get(i, 0)); voisin++) {
                int j = CoresVoiLin.get(voisin, 0);
                int i2 = i + Nbus;
                int j2 = j + Nbus;
                if (i == j) {
                    Jac.set(i, i, -W.get(i2, 0) - BgridLin.get(voisin, 0) * E.get(i2, 0) * E.get(i2, 0));
                    Jac.set(i, i2, W.get(i, 0) / E.get(i2, 0) + GgridLin.get(voisin, 0) * E.get(i2, 0));
                    Jac.set(i2, i, W.get(i, 0) - GgridLin.get(voisin, 0) * E.get(i2, 0) * E.get(i2, 0));
                    Jac.set(i2, i2, W.get(i2, 0) / E.get(i2, 0) - BgridLin.get(voisin, 0) * E.get(i2, 0));

                }
                else {
                    double dTheta_ij = E.get(i, 0) - E.get(j, 0);
                    double inter = E.get(i2, 0) * E.get(j2, 0) * (GgridLin.get(voisin, 0) * sin(dTheta_ij) - BgridLin.get(voisin, 0) * cos(dTheta_ij));
                    Jac.set(i, j, inter);

                    inter = E.get(i2, 0) * (GgridLin.get(voisin, 0) * cos(dTheta_ij) + BgridLin.get(voisin, 0) * sin(dTheta_ij));
                    Jac.set(i, j2, inter);
                    inter = -E.get(i2, 0) * E.get(j2, 0) * (GgridLin.get(voisin, 0) * cos(dTheta_ij) + BgridLin.get(voisin, 0) * sin(dTheta_ij));
                    Jac.set(i2, j, inter);
                    inter = E.get(i2, 0) * (GgridLin.get(voisin, 0) * sin(dTheta_ij) - BgridLin.get(voisin, 0) * cos(dTheta_ij));
                    Jac.set(i2, j2, inter);

                }


            }
            /*for (int j = 0; j < Nbus; j++) {
                int i2 = i + Nbus;
                int j2 = j + Nbus;
                if (i == j) {
                    Jac.set(i, i, -W.get(i2, 0) - Bgrid.get(i, i) * E.get(i2, 0) * E.get(i2, 0));
                    Jac.set(i, i2, W.get(i, 0) / E.get(i2, 0) + Ggrid.get(i, i) * E.get(i2, 0));
                    Jac.set(i2, i, W.get(i, 0) - Ggrid.get(i, i) * E.get(i2, 0) * E.get(i2, 0));
                    Jac.set(i2, i2, W.get(i2, 0) / E.get(i2, 0) - Bgrid.get(i, i) * E.get(i2, 0));

                }
                else {
                    double dTheta_ij = E.get(i, 0) - E.get(j, 0);
                    double inter = (double)E.get(i2, 0) * E.get(j2, 0) * (Ggrid.get(i, j) * sin(dTheta_ij) - Bgrid.get(i, j) * cos(dTheta_ij));
                    Jac.set(i, j, inter);

                    inter = (double)E.get(i2, 0) * (Ggrid.get(i, j) * cos(dTheta_ij) + Bgrid.get(i, j) * sin(dTheta_ij));
                    Jac.set(i, j2, inter);
                    inter = (double)-E.get(i2, 0) * E.get(j2, 0) * (Ggrid.get(i, j) * cos(dTheta_ij) + Bgrid.get(i, j) * sin(dTheta_ij));
                    Jac.set(i2, j, inter);
                    inter = (double)E.get(i2, 0) * (Ggrid.get(i, j) * sin(dTheta_ij) - Bgrid.get(i, j) * cos(dTheta_ij));
                    Jac.set(i2, j2, inter);

                }
            }*/
        }
        // pas de noeud PV

        // noeud de ref
        Jac.set(0, 0, 1);
        Jac.set(Nbus, Nbus, 1);
        
    }
 
 
}


int CPUPF::calcVoltage()
{

#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    calcJac();
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 3, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 3, 1);
#endif // INSTRUMENTATION
    //std::cout << " Jac : " << std::endl;
    
    if (_useDouble) {
#ifdef INSTRUMENTATION
        t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION

#ifdef Eigen
        dED.solveSysEigen(&JacD, &dWD);
#else
    try
        {
            JacD.LUPFactorization(&AD, &PD);
        }
        catch (const std::exception&)
        {
            return 1;
        }
        dED.solveSys(&AD, &PD, &dWD);
#endif
        

        
#ifdef INSTRUMENTATION
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 4, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 4, 1);
#endif // INSTRUMENTATION

        /*JacInvD.invertGaussJordan(&JacD); //inversion matrice Jac
        dED.MultiplyMatVec(&JacInvD, &dWD);// dE = Jac_inv * dW; */
        ED.add(&ED, &dED);// E = E + dE;

    }
    else {
#ifdef INSTRUMENTATION
        t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
        
        #ifdef EIGEN
            dE.solveSysEigen(&Jac, &dW);
        #else
            try
            {
                Jac.LUPFactorization(&A, &P);
            }
            catch (const std::exception&)
            {
                return 1;
            }
            dE.solveSys(&A, &P, &dW);
        #endif
      
#ifdef INSTRUMENTATION
        t2 = std::chrono::high_resolution_clock::now();
        timePerBlock.increment(0, 4, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        occurencePerBlock.increment(0, 4, 1);
#endif // INSTRUMENTATION
        /*JacInv.invertEigen(&Jac); //inversion matrice Jac
        dE.MultiplyMatVec(&JacInv, &dW);// dE = Jac_inv * dW;*/
        E.add(&E, &dE);// E = E + dE;
    }
    return 0;
}



void CPUPF::setE(MatrixCPU* Enew)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    E = *Enew;
    if (_useDouble) {
        E.toMatCPUD(ED);
    }
#ifdef INSTRUMENTATION
    t2 = std::chrono::high_resolution_clock::now();
    timePerBlock.increment(0, 8, std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
    occurencePerBlock.increment(0, 8, 1);
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    calcW();
   // W.display();
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

void CPUPF::setE(MatrixCPUD* Enew)
{
#ifdef INSTRUMENTATION
    t1 = std::chrono::high_resolution_clock::now();
#endif // INSTRUMENTATION
    ED = *Enew;
    if (!_useDouble) {
        E = ED;
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

void CPUPF::setW(MatrixCPU* Wnew)
{
    W = *Wnew;
}



void CPUPF::calcE()
{
    // nothing to do
}

void CPUPF::calcPhi()
{
    calcE();
    for (int l = 0; l < Nline; l++) {
        int i = CoresLineBus.get(l, 0); //from 
        int i2 = i + Nbus;
        int j = CoresLineBus.get(l, 1); // to
        int j2 = j + Nbus;
        //float dTheta_ij = E.get(i, 0) - E.get(j, 0);

       //float p = E.get(i2, 0) * E.get(j2, 0) * (Ggrid.get(i, j) * cos(dTheta_ij) + Bgrid.get(i, j) * sin(dTheta_ij));
       //float q = E.get(i2, 0) * E.get(j2, 0) * (Ggrid.get(i, j) * sin(dTheta_ij) - Bgrid.get(i, j) * cos(dTheta_ij));
       
       float s = E.get(i2, 0) * E.get(j2, 0) * Ggrid2Bgrid2.get(i,j); // s = p + q*j
       
       if (i < j) { // est ce que c'est important le signe ?
           s = abs(s);
       }
       else {
           s = -abs(s);
       }
        
       Phi.set(l, 0, s);
        
    }


}

void CPUPF::calcJacPhiE()
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

MatrixCPU* CPUPF::calcG()
{
    calcJacPhiE();
   
    tempB2N2.MultiplyMatMat(&JacInv, &I_aug);
    
    tempLN2.MultiplyMatMat(&JacPhiE, &tempB2N2);
    
    G.setBloc(0, B2, 0, N2, &tempB2N2);
    
    G.setBloc(B2, Nconstraint, 0, N2, &tempLN2);

    return &G;
}

MatrixCPU* CPUPF::calcY()
{
    for (int i = 0; i < B2; i++) {
        Y.set(i, 0, E.get(i, 0));
    }
    for (int i = 0; i < Nline; i++) {
        Y.set(i + B2, 0, Phi.get(i, 0)); // puissance apparente !
    }
    return &Y;
}
float CPUPF::getPloss()
{
    float s = 0;
    if (_useDouble) {
        for (int i = 0; i < Nbus; i++) {
            s += WD.get(i, 0);
        }
    }
    else {
        for (int i = 0; i < Nbus; i++) {
            s += W.get(i, 0);
        }
    }
   
    return s; // desequilibre !
}

float CPUPF::getQloss()
{
    float s = 0;
    if (_useDouble) {
        for (int i = Nbus; i < B2; i++) {
            s += WD.get(i, 0);
        }
    }
    else {
        for (int i = Nbus; i < B2; i++) {
            s += W.get(i, 0);
        }
    }
    
    return s; // desequilibre !
}

float CPUPF::getRes()
{
    return err;
}

int CPUPF::getIter()
{
    return iter;
}

float CPUPF::getP0()
{
    if (_useDouble) {
        return WD.get(0, 0);
    }
    else {
        return W.get(0 ,0);
    }
}

float CPUPF::getQ0()
{
    if (_useDouble) {
        return WD.get(Nbus, 0);
    }
    else {
        return W.get(Nbus, 0);
    }
}

float CPUPF::getTime()
{
    return  (float) time / CLOCKS_PER_SEC;
}

int CPUPF::getConv()
{
    return status;
}

MatrixCPU CPUPF::getE()
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
    
    return E;
}
MatrixCPU CPUPF::getW()
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
    return W;
}

MatrixCPU CPUPF::getY()
{
    MatrixCPU Y(2 * Nbus + Nline, 1);
    for (int b = 0; b < B2; b++) {
        Y.set(b, 0, E.get(b, 0));
    }
    int line = 0;
    for (int b = 0; b < Nbus; b++) {
        int k = CoresBusLin.get(b, 0);
        float thetai = E.get(b, 0);
        float Vi = E.get(b + Nbus, 0);
        float ei = Vi * cos(thetai);
        float fi = Vi * sin(thetai);
       
        for (int voisin = k + 1; voisin < (k + nLines.get(b, 0)); voisin++) {
            int j = CoresVoiLin.get(voisin, 0);
            if (j > b) {
                float thetaj = E.get(j, 0);
                float Vj = E.get(j + Nbus, 0);
                float B = BgridLin.get(voisin, 0);
                float G = GgridLin.get(voisin, 0);
                float ej = Vj * cos(thetaj);
                float fj = Vj * sin(thetaj);
                
                float Pij = (ei * ei + fi * fi - ei * ej - fi * fj) * G + (ei * fj - ej * fi) * B;

                Y.set(2 * Nbus + line, 0, Pij);
             
                line++;
            }
        }
    }

    return Y;
}

void CPUPF::display()
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

void CPUPF::display2(bool all)
{
    std::cout.precision(3);
    
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
        std::cout << "method " << _name << " converged in " << iter << " iterations." << std::endl;
        std::cout << "Converged in " << (float)time / CLOCKS_PER_SEC << " seconds" << std::endl;
        if (_useDouble) {
            std::cout << " Computation with double precision" << std::endl;
            /*double temp = WD.get(0, 0);
            double temp2 = WD.get(Nbus, 0);
            WD.set(0, 0, W0D.get(0, 0));
            WD.set(Nbus, 0, W0D.get(Nbus, 0));
            dWD.subtract(&W0D, &WD); // dW = W0 - W
            err = dWD.max2(); //err = ||dW||
            for (int b = 0; b < Nbus; b++) {
                std::cout << b << " " << dW.get(b, 0) << " " << dW.get(Nbus + b, 0) << std::endl;
            }
            WD.set(0, 0, temp);
            WD.set(Nbus, 0, temp2);*/
        }
        else {
            std::cout << " Computation with float simple precision" << std::endl;
            /*float temp = W.get(0, 0);
            float temp2 = W.get(Nbus, 0);
            W.set(0, 0, W0.get(0, 0));
            W.set(Nbus, 0, W0.get(Nbus, 0));
            dW.subtract(&W0, &W); // dW = W0 - W
            err = dW.max2(); //err = ||dW||
            for (int b = 0; b < Nbus; b++) {
                std::cout << b << " " << dW.get(b, 0) << " " << dW.get(Nbus + b, 0) << std::endl;
            }
            W.set(0, 0, temp);
            W.set(Nbus, 0, temp2);*/
        }
        
    }
    else {
        std::cout << "method " << _name << " not converged in " << iter << " iterations." << std::endl;
        std::cout << "time taken " << (float)time / CLOCKS_PER_SEC << " seconds" << std::endl;
    }
    std::cout << "The error of this state is " << err << std::endl;
    std::cout << "===============================================================|" << std::endl;
    std::cout << "      System Summary                                           |" << std::endl;
    std::cout << "===============================================================|" << std::endl;
    std::cout << "Buses            " << Nbus << std::endl;
    std::cout << "Branches         " << Nline << std::endl;
    std::cout << "Ploss            " << getPloss() << std::endl;
    std::cout << "Qloss            " << getQloss() << std::endl;


    std::cout <<std::endl << std::endl;
    if (all) {

        std::cout << "===============================================================================================|" << std::endl;
        std::cout << "      Bus Data                                                                                 |" << std::endl;
        std::cout << "===============================================================================================|" << std::endl;
        std::cout << " Bus |          Voltage        |  Power = Generation  + Load   |  Init = Generation  + Load    |" << std::endl;
        std::cout << "  #  |    Mag(pu) |  Ang(deg)  |    P (pu)     |     Q (pu)    |    P (pu)     |     Q (pu)    |" << std::endl;
        std::cout << "-----|------------|------------|---------------|---------------|---------------|---------------|" << std::endl;

        //std::cout << 0 << "      " << E.get(Nbus, 0) << "             " << E.get(0, 0) * (abs(E.get(0, 0)) > 0.0001) * 180 / 3.1415 << "              " << (abs(W.get(0, 0)) > 0.0001) * W.get(0, 0) << "         " << (abs(W.get(Nbus, 0)) > 0.0001) * W.get(Nbus, 0) << std::endl;

        float seuil = 0.0001;
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
        float seuil = 0.0001;
        std::cout << " Bus |          Voltage        |  Power = Generation  + Load   |  Init = Generation  + Load    |" << std::endl;
        std::cout << "  #  |    Mag(pu) |  Ang(deg)  |    P (pu)     |     Q (pu)    |    P (pu)     |     Q (pu)    |" << std::endl;
        std::cout << "-----|------------|------------|---------------|---------------|---------------|---------------|" << std::endl;
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

void CPUPF::saveTimeBlock(std::string fileName)
{
    std::ios_base::openmode mode = std::fstream::in | std::fstream::out | std::fstream::app;
    float factor = 1000000; // go from ns to ms fot the time


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
