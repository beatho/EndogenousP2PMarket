#include "../head/StudyCaseInterface.h"



// - taille 1 : Sbase, Vbase, nAgent, nCons, nGen, nBus, nLine, V0, theta0
//  - taille N : PosBus, a, b, a^q, b^q, Pobj, Pmin, Pmax, Qobj, Qmin, Qmax, zone 
//  - taille B : Gs, Bs, Vmin, Vmax, V0, theta0
//  - taille L : from, to, Ys Real, Ys Im, Yp, tau, theta, Limit=0, zs Real, zs Imag;

StudyCaseInterface::StudyCaseInterface(int N, int B, int L){
    infoCase   = MatrixCPU(1, finInfo_ind);
    agentCase  = MatrixCPU(N, finAgent_ind);
    busCase    = MatrixCPU(B, finBuses_ind);
    branchCase = MatrixCPU(L, finBranch_ind);

    infoCase.set(0, Sbase_ind, 1);
    infoCase.set(0, Vbase_ind, 1);
    infoCase.set(0, nAgent_ind, N);
    // 3 4 
    infoCase.set(0, nBus_ind, B);
    infoCase.set(0, nLine_ind, L);
    infoCase.set(0, V0_ind, 1);
    infoCase.set(0, theta0_ind, 0);

    _N = N;
    _B = B;
    _L = L;

    for (int i = 0; i < _B; i++)
    {
        busCase.set(i, Vinit_ind, 1);
        busCase.set(i, thetainit_ind, 0);
    }
    
}

void StudyCaseInterface::setSbase(float Sbase){
    infoCase.set(0, Sbase_ind, Sbase);
}
void StudyCaseInterface::setVbase(float Vbase){
    infoCase.set(0, Vbase_ind, Vbase);
}
void StudyCaseInterface::setV0(float V0){
    infoCase.set(0, V0_ind, V0);
}
void StudyCaseInterface::setTheta(float theta0){
    infoCase.set(0, theta0_ind, theta0);
}
void StudyCaseInterface::setName(std::string name){
    _name = name;
}
// taille N : PosBus, a, b, a^q, b^q, Pobj, Pmin, Pmax, Qobj, Qmin, Qmax, zone 
void StudyCaseInterface::setPosBus(MatrixCPU PosBus, MatrixCPU zone){
    for (int i = 0; i < _N; i++)
    {
        agentCase.set(i, PosBus_ind, PosBus.get(i,0));
        agentCase.set(i, zone_ind, zone.get(i,0));
    }
}
void StudyCaseInterface::setCostFunction(MatrixCPU a, MatrixCPU b){
    for (int i = 0; i < _N; i++)
    {
        agentCase.set(i, a_ind, a.get(i,0));
        agentCase.set(i, b_ind, b.get(i,0));
        agentCase.set(i, aq_ind, a.get(i + _N,0));
        agentCase.set(i, bq_ind, b.get(i + _N,0));
    }  
}
void StudyCaseInterface::setPower(MatrixCPU Pmin, MatrixCPU Pmax, MatrixCPU Pobj){
    for (int i = 0; i < _N; i++)
    {
        agentCase.set(i, Pobj_ind, Pobj.get(i,0));
        agentCase.set(i, Pmin_ind, Pmin.get(i,0));
        agentCase.set(i, Pmax_ind, Pmax.get(i,0));
        agentCase.set(i, Qobj_ind, Pobj.get(i + _N,0));
        agentCase.set(i, Qmin_ind, Pmin.get(i + _N,0));
        agentCase.set(i, Qmax_ind, Pmax.get(i + _N,0));
    }
    
}
// taille B : Gs, Bs, Vmin, Vmax, thetaMin, thetaMax, V0, theta0
void StudyCaseInterface::setImpedanceBus(MatrixCPU Gs, MatrixCPU Bs){
    for (int i = 0; i < _B; i++)
    {
        busCase.set(i, Gs_ind, Gs.get(i,0));
        busCase.set(i, Bs_ind, Bs.get(i,0));
    }
}
void StudyCaseInterface::setVoltageBound(MatrixCPU Vmin, MatrixCPU Vmax)
{
    for(int i=0; i<_B; i++){
        busCase.set(i, Vmin_ind, Vmin.get(i,0));
        busCase.set(i, Vmax_ind, Vmax.get(i,0));
    }
}
void StudyCaseInterface::setAngleBound(MatrixCPU thetaMin, MatrixCPU thetaMax)
{
    for (int i = 0; i < _B; i++)
    {
        busCase.set(i, thetamin_ind, thetaMin.get(i,0));
        busCase.set(i, thetamax_ind, thetaMax.get(i,0));
    }
}
void StudyCaseInterface::setVoltageInit(MatrixCPU V0, MatrixCPU theta0){
    for (int i = 0; i < _B; i++)
    {
        busCase.set(i, Vinit_ind, V0.get(i,0));
        busCase.set(i, thetainit_ind, theta0.get(i,0));
    }
}
//taille L : from, to, Ys Real, Ys Im, Yp, tau, theta, Limit=0, zs Real, zs Imag;
void StudyCaseInterface::setLink(MatrixCPU from, MatrixCPU to){
    for (int i = 0; i < _L; i++)
    {
        branchCase.set(i, From_ind, from.get(i,0));
        branchCase.set(i, To_ind, to.get(i,0));
    }
}
void StudyCaseInterface::setAdmitance(MatrixCPU YsRe, MatrixCPU YsIm, MatrixCPU Yb){
    for (int i = 0; i < _L; i++)
    {
        branchCase.set(i, YsRe_ind, YsRe.get(i,0));
        branchCase.set(i, YsIm_ind, YsIm.get(i,0));
        branchCase.set(i, Yp_ind, Yb.get(i,0));
    }
}
void StudyCaseInterface::setImpedance(MatrixCPU zsRe, MatrixCPU zsIm){
    for (int i = 0; i < _L; i++)
    {
        branchCase.set(i, ZsRe_ind, zsRe.get(i,0));
        branchCase.set(i, ZsIm_ind, zsIm.get(i,0));
    }
}
void StudyCaseInterface::setTransfo(MatrixCPU tau, MatrixCPU theta){
    for (int i = 0; i < _L; i++)
    {
        branchCase.set(i, tau_ind, tau.get(i,0));
        branchCase.set(i, theta_ind, theta.get(i,0));
    }
}
void StudyCaseInterface::setLineLimit(MatrixCPU limit){
    for (int i = 0; i < _L; i++)
    {
        branchCase.set(i, lim_ind, limit.get(i,0));
    }
    
}
// taille N*N (ou reduit) : matrice de connexion
void StudyCaseInterface::setConnexion(MatrixCPU Newconnexion){
    connexionDefined = true;
    connexion = Newconnexion;
}

// taille N*N (ou reduit) : matrices trade min et max
void StudyCaseInterface::setTradeLim(MatrixCPU min, MatrixCPU max){
    tradeBoundDefined = true;
    
    Lbmat = min;
    Ubmat = max;
    
}
// taille B*B (ou reduit) : matrices impedance
void StudyCaseInterface::setMatImpedance(MatrixCPU Gs, MatrixCPU Bs){
    impendanceDefined = true;
    Gmat = Gs;
    Bmat = Bs;
}

void StudyCaseInterface::setBeta(MatrixCPU beta){
    betaDefined = true;
    Beta = beta;
}

int StudyCaseInterface::getN(){return _N;}
int StudyCaseInterface::getB(){return _B;}
int StudyCaseInterface::getL(){return _L;};

bool StudyCaseInterface::isConnexionDefined(){return connexionDefined;}
bool StudyCaseInterface::isTradeBoundDefined(){return tradeBoundDefined;}
bool StudyCaseInterface::isImpedanceDefined(){return impendanceDefined;}
bool StudyCaseInterface::isBetaDefined(){return betaDefined;}  

MatrixCPU StudyCaseInterface::getInfoCase()
{
    return infoCase;
}
MatrixCPU StudyCaseInterface::getAgentCase()
{
    return agentCase;
}
MatrixCPU StudyCaseInterface::getBranchCase()
{
    return branchCase;
}
MatrixCPU StudyCaseInterface::getBusCase()
{
    return busCase;
}
MatrixCPU StudyCaseInterface::getConnexion()
{
    return connexion;
}
MatrixCPU StudyCaseInterface::getLbMat()
{
    return Lbmat;
}
MatrixCPU StudyCaseInterface::getUbMat()
{
    return Ubmat;
}
MatrixCPU StudyCaseInterface::getGmat()
{
    return Gmat;
}
MatrixCPU StudyCaseInterface::getBmat()
{
    return Bmat;
}
MatrixCPU StudyCaseInterface::getBeta(){return Beta;}

std::string StudyCaseInterface::getName() { return _name; }
void StudyCaseInterface::display(int type)
{
    std::cout << "Study Case of " << _N << " agent " << _B << " buses and " << _L << " lines" << std::endl;
    checkCase(_L);
    std::cout << "**************************************************"<<std::endl << std::endl; 
    if(type == 0 || type == 1){
        std::cout << " Info StudyCase : Sbase, Vbase, nAgent, nCons, nGen, nBus, nLine, V0, theta0 " << std::endl;
        infoCase.display();
        //std::cout << std::endl;
    }
    if(type == 0 || type == 2){
        std::cout << " Agent StudyCase : PosBus, a, b, a^q, b^q, Pobj, Pmin, Pmax, Qobj, Qmin, Qmax, zone  " << std::endl;
        agentCase.display();
        //std::cout << std::endl;
    }
    if(type == 0 || type == 3){
        std::cout << " Bus StudyCase : Gs, Bs, Vmin, Vmax, V0, theta0 " << std::endl;
        busCase.display();
        //std::cout << std::endl;
    }
    if(type == 0 || type == 4){
        std::cout << "Branches StudyCase : from, to, Ys Real, Ys Im, Yp, tau, theta, Limit=0, zs Real, zs Imag " << std::endl;
        branchCase.display();
        //std::cout << std::endl;
    }
    if(tradeBoundDefined && type == 0){
        std::cout << "bound defined min, max" << std::endl;
        Lbmat.display();
        Ubmat.display();
        std::cout << std::endl;
    }
    if(connexionDefined && type == 0){
        std::cout << "connexion defined" << std::endl;
        connexion.display();
        std::cout << std::endl; 
    }
    if(impendanceDefined && type == 0){
        std::cout << "Impedance defined Gmat, Bmat" <<std::endl;
        Gmat.display();
        Bmat.display();
    }


}

void StudyCaseInterface::checkCase(int nLineCons)
{
    int nCons = 0;
    int nGen  = 0; 

    int n=0;
    
    indInfo indnCons = nCons_ind;
    indInfo indnGen = nGen_ind;
    
    while(n< _N && agentCase.get(n, Pmax_ind) <= 0 ){ // consumers
        nCons++;
        n++;
        
    }
    
    while(n <_N && agentCase.get(n,Pmin_ind) >= 0){ // generators
        nGen++;
        n++;
        
    }
    //std::cout << "Checkcase nCons " << nCons << " nGen " << nGen <<std::endl;
    infoCase.set(0, nCons_ind, nCons);
    infoCase.set(0, nGen_ind, nGen);
    if(n!=_N){ // error or prosumers
        // prosumers => Pmin<0, Pmax >0
        // community => Pmin=Pmax=0
        // else erreur
        for(int i=n; i<_N; i++){
            if(agentCase.get(i,Pmin_ind)>0 || agentCase.get(i, Pmax_ind)<0){ // consumer or generator
                std::cout << "[WARNING] : the agents doest not respect a valide order, unexpected result can occur" << std::endl;
                std::cout << "Order must be consumers, generator, other" << std::endl;
                return;
            }
        }
    }


    if(nLineCons != _L){
        int counter = 0;
        for(int i=0; i<_L; i++){
            float lim = branchCase.get(i, lim_ind);
            if(lim>0){
                counter++;
            }
        }
        if(counter != nLineCons){
            std::cout << "[WARNING] : the number of line constraint gived is not coherent with the line data" <<std::endl;
        }
    }

}

//  - taille 1 : Sbase, Vbase, nAgent, nCons, nGen, nBus, nLine, V0, theta0
//  - taille N : PosBus, a, b, a^q, b^q, Pobj, Pmin, Pmax, Qobj, Qmin, Qmax, zone 
//  - taille B : Gs, Bs, Vmin, Vmax, V0, theta0
//  - taille L : from, to, Ys Real, Ys Im, Yp, tau, theta, Limit=0, zs Real, zs Imag;