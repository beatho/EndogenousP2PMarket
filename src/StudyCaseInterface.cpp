#include "../head/StudyCaseInterface.h"


// - taille 1 : Sbase, Vbase, nAgent, nCons, nGenSup, nBus, nLine, V0, theta0
//  - taille N : PosBus, a, b, a^q, b^q, Pobj, Pmin, Pmax, Qobj, Qmin, Qmax, zone 
//  - taille B : Gs, Bs, Vmin, Vmax, V0, theta0
//  - taille L : from, to, Ys Real, Ys Im, Yp, tau, theta, Limit=0, zs Real, zs Imag;

StudyCaseInterface::StudyCaseInterface(int N, int B, int L){
    infoCase   = MatrixCPU(1, 9);
    agentCase  = MatrixCPU(N, 12);
    busCase    = MatrixCPU(B, 6);
    branchCase = MatrixCPU(L, 10);

    infoCase.set(0, 0, 1);
    infoCase.set(0, 1, 1);
    infoCase.set(0, 2, N);
    // 3 4 
    infoCase.set(0, 5, B);
    infoCase.set(0, 6, L);
    infoCase.set(0, 7, 1);
    infoCase.set(0, 8, 0);

    _N = N;
    _B = B;
    _L = L;

    for (int i = 0; i < _B; i++)
    {
        busCase.set(i, 4, 1);
        busCase.set(i, 5, 0);
    }

}

void StudyCaseInterface::setSbase(float Sbase){
    infoCase.set(0, 0, Sbase);
}
void StudyCaseInterface::setVbase(float Vbase){
    infoCase.set(0, 1, Vbase);
}
void StudyCaseInterface::setV0(float V0){
    infoCase.set(0, 7, V0);
}
void StudyCaseInterface::setTheta(float theta0){
    infoCase.set(0, 8, theta0);
}
void StudyCaseInterface::setName(std::string name){
    _name = name;
}
// taille N : PosBus, a, b, a^q, b^q, Pobj, Pmin, Pmax, Qobj, Qmin, Qmax, zone 
void StudyCaseInterface::setPosBus(MatrixCPU PosBus, MatrixCPU zone){
    for (int i = 0; i < _N; i++)
    {
        agentCase.set(i, 0, PosBus.get(i,0));
        agentCase.set(i, 11, zone.get(i,0));
    }
}
void StudyCaseInterface::setCostFunction(MatrixCPU a, MatrixCPU b){
    for (int i = 0; i < _N; i++)
    {
        agentCase.set(i, 1, a.get(i,0));
        agentCase.set(i, 2, b.get(i,0));
        agentCase.set(i, 3, a.get(i + _N,0));
        agentCase.set(i, 4, b.get(i + _N,0));
    }  
}
void StudyCaseInterface::setPower(MatrixCPU Pmin, MatrixCPU Pmax, MatrixCPU Pobj){
    for (int i = 0; i < _N; i++)
    {
        agentCase.set(i, 5, Pobj.get(i,0));
        agentCase.set(i, 6, Pmin.get(i,0));
        agentCase.set(i, 7, Pmax.get(i,0));
        agentCase.set(i, 8, Pobj.get(i + _N,0));
        agentCase.set(i, 9, Pmin.get(i + _N,0));
        agentCase.set(i, 10, Pmax.get(i + _N,0));
    }
    
}
// taille B : Gs, Bs, Vmin, Vmax, thetaMin, thetaMax, V0, theta0
void StudyCaseInterface::setImpedanceBus(MatrixCPU Gs, MatrixCPU Bs){
    for (int i = 0; i < _B; i++)
    {
        busCase.set(i, 0, Gs.get(i,0));
        busCase.set(i, 1, Bs.get(i,0));
    }
}
void StudyCaseInterface::setVoltageBound(MatrixCPU Vmin, MatrixCPU Vmax)
{
    for(int i=0; i<_B; i++){
        busCase.set(i, 2, Vmin.get(i,0));
        busCase.set(i, 3, Vmax.get(i,0));
    }
}
void StudyCaseInterface::setAngleBound(MatrixCPU thetaMin, MatrixCPU thetaMax)
{
    for (int i = 0; i < _B; i++)
    {
        busCase.set(i, 4, thetaMin.get(i,0));
        busCase.set(i, 5, thetaMax.get(i,0));
    }
}
void StudyCaseInterface::setVoltageInit(MatrixCPU V0, MatrixCPU theta0){
    for (int i = 0; i < _B; i++)
    {
        busCase.set(i, 6, V0.get(i,0));
        busCase.set(i, 7, theta0.get(i,0));
    }
}
//taille L : from, to, Ys Real, Ys Im, Yp, tau, theta, Limit=0, zs Real, zs Imag;
void StudyCaseInterface::setLink(MatrixCPU from, MatrixCPU to){
    for (int i = 0; i < _L; i++)
    {
        branchCase.set(i, 0, from.get(i,0));
        branchCase.set(i, 1, to.get(i,0));
    }
}
void StudyCaseInterface::setAdmitance(MatrixCPU YsRe, MatrixCPU YsIm, MatrixCPU Yb){
    for (int i = 0; i < _L; i++)
    {
        branchCase.set(i, 2, YsRe.get(i,0));
        branchCase.set(i, 3, YsIm.get(i,0));
        branchCase.set(i, 4, Yb.get(i,0));
    }
}
void StudyCaseInterface::setImpedance(MatrixCPU zsRe, MatrixCPU zsIm){
    for (int i = 0; i < _L; i++)
    {
        branchCase.set(i, 8, zsRe.get(i,0));
        branchCase.set(i, 9, zsIm.get(i,0));
    }
}
void StudyCaseInterface::setTransfo(MatrixCPU tau, MatrixCPU theta){
    for (int i = 0; i < _L; i++)
    {
        branchCase.set(i, 5, tau.get(i,0));
        branchCase.set(i, 6, theta.get(i,0));
    }
}
void StudyCaseInterface::setLineLimit(MatrixCPU limit){
    for (int i = 0; i < _L; i++)
    {
        branchCase.set(i, 7, limit.get(i,0));
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

int StudyCaseInterface::getN(){return _N;}
int StudyCaseInterface::getB(){return _B;}
int StudyCaseInterface::getL(){return _L;};

bool StudyCaseInterface::isConnexionDefined(){return connexionDefined;}
bool StudyCaseInterface::isTradeBoundDefined(){return tradeBoundDefined;}
bool StudyCaseInterface::isImpedanceDefined(){return impendanceDefined;}

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
std::string StudyCaseInterface::getName(){return _name;}