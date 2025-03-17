#pragma once
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "MatrixCPU.h"

/*
Pour créer un cas :
  - taille 1 : Sbase, Vbase, nAgent, nCons, nGenSup, nBus, nLine, V0, theta0
  - taille N : PosBus, a, b, a^q, b^q, Pobj, Pmin, Pmax, Qobj, Qmin, Qmax, zone 
  - taille B : Gs, Bs, Vmin, Vmax, thetamin, thetamax V0, theta0
  - taille L : from, to, Ys Real, Ys Im, Yp, tau, theta, Limit=0, zs Real, zs Imag;

    facultatif :
  - taille N*N (ou reduit) : matrice de connexion
  - taille N*N (ou reduit) : matrices trade min et max
  - taille B*B (ou reduit) : matrices impedance 
*/
enum indInfo   { Sbase_ind, Vbase_ind, nAgent_ind, nCons_ind, nGen_ind, nBus_ind, nLine_ind, V0_ind, theta0_ind, finInfo_ind };
enum indAgent  { PosBus_ind, a_ind, b_ind, aq_ind, bq_ind, Pobj_ind, Pmin_ind, Pmax_ind, Qobj_ind, Qmin_ind, Qmax_ind, zone_ind, finAgent_ind };
enum indBuses  { Gs_ind, Bs_ind, Vmin_ind, Vmax_ind, thetamin_ind, thetamax_ind, Vinit_ind, thetainit_ind, finBuses_ind};
enum indBranch { From_ind, To_ind, YsRe_ind, YsIm_ind, Yp_ind, tau_ind, theta_ind, lim_ind, ZsRe_ind, ZsIm_ind, finBranch_ind};
const char* const indInfoStr[] = {"Sbase_ind","Vbase_ind","nAgent_ind","nCons_ind","nGen_ind","nBus_ind","nLine_ind","V0_ind","theta0_ind", "finInfo_ind"};



class StudyCaseInterface{
    private :
        MatrixCPU infoCase;
        MatrixCPU agentCase;
        MatrixCPU branchCase;
        MatrixCPU busCase;
        MatrixCPU connexion;
        MatrixCPU Lbmat;
        MatrixCPU Ubmat;
        MatrixCPU Gmat;
        MatrixCPU Bmat;
        MatrixCPU Beta;

        std::string _name = "default";

        int _N = 0;
        int _B = 0;
        int _L = 0;
        bool connexionDefined  = false;
        bool tradeBoundDefined = false;
        bool impendanceDefined = false;
        bool betaDefined       = false;
    public:
        StudyCaseInterface(int N, int B, int L);
        // taille 1 : Sbase, Vbase, nAgent, nCons, nGenSup, nBus, nLine, V0, theta0 
        void setSbase(float Sbase);
        void setVbase(float Vbase);
        void setV0(float V0);
        void setTheta(float theta0);
        void setName(std::string name);
        // taille N : PosBus, a, b, a^q, b^q Pobj, Pmin, Pmax, Qobj, Qmin, Qmax, zone 
        void setPosBus(MatrixCPU PosBus, MatrixCPU zone);
        void setCostFunction(MatrixCPU a, MatrixCPU b);
        void setPower(MatrixCPU Pmin, MatrixCPU Pmax, MatrixCPU Pobj); 
        // taille B : Gs, Bs, Vmin, Vmax, thetamin, thetamax, V0, theta0
        void setImpedanceBus(MatrixCPU Gs, MatrixCPU Bs);
        void setVoltageBound(MatrixCPU Vmin, MatrixCPU Vmax);
        void setAngleBound(MatrixCPU thetaMin, MatrixCPU thetaMax);
        void setVoltageInit(MatrixCPU V0, MatrixCPU theta0);
        //taille L : from, to, Ys Real, Ys Im, Yp, tau, theta, Limit=0, zs Real, zs Imag;
        void setLink(MatrixCPU from, MatrixCPU to);
        void setAdmitance(MatrixCPU YsRe, MatrixCPU YsIm, MatrixCPU Yb);
        void setImpedance(MatrixCPU zsRe, MatrixCPU zsIm);
        void setTransfo(MatrixCPU tau, MatrixCPU theta);
        void setLineLimit(MatrixCPU limit);
        // taille N*N : matrice de connexion
        void setConnexion(MatrixCPU connexion);
        // taille N*N : matrices trade min et max
        void setTradeLim(MatrixCPU min, MatrixCPU max);
        // taille B*B : matrices impedance
        void setMatImpedance(MatrixCPU Gs, MatrixCPU Bs);
        // taille N*N : matrix heterogenous preferences
        void setBeta(MatrixCPU beta);

        // getter 
        int getN();
        int getB();
        int getL();
        bool isConnexionDefined();
        bool isTradeBoundDefined();
        bool isImpedanceDefined();
        bool isBetaDefined();
        
        MatrixCPU getInfoCase();
        MatrixCPU getAgentCase();
        MatrixCPU getBranchCase();
        MatrixCPU getBusCase();
        MatrixCPU getConnexion();
        MatrixCPU getLbMat();
        MatrixCPU getUbMat();
        MatrixCPU getGmat();
        MatrixCPU getBmat();
        MatrixCPU getBeta();
        std::string getName();

        void display(int type=0);
        void checkCase(int nLineCons);

};



    /*static PyModuleDef custommodule = {
        .m_base = PyModuleDef_HEAD_INIT,
        .m_name = "custom",
        .m_doc = "Example module that creates an extension type.",
        .m_size = -1,
    };*/


    /*PyMODINIT_FUNC
    PyInit_beamodule(void)
    {
        PyObject *m;
        if (PyType_Ready(&CustomType) < 0)
            return NULL;

        m = PyModule_Create(&custommodule);
        if (m == NULL)
            return NULL;

        if (PyModule_AddObjectRef(m, "bea", (PyObject *) &CustomType) < 0) {
            Py_DECREF(m);
            return NULL;
        }

        return m;
    }*/
