#pragma once
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "MatrixCPU.h"

/*
Pour créer un cas :
  - taille 1 : Sbase, Vbase, nAgent, nCons, nGenSup, nBus, nLine, V0, theta0
  - taille N : PosBus, a, b, Pobj, Pmin, Pmax, Qobj, Qmin, Qmax, zone 
  - taille B : Gs, Bs, Vmin, Vmax, V0, theta0
  - taille L : from, to, Ys Real, Ys Im, Yp, tau, theta, Limit=0, zs Real, zs Imag;

    facultatif :
  - taille N*N (ou reduit) : matrice de connexion
  - taille N*N (ou reduit) : matrices trade min et max
  - taille B*B (ou reduit) : matrices impedance 
*/

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

        std::string _name = "default";

        int _N = 0;
        int _B = 0;
        int _L = 0;
        bool connexionDefined  = false;
        bool tradeBoundDefined = false;
        bool impendanceDefined = false;
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
        // taille N*N (ou reduit) : matrice de connexion
        void setConnexion(MatrixCPU connexion);
        // taille N*N (ou reduit) : matrices trade min et max
        void setTradeLim(MatrixCPU min, MatrixCPU max);
        // taille B*B (ou reduit) : matrices impedance
        void setMatImpedance(MatrixCPU Gs, MatrixCPU Bs);

        // getter 
        int getN();
        int getB();
        int getL();
        bool isConnexionDefined();
        bool isTradeBoundDefined();
        bool isImpedanceDefined();
        
        MatrixCPU getInfoCase();
        MatrixCPU getAgentCase();
        MatrixCPU getBranchCase();
        MatrixCPU getBusCase();
        MatrixCPU getConnexion();
        MatrixCPU getLbMat();
        MatrixCPU getUbMat();
        MatrixCPU getGmat();
        MatrixCPU getBmat();
        std::string getName();

};

static MatrixCPU convertListToVectorCPUf(PyObject* list){
    MatrixCPU toReturn(PyList_GET_SIZE(list),1);
    for (int i=0; i <PyList_GET_SIZE(list); i++ ){
        PyObject* pItem = PyList_GetItem(list, i);
        double Gs = PyFloat_AsDouble(pItem);
        toReturn.set(i, 0, (float) Gs);
    }
    return toReturn;
}

static MatrixCPU convertListToVectorCPUi(PyObject* list){
    MatrixCPU toReturn(PyList_GET_SIZE(list),1);
    for (int i=0; i <PyList_GET_SIZE(list); i++ ){
        PyObject* pItem = PyList_GetItem(list, i);
        long Gs = PyLong_AsDouble(pItem);
        toReturn.set(i, 0, (float) Gs);
    }
    return toReturn;
}

static MatrixCPU convertListToMatrixCPUf(PyObject* list, int collum, int row){
    MatrixCPU toReturn(row,collum);
    int N = PyList_GET_SIZE(list);
    if(N != collum*row){
        PyErr_SetString(PyExc_TypeError, "Wrong size for the matrix from a list");
        return MatrixCPU(0,0);
    }
    int k = 0;
    for (int i=0; i <row; i++ ){
        for (int j = 0; i < collum; j++)
        {
            PyObject* pItem = PyList_GetItem(list, k);
            double Gs = PyFloat_AsDouble(pItem);
            toReturn.set(i, j, (float) Gs);
            k++;
        }
    }
    return toReturn;
}
static MatrixCPU convertListToMatrixCPUi(PyObject* list, int collum, int row){
    MatrixCPU toReturn(row,collum);
    int N = PyList_GET_SIZE(list);
    if(N != collum*row){
        PyErr_SetString(PyExc_TypeError, "Wrong size for the matrix from a list");
        return MatrixCPU(0,0);
    }
    int k = 0;
    for (int i=0; i <row; i++ ){
        for (int j = 0; i < collum; j++)
        {
            PyObject* pItem = PyList_GetItem(list, k);
            long Gs = PyLong_AsDouble(pItem);
            toReturn.set(i, j, (float) Gs);
            k++;
        }
    }
    return toReturn;
}

extern "C"{

    typedef struct {
        PyObject_HEAD
        StudyCaseInterface* interfaceCase;
        /* Type-specific fields go here. */
    } CustomObject;

    static PyTypeObject CustomStudyCase;

    static int MonObject_init(CustomObject* self, PyObject* args, PyObject* kwds){
        std::cout << "creation d'un interface cas d'etude" << std::endl;
        int N, B, L;
        if(!PyArg_ParseTuple(args, "iii", &N, &B, &L)){
            return -1;
        }
        if(N<0 || B<0 || L<0){
            std::cout << "[Error] : N, B, L are sizes, must be positive" << std::endl;
            return -1;
        }
        self->interfaceCase = new StudyCaseInterface(N, B, L);
        return 0;
    }
    static void MonObject_dealloc(CustomObject* self){
        delete self->interfaceCase;
    }

    static PyObject* StudyCase_setSbase(CustomObject* self, PyObject* args){
        float Sbase;
        if(!PyArg_ParseTuple(args, "f", &Sbase)){
            return NULL;
        }
        if(Sbase<=0){
            std::cout << "[ERROR] : Sbase must be positive" <<std::endl;
            return NULL;
        }
        self->interfaceCase->setSbase(Sbase);
        return Py_None;
    }
    static PyObject* StudyCase_setVbase(CustomObject* self, PyObject* args){
        float Vbase;
        if(!PyArg_ParseTuple(args, "f", &Vbase)){
            return NULL;
        }
        if(Vbase<=0){
            std::cout << "[ERROR] : Vbase must be positive" <<std::endl;
            return NULL;
        }
        self->interfaceCase->setVbase(Vbase);
        return Py_None;
    }
    static PyObject* StudyCase_setV0(CustomObject* self, PyObject* args){
        float V0;
        if(!PyArg_ParseTuple(args, "f", &V0)){
            return NULL;
        }
        if(V0<=0){
            std::cout << "[ERROR] : V0 must be positive" <<std::endl;
            return NULL;
        }
        self->interfaceCase->setV0(V0);
        return Py_None;
    }
    static PyObject* StudyCase_setTheta(CustomObject* self, PyObject* args){
        float theta0;
        if(!PyArg_ParseTuple(args, "f", &theta0)){
            return NULL;
        }
        if(theta0<-3.14159265359 || theta0>3.14159265359){
            std::cout << "[ERROR] : theta0 must be in radian" <<std::endl;
            return NULL;
        }
        self->interfaceCase->setTheta(theta0);
        return Py_None;
    }
    static PyObject* StudyCase_setName(CustomObject* self, PyObject* args){
        std::string name; 
        char* buffer;
        if(!PyArg_ParseTuple(args, "s", &buffer)){
            PyErr_SetString(PyExc_TypeError, "parameter must be a list.");
            return NULL;
        }

        name = buffer;
        self->interfaceCase->setName(name);
        
    }
    // taille N : PosBus, a, b, a^q, b^q Pobj, Pmin, Pmax, Qobj, Qmin, Qmax, zone 
    static PyObject* StudyCase_setPosBus(CustomObject* self, PyObject* args){
        PyObject *pList;
        PyObject *pList2;
        Py_ssize_t n;

        if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &pList, &PyList_Type, &pList2)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be a list.");
            return NULL;
        }

        n = PyList_Size(pList);
        if(n != self->interfaceCase->getN()){
            PyErr_SetString(PyExc_TypeError, "The size must be equal to N");
            return NULL;
        }
        if(n != PyList_Size(pList2)){
            PyErr_SetString(PyExc_TypeError, "The size of the second list must be the same as the first");
            return NULL;
        }
        MatrixCPU PosBus = convertListToVectorCPUi(pList);
        MatrixCPU zone = convertListToVectorCPUi(pList2);   
        
        for (int i=0; i<n; i++) {
            int bus = PosBus.get(i,0);
            if(bus<-1 || bus >self->interfaceCase->getB()){
                PyErr_SetString(PyExc_TypeError, "Position must be a int between -1 (not on the grid) and B (exclude)");
                return NULL;
            }
        }
        self->interfaceCase->setPosBus(PosBus, zone);
    }
    static PyObject* StudyCase_setCostFunction(CustomObject* self, PyObject* args){
        PyObject *pList;
        PyObject *pList2;
        Py_ssize_t n;

        if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &pList, &PyList_Type, &pList2)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be two lists.");
            return NULL;
        }

        n = PyList_Size(pList);
        int N = self->interfaceCase->getN();
        if(n != N && n!= 2*N){
            PyErr_SetString(PyExc_TypeError, "The size must be equal to N or to 2*N");
            return NULL;
        }
        if(n!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_TypeError, "The size of the second list must be the same as the first");
            return NULL;
        }

        MatrixCPU temp_a = convertListToVectorCPUf(pList);
        MatrixCPU temp_b = convertListToVectorCPUf(pList2);

        MatrixCPU a_vect(2*N, 1);
        MatrixCPU b_vect(2*N, 1);

        for (int i=0; i<n; i++) {
            float a = temp_a.get(i,0);
            if(a<=0){
                PyErr_SetString(PyExc_TypeError, "a must be positive");
                return NULL;
            }
            a_vect.set(i, 0, a);
            b_vect.set(i, 0, temp_b.get(i,0));
        }
        self->interfaceCase->setCostFunction(a_vect, b_vect);
    }
    static PyObject* StudyCase_setPower(CustomObject* self, PyObject* args){
        PyObject *pList;
        PyObject *pList2;
        PyObject *pList3 = nullptr;
        Py_ssize_t n;

        if (!PyArg_ParseTuple(args, "O!O!|O!", &PyList_Type, &pList, &PyList_Type, &pList2, &PyList_Type, &pList3)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be two or three lists.");
            return NULL;
        }

        n = PyList_Size(pList);
        int N = self->interfaceCase->getN();

        if(n != N && n!=2*N){
            PyErr_SetString(PyExc_TypeError, "The size must be equal to N or to 2*N");
            return NULL;
        } if(n!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_TypeError, "The size of the second list must be the same as the first");
            return NULL;
        } if(pList3!=nullptr){
            if(n!=PyList_Size(pList3)){
                PyErr_SetString(PyExc_TypeError, "The size of the thrid list must be the same as the first");
                return NULL;
            }
        }

        MatrixCPU Pmin_temp = convertListToVectorCPUf(pList);
        MatrixCPU Pmax_temp = convertListToVectorCPUf(pList2);
        MatrixCPU Pobj_temp;
        if(pList3){
            Pobj_temp = convertListToVectorCPUf(pList3);
        }

        MatrixCPU Pmin_vect(2*N,1);
        MatrixCPU Pmax_vect(2*N,1);
        MatrixCPU Pobj_vect(2*N,1);
        for (int i=0; i<n; i++) {
            float Pmin = Pmin_temp.get(i,0);
            Pmin_vect.set(i, 0, Pmin);
            
           
            float Pmax = Pmax_temp.get(i,0);
            Pmax_vect.set(i, 0, Pmax);

            if(Pmax<Pmin){
                PyErr_SetString(PyExc_TypeError, "Pmax must be greater than Pmin");
                return NULL;
            }
            float Pobj = (Pmin + Pmax)/2;
            if(pList3){
               Pobj = Pobj_temp.get(i,0);
            }
            if(Pobj<Pmin || Pobj>Pmax){
                std::cout << "[WARNING] : Pobj is not within the bound !!!!" << std::endl;
            }
            Pobj_vect.set(i, 0, Pobj);
        }
        self->interfaceCase->setPower(Pmin_vect, Pmax_vect, Pobj_vect);
    }
    
    // taille B : Gs, Bs, Vmin, Vmax, V0, theta0
    static PyObject* StudyCase_setImpedanceBus(CustomObject* self, PyObject* args){
        PyObject *pList;
        PyObject *pList2;
        Py_ssize_t b;

        if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &pList, &PyList_Type, &pList2)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be two lists.");
            return NULL;
        }

        
        b = PyList_Size(pList);
        int B = self->interfaceCase->getB();
        if(b != B){
            PyErr_SetString(PyExc_TypeError, "The size must be equal to B");
            return NULL;
        }
        if(b!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_TypeError, "The size of the second list must be the same as the first");
            return NULL;
        }

        MatrixCPU Gs_vect = convertListToVectorCPUf(pList);
        MatrixCPU Bs_vect = convertListToVectorCPUf(pList2);

        self->interfaceCase->setImpedanceBus(Gs_vect, Bs_vect);
    }
    static PyObject* StudyCase_setVoltageBound(CustomObject* self, PyObject* args){
        PyObject *pList;
        PyObject *pList2;
        Py_ssize_t b;

        if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &pList, &PyList_Type, &pList2)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be two lists.");
            return NULL;
        }

        
        b = PyList_Size(pList);
        int B = self->interfaceCase->getB();
        if(b != B){
            PyErr_SetString(PyExc_TypeError, "The size must be equal to B");
            return NULL;
        }
        if(b!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_TypeError, "The size of the second list must be the same as the first");
            return NULL;
        }

        MatrixCPU Vmin_vect = convertListToVectorCPUf(pList);
        MatrixCPU Vmax_vect = convertListToVectorCPUf(pList2);

        for(int i=0; i<B; i++){
            float Vmin = Vmin_vect.get(i,0);
            float Vmax = Vmax_vect.get(i,0);

            if(Vmin>Vmax){
                PyErr_SetString(PyExc_TypeError, "Vmax must be greater than Vmin");
            }

        }
        self->interfaceCase->setVoltageBound(Vmin_vect, Vmax_vect);
    }
    static PyObject* StudyCase_setThetaBound(CustomObject* self, PyObject* args){
        PyObject *pList;
        PyObject *pList2;
        Py_ssize_t b;

        if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &pList, &PyList_Type, &pList2)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be two lists.");
            return NULL;
        }

        
        b = PyList_Size(pList);
        int B = self->interfaceCase->getB();
        if(b != B){
            PyErr_SetString(PyExc_TypeError, "The size must be equal to B");
            return NULL;
        }
        if(b!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_TypeError, "The size of the second list must be the same as the first");
            return NULL;
        }

        MatrixCPU thetamin_vect = convertListToVectorCPUf(pList);
        MatrixCPU thetamax_vect = convertListToVectorCPUf(pList2);

        for(int i=0; i<B; i++){
            float thetamin = thetamin_vect.get(i,0);
            float thetamax = thetamax_vect.get(i,0);

            if(thetamin>thetamax){
                PyErr_SetString(PyExc_TypeError, "thetamax must be greater than thetamin");
                return NULL;
            }
            if(thetamin<-3.14159265359 || thetamax>3.14159265359){
                std::cout << "[ERROR] : theta must be in radian" <<std::endl;
                return NULL;
            }
        }
        self->interfaceCase->setAngleBound(thetamin_vect, thetamax_vect);
    }
    static PyObject* StudyCase_setVoltageInit(CustomObject* self, PyObject* args){
        PyObject *pList;
        PyObject *pList2;
        Py_ssize_t b;

        if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &pList, &PyList_Type, &pList2)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be two lists.");
            return NULL;
        }

        
        b = PyList_Size(pList);
        int B = self->interfaceCase->getB();
        if(b != B){
            PyErr_SetString(PyExc_TypeError, "The size must be equal to B");
            return NULL;
        }
        if(b!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_TypeError, "The size of the second list must be the same as the first");
            return NULL;
        }

        MatrixCPU V0_vect = convertListToVectorCPUf(pList);
        MatrixCPU theta0_vect = convertListToVectorCPUf(pList2);

        for(int i=0; i<B; i++){
            float V0 = V0_vect.get(i,0);
            float theta0 = theta0_vect.get(i,0);

            if(V0< 0 || V0 > 2){
                PyErr_SetString(PyExc_TypeError, "V0 must be in p.u (between 0 and 2)");
                return NULL;
            }
            if(theta0<-3.14159265359 || theta0>3.14159265359){
                std::cout << "[ERROR] : theta must be in radian" <<std::endl;
                return NULL;
            }
        }
        self->interfaceCase->setVoltageInit(V0_vect, theta0_vect);
    }
    
    //taille L : from, to, Ys Real, Ys Im, Yp, tau, theta, Limit=0, zs Real, zs Imag;
    static PyObject* StudyCase_setLink(CustomObject* self, PyObject* args){
        PyObject *pList;
        PyObject *pList2;
        Py_ssize_t l;

        if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &pList, &PyList_Type, &pList2)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be two lists.");
            return NULL;
        }

        
        l = PyList_Size(pList);
        int L = self->interfaceCase->getL();
        if(l != L){
            PyErr_SetString(PyExc_TypeError, "The size must be equal to L");
            return NULL;
        }
        if(l !=PyList_Size(pList2)){
            PyErr_SetString(PyExc_TypeError, "The size of the second list must be the same as the first");
            return NULL;
        }

        MatrixCPU From_vect = convertListToVectorCPUi(pList);
        MatrixCPU To_vect = convertListToVectorCPUi(pList2);
        int B = self->interfaceCase->getB();

        for(int i=0; i<L; i++){
            int from = From_vect.get(i,0);
            int to = To_vect.get(i,0);

            if(from < 0 || from > (B-1)){
                PyErr_SetString(PyExc_TypeError, "From must be a bus id, between 0 and B (exclude)");
                return NULL;
            }
            if(from < 0 || from > (B-1)){
                PyErr_SetString(PyExc_TypeError, "To must be a bus id, between 0 and B (exclude)");
                return NULL;
            }
        }
        self->interfaceCase->setLink(From_vect, To_vect);
    }
    static PyObject* StudyCase_setAdmitance(CustomObject* self, PyObject* args) {
        PyObject *pList;
        PyObject *pList2;
        PyObject *pList3;
        Py_ssize_t l;

        if (!PyArg_ParseTuple(args, "O!O!O!", &PyList_Type, &pList, &PyList_Type, &pList2, &PyList_Type, &pList3)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be three lists.");
            return NULL;
        }

        
        l = PyList_Size(pList);
        int L = self->interfaceCase->getL();
        if(l != L){
            PyErr_SetString(PyExc_TypeError, "The size must be equal to L");
            return NULL;
        }
        if(l !=PyList_Size(pList2)){
            PyErr_SetString(PyExc_TypeError, "The size of the second list must be the same as the first");
            return NULL;
        }
        if(l !=PyList_Size(pList3)){
            PyErr_SetString(PyExc_TypeError, "The size of the third list must be the same as the first");
            return NULL;
        }

        MatrixCPU YsRe_vect = convertListToVectorCPUf(pList);
        MatrixCPU YsIm_vect = convertListToVectorCPUf(pList2);
        MatrixCPU Yp_vect = convertListToVectorCPUf(pList3);
       
        self->interfaceCase->setAdmitance(YsRe_vect, YsIm_vect, Yp_vect);
    }
    static PyObject* StudyCase_setImpedance(CustomObject* self, PyObject* args){
        PyObject *pList;
        PyObject *pList2;
        Py_ssize_t l;

        if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &pList, &PyList_Type, &pList2)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be two lists.");
            return NULL;
        }

        
        l = PyList_Size(pList);
        int L = self->interfaceCase->getL();
        if(l != L){
            PyErr_SetString(PyExc_TypeError, "The size must be equal to L");
            return NULL;
        }
        if(l !=PyList_Size(pList2)){
            PyErr_SetString(PyExc_TypeError, "The size of the second list must be the same as the first");
            return NULL;
        }
        
        MatrixCPU ZsRe_vect = convertListToVectorCPUf(pList);
        MatrixCPU ZsIm_vect = convertListToVectorCPUf(pList2);
        
        self->interfaceCase->setImpedance(ZsRe_vect, ZsIm_vect);
    }
    static PyObject* StudyCase_setTransfo(CustomObject* self, PyObject* args){
        PyObject *pList;
        PyObject *pList2;
        Py_ssize_t l;

        if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &pList, &PyList_Type, &pList2)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be two lists.");
            return NULL;
        }

        
        l = PyList_Size(pList);
        int L = self->interfaceCase->getL();
        if(l != L){
            PyErr_SetString(PyExc_TypeError, "The size must be equal to L");
            return NULL;
        }
        if(l !=PyList_Size(pList2)){
            PyErr_SetString(PyExc_TypeError, "The size of the second list must be the same as the first");
            return NULL;
        }

        MatrixCPU Tau_vect = convertListToVectorCPUf(pList);
        MatrixCPU theta_vect = convertListToVectorCPUf(pList2);
        std::cout << "WIP : pas de verification des paramètres" <<std::endl;
               
        self->interfaceCase->setTransfo(Tau_vect, theta_vect);
    }
    static PyObject* StudyCase_setLineLimit(CustomObject* self, PyObject* args){
        PyObject *pList;
        Py_ssize_t l;

        if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &pList)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be one lists.");
            return NULL;
        }

        
        l = PyList_Size(pList);
        int L = self->interfaceCase->getL();
        if(l != L){
            PyErr_SetString(PyExc_TypeError, "The size must be equal to L");
            return NULL;
        }
    
        MatrixCPU limit_vect = convertListToVectorCPUf(pList);
        for(int i=0; i<L;i++){
            float limit = limit_vect.get(i,0);
            if(limit<0){
                PyErr_SetString(PyExc_TypeError, "Limits must be positve or null");
                return NULL;
            }
        }
          
        self->interfaceCase->setLineLimit(limit_vect);
    }
    // taille N*N : matrice de connexion
    static PyObject* StudyCase_setConnexion(CustomObject* self, PyObject* args){
        PyObject *pList;
        Py_ssize_t m;

        if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &pList)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be one lists.");
            return NULL;
        }

        
        m = PyList_Size(pList);
        int N = self->interfaceCase->getN();
        int M = N*N;
        if(m != M){
            PyErr_SetString(PyExc_TypeError, "The size must be equal to N*N");
            return NULL;
        }
    
        MatrixCPU connex_mat = convertListToMatrixCPUi(pList, N, N);
        for(int i=0; i<N;i++){
            for (int j = 0; j < N; j++)
            {
                float connexion = connex_mat.get(i,0);
                if(connexion!=1 && connexion!=0){
                    PyErr_SetString(PyExc_TypeError, "connexion is equal to 1 if connected or 0 if not ");
                    return NULL;
                }
            }
        } 
        self->interfaceCase->setConnexion(connex_mat);
    }
    // taille N*N (ou reduit) : matrices trade min et max
    static PyObject* StudyCase_setTradeLim(CustomObject* self, PyObject* args){
        PyObject *pList;
        PyObject *pList2;
        Py_ssize_t m;

        if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &pList, &PyList_Type, &pList2)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be two lists.");
            return NULL;
        }

        
        m = PyList_Size(pList);
        int N = self->interfaceCase->getN();
        int M = N*N;
        if(m != M){
            PyErr_SetString(PyExc_TypeError, "The size must be equal to N*N");
            return NULL;
        }
        if(m!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_TypeError, "The size of the second list must be the same as the first");
            return NULL;
        }
    
        MatrixCPU tradeMin_mat = convertListToMatrixCPUf(pList, N, N);
        MatrixCPU tradeMax_mat = convertListToMatrixCPUf(pList2, N, N);
        for(int i=0; i<N;i++){
            for (int j = 0; j < N; j++)
            {
                float tradeMin = tradeMin_mat.get(i,j);
                float tradeMax = tradeMax_mat.get(i,j);
                if(tradeMin > tradeMax){
                    PyErr_SetString(PyExc_TypeError, "tradeMax must be greater than tradeMin");
                    return NULL;
                }
            }
        } 
          
        self->interfaceCase->setTradeLim(tradeMin_mat, tradeMax_mat);
    }
    // taille B*B (ou reduit) : matrices impedance
    static PyObject* StudyCase_setMatImpedance(CustomObject* self, PyObject* args){
        PyObject *pList;
        PyObject *pList2;
        Py_ssize_t b2;

        if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &pList, &PyList_Type, &pList2)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be two lists.");
            return NULL;
        }

        
        b2 = PyList_Size(pList);
        int B = self->interfaceCase->getB();
        int B2 = B*B;
        if(b2 != B2){
            PyErr_SetString(PyExc_TypeError, "The size must be equal to B*B");
            return NULL;
        }
        if(b2!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_TypeError, "The size of the second list must be the same as the first");
            return NULL;
        }
    
        MatrixCPU Gs_mat = convertListToMatrixCPUf(pList, B, B);
        MatrixCPU Bs_mat = convertListToMatrixCPUf(pList2, B, B);
                  
        self->interfaceCase->setMatImpedance(Gs_mat, Bs_mat);
    }

    /*static PyObject* MonObject_getX(CustomObject* self, PyObject* args){
        int x = self->mavar->getX();
        return PyLong_FromLong((long) x);
    }

    static PyObject* MonObject_getY(CustomObject* self, PyObject* args){
        int y = self->mavar->getY();
        return PyLong_FromLong((long) y);
    }*/
}   



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
