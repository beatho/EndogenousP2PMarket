#pragma once

#include <Python.h>
#include "MatrixCPU.h"
#include "StudyCaseInterface.h"
#include "ParamInterface.h"


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
    bool warning = false;
    for (int i=0; i <PyList_GET_SIZE(list); i++ ){
        PyObject* pItem = PyList_GetItem(list, i);
        if(PyLong_Check(pItem)){
            long Gs = PyLong_AsLong(pItem);
            toReturn.set(i, 0, (float) Gs);
        } 
        else if(PyFloat_Check(pItem)){
            double Gs = PyFloat_AsDouble(pItem);
            Gs = (int) Gs;
            warning = true;
            toReturn.set(i, 0, (float) Gs);
        }
        else{
            std::cout << "[WARNING] : unknow type, nothing to do ??? " << std::endl;
            //throw std::invalid_argument("Wrong type for input, must be integer");
        }
    }
    if(warning){
        std::cout << "[WARNING] : implicit conversion from float to int " << std::endl;
    }
    return toReturn;
}
static MatrixCPU convertListToMatrixCPUf(PyObject* list, int row, int collum){
    MatrixCPU toReturn(row,collum);
    int N = PyList_GET_SIZE(list);
    if(N != collum*row){
        PyErr_SetString(PyExc_ValueError, "Wrong size for the matrix from a list");
        return MatrixCPU(0,0);
    }
    int k = 0;
    for (int i=0; i <row; i++ ){
        for (int j = 0; j < collum; j++)
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
        PyErr_SetString(PyExc_ValueError, "Wrong size for the matrix from a list");
        return MatrixCPU(0,0);
    }
    int k = 0;
    bool warning = false;
    for (int i=0; i <row; i++ ){
        for (int j = 0; j < collum; j++)
        {
            PyObject* pItem = PyList_GetItem(list, k);
            if(PyLong_Check(pItem)){
                long Gs = PyLong_AsLong(pItem);
                toReturn.set(i, j, (float) Gs);
            } 
            else if(PyFloat_Check(pItem)){
                double Gs = PyFloat_AsDouble(pItem);
                Gs = (int) Gs;
                warning = true;
                toReturn.set(i, j, (float) Gs);
            }
            else{
                //throw std::invalid_argument("Wrong type for input, must be integer");
            }
            k++;
        }
    }
    if(warning){
        std::cout << "[WARNING] : implicit conversion from float to int " << std::endl;
    }
    return toReturn;
}

static PyObject* convertMatrixCPUtoList(MatrixCPU mat){
    int row = mat.getNLin();
    int collumn = mat.getNCol();
    PyObject* lst = PyList_New(row*collumn);
    int k = 0;
    for(int i=0; i<row; i++){
        for(int j=0; j<collumn; j++){
            double value = mat.get(i,j);
            PyList_SET_ITEM(lst, k, PyFloat_FromDouble(value));
            k++;
        }
    }
    return lst;
}
//static PyObject* convertVectorCPUtoList(MatrixCPU mat);

extern "C"{

    typedef struct {
        PyObject_HEAD
        StudyCaseInterface* interfaceCase = nullptr;
        ParamInterface*    paramInterface = nullptr;
        ResultInterface*  resultInterface = nullptr;
        /* Type-specific fields go here. */
    } CustomObject;

    static PyTypeObject CustomInterface;

    static int MonObject_init(CustomObject* self, PyObject* args, PyObject* kwds){
        //std::cout << "creation d'un interface cas d'etude" << std::endl;
        DELETEB(self->interfaceCase)
        DELETEB(self->paramInterface)
        DELETEB(self->resultInterface)
        int N, B, L;
        int Lconst = 0;
        if(!PyArg_ParseTuple(args, "iii|i", &N, &B, &L, &Lconst)){
            PyErr_SetString(PyExc_TypeError, "[Error] : N, B, L are sizes, must be 3 or 4 integers");
            return -1;
        }
        if(N<0 || B<0 || L<0 || Lconst < 0){
            PyErr_SetString(PyExc_ValueError, "[Error] : N, B, L are sizes, must be positive");
            return -1;
        }
        self->interfaceCase = new StudyCaseInterface(N, B, L);
        int Lparam = L;
        if(Lconst){
            Lparam = Lconst;
        }
        self->paramInterface  = new ParamInterface(N, B, L, Lparam);
        self->resultInterface = new ResultInterface(N, B, L, Lparam);
        return 0;
    }
    static void MonObject_dealloc(CustomObject* self){
        DELETEB(self->interfaceCase)
        DELETEB(self->paramInterface)
        DELETEB(self->resultInterface)
    }

    static PyObject* Interface_setSbase(CustomObject* self, PyObject* args){
        float Sbase;
        if(!PyArg_ParseTuple(args, "f", &Sbase)){
            PyErr_SetString(PyExc_TypeError, "parameter must be a float.");
            return NULL;
        }
        if(Sbase<=0){
            PyErr_SetString(PyExc_ValueError, "Sbase must be positive");
            return NULL;
        }
        self->interfaceCase->setSbase(Sbase);
        self->paramInterface->setSbase(Sbase);
        Py_IncRef(Py_None);
        return Py_None;
    }
    static PyObject* Interface_setVbase(CustomObject* self, PyObject* args){
        float Vbase;
        if(!PyArg_ParseTuple(args, "f", &Vbase)){
            PyErr_SetString(PyExc_TypeError, "Vbase must be a float");
            return NULL;
        }
        if(Vbase<=0){
            PyErr_SetString(PyExc_ValueError, "Vbase must be positive");
            return NULL;
        }
        self->interfaceCase->setVbase(Vbase);
        self->paramInterface->setVbase(Vbase);
        Py_IncRef(Py_None);
        return Py_None;
    }
    static PyObject* Interface_display(CustomObject* self, PyObject* args){
        int type = 0;
        if(!PyArg_ParseTuple(args, "|i", &type)){
            PyErr_SetString(PyExc_TypeError, "parameter must be none or an integer");
            return NULL;
        }
        
        int offsetCase  = 1;
        int offsetParam = 6;
        int offsetRes = 9;
        int offsetAll = 13;
        std::cout << std::endl;
        if(type==0){
            self->interfaceCase->display(type);
            self->paramInterface->display(type);
            self->resultInterface->display(self->interfaceCase,type);
        }if(type < offsetParam){
            self->interfaceCase->display(type - offsetCase);
        } else if(type < offsetRes){
            self->paramInterface->display(type - offsetParam);
        } else if(type < offsetAll){
            self->resultInterface->display(self->interfaceCase, type-offsetRes);
        } else{
            std::cout << " command for display : " << std::endl;
            std::cout << " 0  : all without details" << std::endl;
            std::cout << " 1  : Study Case, All    info" << std::endl;
            std::cout << " 2  : Study Case, Case   info" << std::endl;
            std::cout << " 3  : Study Case, Agent  info" << std::endl;
            std::cout << " 4  : Study Case, Bus    info" << std::endl;
            std::cout << " 5  : Study Case, Branch info" << std::endl;
            std::cout << " 6  : Param, All info" << std::endl;
            std::cout << " 7  : Param, parameters details" << std::endl;
            std::cout << " 8  : Param, using matrix" << std::endl;
            std::cout << " 9  : Results, using matrix" << std::endl;
            std::cout << " 10 : Results, for market" << std::endl;
            std::cout << " 11 : Results, for Power Flow " << std::endl;
            std::cout << " 12 : Results, for EndoMarket or OPF" << std::endl;
            std::cout << ">12 : This print" << std::endl;
        }
        
        Py_IncRef(Py_None);
        return Py_None;
    }
      
    static PyObject* Interface_initProblem(CustomObject* self, PyObject* args){
        PyObject *pList;
        PyObject *pList2;
        Py_ssize_t n2;

        if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &pList, &PyList_Type, &pList2)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be two lists.");
            return NULL;
        }
        n2 = PyList_Size(pList);
        int N = self->interfaceCase->getN();
        int N2 = N*N;
        int N3 = 2*N*N;
        if(n2 != N2 && n2 != N3 ){
            PyErr_SetString(PyExc_ValueError, "The size must be equal to N*N or 2N*N");
            return NULL;
        }
        if(n2/N!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_ValueError, "The size of the second list must be the first divided by N");
            return NULL;
        }
        if(n2==N2){
            MatrixCPU trade_mat = MatrixCPU(2*(N + 1), N + 1);
            MatrixCPU pn_mat = MatrixCPU(2*(N + 1), 1);
            MatrixCPU trade_temp = convertListToMatrixCPUf(pList, N, N);
            MatrixCPU pn_temp = convertListToVectorCPUf(pList2);
            for(int i=0; i <N; i++){
                for(int j=0; j<N; j++){
                    trade_mat.set(i + 1, j + 1, trade_temp.get(i,j));
                }
                pn_mat.set(i + 1, 0, pn_temp.get(i,0));
            }

            self->paramInterface->initProbleme(trade_mat, pn_mat);
            self->resultInterface->setProbleme(trade_mat, pn_mat);
        } else{
            MatrixCPU trade_mat = MatrixCPU(2*(N + 1), N + 1);
            MatrixCPU pn_mat = MatrixCPU(2*(N + 1), 1);
            MatrixCPU trade_temp = convertListToMatrixCPUf(pList, 2*N, N);
            MatrixCPU pn_temp = convertListToVectorCPUf(pList2);
            for(int i=0; i <N; i++){
                for(int j=0; j<N; j++){
                    trade_mat.set(    i + 1, j + 1, trade_temp.get(    i, j));
                    trade_mat.set(N + i + 2, j + 1, trade_temp.get(N + i, j));
                }
                pn_mat.set(    i + 1, 0, pn_temp.get(i,0));
                pn_mat.set(N + i + 2, 0, pn_temp.get(N + i,0));
            }
            self->paramInterface->initProbleme(trade_mat, pn_mat);
            self->resultInterface->setProbleme(trade_mat, pn_mat);
        }

        Py_IncRef(Py_None);
        return Py_None;
    }
    static PyObject* Interface_initDual(CustomObject* self, PyObject* args){
        PyObject *pList;
        PyObject *pList2;
        Py_ssize_t n2;

        if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &pList, &PyList_Type, &pList2)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be two lists.");
            return NULL;
        }
        n2 = PyList_Size(pList);
        int N = self->interfaceCase->getN();
        int N2 = N*N;
        int N3 = 2*N*N;
        if(n2 != N2 && n2 != N3 ){
            PyErr_SetString(PyExc_ValueError, "The size must be equal to N*N or 2N*N");
            return NULL;
        }
        if(n2/N!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_ValueError, "The size of the second list must be the first divided by N");
            return NULL;
        }
        if(n2==N2){
            MatrixCPU lambda_mat = MatrixCPU(2*(N+1), N + 1);
            MatrixCPU mu_mat = MatrixCPU(2* (N + 1), 1);
            MatrixCPU lambda_temp = convertListToMatrixCPUf(pList, N, N);
            MatrixCPU mu_temp = convertListToVectorCPUf(pList2);
            for(int i=0; i <N; i++){
                for(int j=0; j<N; j++){
                    lambda_mat.set(i + 1,j + 1, lambda_temp.get(i,j));
                }
                mu_mat.set(i + 1, 0, mu_temp.get(i,0));
            }

            self->paramInterface->initDual(lambda_mat, mu_mat);
            self->resultInterface->setDual(lambda_mat, mu_mat);
        } else{
            MatrixCPU lambda_mat = MatrixCPU(2*(N + 1), N + 1);
            MatrixCPU mu_mat = MatrixCPU(2*(N + 1), 1);
            MatrixCPU lambda_temp = convertListToMatrixCPUf(pList, 2*N, N);
            MatrixCPU mu_temp     = convertListToVectorCPUf(pList2);
            for(int i=0; i <N; i++){
                for(int j=0; j<N; j++){
                    lambda_mat.set(    i + 1, j + 1, lambda_temp.get(    i, j));
                    lambda_mat.set(N + i + 2, j + 1, lambda_temp.get(N + i, j));
                }
                mu_mat.set(    i + 1, 0, mu_temp.get(i,0));
                mu_mat.set(N + i + 2, 0, mu_temp.get(N + i,0));
            }
            //std::cout<< "affichage dans interface taille 2*N" <<std::endl;
            //lambda_mat.display();
            self->paramInterface->initDual(lambda_mat, mu_mat);
            self->resultInterface->setDual(lambda_mat, mu_mat);
        }
        
        Py_IncRef(Py_None);
        return Py_None;
    }
    static PyObject* Interface_initDelta(CustomObject* self, PyObject* args){
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
            PyErr_SetString(PyExc_ValueError, "The size must be equal to L");
            return NULL;
        }
        if(l !=PyList_Size(pList2)){
            PyErr_SetString(PyExc_ValueError, "The size of the second list must be the first divided");
            return NULL;
        }
        

        MatrixCPU delta1 = convertListToVectorCPUf(pList);
        MatrixCPU delta2 = convertListToVectorCPUf(pList2);
        self->paramInterface->initDelta(delta1, delta2);
        self->resultInterface->setDelta(delta1, delta2);
        
        
        Py_IncRef(Py_None);
        return Py_None;
    }


    
    /* ******************** Study Case*****************************     */
    static PyObject* StudyCase_setV0(CustomObject* self, PyObject* args){
        float V0;
        if(!PyArg_ParseTuple(args, "f", &V0)){
            PyErr_SetString(PyExc_TypeError, "V0 must be a float");
            return NULL;
        }
        if(V0<=0){
            PyErr_SetString(PyExc_ValueError, "V0 must be positive");
            return NULL;
        }
        self->interfaceCase->setV0(V0);
        Py_IncRef(Py_None);
        return Py_None;
    }
    static PyObject* StudyCase_setTheta(CustomObject* self, PyObject* args){
        float theta0;
        if(!PyArg_ParseTuple(args, "f", &theta0)){
            PyErr_SetString(PyExc_TypeError, "theta must be a float");
            return NULL;
        }
        if(theta0<-3.14159265359 || theta0>3.14159265359){
            PyErr_SetString(PyExc_ValueError, "theta0 must be in radian");
            return NULL;
        }
        self->interfaceCase->setTheta(theta0);
        Py_IncRef(Py_None);
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
        Py_IncRef(Py_None);
        return Py_None;
    }
    // taille N : PosBus, a, b, a^q, b^q Pobj, Pmin, Pmax, Qobj, Qmin, Qmax, zone 
    static PyObject* StudyCase_setPosBus(CustomObject* self, PyObject* args){
        PyObject *pList;
        PyObject *pList2 = nullptr;
        Py_ssize_t n;

        if (!PyArg_ParseTuple(args, "O!|O!", &PyList_Type, &pList, &PyList_Type, &pList2)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be one or two list.");
            return NULL;
        }

        n = PyList_Size(pList);
        if(n != self->interfaceCase->getN()){
            PyErr_SetString(PyExc_ValueError, "The size must be equal to N");
            return NULL;
        }
        if(pList2!=nullptr){
            if(n != PyList_Size(pList2)){
                PyErr_SetString(PyExc_ValueError, "The size of the second list must be the same as the first");
                return NULL;
            }
        }
        
        MatrixCPU PosBus = convertListToVectorCPUi(pList);
        MatrixCPU zone(n, 1);
        if(pList2!=nullptr){
            zone = convertListToVectorCPUi(pList2);   
        } 
        for (int i=0; i<n; i++) {
            int bus = PosBus.get(i,0);
            if(bus<-1 || bus >self->interfaceCase->getB()){
                PyErr_SetString(PyExc_ValueError, "Position must be a int between -1 (not on the grid) and B (exclude)");
                return NULL;
            }
        }
        self->interfaceCase->setPosBus(PosBus, zone);
        Py_IncRef(Py_None);
        return Py_None;
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
            PyErr_SetString(PyExc_ValueError, "The size must be equal to N or to 2*N");
            return NULL;
        }
        if(n!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_ValueError, "The size of the second list must be the same as the first");
            return NULL;
        }

        MatrixCPU temp_a = convertListToVectorCPUf(pList);
        MatrixCPU temp_b = convertListToVectorCPUf(pList2);

        MatrixCPU a_vect(2*N, 1);
        MatrixCPU b_vect(2*N, 1);

        for (int i=0; i<n; i++) {
            float a = temp_a.get(i,0);
            if(a<=0){
                PyErr_SetString(PyExc_ValueError, "a must be positive");
                return NULL;
            }
            a_vect.set(i, 0, a);
            b_vect.set(i, 0, temp_b.get(i,0));
        }
        self->interfaceCase->setCostFunction(a_vect, b_vect);
        
        Py_IncRef(Py_None);
        return Py_None;
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
            PyErr_SetString(PyExc_ValueError, "The size must be equal to N or to 2*N");
            return NULL;
        } if(n!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_ValueError, "The size of the second list must be the same as the first");
            return NULL;
        } if(pList3!=nullptr){
            if(n!=PyList_Size(pList3)){
                PyErr_SetString(PyExc_ValueError, "The size of the thrid list must be the same as the first");
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
                PyErr_SetString(PyExc_ValueError, "Pmax must be greater than Pmin");
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
        Py_IncRef(Py_None);
        return Py_None;
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
            PyErr_SetString(PyExc_ValueError, "The size must be equal to B");
            return NULL;
        }
        if(b!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_ValueError, "The size of the second list must be the same as the first");
            return NULL;
        }

        MatrixCPU Gs_vect = convertListToVectorCPUf(pList);
        MatrixCPU Bs_vect = convertListToVectorCPUf(pList2);

        self->interfaceCase->setImpedanceBus(Gs_vect, Bs_vect);
        Py_IncRef(Py_None);
        return Py_None;
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
            PyErr_SetString(PyExc_ValueError, "The size must be equal to B");
            return NULL;
        }
        if(b!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_ValueError, "The size of the second list must be the same as the first");
            return NULL;
        }

        MatrixCPU Vmin_vect = convertListToVectorCPUf(pList);
        MatrixCPU Vmax_vect = convertListToVectorCPUf(pList2);

        for(int i=0; i<B; i++){
            float Vmin = Vmin_vect.get(i,0);
            float Vmax = Vmax_vect.get(i,0);

            if(Vmin>Vmax){
                PyErr_SetString(PyExc_ValueError, "Vmax must be greater than Vmin");
            }

        }
        self->interfaceCase->setVoltageBound(Vmin_vect, Vmax_vect);
        Py_IncRef(Py_None);
        return Py_None;
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
            PyErr_SetString(PyExc_ValueError, "The size must be equal to B");
            return NULL;
        }
        if(b!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_ValueError, "The size of the second list must be the same as the first");
            return NULL;
        }

        MatrixCPU thetamin_vect = convertListToVectorCPUf(pList);
        MatrixCPU thetamax_vect = convertListToVectorCPUf(pList2);

        for(int i=0; i<B; i++){
            float thetamin = thetamin_vect.get(i,0);
            float thetamax = thetamax_vect.get(i,0);

            if(thetamin>thetamax){
                PyErr_SetString(PyExc_ValueError, "thetamax must be greater than thetamin");
                return NULL;
            }
            if(thetamin<-3.14159265359 || thetamax>3.14159265359){
                PyErr_SetString(PyExc_ValueError, " theta must be in radian");
                return NULL;
            }
        }
        self->interfaceCase->setAngleBound(thetamin_vect, thetamax_vect);
        Py_IncRef(Py_None);
        return Py_None;
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
            PyErr_SetString(PyExc_ValueError, "The size must be equal to B");
            return NULL;
        }
        if(b!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_ValueError, "The size of the second list must be the same as the first");
            return NULL;
        }

        MatrixCPU V0_vect = convertListToVectorCPUf(pList);
        MatrixCPU theta0_vect = convertListToVectorCPUf(pList2);

        for(int i=0; i<B; i++){
            float V0 = V0_vect.get(i,0);
            float theta0 = theta0_vect.get(i,0);

            if(V0< 0 || V0 > 2){
                PyErr_SetString(PyExc_ValueError, "V0 must be in p.u (between 0 and 2)");
                return NULL;
            }
            if(theta0<-3.14159265359 || theta0>3.14159265359){
                PyErr_SetString(PyExc_ValueError, "theta must be in radian");
                return NULL;
            }
        }
        self->interfaceCase->setVoltageInit(V0_vect, theta0_vect);
        Py_IncRef(Py_None);
        return Py_None;
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
            PyErr_SetString(PyExc_ValueError, "The size must be equal to L");
            return NULL;
        }
        if(l !=PyList_Size(pList2)){
            PyErr_SetString(PyExc_ValueError, "The size of the second list must be the same as the first");
            return NULL;
        }

        MatrixCPU From_vect = convertListToVectorCPUi(pList);
        MatrixCPU To_vect = convertListToVectorCPUi(pList2);
        int B = self->interfaceCase->getB();

        for(int i=0; i<L; i++){
            int from = From_vect.get(i,0);
            int to = To_vect.get(i,0);

            if(from < 0 || from > (B-1)){
                PyErr_SetString(PyExc_ValueError, "From must be a bus id, between 0 and B (exclude)");
                return NULL;
            }
            if(from < 0 || from > (B-1)){
                PyErr_SetString(PyExc_ValueError, "To must be a bus id, between 0 and B (exclude)");
                return NULL;
            }
        }
        self->interfaceCase->setLink(From_vect, To_vect);
        Py_IncRef(Py_None);
        return Py_None;
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
            PyErr_SetString(PyExc_ValueError, "The size must be equal to L");
            return NULL;
        }
        if(l !=PyList_Size(pList2)){
            PyErr_SetString(PyExc_ValueError, "The size of the second list must be the same as the first");
            return NULL;
        }
        if(l !=PyList_Size(pList3)){
            PyErr_SetString(PyExc_ValueError, "The size of the third list must be the same as the first");
            return NULL;
        }

        MatrixCPU YsRe_vect = convertListToVectorCPUf(pList);
        MatrixCPU YsIm_vect = convertListToVectorCPUf(pList2);
        MatrixCPU Yp_vect = convertListToVectorCPUf(pList3);
       
        self->interfaceCase->setAdmitance(YsRe_vect, YsIm_vect, Yp_vect);
        Py_IncRef(Py_None);
        return Py_None;
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
            PyErr_SetString(PyExc_ValueError, "The size must be equal to L");
            return NULL;
        }
        if(l !=PyList_Size(pList2)){
            PyErr_SetString(PyExc_ValueError, "The size of the second list must be the same as the first");
            return NULL;
        }
        
        MatrixCPU ZsRe_vect = convertListToVectorCPUf(pList);
        MatrixCPU ZsIm_vect = convertListToVectorCPUf(pList2);
        
        self->interfaceCase->setImpedance(ZsRe_vect, ZsIm_vect);
        Py_IncRef(Py_None);
        return Py_None;
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
            PyErr_SetString(PyExc_ValueError, "The size must be equal to L");
            return NULL;
        }
        if(l !=PyList_Size(pList2)){
            PyErr_SetString(PyExc_ValueError, "The size of the second list must be the same as the first");
            return NULL;
        }

        MatrixCPU Tau_vect = convertListToVectorCPUf(pList);
        MatrixCPU theta_vect = convertListToVectorCPUf(pList2);
        std::cout << "WIP : pas de verification des paramètres" <<std::endl;
               
        self->interfaceCase->setTransfo(Tau_vect, theta_vect);
        Py_IncRef(Py_None);
        return Py_None;
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
            PyErr_SetString(PyExc_ValueError, "The size must be equal to L");
            return NULL;
        }
    
        MatrixCPU limit_vect = convertListToVectorCPUf(pList);
        for(int i=0; i<L;i++){
            float limit = limit_vect.get(i,0);
            if(limit<0){
                PyErr_SetString(PyExc_ValueError, "Limits must be positve or null");
                return NULL;
            }
        }
          
        self->interfaceCase->setLineLimit(limit_vect);
        Py_IncRef(Py_None);
        return Py_None;
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
            PyErr_SetString(PyExc_ValueError, "The size must be equal to N*N");
            return NULL;
        }
    
        MatrixCPU connex_mat = convertListToMatrixCPUi(pList, N, N);
        for(int i=0; i<N;i++){
            for (int j = 0; j < N; j++)
            {
                float connexion = connex_mat.get(i,0);
                if(connexion!=1 && connexion!=0){
                    PyErr_SetString(PyExc_ValueError, "connexion is equal to 1 if connected or 0 if not ");
                    return NULL;
                }
            }
        } 
        self->interfaceCase->setConnexion(connex_mat);
        Py_IncRef(Py_None);
        return Py_None;
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
            PyErr_SetString(PyExc_ValueError, "The size must be equal to N*N");
            return NULL;
        }
        if(m!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_ValueError, "The size of the second list must be the same as the first");
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
                    PyErr_SetString(PyExc_ValueError, "tradeMax must be greater than tradeMin");
                    return NULL;
                }
            }
        } 
          
        self->interfaceCase->setTradeLim(tradeMin_mat, tradeMax_mat);
        Py_IncRef(Py_None);
        return Py_None;
    }
    
    static PyObject* StudyCase_setBeta(CustomObject* self, PyObject* args){
        PyObject *pList;
        Py_ssize_t m;

        if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &pList)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be one list.");
            return NULL;
        }

        
        m = PyList_Size(pList);
        int N = self->interfaceCase->getN();
        int M = N*N;
        if(m != M){
            PyErr_SetString(PyExc_ValueError, "The size must be equal to N*N");
            return NULL;
        }
    
        MatrixCPU beta_mat = convertListToMatrixCPUf(pList, N, N);
        self->interfaceCase->setBeta(beta_mat);
        Py_IncRef(Py_None);
        return Py_None;
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
            PyErr_SetString(PyExc_ValueError, "The size must be equal to B*B");
            return NULL;
        }
        if(b2!=PyList_Size(pList2)){
            PyErr_SetString(PyExc_ValueError, "The size of the second list must be the same as the first");
            return NULL;
        }
    
        MatrixCPU Gs_mat = convertListToMatrixCPUf(pList, B, B);
        MatrixCPU Bs_mat = convertListToMatrixCPUf(pList2, B, B);
                  
        self->interfaceCase->setMatImpedance(Gs_mat, Bs_mat);
        Py_IncRef(Py_None);
        return Py_None;
    }



    /* MatrixCPU infoCase;
        MatrixCPU agentCase;
        MatrixCPU branchCase;
        MatrixCPU busCase;*/
    static PyObject* StudyCase_getInfo(CustomObject* self, PyObject* arg){
        MatrixCPU info = self->interfaceCase->getInfoCase();
        return convertMatrixCPUtoList(info);
    }
    static PyObject* StudyCase_getAgent(CustomObject* self, PyObject* arg){
        MatrixCPU Agent = self->interfaceCase->getAgentCase();
        return convertMatrixCPUtoList(Agent);
    }
    static PyObject* StudyCase_getBranch(CustomObject* self, PyObject* arg){
        MatrixCPU Branch = self->interfaceCase->getBranchCase();
        return convertMatrixCPUtoList(Branch);
    }
    static PyObject* StudyCase_getBus(CustomObject* self, PyObject* arg){
        MatrixCPU Bus = self->interfaceCase->getBusCase();
        return convertMatrixCPUtoList(Bus);
    }

    
    static PyObject* StudyCase_checkCase(CustomObject* self, PyObject* args){
        self->interfaceCase->checkCase();
        Py_IncRef(Py_None);
        return Py_None;
    }
    /*static PyObject* MonObject_getX(CustomObject* self, PyObject* args){
        int x = self->mavar->getX();
        return PyLong_FromLong((long) x);
    }

    static PyObject* MonObject_getY(CustomObject* self, PyObject* args){
        int y = self->mavar->getY();
        return PyLong_FromLong((long) y);
    }*/

     /* ******************** Param *****************************     */
    /*
Pour créer les paramètres :
  - taille 1 : - IterG, iterL, iterIntern, ,
               - stepG, stepL, stepIntern,
               - epsG, epsL, epsX, epsInter, 
               - rho, rhoL, rho1
               - Sbase, Vbase
               - N, B, L, 
               - iterFinal, temps, fc
  - taille N*N : lambda, trade,
  - taille N : Pn, MU
  - taille L : _delta1, _delta2
  - taille iter : _resF
    */
    static PyObject* ParamInterface_setIter(CustomObject* self, PyObject* args){
        int iterG;
        int iterL = 0;
        int iterIntern = 0;
        if(!PyArg_ParseTuple(args, "i|ii", &iterG, &iterL, &iterIntern)){
            PyErr_SetString(PyExc_TypeError, "parameter must be one, two or three int.");
            return NULL;
        }
        if(iterG<0 || iterL < 0 || iterIntern < 0){
            PyErr_SetString(PyExc_ValueError, "iter must be positive or null for default value");
            return NULL;
        }
        self->paramInterface->setIter(iterG, iterL, iterIntern);
        self->resultInterface->changeIterStep(iterG, 0);
        Py_IncRef(Py_None);
        return Py_None;
    }
    static PyObject* ParamInterface_setStep(CustomObject* self, PyObject* args){
        int stepG;
        int stepL = 0;
        int stepIntern = 0;
        if(!PyArg_ParseTuple(args, "i|ii", &stepG, &stepL, &stepIntern)){
            PyErr_SetString(PyExc_TypeError, "parameter must be one, two or three int.");
            return NULL;
        }
        if(stepG<0 || stepL < 0 || stepIntern < 0){
            PyErr_SetString(PyExc_ValueError, "step must be positive or null for default value");
            return NULL;
        }
        self->paramInterface->setStep(stepG, stepL, stepIntern);
        self->resultInterface->changeIterStep(0, stepG);
        Py_IncRef(Py_None);
        return Py_None;
    }
    static PyObject* ParamInterface_setEps(CustomObject* self, PyObject* args){
        float epsG;
        float epsL;
        float epsX = 0;
        float epsIntern = 0;
        if(!PyArg_ParseTuple(args, "ff|ff", &epsG, &epsL, &epsX, &epsIntern)){
            PyErr_SetString(PyExc_TypeError, "parameter must be two, three, four float.");
            return NULL;
        }
        if(epsG<0 || epsL < 0 || epsX < 0 || epsIntern){
            PyErr_SetString(PyExc_ValueError, "eps must be positive or null for default value");
            return NULL;
        }
        self->paramInterface->setEps(epsG, epsL, epsX, epsIntern);

        Py_IncRef(Py_None);
        return Py_None;
    }
    static PyObject* ParamInterface_setRho(CustomObject* self, PyObject* args){
        float rhoG;
        float rhoL = 0;
        float rhoX = 0;
        
        if(!PyArg_ParseTuple(args, "ff|f", &rhoG, &rhoL, &rhoX)){
            PyErr_SetString(PyExc_TypeError, "parameter must be two, three, four float.");
            return NULL;
        }
        if(rhoG < 0 || rhoL < 0 || rhoX < 0){
            PyErr_SetString(PyExc_ValueError, "rho must be positive or null for default value");
            return NULL;
        }
        self->paramInterface->setRho(rhoG, rhoL, rhoX);

        Py_IncRef(Py_None);
        return Py_None;
    }
   
   
    static PyObject* ResultInterface_getResults(CustomObject* self, PyObject* args){
        MatrixCPU res = self->resultInterface->getResults();
        return convertMatrixCPUtoList(res);
    }
    static PyObject* ResultInterface_getPn(CustomObject* self, PyObject* args){
        MatrixCPU pn = self->resultInterface->getPn();
        return convertMatrixCPUtoList(pn);
    }
    static PyObject* ResultInterface_getLambda(CustomObject* self, PyObject* args){
        MatrixCPU lambda = self->resultInterface->getLambda();
        return convertMatrixCPUtoList(lambda);
    }
    static PyObject* ResultInterface_getTrade(CustomObject* self, PyObject* args){
        MatrixCPU trade = self->resultInterface->getTrade();
        return convertMatrixCPUtoList(trade);
    }
    static PyObject* ResultInterface_getDelta(CustomObject* self, PyObject* args){
        MatrixCPU delta = self->resultInterface->getDelta();
        return convertMatrixCPUtoList(delta);
    }
    static PyObject* ResultInterface_getMu(CustomObject* self, PyObject* args){
        MatrixCPU mu = self->resultInterface->getMU();
        return convertMatrixCPUtoList(mu);
    }
    static PyObject* ResultInterface_getPb(CustomObject* self, PyObject* args){
        MatrixCPU Pb = self->resultInterface->getPb();
        return convertMatrixCPUtoList(Pb);
    }
    static PyObject* ResultInterface_getPhi(CustomObject* self, PyObject* args){
        MatrixCPU Phi = self->resultInterface->getPhi();
        return convertMatrixCPUtoList(Phi);
    }
    static PyObject* ResultInterface_getE(CustomObject* self, PyObject* args){
        MatrixCPU E = self->resultInterface->getE();
        return convertMatrixCPUtoList(E);
    }

}   
