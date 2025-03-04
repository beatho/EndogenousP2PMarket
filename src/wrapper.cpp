/*
Objectif realiser un wrapper entre un code python et un mon code en c++/cuda
build module : python .\setup.py build


entree python :
    - cas d'étude : nom de fichier
    - methode     : nom de la methode
    - paramètre   : iter^3, eps^3, step^3, rho^3

sortie python :
    - trade
    - puissance
    - residu
    - iter / temps

*/
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include <string>
#include "../head/System.h"
#include "../head/StudyCaseInterface.h"

PyObject* setStudyCaseFromFile(PyObject* self, PyObject* args){

    std::string filename; 
    char* buffer;
    PyArg_ParseTuple(args, "s", &buffer);

    filename = buffer;
    std::cout << "creation du cas " << filename <<std::endl; 

    System sys;
    sys.setStudyCase(filename);


    return Py_None;
}

PyObject* solveACFromFile(PyObject* self, PyObject* args){

    std::string filename; 
    std::string MethodName;
    char* bufferCas;
    char* bufferMethod;
    PyArg_ParseTuple(args, "ss", &bufferCas, &bufferMethod);

    filename = bufferCas;
    MethodName= bufferMethod;
    std::cout << "creation du cas " << filename <<std::endl; 
    std::cout << "Resolution avec la methode" << MethodName << std::endl;

    System sys;
    try
    {
        sys.setStudyCase(filename);
        try
        {
            sys.setMethod(MethodName);
            try
            {
                Simparam res = sys.solve();
                res.display(1); 
                res.display(2);
            }
            catch(const std::exception& e)
            {
                std::cout << "problème lors de la resolution" <<std::endl;
                std::cerr << e.what() << '\n';
            }
            
        }
        catch(const std::exception& e)
        {
            std::cout << "problème lors de la selection de la methode" <<std::endl;
            std::cerr << e.what() << '\n';
        }
        
    }
    catch(const std::exception& e)
    {
        std::cout << "problème lors de la selection du cas" <<std::endl;
        std::cerr << e.what() << '\n';
    }
    
    return Py_None;
}


PyObject* solvePFFromFile(PyObject* self, PyObject* args){
    std::string filename; 
    std::string MethodName;
    char* bufferCas;
    char* bufferMethod;
    bool computeInDouble;
    PyArg_ParseTuple(args, "ssp", &bufferCas, &bufferMethod, &computeInDouble);

    filename = bufferCas;
    MethodName= bufferMethod;
    std::cout << "creation du cas " << filename <<std::endl; 
    std::cout << "Resolution avec la methode" << MethodName << std::endl;
    std::cout << "Calcul en " << (computeInDouble? "double" : "simple") << " precision" << std::endl;

  
    System sys;
    try
    {
        sys.setStudyCase(filename);
        try
        {
            sys.setMethodPF(MethodName, computeInDouble);
            try
            {
                Simparam res = sys.solvePF();
                //res.display(1); 
                //res.display(2);
            }
            catch(const std::exception& e)
            {
                std::cout << "problème lors de la resolution" <<std::endl;
                std::cerr << e.what() << '\n';
            }
            
        }
        catch(const std::exception& e)
        {
            std::cout << "problème lors de la selection de la methode" <<std::endl;
            std::cerr << e.what() << '\n';
        }
        
    }
    catch(const std::exception& e)
    {
        std::cout << "problème lors de la selection du cas" <<std::endl;
        std::cerr << e.what() << '\n';
    }
    
    return Py_None;
}



PyMethodDef MonObject_Methods[] = {
    {"setSbase", (PyCFunction) StudyCase_setSbase, METH_VARARGS, "StudyCase set Sbase"},
    {"setVoltageInit", (PyCFunction) StudyCase_setVoltageInit, METH_VARARGS, "sets Voltage init"},
    {"setVbase", (PyCFunction) StudyCase_setVbase, METH_VARARGS, "sets Vbase"},
    {"setV0", (PyCFunction) StudyCase_setV0, METH_VARARGS, "sets V0"},
    {"setTheta", (PyCFunction) StudyCase_setTheta, METH_VARARGS, "sets Theta"},
    {"setName", (PyCFunction) StudyCase_setName, METH_VARARGS, "sets Name"},
    {"setPosBus", (PyCFunction) StudyCase_setPosBus, METH_VARARGS, "sets PosBus"},
    {"setCostFunction", (PyCFunction) StudyCase_setCostFunction, METH_VARARGS, "sets Cost Function"},
    {"setPower", (PyCFunction) StudyCase_setPower, METH_VARARGS, "set Power Pmin, Pmax, Pobj (opt)"},
    {"setImpedanceBus", (PyCFunction) StudyCase_setImpedanceBus, METH_VARARGS, "sets bus Impedance"},
    {"setVoltageBound", (PyCFunction) StudyCase_setVoltageBound, METH_VARARGS, "sets voltage bound"},
    {"setThetaBound", (PyCFunction) StudyCase_setThetaBound, METH_VARARGS, "set theta Bound"},
    {"setLink", (PyCFunction) StudyCase_setLink, METH_VARARGS, "sets Link"},
    {"setAdmitance", (PyCFunction) StudyCase_setAdmitance, METH_VARARGS, "sets Admitance"},
    {"setImpedance", (PyCFunction) StudyCase_setImpedance, METH_VARARGS, "sets Impedance"},
    {"setTransfo", (PyCFunction) StudyCase_setTransfo, METH_VARARGS, "sets Transfo"},
    {"setLineLimit", (PyCFunction) StudyCase_setLineLimit, METH_VARARGS, "sets LineLimit"},
    {"setConnexion", (PyCFunction) StudyCase_setConnexion, METH_VARARGS, "gets Connexion"},
    {"setTradeLim", (PyCFunction) StudyCase_setTradeLim, METH_VARARGS, "sets TradeLim"},
    {"setMatImpedance", (PyCFunction) StudyCase_setMatImpedance, METH_VARARGS, "sets Impedance Matrix"},
    {NULL, NULL, 0, NULL}
};


PyMethodDef EndoCudaFunction[] = {
    {
    "setStudyCaseByFile",
    setStudyCaseFromFile,
    METH_VARARGS,
    "this function create a study case"
    },
     {
    "solveACFromFile",
    solveACFromFile,
    METH_VARARGS,
    "this function solve a study case"
    },
    {
    "solvePFFromFile",
    solvePFFromFile,
    METH_VARARGS,
    "this function solve a Power flow on a study case"
    },
    {NULL, NULL, 0, NULL}
};

PyModuleDef EndoCudaModule = {
    PyModuleDef_HEAD_INIT,
    "EndoCuda",
    "le super module de bea",
    -1,
    EndoCudaFunction
};

PyMODINIT_FUNC PyInit_EndoCuda(){
    
    /* creation de la classe custom */
    CustomStudyCase.tp_name = "EndoCuda.interfaceStudyCase";
    CustomStudyCase.tp_basicsize = sizeof(CustomObject);
    CustomStudyCase.tp_itemsize = 0;
    CustomStudyCase.tp_dealloc = (destructor)MonObject_dealloc;
    CustomStudyCase.tp_flags = Py_TPFLAGS_DEFAULT;
    CustomStudyCase.tp_doc = PyDoc_STR("Custom objects");
    CustomStudyCase.tp_methods = MonObject_Methods;
    CustomStudyCase.tp_init = (initproc) MonObject_init;
    CustomStudyCase.tp_new = PyType_GenericNew;

    PyObject *m;
    if (PyType_Ready(&CustomStudyCase) < 0)
        return NULL;

    m = PyModule_Create(&EndoCudaModule);
    if (m == NULL)
        return NULL;

    if (PyModule_AddObjectRef(m, "interfaceStudyCase", (PyObject *) &CustomStudyCase) < 0) {
        Py_DECREF(m);
        return NULL;
    }

    return m;
}


