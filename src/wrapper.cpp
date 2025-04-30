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
#include "../head/interface.h"
#include "../head/main.h"


PyObject* callMainFunction(PyObject* self, PyObject* arg){
    main2();
    Py_INCREF(Py_None);
    return Py_None;
}

PyObject* setStudyCaseFromFile(PyObject* self, PyObject* args){

    std::string filename; 
    char* buffer;
    PyArg_ParseTuple(args, "s", &buffer);

    filename = buffer;
    std::cout << "creation du cas " << filename <<std::endl; 

    System sys;
    sys.setStudyCase(filename);

    Py_INCREF(Py_None);
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
    Py_INCREF(Py_None);
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
    
    Py_INCREF(Py_None);
    return Py_None;
}

PyObject* testAffichage(PyObject* self, PyObject* args){
    CustomObject* interface;

    if (!PyArg_ParseTuple(args, "O!", &CustomInterface, &interface)) {
        PyErr_SetString(PyExc_TypeError, "must be the interface");
        return NULL;
    } 
    interface->interfaceCase->display();
    
    Py_INCREF(Py_None);
    return Py_None;
}

PyObject* solveMarketFromInterface(PyObject* self, PyObject* args){

    CustomObject* interface;
    std::string MethodName;
    char* bufferMethod = nullptr;
    if (!PyArg_ParseTuple(args, "O!|s", &CustomInterface, &interface, &bufferMethod)) {
        PyErr_SetString(PyExc_TypeError, "must be the interface and the methode Name");
        return NULL;
    } 

    if(bufferMethod){
        MethodName= bufferMethod;
    } else {
        MethodName = "ADMM";
    }
    
    //std::cout << "Resolution avec la methode " << MethodName << std::endl;

    System sys;
    try
    {
        
        sys.setMethod(MethodName);
        try
        {
            interface->resultInterface = sys.solve(interface->resultInterface, interface->paramInterface, interface->interfaceCase, false);
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
        
    
    Py_INCREF(Py_None);
    return Py_None;
}

PyObject* solveACFromInterface(PyObject* self, PyObject* args){

    CustomObject* interface;
    std::string MethodName;
    char* bufferMethod = nullptr;
    if (!PyArg_ParseTuple(args, "O!|s", &CustomInterface, &interface, &bufferMethod)) {
        PyErr_SetString(PyExc_TypeError, "must be the interface and the methode Name");
        return NULL;
    } 

    if(bufferMethod){
        MethodName= bufferMethod;
    } else {
        MethodName = "ADMM";
    }
    
    std::cout << "Resolution avec la methode" << MethodName << std::endl;

    System sys;
    try
    {
        
        sys.setMethod(MethodName);
        try
        {
            interface->resultInterface = sys.solve(interface->resultInterface, interface->paramInterface, interface->interfaceCase, true);
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
        
    
    Py_INCREF(Py_None);
    return Py_None;
}

PyObject* solvePFFromInterface(PyObject* self, PyObject* args){

    CustomObject* interface;
    std::string MethodName;
    char* bufferMethod = nullptr;
    if (!PyArg_ParseTuple(args, "O!|s", &CustomInterface, &interface, &bufferMethod)) {
        PyErr_SetString(PyExc_TypeError, "must be the interface and the methode Name");
        return NULL;
    } 

    if(bufferMethod){
        MethodName= bufferMethod;
    } else {
        MethodName = "NR";
    }
    
    std::cout << "Resolution avec la methode" << MethodName << std::endl;

    System sys;
    try
    {
        
        sys.setMethodPF(MethodName, false);
        try
        {
            interface->resultInterface = sys.solvePF(interface->resultInterface, interface->paramInterface, interface->interfaceCase);
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
        
    
    Py_INCREF(Py_None);
    return Py_None;
}

PyMethodDef MonObject_Methods[] = {
    {"setSbase", (PyCFunction) Interface_setSbase, METH_VARARGS, " set Sbase"},
    {"setVbase", (PyCFunction) Interface_setVbase, METH_VARARGS, "sets Vbase"},
    {"display", (PyCFunction) Interface_display, METH_VARARGS, "display all data, >13 for help"},
    {"initProblem", (PyCFunction) Interface_initProblem, METH_VARARGS, "set trade, pn"},
    {"initDual", (PyCFunction) Interface_initDual, METH_VARARGS, "set lambda mu"},
    {"initDelta", (PyCFunction) Interface_initDelta, METH_VARARGS, "set delta"},
    {"setVoltageInit", (PyCFunction) StudyCase_setVoltageInit, METH_VARARGS, "sets Voltage init"},
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
    {"setConnexion", (PyCFunction) StudyCase_setConnexion, METH_VARARGS, "sets Connexion"},
    {"setBeta", (PyCFunction) StudyCase_setBeta, METH_VARARGS, "sets Beta"},
    {"setTradeLim", (PyCFunction) StudyCase_setTradeLim, METH_VARARGS, "sets TradeLim"},
    {"setMatImpedance", (PyCFunction) StudyCase_setMatImpedance, METH_VARARGS, "sets Impedance Matrix"},
    {"getInfo", (PyCFunction) StudyCase_getInfo, METH_VARARGS, "get Info"},
    {"getAgent", (PyCFunction) StudyCase_getAgent, METH_VARARGS, "gets Agent"},
    {"getBus", (PyCFunction) StudyCase_getBus, METH_VARARGS, "gets Bus"},
    {"getBranch", (PyCFunction) StudyCase_getBranch, METH_VARARGS, "gets Branch"},
    {"chekcase", (PyCFunction) StudyCase_checkCase, METH_VARARGS, "compute nCons nGen"},
    {"setIter", (PyCFunction) ParamInterface_setIter, METH_VARARGS, "sets iter"},
    {"setStep", (PyCFunction) ParamInterface_setStep, METH_VARARGS, "sets step"},
    {"setEps", (PyCFunction) ParamInterface_setEps, METH_VARARGS, "set eps"},
    {"setRho", (PyCFunction) ParamInterface_setRho, METH_VARARGS, "set rho"},
    {"getResults", (PyCFunction) ResultInterface_getResults, METH_VARARGS, "set eps"},
    {"getPn", (PyCFunction) ResultInterface_getPn, METH_VARARGS, "get Pn"},
    {"getLambda", (PyCFunction) ResultInterface_getLambda, METH_VARARGS, "get lambda"},
    {"getTrade", (PyCFunction) ResultInterface_getTrade, METH_VARARGS, "get trade"},
    {"getDelta", (PyCFunction) ResultInterface_getDelta, METH_VARARGS, "get delta"},
    {"getMu", (PyCFunction) ResultInterface_getMu, METH_VARARGS, "get Mu"},
    {"getPb", (PyCFunction) ResultInterface_getPb, METH_VARARGS, "get Pb"},
    {"getPhi", (PyCFunction) ResultInterface_getPhi, METH_VARARGS, "get Phi"},
    {"getE", (PyCFunction) ResultInterface_getE, METH_VARARGS, "get E"},
    {NULL, NULL, 0, NULL}
};

/**/

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
    {
        "solveMarketFromInterface",
        solveMarketFromInterface,
        METH_VARARGS,
        "this function test the interaction with custom Python object in c"
    },
    {
        "solveACFromInterface",
        solveACFromInterface,
        METH_VARARGS,
        "this function test the interaction with custom Python object in c"
    },
    {
        "solvePFFromInterface",
        solvePFFromInterface,
        METH_VARARGS,
        "this function test the interaction with custom Python object in c"
    },
    {
        "testAffichage",
        testAffichage,
        METH_VARARGS,
        "this function test the interaction with custom Python object in c"
    },
    {
        "callMainFunction",
        callMainFunction,
        METH_VARARGS,
        "this function call the main of the c program"
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
    CustomInterface.tp_name = "EndoCuda.interface";
    CustomInterface.tp_basicsize = sizeof(CustomObject);
    CustomInterface.tp_itemsize = 0;
    CustomInterface.tp_dealloc = (destructor)MonObject_dealloc;
    CustomInterface.tp_flags = Py_TPFLAGS_DEFAULT;
    CustomInterface.tp_doc = PyDoc_STR("Custom objects");
    CustomInterface.tp_methods = MonObject_Methods;
    CustomInterface.tp_init = (initproc) MonObject_init;
    CustomInterface.tp_new = PyType_GenericNew;

    PyObject *m;
    if (PyType_Ready(&CustomInterface) < 0)
        return NULL;

    m = PyModule_Create(&EndoCudaModule);
    if (m == NULL)
        return NULL;

    if (PyModule_AddObjectRef(m, "interface", (PyObject *) &CustomInterface) < 0) {
        Py_DECREF(m);
        return NULL;
    }

    return m;
}


