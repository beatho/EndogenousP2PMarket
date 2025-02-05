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
    PyObject* module;
    
    module = PyModule_Create(&EndoCudaModule);
    return module;
}


