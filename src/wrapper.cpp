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

PyObject* setStudyCase(PyObject* self, PyObject* args){

    std::string filename; 
    char* buffer;
    PyArg_ParseTuple(args, "s", &buffer);




    return Py_None;
}








PyMethodDef EndoCudaFunction[] = {
    {
    "setStudyCaseByFile",
    setStudyCase,
    METH_VARARGS,
    "this function create a study case"
    },
    {NULL, NULL, 0, NULL}
};

PyModuleDef EndoCudaModule = {
    PyModuleDef_HEAD_INIT,
    "verlet",
    "le super module du cours de C",
    -1,
    verletFunctions
};

PyMODINIT_FUNC PyInit_EndoCuda(){
    PyObject* module;
    
    module = PyModule_Create(&EndoCudaModule);
    return module;
}


