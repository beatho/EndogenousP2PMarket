#pragma once


#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "MatrixCPU.h"
#include "MatrixCPUD.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <string>
#include "StudyCaseInterface.h"



#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}

/* existe mais avec impedance infini
_nBus = 2;
    _nLine = 1;
    _CoresBusAgent = MatrixCPU(_nBus, _nAgent);

    _LineImpedance = MatrixCPU(_nLine, _nLine);
    _CoresBusLine = MatrixCPU(_nBus, _nLine);
    _SensiPower = MatrixCPU(_nLine, _nAgent); // G = 0 !!!!
    _lineLimits = MatrixCPU(_nLine, 1, 0);
    _nLineConstraint = _nLine;
    _SensiPowerReduce = MatrixCPU(_SensiPower);
    _lineLimitsReduce = MatrixCPU(_lineLimits);

*/

enum RadialType { Normal, Line, Balance, OneStep };

class StudyCaseDCGrid
{
protected :
    static constexpr float LINELIMITMAX = 10000000;
   

    int _nBus = 0;
    int _nLine = 0;
    int _nLineConstraint = 0;
   
    std::string _name = "None";
    
    
    float _timeInit;
    

    MatrixCPU _LineImpedance; // B
    MatrixCPU _CoresBusLine; // C
    MatrixCPU _CoresLineBus; // Perso from, to
    MatrixCPU _CoresLineBusReduce; // Perso

    MatrixCPU _SensiBusLine; // A
    MatrixCPU _SensiBusLineReduce; // Areduce

    MatrixCPU _lineLimits; // l
    MatrixCPU _lineLimitsReduce; // l
    MatrixCPU _lineLimitsChange; // l
   
    MatrixCPU fileCoresBus;
   
    /*MatrixCPU _CoresVoiLin;
    MatrixCPU _CoresBusLin;
    MatrixCPU _nLines;*/


    std::vector<int> _indiceLineConstraint;// numero de ligne contrainte
    std::vector<int> _indiceLineNonConstraint;// numero de ligne non contrainte
    
    MatrixCPU _zoneBus; // taille B*1 indique pour chaque bus la zone o� il est 
    std::vector<std::string> _nameZone;

    void setGridFromFile(const std::string& path, MatrixCPU* fileBusAgent);
    virtual void setBusFromFile(const std::string& path, MatrixCPU* fileCoresBus);
    
    void CalcGridSensi();
    void ReduceSensi();

    float rand1() const; 
    int randab(int a, int b) const; // return a random number between a and b
    int getNFileline(std::string nameFile);
    

public:
   
    int _invertMethod = 0;
    bool toReduce = false;
    float lineMin = 0;
    float lineoffset = 0;
    
    float _Zbase = 1;
    float _Sbase = 1;


    StudyCaseDCGrid();
    StudyCaseDCGrid(const StudyCaseDCGrid& s);
    StudyCaseDCGrid& operator= (const StudyCaseDCGrid& s);
    ~StudyCaseDCGrid();

    // generator
    
    void genGridFromFile(std::string path, bool alreadyDefine=true); 
    void genLineLimit(int nLine, float limit, float dlimit );

    //Set case
    void Set39Bus(std::string path = "data/", bool alreadyDefine=false);
    void Set3Bus(std::string path = "data/");
    void Set4nodeBis(std::string path);
    void Set2nodeConstraint(float lim = 0.8);
    void SetEuropeP0(const std::string& path, bool alreadyDefine=0);
    void SetStudyCaseDCGrid(std::string path, std::string name, int nBus, bool alreadyDefine);
    virtual void setFromInterface(StudyCaseInterface* interface);
    
    // set Parametters
    void setLineLimitMin(float lineMin);
    void setLineLimitRelaxation(float eps);
    void setLineLimit(int line, float limit);


    // getter 
 

    MatrixCPU getPowerSensiBus(bool force=false) const;
    MatrixCPU getPowerSensiBusReduce() const;
    MatrixCPU getLineLimit() const;
    MatrixCPU getCoresLineBus(bool force=false) const;
    MatrixCPU getfileCoresBus() const;
    MatrixCPU getZones() const;


    float getTimeInit() const;
    int getNLine(bool force=false) const;
    int getNLineConstraint() const;
    int getNBus() const;

    std::string getName() const;



    // modification
 
    void saveCSV(const std::string& fileName);
    void display(int type=0) const;
    void displayLineCores(MatrixCPU* g, bool all = true);
 

};
