#pragma once


#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include "MatrixCPU.h"
#include "MatrixGPU.cuh"
#include "MatrixCPUD.h"
#include <math.h>
#include <iostream>
#include <vector>


// ne gère plus les zones dans setACfile
// displayLineCores doit être modifier pour montrer toutes les contraintes

#define DELETEB(x) if (x!=nullptr) {delete x; x = nullptr;}
#define DELETEA(x) if (x!=nullptr) {delete[] x; x = nullptr;}


class StudyCaseACGrid
{
    static constexpr float LINELIMITMAX = 100000;
    std::string _name = "None";
    float _timeInit;
  

    // AC et DC
    MatrixCPU _zoneBus; // taille B*1 indique pour chaque agent la zone où il est 
    std::vector<std::string> _nameZone;
    int _nBus = 0;
    int _nLine = 0;
    int _nLineConstraint = 0;
    MatrixCPU _CoresLineBus; // Perso from, to
    MatrixCPU _lineLimits; // l
    //MatrixCPU _lineLimitsReduce; // l  //MatrixCPU _CoresLineBusReduce; // Perso

   

    // AC
    int _nConstraint = 0;

  
    MatrixCPU _lineSuceptance; // Bij
    MatrixCPU _lineReactance; // Gij
    MatrixCPUD _lineSuceptanceD; // Bij
    MatrixCPUD _lineReactanceD; // Gij

    MatrixCPU _lineSuceptanceLin;// Bij
    MatrixCPU _lineReactanceLin;// Gij
    MatrixCPUD _lineSuceptanceLinD;// Bij
    MatrixCPUD _lineReactanceLinD;// Gij

    MatrixCPU _lineSuceptanceLin2; // only line impedance without self impedance
    MatrixCPU _lineReactanceLin2;

    
    MatrixCPU _CoresVoiLin;
    MatrixCPU _CoresBusLin;
    MatrixCPU _nLines;
    MatrixCPU _nLinesBegin; 



    MatrixCPU _busSuceptance; // Yd for distribution network
    MatrixCPU _lineImpedanceReal; // real(zs) for distribution network
    MatrixCPU _lineImpedanceImag; // imag(zs) for distribution network

    MatrixCPU _upperBound; // overline(Y) : angle, voltage, powerFlow
    MatrixCPU _lowerBound; // underline(Y)  : angle, voltage, powerFlow

    MatrixCPUD _SolutionPF; // Matpower
    MatrixCPUD _VoltageInitD; // warmstart for PF
    MatrixCPU _VoltageInit; // warmstart for PF


   
    void setGridACFromFile(const std::string& path, MatrixCPU* fileBusAgent);
    void setBusFromFile(const std::string& path, MatrixCPU* fileCoresBus);
    
   
    float rand1() const; 
    int randab(int a, int b) const; // return a random number between a and b
    int getNFileline(std::string nameFile);
    
    void LinearizeImp();

public:
   
   
    float _Zbase = 1;
    float _Sbase = 1;
    float _Vbase = 1;
    double _V0 = 1; // tension au bus 0
    double _theta0 = 0; // angle au bus 0
    bool radial = false;

    StudyCaseACGrid();
    
    StudyCaseACGrid(const StudyCaseACGrid& s);
    StudyCaseACGrid& operator= (const StudyCaseACGrid& s);

    // generator
   
    
    void genGrid(int _nBus, int _nMajorLine, int _minorLine, float ReacMajor, float DeltaReacMajor, float ReacMinor, float DeltaReacMinor, float LlimitMajor, float dLlimitMajor, float LlimitMinor, float dLlimitMinor);
    void genGridBT(int _nBus, int Nbranch, int Ndeep, float length, float dlength);
    void genGridHTB(int nBus, int nLine, int dnLine, float length, float dlength);

       
    //Set case
   
    
    void SetAC39Bus(std::string path = "data/", bool alreadyDefine = false);
    void SetAC3Bus(std::string path = "data/");
    void SetACFromFile(std::string name, std::string path = "data/ACGrid/");
    void SetEuropeTestFeeder(std::string path = "data/ACGrid/");
    void setLineLimit(int line, float limit);
    
    // getter 

    MatrixCPU getLineLimit() const;
    MatrixCPU getCoresLineBus() const;
    MatrixCPU getLineSuceptance() const;
    MatrixCPU getLineReactance() const;
    MatrixCPUD getLineSuceptanceD() const;
    MatrixCPUD getLineReactanceD() const;
    MatrixCPU getUpperBound() const;
    MatrixCPU getLowerBound() const;
    
   
    MatrixCPUD getSolPF() const;
    double getV0() const;
    double gettheta0() const;
    MatrixCPU getZsRe() const;
    MatrixCPU getZsImag() const;
    MatrixCPU getYd() const;
    MatrixCPU getGlin() const;
    MatrixCPU getBlin() const;
    MatrixCPU getGlin2() const;
    MatrixCPU getBlin2() const;
    MatrixCPU getVoltageInit() const;
    MatrixCPUD getVoltageInitD() const;
    MatrixCPUD getGlinD() const;
    MatrixCPUD getBlinD() const;
    MatrixCPU getNLines() const;
    MatrixCPU getCoresBusLin() const;
    MatrixCPU getCoresVoiLin() const;
    MatrixCPU getZones() const;
    MatrixCPU getNLinesBegin() const;
    
    float getTimeInit() const;
    
    int getNLine() const;
    int getNBus() const;
   
    std::string getName() const;


    void saveCSV(const std::string& fileName, bool all = true);
    void display(int type=0);
    void displayLineCores(MatrixCPU* g, bool all = true);
    ~StudyCaseACGrid();

};
