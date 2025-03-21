# Endogenous Peer to Peer market 

## Dependencies
### Mandatory
- cl.exe
- Cuda
- SetupTool

### Optional 
- OSQP https://osqp.org/ 
- Eigen : https://eigen.tuxfamily.org/index.php?title=Main_Page

## How to install :
- download or clone the git repository
- download and install Visual Studio OR only "C++ build Tools" From "Build Tools for Visual Studio" (https://visualstudio.microsoft.com/fr/downloads/)
- download Cuda (https://developer.nvidia.com/cuda-downloads)
- (optional :) if cl.exe is not found it may be useful to add it on the PATH (year or version may be different):
    - ==> C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.33.31629\bin\Hostx64\x64
- install the python module setupTools (**pip install setuptools**)
- in the folder EndogenousP2PMarket, compile by using the commande **python setup.py build --force**
- the module **EndoCuda.cp311-win_amd64.pyd** located in *build/lib.win-amd64-cpython-311/* can be moved on the folder where you want to use the code. 

## How to use :
### Optionnal
- the file "test_interface.py" is composed by all unit testings and show examples
- it must be at the same place that *EndoCuda.cp311-win_amd64.pyd*
- it can be used by using **python3 -m unittest -v** (unittest may be install with **pip install unittest** if not already done)

### Mandatory
- import the module by using : import EndoCuda
- an interface (*here named resolution*) that store all data for N agents, B buses, L branches, and Lconst csontrainted branches (optional) : study case, computation's parameters and solver is created by **resolution = EndoCuda.interface(N, B, L, Lconst)**
- each data must be set by methods (ex : setPower, setCostFunction ....)
- solve by using :
    - **solveMarketFromInterface(interface resolution, string MethodName)** : for a DC market, DC Endogenous Market
    - **solveACFromInterface(interface resolution, string MethodName)** : for a AC market, AC-OPF, AC Endogenous Market
    - **solvePFFromInterface(interface resolution, string MethodName)** : for a Power Flow
- display informations by using **interface.display(Id)**, (Id = 0 shows all, Id > 12 shows effect for each Id)
- get result by using methods (ex : getPn, getTrade)

## Data 
Data come from serveral open sources data bases, scripts to convert these data into usable files are in the data folder :
- Matpower cases (https://matpower.org/) can be treated with the matlab files
- European transport grid (https://doi.org/10.1038/sdata.2017.175) can be traited with python files
- European testFeeder grid (http ://ewh.ieee.org/soc/pes/dsacom/testfeeders.html)



