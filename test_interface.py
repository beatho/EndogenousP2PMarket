import unittest
import EndoCuda

#python ou python3 !
#Pour build
# python setup.py build --force
# python setup.py build_ext --force       pour que cpp et linkage

# Pour lancer le test
# python3 -m unittest -v

''' Param et Res
enum indParam {iterG_ind, iterL_ind, iterItern_ind, stepG_ind, stepL_ind, stepIntern_ind, epsG_ind,
     epsL_ind, epsX_ind, epsIntern_ind, rho_ind, rhoL_ind, rho1_ind, SbaseP_ind, VbaseP_ind, finParam_ind};
enum indSizes {nAgentP_ind, nBusP_ind, nLineP_ind, nLineCons_ind, finSizes_ind};
enum indRes   {iterF_ind, temps_ind, fc_ind, iterMax_ind, stepG2_ind, resR_ind, resS_ind, resX_ind, finRes_ind};

StudyCase
enum indInfo   { Sbase_ind, Vbase_ind, nAgent_ind, nCons_ind, nGen_ind, nBus_ind, nLine_ind, V0_ind, theta0_ind, finInfo_ind };
enum indAgent  { PosBus_ind, a_ind, b_ind, aq_ind, bq_ind, Pobj_ind, Pmin_ind, Pmax_ind, Qobj_ind, Qmin_ind, Qmax_ind, zone_ind, finAgent_ind };
enum indBuses  { Gs_ind, Bs_ind, Vmin_ind, Vmax_ind, thetamin_ind, thetamax_ind, Vinit_ind, thetainit_ind, finBuses_ind};
enum indBranch { From_ind, To_ind, YsRe_ind, YsIm_ind, Yp_ind, tau_ind, theta_ind, lim_ind, ZsRe_ind, ZsIm_ind, finBranch_ind};

'''
class TestInterfaceMethod(unittest.TestCase):

    def test_Creation(self):
        test = EndoCuda.interface(1,1,1)
        R = test.getInfo()
        self.assertEqual(R[ind["nAgent_ind"]], 1)
        self.assertEqual(R[ind["nBus_ind"]], 1)
        self.assertEqual(R[ind["nLine_ind"]], 1)
        test = EndoCuda.interface(4,3, 2, 1)
        R = test.getInfo()
        self.assertEqual(R[ind["nAgent_ind"]], 4)
        self.assertEqual(R[ind["nBus_ind"]], 3)
        self.assertEqual(R[ind["nLine_ind"]], 2)

        with self.assertRaises(TypeError):
            test = EndoCuda.interface(3,2)

        with self.assertRaises(ValueError):
            test = EndoCuda.interface(-3, 2, 1)

        with self.assertRaises(TypeError):
            test = EndoCuda.interface(3.3, 2, 1)
    def test_SetSbase(self):
        test = EndoCuda.interface(4,4,3)
        Sbase = 10
        R1 = test.getInfo()
        S1 = R1[ind["Sbase_ind"]]
        self.assertEqual(S1, 1)
        test.setSbase(Sbase)
        R2 = test.getInfo()
        S2 = R2[ind["Sbase_ind"]]
        self.assertEqual(S2, Sbase)

        with self.assertRaises(TypeError):
            test.setSbase("10")
        with self.assertRaises(ValueError):
            test.setSbase(0)
        with self.assertRaises(ValueError):
            test.setSbase(-1)
    def test_SetVbase(self):
        test = EndoCuda.interface(4,4,3)
        Vbase = 10
        R1 = test.getInfo()
        V1 = R1[ind["Vbase_ind"]]
        self.assertEqual(V1, 1)
        test.setVbase(Vbase)
        R2 = test.getInfo()
        V2 = R2[ind["Vbase_ind"]]
        self.assertEqual(V2, Vbase)

        with self.assertRaises(TypeError):
            test.setVbase("10")
        with self.assertRaises(ValueError):
            test.setVbase(0)
        with self.assertRaises(ValueError):
            test.setVbase(-1)
    
    def test_SetV0(self):
        test = EndoCuda.interface(4,4,3)
        V0 = 10
        R1 = test.getInfo()
        V1 = R1[ind["V0_ind"]]
        self.assertEqual(V1, 1)
        test.setV0(V0)
        R2 = test.getInfo()
        V2 = R2[ind["V0_ind"]]
        self.assertEqual(V2, V0)

        with self.assertRaises(TypeError):
            test.setV0("10")
        with self.assertRaises(ValueError):
            test.setV0(0)
        with self.assertRaises(ValueError):
            test.setV0(-1)
    def test_SetTheta(self):
        test = EndoCuda.interface(4,4,3)
        theta = 0.2
        R1 = test.getInfo()
        T1 = R1[ind["theta0_ind"]]
        self.assertEqual(T1, 0)
        test.setTheta(theta)
        R2 = test.getInfo()
        V2 = R2[ind["theta0_ind"]]
        self.assertAlmostEqual(V2, theta)

        with self.assertRaises(TypeError):
            test.setTheta("10")
        with self.assertRaises(ValueError):
            test.setTheta(10)
        with self.assertRaises(ValueError):
            test.setTheta(-10)

    def test_PosBus(self):
        test = EndoCuda.interface(3,2,1)
        L  = [0, 1, -1]
        L2 = [3, 1, 3]
        L3 = [1, 0]
        L4 = [-4, 0, 0]
        L5 = [3.3, 0, 0]
        test.setPosBus(L)
        R = getColumn(test.getAgent(),ind["PosBus_ind"], 3, ind["finAgent_ind"])
        self.assertEqual(R, L)
        
        test.setPosBus(L,L2)
        R = test.getAgent()
        R1 = getColumn(test.getAgent(),ind["PosBus_ind"], 3, ind["finAgent_ind"])
        R2 = getColumn(test.getAgent(),ind["zone_ind"], 3, ind["finAgent_ind"])
        self.assertEqual(R1, L)
        self.assertEqual(R2, L2)

        with self.assertRaises(ValueError):
            test.setPosBus(L2, L) #valeur>B
        with self.assertRaises(TypeError):
            test.setPosBus(L, L3, L2) #nb de listes
        with self.assertRaises(ValueError):
            test.setPosBus(L3) # taille de liste
        with self.assertRaises(ValueError):
            test.setPosBus(L4) # valeur <-1
        with self.assertRaises(ValueError):
            test.setPosBus(L5) # valeur float
    def test_CostFunction(self):
        test = EndoCuda.interface(3,2,1)
        L = [1, 1, 1]
        L2 = [2, 4, 3]
        L3 = [4, 5]
        L4 = [-1, 1, 1]
        test.setCostFunction(L2, L)
        test.setCostFunction(L, L2)
        R = test.getAgent()
        R1 = getColumn(R,ind["a_ind"], 3, ind["finAgent_ind"])
        R2 = getColumn(R,ind["b_ind"], 3, ind["finAgent_ind"])
        self.assertEqual(R1, L)
        self.assertEqual(R2, L2)
        with self.assertRaises(TypeError):
            test.setCostFunction(L)
        with self.assertRaises(ValueError):
            test.setCostFunction(L, L3)
        with self.assertRaises(ValueError):
            test.setCostFunction(L4, L3)
    def test_SetPower(self):
        test = EndoCuda.interface(3,2,1)
        L = [1, 1, 1]
        L2 = [2, 4, 3]
        L3 = [-1, 1, 1]
        L4 = [4, 5]
        L0 = [0, 0, 0]
        LL3 = [0, 1, 1]
        L11 = [-1, 3, -4, 2, 5, 6]
        L22 = [2, 5, 2, 4, 7, 8]
        L33 = [1, 4, 1, 4, 6, 5]
        L1122 = [0.5, 4, -1, 3, 6, 7]

        test.setPower(L3, L)
        R = test.getAgent()
        self.CheckAllPower(R, L3, L, LL3, L0, L0, L0, 3)
       
        test.setPower(L, L2, L3)
        R = test.getAgent()
        self.CheckAllPower(R, L, L2, L3, L0, L0, L0, 3)

        test.setPower(L11, L22)
        R = test.getAgent()
        self.CheckAllPower(R, L11, L22, L1122, L0, L0, L0, 3)

        test.setPower(L11, L22, L33)
        R = test.getAgent()
        self.CheckAllPower(R, L11, L22, L33, L0, L0, L0, 3)

        with self.assertRaises(TypeError):
            test.setPower(L) # pas assez
        with self.assertRaises(TypeError):
            test.setPower(L, L2, L3, L0) # trop
        with self.assertRaises(ValueError):
            test.setPower(L, L4) # taille 2
        with self.assertRaises(ValueError):
            test.setPower(L4, L) # taille 1
        with self.assertRaises(ValueError):
            test.setPower(L3, L, L11) # taille 3
        with self.assertRaises(ValueError):
            test.setPower(L2, L) # max < min
    def CheckAllPower(self, Agents, Pmin, Pmax, Pobj, Qmin, Qmax, Qobj, N):
        R1 = getColumn(Agents, ind["Pmin_ind"], 3, ind["finAgent_ind"])
        R2 = getColumn(Agents, ind["Pmax_ind"], 3, ind["finAgent_ind"])
        R3 = getColumn(Agents, ind["Pobj_ind"], 3, ind["finAgent_ind"])
        R4 = getColumn(Agents, ind["Qmin_ind"], 3, ind["finAgent_ind"])
        R5 = getColumn(Agents, ind["Qmax_ind"], 3, ind["finAgent_ind"])
        R6 = getColumn(Agents, ind["Qobj_ind"], 3, ind["finAgent_ind"])

        if(len(Pmin) == 2*N):
            Qmin = []
            Qmax = []
            Qobj = []
            for i in range(N):
                Qmin.append(Pmin[N + i])
                Qmax.append(Pmax[N + i])
                Qobj.append(Pobj[N + i])
        elif(len(Pmin) != N):
            print("[ERROR] : wrong size of Vector used ")
            self.assertEqual(0, 1)

        for i in range(N):
            self.assertAlmostEqual(R1[i], Pmin[i])
            self.assertAlmostEqual(R2[i], Pmax[i])
            self.assertAlmostEqual(R3[i], Pobj[i])
            self.assertAlmostEqual(R4[i], Qmin[i])
            self.assertAlmostEqual(R5[i], Qmax[i])
            self.assertAlmostEqual(R6[i], Qobj[i])
    def test_setImpedanceBus(self):
        test = EndoCuda.interface(3,2,1)
        L = [1, 1]
        L2 = [2, 4]
        L3 = [4, 5, 4]
        test.setImpedanceBus(L2, L)
        test.setImpedanceBus(L, L2)
        R = test.getBus()
        R1 = getColumn(R,ind["Gs_ind"], 2, ind["finBuses_ind"])
        R2 = getColumn(R,ind["Bs_ind"], 2, ind["finBuses_ind"])
        self.assertEqual(R1, L)
        self.assertEqual(R2, L2)
        with self.assertRaises(TypeError):
            test.setImpedanceBus(L)
        with self.assertRaises(ValueError):
            test.setImpedanceBus(L, L3)
        with self.assertRaises(ValueError):
            test.setImpedanceBus(L3, L)
    def test_setVoltageBound(self):
        test = EndoCuda.interface(3,2,1)
        L = [1, 1]
        L2 = [2, 4]
        L3 = [4, 5, 4]
        L4 = [0.5, -3]
        test.setVoltageBound(L, L2)
        test.setVoltageBound(L, L2)
        R = test.getBus()
        R1 = getColumn(R,ind["Vmin_ind"], 2, ind["finBuses_ind"])
        R2 = getColumn(R,ind["Vmax_ind"], 2, ind["finBuses_ind"])
        self.assertEqual(R1, L)
        self.assertEqual(R2, L2)
        with self.assertRaises(TypeError):
            test.setVoltageBound(L)
        with self.assertRaises(ValueError):
            test.setVoltageBound(L, L3)
        with self.assertRaises(ValueError):
            test.setVoltageBound(L3, L)
        with self.assertRaises(ValueError):
            test.setVoltageBound(L4, L)
        with self.assertRaises(ValueError):
            test.setVoltageBound(L2, L)

    def test_setProblem(self): # et conversion
        test = EndoCuda.interface(2,1,0)
        L1 = [0, -1, 1, 0]
        M1 = [-1, 1]
        test.initProblem(L1, M1)
        R1  = test.getTrade()
        RM1 = test.getPn()
        self.assertEqual(convertNTo2N1Matrix(L1, 2), R1)
        self.assertEqual(convertNTo2N1Vector(M1), RM1)
        self.assertEqual(L1, convert2N1toNMatrix(R1, 2))
        self.assertEqual(M1, convert2N1toNVector(RM1))

        L2 = [0, -1, 1, 0, 0, 2, -2, 0]
        M2 = [-1, 1, 2, -2]
        test.initProblem(L2, M2)
        R2 = test.getTrade()
        RM2 = test.getPn()
        self.assertEqual(addLossAgentMatrix(L2, 2), R2)
        self.assertEqual(addLossAgentVector(M2), RM2)
        
        L3 = [ 1, 1, 1, 1, 1]
        


        with self.assertRaises(TypeError):
            test.initProblem(L1)
        with self.assertRaises(ValueError):
            test.initProblem(L3, M1)
        with self.assertRaises(ValueError):
            test.initProblem(L1, M2)
    def test_setDual(self): # et conversion
        test = EndoCuda.interface(2,1,0)
        L1 = [0, -6, 6, 0]
        M1 = [3, 4]
        test.initDual(L1, M1)
        R1  = test.getLambda()
        RM1 = test.getMu()
        self.assertEqual(convertNTo2N1Matrix(L1, 2), R1)
        self.assertEqual(convertNTo2N1Vector(M1), RM1)
        self.assertEqual(L1, convert2N1toNMatrix(R1, 2))
        self.assertEqual(M1, convert2N1toNVector(RM1))

        L2 = [0, -6, 6, 0, 0, 2, 2, 0]
        M2 = [1, 2, 3, 4]
        test.initDual(L2, M2)
        R2 = test.getLambda()
        RM2 = test.getMu()
        self.assertEqual(addLossAgentMatrix(L2, 2), R2)
        self.assertEqual(addLossAgentVector(M2), RM2)
        
        L3 = [ 1, 1, 1, 1, 1]
        


        with self.assertRaises(TypeError):
            test.initDual(L1)
        with self.assertRaises(ValueError):
            test.initDual(L3, M1)
        with self.assertRaises(ValueError):
            test.initDual(L1, M2)
    def test_setDelta(self):
        test = EndoCuda.interface(4, 3, 2)
        d1 = [2, 0]
        d2 = [0, 3]
        test.initDelta(d1, d2)
        R = test.getDelta()
        R1 = getColumn(R, 0, 2, 2)
        R2 = getColumn(R, 1, 2, 2)
        self.assertEqual(d1, R1)
        self.assertEqual(d2, R2)
        d3 = [2, 1, 3]

        with self.assertRaises(TypeError):
            test.initDelta(d1)
        with self.assertRaises(ValueError):
            test.initDelta(d3, d1)
        with self.assertRaises(ValueError):
            test.initDelta(d1, d3)


    def test_display(self):
        test = EndoCuda.interface(2,2,1)
        #test.display(20)


    
    '''float Plim1[2] = { -30, 0 };
	float Plim2[2] = { 0, 60 };
	float Cost1[2] = { 1, 1 };
	float Cost2[2] = { 8, 4 };
	float Qobj[2] = { -1, 0 };'''
    def test_solveMarket(self):
        test = EndoCuda.interface(2,1,0)
        Pmin = [-30, 0]
        Pmax = [0, 60]
        Cost1 = [1, 1]
        Cost2 = [8, 4]
        Beta = [0, -1, 1, 0]
        methodes = ["ADMM", "ADMMMP", "ADMMGPU"]
        test.setCostFunction(Cost1, Cost2)
        test.setPower(Pmin, Pmax)
        test.setBeta(Beta)
        PnExpected = [-1, 1]
        TradeExpected = [0, -1, 1, 0]
        PnExpected = convertNTo2N1Vector(PnExpected)
        TradeExpected = convertNTo2N1Matrix(TradeExpected, 2)
        
        #EndoCuda.testAffichage(test)
        for methode in methodes :
            EndoCuda.solveMarketFromInterface(test, methode)
            #test.display(10)
            
            Pn = test.getPn()
            Trade = test.getTrade()
          
            for i in range(2*(2+1)):
                self.assertAlmostEqual(Pn[i], PnExpected[i], 4)
            for i in range(2*(2+1)*(2+1)):
                self.assertAlmostEqual(Trade[i], TradeExpected[i], 4)

    def test_solveMarketinit(self): #initialisation au bon resultat, iter == 1
        test = EndoCuda.interface(2,1,0)
        Pmin = [-30, 0]
        Pmax = [0, 60]
        Cost1 = [1, 1]
        Cost2 = [8, 4]
        Beta = [0, -1, 1, 0]
        PnExpected = [-1, 1]
        TradeExpected = [0, -1, 1, 0]
        LambdaExpected = [0, -6, -6, 0]
        Mu = [0, 0]
        methodes = ["ADMM", "ADMMMP", "ADMMGPU"]
        test.setCostFunction(Cost1, Cost2)
        test.setPower(Pmin, Pmax)
        test.setBeta(Beta)
        test.initProblem(TradeExpected, PnExpected)
        test.initDual(LambdaExpected, Mu)

        
        PnExpected = convertNTo2N1Vector(PnExpected)
        TradeExpected = convertNTo2N1Matrix(TradeExpected, 2)
        
        #EndoCuda.testAffichage(test)
        for methode in methodes :
            EndoCuda.solveMarketFromInterface(test, methode)
            #test.display(10)
            res = test.getResults()
            Pn = test.getPn()
            Trade = test.getTrade()

            self.assertEqual(res[ind["iterF_ind"]], 1)
            for i in range(2*(2+1)):
                self.assertAlmostEqual(Pn[i], PnExpected[i], 4)
            for i in range(2*(2+1)*(2+1)):
                self.assertAlmostEqual(Trade[i], TradeExpected[i], 4)

    def test_solveDCEndoMarket(self):
        test = EndoCuda.interface(2,2,1)
        Pmin = [-30, 0]
        Pmax = [0, 60]
        Cost1 = [1, 1]
        Cost2 = [8, 4]
        Beta = [0, -1, 1, 0]
        PosBus = [0, 1]
        reactance = [-0.001]
        lim = [0.8]
        resistance = [ 0 ] # la valeur est inutile
        methodes = ["DCEndoMarket", "DCEndoMarketGPU"]
        test.setCostFunction(Cost1, Cost2)
        test.setPower(Pmin, Pmax)
        test.setBeta(Beta)
        test.setPosBus(PosBus)
        test.setImpedance(resistance, reactance)
        test.setLink([0], [1])
        test.setLineLimit(lim)
        test.setEps(0,0,0.00001)

        PnExpected = [-lim[0], lim[0]]
        TradeExpected = [0, -lim[0], lim[0], 0]
        PnExpected = convertNTo2N1Vector(PnExpected)
        TradeExpected = convertNTo2N1Matrix(TradeExpected, 2)
        
        #EndoCuda.testAffichage(test)
        for methode in methodes :
            EndoCuda.solveMarketFromInterface(test, methode)
            #test.display(10)
            
            Pn = test.getPn()
            Trade = test.getTrade()
          
            for i in range(2*(2+1)):
                self.assertAlmostEqual(Pn[i], PnExpected[i], 4)
            for i in range(2*(2+1)*(2+1)):
                self.assertAlmostEqual(Trade[i], TradeExpected[i], 4)

    def test_solveOPF(self):
        N = 11
        test = EndoCuda.interface(N, 10, 9)
        self.setAgentCase10ba(test)
        self.setGridCase10ba(test)
        test.setSbase(10)
        test.setVbase(23)
        #test.display(1)
        test.setIter(1000,1000)
        test.setStep(1, 5)
        test.setEps(0.001, 0.0001)
        test.setRho(1)
        PnExpected = [0, -0.177087, -0.0911295, -0.172175, -0.153039, -0.154113, -0.0711774, -0.108034, -0.0906724, -0.156386, 1.24373, 0, 0, -0.0470297, -0.0350324, -0.0456224, -0.185028, 
-0.0609725, -0.01155, -0.0063, -0.01365, -0.0208353, 0.417576, 0.102623]
        #test.display(7)
        methodes = ["OPFADMM", "OPFADMM2", "OPFADMMGPU", "OPFADMMGPU2"]
        for methode in methodes :
            EndoCuda.solveACFromInterface(test, methode)
            #test.display(13)
            res = test.getResults()
            Pn = test.getPn()
            for i in range(len(Pn)):
                self.assertAlmostEqual(Pn[i], PnExpected[i], 2)
                #print("{0:g} ".format(Pn[i]), end="")
            #print("\n perte = ", sum(Pn))
            #print("\n", res)
            #print("***************************************")
       
        #EndoCuda.solveACFromFile("case10ba","OPFADMM")

    def test_solveACEndoMarket(self):
        N = 11
        Sbase = 10
        test = EndoCuda.interface(N, 10, 9)
        self.setAgentCase10ba(test)
        self.setGridCase10ba(test)
        test.setSbase(Sbase)
        test.setVbase(23)
        #test.display(1)
        test.setIter(10000,5000, 5000)
        test.setStep(1, 1, 1)
        test.setEps(0.001, 0.0001, 0.001, 0.0001)
        #test.setRho(10)
        PnExpected = [-0.0682349, -0.178785, -0.0926228, -0.172995, -0.153318, -0.153767, -0.0704824, -0.107069, -0.089416, -0.15496,
             1.24178, 0, -0.0913371, -0.0479741, -0.0357, -0.0460836, -0.185217, -0.0609629, -0.01155, -0.0063, -0.0136173, -0.0205507, 0.4165, 0.101547]
        #test.display(7)
        methodes = ["EndoMarketCons", "EndoMarketDirect", "EndoMarketConsGPU", "EndoMarketDirectGPU"]
        for methode in methodes :
            EndoCuda.solveACFromInterface(test, methode)
            #test.display(13)
            res = test.getResults()
            Pn = test.getPn()
            for i in range(len(Pn)):
                self.assertAlmostEqual(Pn[i], PnExpected[i], 2)
                #print("{0:g} ".format(Pn[i]), end="")
            #print("\n déséquilibre = ", sum(Pn))
            #print("\n", res)
            #print("***************************************")
       
        #EndoCuda.solveACFromFile("case10ba","OPFADMM")
    
    def test_solvePF(self):
        N = 11
        Sbase = 10
        test = EndoCuda.interface(N, 10, 9)
        self.setAgentCase10ba(test)
        self.setGridCase10ba(test)
        test.setSbase(Sbase)
        test.setVbase(23)
        #test.display(1)
        test.setIter(1000,1000)
        test.setStep(1, 5)
        test.setEps(0.001, 0.0001)
        test.setRho(1)
        	
        PbExpected  = [13.1517784516434, -1.84, -0.98, -1.79, -1.598, -1.61, -0.78, -1.15, -0.98, -1.64,
                       5.22247441139584, -0.46, -0.34, -0.446, -1.84, -0.6, -0.11, -0.06, -0.13, -0.2] 
        #test.display(7)
        methodes = ["NR", "GS", "DistPQ", "NRGPU", "GSGPU", "DistPQGPU"]
        for methode in methodes :
            EndoCuda.solvePFFromInterface(test, methode)
            #test.display(13)
            res = test.getResults()
            Pb = test.getPb()
            for i in range(len(Pb)):
                self.assertAlmostEqual(Pb[i], PbExpected[i]/ Sbase, 2)
                print("{0:g} ".format(Pb[i]), end="")
            #print("\n perte = ", sum(Pn))
            #print("\n", res)
            #print("***************************************")
       
        #EndoCuda.solveACFromFile("case10ba","OPFADMM")
    
    def setAgentCase10ba(self, interface):
        N = 11
        PosBus = [1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 0]
        Cost_a = [0.1 for i in range(2*N)]
        Cost_b = [0.184, 0.098, 0.179, 0.1598, 0.161, 0.078, 0.115, 0.098, 0.164, -1.2368, 20]
        Pobj   = [-1.84, -0.98, -1.79, -1.598, -1.61, -0.78, -1.15, -0.98, -1.64, 12.368, 0]
        Pmin   = [-2.208, -1.176, -2.148, -1.9176, -1.932, -0.936, -1.38, -1.176, -1.968, 9.8944, 0]
        Pmax   = [-1.472, -0.784, -1.432, -1.2784, -1.288, -0.624, -0.92, -0.784, -1.312, 14.84, 10]
        Qobj   = [-0.46, -0.34, -0.446, -1.84, -0.6, -0.11, -0.06, -0.13, -0.2, 4.186, 1.03647441139584]
        Qmin   = [-0.483, -0.357, -0.4683, -1.932, -0.63, -0.1155, -0.063, -0.1365, -0.21, 3.9767, -10]
        Qmax   = [-0.437, -0.323, -0.4237, -1.748, -0.57, -0.1045, -0.057, -0.1235, -0.19, 4.3953, 10]
        for i in range(N):
            Cost_b.append(-0.1*Qobj[i])
            Pobj.append(Qobj[i])
            Pmin.append(Qmin[i])
            Pmax.append(Qmax[i])
        interface.setPosBus(PosBus)
        interface.setCostFunction(Cost_a, Cost_b)
        interface.setPower(Pmin, Pmax, Pobj)
    def setGridCase10ba(self, interface):
        B = 10
        Vmin = [1, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.8, 0.8, 0.8]
        Vmax = [1.1 for i in range(B)]
        Vmax[0] = 1

        From = [ i for i in range(B-1)]
        To   = [ (i+1) for i in range(B-1)]
        YRe = [35.1575196289778, 2.02160816882563, 19.651315209743, 43.0642207686485, 15.1657701044738, 33.2235613793964, 19.488275265677, 8.35226725110039, 7.49557705154355 ]
        Yim = [-117.6764667549, -87.3767930683134, -31.7296460240401, -37.5147077829979, -13.2118321983202, -28.940793663749, -11.0375400979214, -4.73062328821735,  -4.24535209581753]
        Yp  = [0 for i in range(B-1)]
        ZRe = [ 0.00233081285444234, 0.000264650283553875,  0.0141077504725898, 0.0132022684310019, 0.0374877126654064, 0.0171134215500945, 0.0388506616257089,  0.090648393194707, 0.101009451795841 ]
        ZsIm = [0.00780151228733459, 0.0114385633270321, 0.0227788279773157, 0.0115009451795841, 0.0326578449905482, 0.0149073724007561, 0.0220037807183365, 0.0513421550094518, 0.0572098298676749]
        interface.setLink(From, To)
        interface.setVoltageBound(Vmin, Vmax)
        interface.setAdmitance(YRe, Yim, Yp)
        interface.setImpedance(ZRe, ZsIm)







if __name__ == '__main__':
    print("[WARNING] : index are defined statically, may be cause of error on the test")
    unittest.main()


indi = []
indi.append("iterG_ind, iterL_ind, iterItern_ind, stepG_ind, stepL_ind, stepIntern_ind, epsG_ind, epsL_ind, epsX_ind, epsIntern_ind, rho_ind, rhoL_ind, rho1_ind, SbaseP_ind, VbaseP_ind, finParam_ind")
indi.append("nAgentP_ind, nBusP_ind, nLineP_ind, nLineCons_ind, finSizes_ind")
indi.append("iterF_ind, temps_ind, fc_ind, iterMax_ind, stepG2_ind, resR_ind, resS_ind, resX_ind, finRes_ind")

indi.append("Sbase_ind, Vbase_ind, nAgent_ind, nCons_ind, nGen_ind, nBus_ind, nLine_ind, V0_ind, theta0_ind, finInfo_ind")
indi.append("PosBus_ind, a_ind, b_ind, aq_ind, bq_ind, Pobj_ind, Pmin_ind, Pmax_ind, Qobj_ind, Qmin_ind, Qmax_ind, zone_ind, finAgent_ind")
indi.append("Gs_ind, Bs_ind, Vmin_ind, Vmax_ind, thetamin_ind, thetamax_ind, Vinit_ind, thetainit_ind, finBuses_ind")
indi.append("From_ind, To_ind, YsRe_ind, YsIm_ind, Yp_ind, tau_ind, theta_ind, lim_ind, ZsRe_ind, ZsIm_ind, finBranch_ind")


ind = {}
for names in indi:
    i=0
    for name in names.split(', '):
        ind[name] = i
        i = i + 1

def getColumn(List, collumn, nbRow, nbColumn):
    L =[]
    for i in range(nbRow):
        L.append(List[i*nbColumn + collumn])
    return L

def printMatrix(List, nbRow, nbColumn):
    k = 0
    for i in range(nbRow):
        for j in range(nbColumn):
            print(List[k], " ")
            k = k + 1
        print("\n")

def convertNTo2N1Vector(List):
    n = len(List)
    L = [0]
    for i in range(n):
        L.append(List[i])
    for i in range(n+1):
        L.append(0)
    return L
def convertNTo2N1Matrix(List, N):
    M = len(List)
    L = []
    if(M != N*N):
        print("[ERROR ] wrong size for the matrix (convertNTo2N1Matrix)")
        return L
    
    for i in range(N + 1):
        L.append(0) #ligne agent des pertes
    k = 0
    for i in range(N):
        L.append(0) # colonne de l'agent des pertes
        for j in range(N):
            L.append(List[k])
            k = k + 1
    for i in range(N + 1):
        for j in range(N + 1):
           L.append(0) # partie des Q
    return L

def addLossAgentVector(List):
    L = [0]
    M = int(len(List)/2)
    if(2*M != len(List)):
        print("[ERROR] : problem on the size (addLossAgentVector) ")
    for i in range(M):
        L.append(List[i])
    L.append(0)
    for i in range(M):
        L.append(List[i + M])
    return L
def addLossAgentMatrix(List, N):
    M = len(List)
    L = []
    if(M != 2*N*N):
        print("[ERROR ] wrong size for the matrix (addLossAgentMatrix)")
    
    for i in range(N + 1):
        L.append(0) #ligne agent des pertes
    k = 0
    for i in range(N):
        L.append(0) # colonne de l'agent des pertes
        for j in range(N):
            L.append(List[k])
            k = k + 1
   
    for i in range(N + 1):
        L.append(0) #ligne agent des pertes
   
    for i in range(N):
        L.append(0) # colonne de l'agent des pertes
        for j in range(N):
            L.append(List[k])
            k = k + 1
    return L

def convert2N1toNVector(List):
    N = int(len(List)/2 - 1)
    if(2*(N + 1) != len(List)):
        print("[ERROR] : problem on the size (convert2N1toNVector) ")
    L = []
    for i in range(N):
        L.append(List[i + 1])
    return L
def convert2N1toNMatrix(List, N):
    if(2*(N+1)* (N+1) != len(List)):
        print("[ERROR] : problem on the size (convert2N1toNMatrix) ")
    L = []
    k = N + 1
    for i in range(N):
        for j in range(N + 1):
            if(j>0):
                L.append(List[k])
            k = k+1

    return L







