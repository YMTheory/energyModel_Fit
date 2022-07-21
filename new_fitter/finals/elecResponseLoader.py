import numpy as np

class elecResponseLoader(object):

    def __init__(self, qmode, cmode, rmode, det):
        """configure the detector and energy response description"""
        self.quench_mode = qmode
        self.quenchNonl_E = []
        self.quenchNonl_val = []
        self.loadQuenchNonlFile_flag = False
        self.Ys = 1400
        self.kB = 6.5e-3
        self.kBflag = False

        self.cerenkov_mode = cmode
        self.loadCerenkovYieldFile_flag = False
        self.gNcer = None
        self.p0 = 0
        self.p1 = 0
        self.p2 = 0
        self.E0 = 0
        self.kC = 1

        self.resolution_mode = rmode
        self.a = 0
        self.b = 0
        self.n = 0
        self.a1 = 0
        self.n1 = 0

        self.Ncer_truth = []
        self.NPE_truth = []
        self.NPE_sigma_truth = []

        self.detconfig = det

        if self.detconfig == "Det1":
            # Det1
            self.Y = 3111.440 / 2.223
            self.m_npeGe68 = 1.30880e+03;
            self.m_sigmaGe68 = 3.80672e+01;

        elif self.detconfig == "Det5" :
            # Det5
            self.m_npeGe68 = 1.64037e+02;
            self.m_sigmaGe68 = 1.28369e+01;
            self.Y =  3.89414e+02/2.223;

        elif self.detconfig == "Det6" :
            # Det 6
            self.m_npeGe68 = 6.54907e+02;
            self.m_sigmaGe68 = 2.64157e+01;
            self.Y = 1.55597e+03/2.223;

        elif self.detconfig == "Det7" :
            # Det 7 
            self.m_npeGe68 = 3.27754e+02
            self.m_sigmaGe68 = 1.82715e+01
            self.Y = 7.78459e+02/2.223;

        elif self.detconfig == "Det8" :
            # Det 8
            self.Y = 6.22873e3 / 2.223
            self.m_npeGe68 = 2.62131e2
            self.m_sigmaGe68 = 1.64129e1

        elif self.detconfig == "DYBnonl" :
            self.Y = 2.94122e3 / 2.223
            self.m_npeGe68 = 1.19868e3
            self.m_sigmaGe68 = 4.02389e1



    def getY(self):
        """return energy scale of current detector configuration"""
        return self.Y
    
    def getNPE_Ge68(self):
        return self.m_npeGe68
    
    def getSigma_Ge68(self):
        return self.m_sigmaGe68

    def setp0(self, p0):
        self.p0 = p0

    def setp1(self, p1):
        self.p1 = p1

    def setp2(self, p2):
        self.p2 = p2

    def setE0(self, E0):
        self.E0 = E0

    def setkC(self, kC):
        self.kC = kC

    def setYs(self, Ys):
        self.Ys = Ys

    def setkB(self, kB):
        self.kB = kB

    def IskBHigh(self, flag):
        self.kBflag = flag

    def seta(self, a):
        self.a = a

    def setb(self, b):
        self.b = b

    def setn(self, n):
        self.n = n

    def seta1(self, a):
        self.a1 = a

    def setn1(self, n):
        self.n1 = n

    def getCerenkovMode(self):
        return self.cerenkov_mode


    def loadQuenchNonlFile(self):
        """load quenching nonlinearity files """
        import ROOT
        if self.quench_mode == "Sim":
            filename = "/junofs/users/miaoyu/energy_model/fitter/energyModel_Fit/new_fitter/data/electron/quench_local.root"
            f = ROOT.TFile(filename, "read")
            kBmin, kBmax = 51, 89
            for kBval in range(kBmin, kBmax+1, 1):
                histname = "kB%d"%kBval
                h = f.Get(histname)
                tmp_E, tmp_val = [], []
                for i in range(h.GetNbinsX()):
                    tmp_E.append(h.GetBinCenter(i+1))
                    tmp_val.append(h.GetBinContent(i+1))
                self.quenchNonl_E.append(tmp_E)
                self.quenchNonl_val.append(tmp_val)
                self.loadQuenchNonlFile_flag = True
            print("+++++++ Load quench file successfully!!! +++++++")
        elif self.quench_mode == "Int":
            if not self.kBflag:
                filename = "/junofs/users/miaoyu/energy_model/fitter/energyModel_Fit/new_fitter/data/electron/Quench_NumInt.root"
            else:
                filename = "/junofs/users/miaoyu/energy_model/fitter/energyModel_Fit/new_fitter/data/electron/Quench_NumInt_DYBkB.root"
            f = ROOT.TFile(filename, "read")
            if not self.kBflag:
                kBmin, kBmax = 45, 95
            else:
                kBmin, kBmax = 120, 179
            for kBval in range(kBmin, kBmax+1, 1):
                h = f.Get("kB%d"%kBval)
                tmp_E, tmp_val = [], []
                for i in range(h.GetNbinsX()):
                    tmp_E.append(h.GetBinCenter(i+1))
                    tmp_val.append(h.GetBinContent(i+1))
                self.quenchNonl_E.append(tmp_E)
                self.quenchNonl_val.append(tmp_val)
                self.loadQuenchNonlFile_flag = True
            print("+++++++ Load quench file successfully!!! +++++++")
        else:
            print("******* Error! No such quench mode!! *******")


    def getQuenchNonl(self, E):
        """ get quenching nonlinearity value for given E and kB """
        if E > 14:
            E = 14
        if not self.loadQuenchNonlFile_flag:
            self.loadQuenchNonlFile()
        kBidx = int(self.kB*1e4)
        if self.quench_mode == "Sim":
            kBmin, kBmax = 51, 89 
            if E < 0.1:
                Eidx = int(E/0.001)
            else:
                Eidx = int((E-0.1)/0.01) + 100
        elif self.quench_mode == "Int":
            if not self.kBflag:
                kBmin, kBmax = 45, 95
            else:
                kBmin, kBmax = 120, 179
            Eidx = int(E/0.001)
        else:
            print("******* Error! No such quench mode!! *******")

        if kBidx < kBmin:
            return self.quenchNonl_val[0][Eidx]
        elif kBidx > kBmax:
            return self.quenchNonl_val[-1][Eidx]
        else:
            kBResid = self.kB*1e4 - kBidx 
            EResid = (E - self.quenchNonl_E[0][Eidx]) / (self.quenchNonl_E[0][Eidx+1]-self.quenchNonl_E[0][Eidx])
            qnl_low = self.quenchNonl_val[kBidx-kBmin]
            qnl_hig = self.quenchNonl_val[kBidx+1-kBmin]
            
            quenchNonl_val_low = qnl_low[Eidx] * (1-EResid) + qnl_hig[Eidx] * EResid
            quenchNonl_val_hig = qnl_hig[Eidx] * (1-EResid) + qnl_hig[Eidx] * EResid
            quenchNonl_val_mid = quenchNonl_val_hig * kBResid + quenchNonl_val_low * (1-kBResid)

            return quenchNonl_val_mid


    def getScintillationNumber(self, E):
        """ return scintillation photno numbers"""
        return self.Ys * self.getQuenchNonl(E) * E


    def loadCerenkovYieldFile(self):
        """ load Ncer file from simulation """
        import ROOT
        filename = "/junofs/users/miaoyu/energy_model/fitter/energyModel_Fit/new_fitter/data/electron/Ncer_local.root"
        f = ROOT.TFile(filename, "read")
        self.gNcer = f.Get("Ncer")

        self.loadCerenkovYieldFile_flag = True
        print("+++++++ Load Cerenkov sim file successfully!!! +++++++")


    def getSimNcer(self, E):
        """ absolute Ncer value from simulation """
        if E > 14:
            E = 14
        if not self.loadCerenkovYieldFile_flag:
            self.loadCerenkovYieldFile()
        return self.gNcer.Eval(E)


    def func1_Ncer(self, x):
        E = x - self.E0
        if E < 0:
            return 0
        else:
            return self.p0 * E * E / (E + self.p1 * np.exp(-self.p2 * E))


    def getCerenkovNumber(self, E):
        """ return Cerenkov photon numbers """
        if self.cerenkov_mode == "Sim":
            return self.getSimNcer(E) * self.kC
        elif self.cerenkov_mode == "Ana1":
            return self.func1_Ncer(E)

        

    def func1_sigma(self, x):
        N = x * self.Y
        return np.sqrt(self.a**2*N + self.b**2*np.power(N, self.n)) / self.Y

    def func2_sigma(self, x):
        N = x * self.Y
        return np.sqrt(self.a1**2*np.power(N, self.n1)) / self.Y


    def getEvisSigma(self, E):
        """ return Evis sigma """
        if self.resolution_mode == "New":
            return self.func1_sigma(E)
        elif self.resolution_mode == "New1":
            return self.func2_sigma(E)






                
