import numpy as np
import elecLoader as eloader
import prmBetaLoader as bloader
import uproot as up
from scipy.stats import norm
import matplotlib.pyplot as plt
import ROOT
from ROOT import TH1D, TF1
import time

ROOT.gROOT.SetBatch(ROOT.kFALSE)

class singleGamma(object):
   
    def __init__(self, name):
        self.name = name
        self.nSamples = 5000
        self.npe = 0
        self.spe = 0
        self.snpe = 0
        self.cnpe = 0
        self.sspe = 0
        self.cspe = 0
        self.npe_sim = 0
        self.npeerr_sim = 0
        self.spe_sim = 0
        self.speerr_sim = 0
        self.prmBetaArr = []
        self.prmAntiBetaArr = []
        self.loadPrmBetaFlag = False
        self.chi2 = 0
        self.sigma_part1 = 0
        self.sigma_part2 = 0
        self.phys = "Livermore"   # Penelope

        self.mode = "Fit"  # "Fit" / "Truth"

        ## Fitting config
        self.es   = 3134.078/2.223

        self.Asct = 1408
        self.kB   = 6.21e-3
        self.kC   = 0.996
        self.m_p0 = 4.20
        self.m_p1 = 3.64
        self.m_p2 = 0.11
        self.m_p3 = 405.02
        self.m_p4 = -2.09
        self.a    = 0.988
        self.b    = 7.89e-3
        self.c    = 0

    def getName(self):
        #print("Current gamma is " + self.name)
        return self.name

    def getNPE(self):
        return self.npe

    def getSPE(self):
        return self.spe

    def getSNPE(self):
        return self.snpe

    def getSSPE(self):
        return self.sspe

    def getCNPE(self):
        return self.cnpe

    def getCSPE(self):
        return self.cspe

    def getNPESim(self):
        return self.npe_sim

    def getNPEErrSim(self):
        return self.npeerr_sim

    def getSPESim(self):
        return self.spe_sim

    def getSPEErrSim(self):
        return self.speerr_sim

    def getSigmaPart1(self):
        return self.sigma_part1

    def getSigmaPart2(self):
        return self.sigma_part2

    def setPhys(self, m_phys):
        self.phys = m_phys


    def setAsct(self, m_Asct):
        self.Asct = m_Asct

    def setkB(self, m_kB):
        self.kB = m_kB

    def setkC(self, m_kC):
        self.kC = m_kC


    def seta(self, m_a):
        self.a = m_a

    def setb(self, m_b):
        self.b = m_b

    def setc(self, m_c):
        self.c = m_c

    
    def getFitNPE(self, Etrue):
        NPE = eloader.getQPE(Etrue, self.kB, self.Asct) + self.kC * eloader.getCerNPE(Etrue)
        return NPE

    def getFitSPE(self, E):
        sigma = np.sqrt(self.c**2 + self.a**2*E/self.es + self.b**2*E**2)*self.es
        return sigma

    def getFitCNPE(self, Etrue):
        #NPE = self.kC * eloader.getCerNPE(Etrue)
        NPE = (self.m_p3*Etrue*Etrue) / (1+self.m_p0*Etrue + self.m_p1*np.exp(-self.m_p2*Etrue))
        return NPE

    def getFitSNPE(self, Etrue):
        NPE = eloader.getQPE(Etrue, self.kB, self.Asct)
        return NPE



    def loadPrmBeta(self):
        filename = "../data/gamma/" + self.name + "_J19.txt"
        prmBeta, prmAntiBeta = bloader.loadPrmBeta(filename)
        return prmBeta, prmAntiBeta



    def loadPrmBeta(self):
        st = time.time()
        filename = ""
        if self.phys == "Penelope" :
            filename = "../data/gamma/" + self.name + "_pene.txt"
        if self.phys == "Livermore":
            filename = "../data/gamma/" + self.name + "_J19.txt"
            if self.name == "nC12":
                filename = "../data/gamma/" + self.name + "_J19_multiple.txt"
        #print("Primary beta distribution from %s" %filename)
        prmBeta, prmAntiBeta = bloader.loadPrmBeta(filename)
        et = time.time()
        self.loadPrmBetaFlag = True
        #print("Primary beta loading time : %.3f s" %(et-st))
        self.prmBetaArr = prmBeta
        self.prmAntiBetaArr = prmAntiBeta


    def calcSingleEvent(self):
        st = time.time()
        if self.loadPrmBetaFlag == False:
            self.loadPrmBeta()
        mu_arr, sigma_arr = [], []
        smu_arr, ssigma_arr = [], []
        cmu_arr, csigma_arr = [], []
        self.nSamples = len(self.prmBetaArr)
        for i in range(len(self.prmBetaArr)):
            tmpnpe, tmpspe = 0, 0
            tmpcnpe, tmpsnpe = 0, 0
            tmpcspe, tmpsspe = 0, 0
            for j in self.prmBetaArr[i]:
                if j == 0:
                    break
                if self.mode == "Truth":
                    tmpnpe += (eloader.getNPE(j))
                    tmpspe += (eloader.getSPE(j)**2)

                if self.mode == "Fit":
                    tmpnpe += (self.getFitNPE(j))
                    tmpspe += (self.getFitSPE(j)**2)
                    tmpcnpe += self.getFitCNPE(j)
                    tmpsnpe += self.getFitSNPE(j)



            for j in self.prmAntiBetaArr[i]:
                if j == 0:
                    break
                if self.mode == "Truth":
                    tmpnpe += (eloader.getNPE(j) + 2*660.8)
                    tmpspe += (eloader.getSPE(j)**2 + 27.07**2*2)
                if self.mode == "Fit":
                    tmpnpe += (self.getFitNPE(j) + 2*660.8)
                    tmpspe += (self.getFitSPE(j)**2 + 27.07**2*2)
                    tmpcnpe += (self.getFitCNPE(j) + 2.49)
                    tmpsnpe += (self.getFitSNPE(j) + 2*660.8-2.49)

            mu_arr.append(tmpnpe)
            sigma_arr.append(np.sqrt(tmpspe))
            cmu_arr.append(tmpcnpe)
            csigma_arr.append(np.sqrt(tmpcspe))
            smu_arr.append(tmpsnpe)
            ssigma_arr.append(np.sqrt(tmpsspe))

        mu_arr = np.array(mu_arr)
        sigma_arr = np.array(sigma_arr)
        cmu_arr = np.array(cmu_arr)
        csigma_arr = np.array(csigma_arr)
        smu_arr = np.array(smu_arr)
        ssigma_arr = np.array(ssigma_arr)

        et = time.time()
        #print("Primary beta dist calculation time : %.3f s" %(et-st))
        return mu_arr, sigma_arr, smu_arr, ssigma_arr, cmu_arr, csigma_arr

    def ModelPrediction(self):
        st = time.time()
        #mu, sigma = self.calcSingleEvent()
        mu, sigma, smu, ssigma, cmu, csigma = self.calcSingleEvent()
        self.npe, self.spe = 0, 0
        self.snpe, self.sspe = 0, 0
        self.cnpe, self.cspe = 0, 0
        self.sigma_part1, self.sigma_part2 = 0, 0
        for i in mu:
            self.npe += i
        self.npe = self.npe/self.nSamples
        for i, j in zip(mu, sigma):
            self.spe += (i-self.npe)**2 + j**2
            self.sigma_part1 += (i-self.npe)**2
            self.sigma_part2 += j**2

        for i in cmu:
            self.cnpe += i
        self.cnpe = self.cnpe/self.nSamples
        for i, j in zip(cmu, csigma):
            self.cspe += (i-self.cnpe)**2 + j**2


        for i in smu:
            self.snpe += i
        self.snpe = self.snpe/self.nSamples
        for i, j in zip(smu, ssigma):
            self.sspe += (i-self.snpe)**2 + j**2


        self.spe = np.sqrt(self.spe/self.nSamples)
        self.sigma_part1 /= self.nSamples
        self.sigma_part2 /= self.nSamples
        self.sspe = np.sqrt(self.sspe/self.nSamples)
        self.cspe = np.sqrt(self.cspe/self.nSamples)
        
        et = time.time()
        #print("total prediction time : %.3f s" %(et-st))
    

    def loadSimTruth(self):
        filename = "../data/gamma/spectrum/" + self.name + "_totpe.root"
        evt = up.open(filename)["evt"]
        totpe = evt["totalPE"].array()
        return totpe


    def calcTruth(self):
        st = time.time()
        totpe = self.loadSimTruth()
        low = np.min(totpe) - 50
        high = np.max(totpe) + 50
        f1 = TF1("f1", "gaus", low, high)
        h1 = TH1D(self.name+"h1", "", 100, low, high)
        for i in totpe:
            h1.Fill(i)
        h1.Fit(f1, "REQ")
        
        self.npe_sim = f1.GetParameter(1)
        self.npeerr_sim = f1.GetParError(1)
        self.spe_sim = f1.GetParameter(2)
        self.speerr_sim = f1.GetParError(2)

        et = time.time()
        #print("Truth loading time : %.3f s" %(et-st))


    def Compare(self):
        totpe = self.loadSimTruth()
        self.ModelPrediction()
        nSimEvt = len(totpe)
        minpe = np.min(totpe)
        maxpe = np.max(totpe)
        X = np.arange(minpe, maxpe, 1)
        calcpe = norm.pdf(X, loc=self.npe, scale=self.spe)
        plt.figure(0)
        plt.hist(totpe, bins=80, density=True, range=(minpe, maxpe), label="simulation")
        plt.plot(X, calcpe , "-", label="prediction")
        #plt.legend()
        plt.xlabel("# P.E.")
        plt.title(self.name)
        plt.savefig(self.name + "cpr.pdf")
        #plt.clf()
        #plt.show()


    def GetChi2(self):
        ModelPrediction()
        self.chi2 = 0
        self.chi2 += (self.npe - self.npe_sim)**2 / self.npeerr_sim**2

        return self.chi2
