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
        self.npe_sim = 0
        self.npeerr_sim = 0
        self.spe_sim = 0
        self.speerr_sim = 0
        self.prmBetaArr = []
        self.prmAntiBetaArr = []
        self.loadPrmBetaFlag = False
        self.chi2 = 0

    def getName(self):
        print("Current gamma is " + self.name)
        return self.name

    def getNPE(self):
        return self.npe

    def getSPE(self):
        return self.spe

    def getNPESim(self):
        return self.npe_sim

    def getNPEErrSim(self):
        return self.npeerr_sim

    def getSPESim(self):
        return self.spe_sim

    def getSPEErrSim(self):
        return self.speerr_sim

    def loadPrmBeta(self):
        filename = "../data/gamma/" + self.name + "_J19.txt"
        prmBeta, prmAntiBeta = bloader.loadPrmBeta(filename)
        return prmBeta, prmAntiBeta



    def loadPrmBeta(self):
        st = time.time()
        filename = "../data/gamma/" + self.name + "_J19.txt"
        prmBeta, prmAntiBeta = bloader.loadPrmBeta(filename)
        et = time.time()
        self.loadPrmBetaFlag = True
        print("Primary beta loading time : %.3f s" %(et-st))
        self.prmBetaArr = prmBeta
        self.prmAntiBetaArr = prmAntiBeta


    def calcSingleEvent(self):
        st = time.time()
        if self.loadPrmBetaFlag == False:
            self.loadPrmBeta()
        mu_arr, sigma_arr = [], []
        self.nSamples = len(self.prmBetaArr)
        for i in range(len(self.prmBetaArr)):
            tmpnpe, tmpspe = 0, 0
            for j in self.prmBetaArr[i]:
                if j == 0:
                    break
                tmpnpe += (eloader.getNPE(j))
                tmpspe += (eloader.getSPE(j)**2)
            for j in self.prmAntiBetaArr[i]:
                if j == 0:
                    break
                tmpnpe += (eloader.getNPE(j) + 2*660.8)
                tmpspe += (eloader.getSPE(j)**2 + 27.07**2*2)

            mu_arr.append(tmpnpe)
            sigma_arr.append(np.sqrt(tmpspe))

        et = time.time()
        print("Primary beta dist calculation time : %.3f s" %(et-st))
        return mu_arr, sigma_arr

    def ModelPrediction(self):
        st = time.time()
        mu, sigma = self.calcSingleEvent()
        self.npe, self.spe = 0, 0
        for i in mu:
            self.npe += i
        self.npe = self.npe/self.nSamples
        for i, j in zip(mu, sigma):
            self.spe += (i-self.npe)**2 + j**2

        self.spe = np.sqrt(self.spe/self.nSamples)
        
        et = time.time()
        print("total prediction time : %.3f s" %(et-st))
    

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
        print("Truth loading time : %.3f s" %(et-st))


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
