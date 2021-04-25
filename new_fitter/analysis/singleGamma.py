import numpy as np
import elecLoader as eloader
import prmBetaLoader as bloader
import uproot as up
from scipy.stats import norm
import matplotlib.pyplot as plt

class singleGamma(object):
   
    def __init__(self, name):
        self.name = name
        self.nSamples = 5000
        self.npe = 0
        self.spe = 0

    def getName(self):
        print("Current gamma is " + self.name)
        return self.name

    def getNPE(self):
        return self.npe

    def getSPE(self):
        return self.spe


    def loadPrmBeta(self):
        filename = "../data/gamma/" + self.name + "_J19.txt"
        prmBeta, prmAntiBeta = bloader.loadPrmBeta(filename)
        return prmBeta, prmAntiBeta


    def calcSingleEvent(self):
        prmBeta, prmAntiBeta = self.loadPrmBeta()
        mu_arr, sigma_arr = [], []
        self.nSamples = len(prmBeta)
        for i in range(len(prmBeta)):
            tmpnpe, tmpspe = 0, 0
            for j in prmBeta[i]:
                if j == 0:
                    break
                tmpnpe += (eloader.getNPE(j))
                tmpspe += (eloader.getSPE(j)**2)
            for j in prmAntiBeta[i]:
                if j == 0:
                    break
                tmpnpe += (eloader.getNPE(j) + 2*660.8)
                tmpspe += (eloader.getSPE(j)**2 + 27.07**2*2)

            mu_arr.append(tmpnpe)
            sigma_arr.append(np.sqrt(tmpspe))

        return mu_arr, sigma_arr

    def ModelPredictino(self):
        mu, sigma = self.calcSingleEvent()
        self.npe, self.spe = 0, 0
        for i, j in zip(mu, sigma):
            self.npe += i
            self.spe += j**2

        self.npe = self.npe/self.nSamples
        self.spe = np.sqrt(self.spe/self.nSamples)
        
    

    def loadSimTruth(self):
        filename = "../data/gamma/spectrum/" + self.name + "_totpe.root"
        evt = up.open(filename)["evt"]
        totpe = evt["totalPE"].array()
        return totpe


    def Compare(self):
        totpe = self.loadSimTruth()
        self.ModelPredictino()
        nSimEvt = len(totpe)
        minpe = np.min(totpe)
        maxpe = np.max(totpe)
        X = np.arange(minpe, maxpe, 1)
        calcpe = norm.pdf(X, loc=self.npe, scale=self.spe)
        nbins = maxpe - minpe + 100
        plt.hist(totpe, bins=nbins, range=(minpe, maxpe))
        plt.plot(X, calcpe, "-")
        plt.show()