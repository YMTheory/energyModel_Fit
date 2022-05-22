import numpy as np
import uproot as up
import ROOT

from electronResponse import electronResponse
from parameters import parameters

class gamma(object):

    gamma_name = "gamma"
    gamma_Edep = 1
    
    nSecMax = 100
    secondaries_elec = np.zeros((parameters.get_Ngamma_samples(), nSecMax))
    secondaries_posi = np.zeros((parameters.get_Ngamma_samples(), nSecMax))

    npe_mean = np.zeros(parameters.get_Ngamma_samples())
    npe_sigma = np.zeros(parameters.get_Ngamma_samples())

    pred_npe_mean = 0
    pred_npe_sigma = 0

    npe_data = None
    load_npe_spectrum_flag = False
    load_gamma_samples_flag = False

    def __init__(self, name, Edep):
        self.gamma_name = name
        self.gamma_Edep = Edep


    def get_npe_data(self):
        if not self.load_npe_spectrum_flag:
            self.load_npe_spectrum()
        return self.npe_data

    def get_pred_npe_mean(self):
        return self.pred_npe_mean

    def get_pred_npe_sigma(self):
        return self.pred_npe_sigma


    def load_npe_spectrum(self):
        filename = parameters.get_gamma_spectrum(self.gamma_name)
        ff = up.open(filename)
        self.npe_data = ff["photon"]["totPE"].array()

        self.load_npe_spectrum_flag = True

    def load_gamma_samples(self):
        filename = parameters.get_gamma_samples(self.gamma_name)
        ff = ROOT.TFile(filename, "read")
        helec = ff.Get(self.gamma_name+"_elec")
        hposi = ff.Get(self.gamma_name+"_posi")
        for i in range(parameters.get_Ngamma_samples()) :
            for j in range(self.nSecMax):
               self.secondaries_elec[i, j] = helec.GetBinContent(i+1, j+1) 
               self.secondaries_posi[i, j] = hposi.GetBinContent(i+1, j+1) 
        self.load_gamma_samples_flag = True

    
    def preCalculation(self):
        if not self.load_gamma_samples_flag:
            self.load_gamma_samples()
        Y = parameters.get_Y()
        Ysct = parameters.get_Ysct()
        kB = parameters.get_kB()
        kC = parameters.get_kC()
        a = parameters.get_res_a()
        b = parameters.get_res_b()
        n = parameters.get_res_n()

        N = parameters.get_Ngamma_samples()
        for i in range(N):
            tmp_pe, tmp_sigma = 0, 0

            for ielec in range(self.nSecMax):
                tmp_E = self.secondaries_elec[i, ielec]
                if tmp_E == 0:
                    break
                # working in the NPE level
                single_npe = electronResponse.get_Ncer(tmp_E, kC) + electronResponse.get_Nsct(tmp_E, kB, Ysct)
                tmp_pe += single_npe
                tmp_sigma += electronResponse.get_Nsigma(single_npe)**2


            for iposi in range(self.nSecMax):
                tmp_E = self.secondaries_posi[i, iposi]
                if tmp_E == 0:
                    break
                # working in the NPE level
                single_npe = electronResponse.get_Ncer(tmp_E, kC) + electronResponse.get_Nsct(tmp_E, kB, Ysct)
                tmp_pe += single_npe
                tmp_sigma += electronResponse.get_Nsigma(single_npe)**2

                tmp_pe += parameters.get_Ge68mean()
                tmp_sigma += parameters.get_Ge68sigma()**2

            self.npe_mean[i] = tmp_pe
            self.npe_sigma[i] = np.sqrt(tmp_sigma)


    def prediction(self):
        nSamples = parameters.get_Ngamma_samples()
        self.preCalculation()

        pred_mean, pred_sigma = 0, 0
        for i in range(nSamples):
            pred_mean += self.npe_mean[i]
        pred_mean /= nSamples

        for i in range(nSamples):
            pred_sigma += (self.npe_mean[i] - pred_mean)**2 + self.npe_sigma[i]**2
        pred_sigma = np.sqrt(pred_sigma / nSamples)

        self.pred_npe_mean = pred_mean
        self.pred_npe_sigma = pred_sigma




    










