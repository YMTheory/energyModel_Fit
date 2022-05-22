import numpy as np
import ROOT
from parameters import parameters

class electronResponse(object):

    init_flag = False
    kBmin, kBmax = 5.1e-3, 7.5e-3
    
    quenchNL_nbins = 1688
    quenchNL_E = np.zeros(quenchNL_nbins)
    quenchNL = np.zeros((int((kBmax - kBmin)*10000)+2, quenchNL_nbins))

    cer_nbins = 101
    Ncer_E = np.zeros(cer_nbins)
    Ncer = np.zeros(cer_nbins)

    # electronResponse parameters :
    Ysct = 1400
    kB = 6.5e-3
    kC = 1
    a, b, n = 0.9, 0.044, 1.3

    @staticmethod
    def initialize():
        print("------ Initialzed electron response -----")

        fquenchNL = ROOT.TFile(parameters.get_quenchNL_file(), "read")
        for ikb in range(51, 76, 1):
            hname = "kB" + str(ikb)
            hquenchNL = fquenchNL.Get(hname)
            for ibin in range(electronResponse.quenchNL_nbins):
                if ikb == 65:
                    electronResponse.quenchNL_E[ibin] = hquenchNL.GetBinCenter(ibin+1)
                electronResponse.quenchNL[ikb-51, ibin] = hquenchNL.GetBinContent(ibin+1)


        fquenchCer = ROOT.TFile(parameters.get_Ncer_file(), "read")
        gname = "Ncer"
        gNcer = fquenchCer.Get(gname)
        for i in range(gNcer.GetN()):
            electronResponse.Ncer_E[i] = gNcer.GetPointX(i)
            electronResponse.Ncer[i] = gNcer.GetPointY(i)

        electronResponse.init_flag = True



    @staticmethod
    def quenchNL_interpolate(E, kB):
        if kB <= electronResponse.kBmin:
            kB = electronResponse.kBmin
        if kB >= electronResponse.kBmax:
            kB = electronResponse.kBmax
        kBlow  = int((kB - electronResponse.kBmin)*10000)
        kBhigh = kBlow + 1 
        r = ((kB - electronResponse.kBmin )*10000 - kBlow) / (kBhigh - kBlow)

        quenchNL_low = np.interp(E, electronResponse.quenchNL_E, electronResponse.quenchNL[kBlow])
        quenchNL_high = np.interp(E, electronResponse.quenchNL_E, electronResponse.quenchNL[kBhigh])

        return r * quenchNL_high + (1-r)*quenchNL_low


    @staticmethod
    def get_Nsct(E, kB, Ysct):
        if not electronResponse.init_flag:
            electronResponse.initialize()
        nl = electronResponse.quenchNL_interpolate(E, kB)
        N = nl * Ysct * E
        return N


    @staticmethod
    def get_Ncer(E, kC):
        if not electronResponse.init_flag:
            electronResponse.initialize()
        N0 = np.interp(E, electronResponse.Ncer_E, electronResponse.Ncer)
        return N0 * kC



    @staticmethod
    def get_Nsigma(E):
        a, b, n = parameters.a, parameters.b, parameters.n
        return np.sqrt(a**2*E + b**2*np.power(E, n))





