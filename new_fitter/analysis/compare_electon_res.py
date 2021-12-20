import numpy as np
import matplotlib.pyplot as plt

def loadlec():
    evis, totpe, totpeerr, sigma, sigmaerr = [], [], [], [], []
    with open("../data/electron/elecResol1.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            evis.append(float(data[0]))
            totpe.append(float(data[1]))
            totpeerr.append(float(data[2]))
            sigma.append(float(data[3]))
            sigmaerr.append(float(data[4]))

    return evis, totpe, totpeerr, sigma, sigmaerr


def readGraph(ff):
    f1 = ROOT.TFile(ff)
    gg = f1.Get("nonl")
    E, totpe = [], []
    for i in range(gg.GetN()):
        E.append(gg.GetPointX(i))
        totpe.append(gg.GetPointY(i))
    return E, totpe


def res_func(p1, p2, x):
    return np.sqrt(p1*x+p2*x**2)

def readFile(filename):
    totpe, sigma, sigmaerr = [], [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            totpe.append(float(data[1]))
            sigma.append(float(data[3]))
            sigmaerr.append(float(data[4]))
    return totpe, sigma, sigmaerr


import ROOT
def main_nonl():
    evis, totpe, totpeerr, sigma, sigmaerr = loadElec()
    eb12, totpeb12 = readGraph("../elec_nonl_onlyB12.root")
    #egam, totpegam = readGraph("../elec_nonl_onlyGam.root")
    eboth, totpeboth = readGraph("../elec_nonl_Both.root")

    
    evis = np.array(evis)
    totpe = np.array(totpe)
    eb12 = np.array(eb12)
    totpeb12 = np.array(totpeb12)
    #egam = np.array(egam)
    #totpegam = np.array(totpegam)
    eboth = np.array(eboth)
    totpeboth = np.array(totpeboth)

    scale = 1409.8

    plt.plot(evis, totpe/scale/evis, "-", label="sim")
    plt.plot(eb12, totpeb12/scale/eb12, '--', label="fitting onlyB12")
    #plt.plot(egam, totpegam/scale/egam, '--', label="fitting onlyGam")
    #plt.plot(eboth, totpeboth/scale/eboth, '--', label="fitting both")
    
    plt.xlim(0.1, 15)
    plt.semilogx()
    plt.semilogy()
    plt.legend()
    plt.show()

import ROOT
def resolFit(npe, resol, resolerr):
    ge = ROOT.TGraphErrors()
    for i in range(len(npe)):
        ge.SetPoint(i, npe[i], resol[i]*resol[i])
        ge.SetPointError(i, 0, resolerr[i]*2*resol[i])
    
    ff = ROOT.TF1("ff", "[0]*x+[1]*x*x", 0, 23000)
    ge.Fit(ff)

    p1 = ge.GetFunction("ff").GetParameter(0)
    p2 = ge.GetFunction("ff").GetParameter(1)

    return p1, p2



def main_res():
    #evis, totpe, totpeerr, sigma, sigmaerr = loadElec()
    #sigma = np.array(sigma)
    #totpe1, sigma1 = [], []
    #for i,j in zip(totpe, sigma):
    #    if i > 12000:
    #        totpe1.append(i)
    #        sigma1.append(j)
    #sigma1 = np.array(sigma1)

    totpe1, sigma1, sigmaerr1 = readFile("resol.txt")
    totpe1 = np.array(totpe1)
    sigma1 = np.array(sigma1)
    sigmaerr1 = np.array(sigmaerr1)

    fitp1, fitp2 = resolFit(totpe1, sigma1, sigmaerr1)

    
    npe_calc = np.arange(1000, 30000, 100)
    sigma_onlygam = res_func(0.979, 6.13e-5, npe_calc)
    sigma_onlyb12 = res_func(8.41665e-01, 6.73931e-05, npe_calc)
    sigma_both = res_func(0.979, 6.038e-5, npe_calc)
    sigma_bestFit = res_func(fitp1, fitp2, npe_calc)
    sigma_highGam = res_func(9.97128e-01, 4.13758e-05, npe_calc)

    # LS
    npe_mic = [7.76149e+04]
    sigma_mic = [8.37856e+02]
    sigmaerr_mic = [2.49857e+02]

    # water
    #npe_mic = [8.07571e+02]
    #sigma_mic = [6.43617e+01]
    #sigmaerr_mic = [1.62765e+01]

    npe_mic = np.array(npe_mic)
    sigma_mic = np.array(sigma_mic)
    sigmaerr_mic = np.array(sigmaerr_mic)

    plt.errorbar(totpe1, sigma1/totpe1, yerr=sigmaerr1/totpe1, fmt="o", ms=4)
    #plt.plot(npe_calc, sigma_onlygam, label="onlygam")
    #plt.plot(npe_calc, sigma_onlyb12/npe_calc, label="onlyB12")
    plt.plot(npe_calc, sigma_both/npe_calc, label="only gammas below 7MeV")
    #plt.plot(npe_calc, sigma_bestFit/npe_calc, label="only gamma below 7MeV")
    #plt.errorbar(npe_mic, sigma_mic/npe_mic/np.sqrt(2), yerr=sigmaerr_mic/npe_mic, fmt='v', label="Michel electron")
    plt.plot(npe_calc, sigma_highGam/npe_calc, "--", label="w/ 15MeV gamma")

    #plt.semilogx()
    plt.xlabel(r"$N_{totpe}$")
    #plt.ylabel(r"$\sigma_{NPE}$")
    plt.ylabel(r"$\sigma(N)/N$")
    plt.legend()
    plt.grid(True)
    plt.savefig("resol_15MeVGam.pdf")
    plt.show()

if __name__ == "__main__":
    #main_nonl()
    main_res()
