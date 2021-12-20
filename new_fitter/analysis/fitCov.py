import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import gridspec
import elecLoader as eloader

def func(x, p0, p1, p2):
    return p0 + p1*x + p2*x**2

from scipy.optimize import curve_fit 
def readSPEfile(name):
    data = pd.read_csv(name, sep=" ")
    ke = data["KE"].to_numpy()
    ntot = data["Ntot"].to_numpy()
    stot = data["Stot"].to_numpy()
    nsct = data["Nsct"].to_numpy()
    ssct = data["Ssct"].to_numpy()
    ncer = data["Ncer"].to_numpy()
    scer = data["Scer"].to_numpy()
    cov = data["Cov"].to_numpy()
    return ke, ntot, stot, nsct, ssct, ncer, scer, cov

def readFile(filename):
    Etrue, totpe, sigma, sigmaErr = [], [], [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Etrue.append(float(data[0]))
            totpe.append(float(data[1]))
            sigma.append(float(data[2]))
            sigmaErr.append(float(data[3]))

    Etrue = np.array(Etrue)
    totpe = np.array(totpe)
    sigma = np.array(sigma)
    sigmaErr = np.array(sigmaErr)

    return Etrue, totpe, sigma, sigmaErr


if __name__ == "__main__":

    es = 3134.078/2.223

    # electron
    covE, covPE, cov, covErr = readFile("../data/electron/elecPECov1.txt")
    sctE, sctPE, sctSigma, sctSigmaErr = readFile("../data/electron/elecSctPEResol1.txt")
    cerE, cerPE, cerSigma, cerSigmaErr = readFile("../data/electron/elecCerPEResol1.txt")
    resolE, totpe, resolData, resolerr = eloader.getResolArray()
    totpe = np.array(totpe)

    E, correlation, correlation_err = [], [], []

    sctNPE, cerNPE, covArr, covErrArr = [], [], [], []
    sctSPE, cerSPE, sctSPEerr, cerSPEerr = [], [], [], []

    for i in range(55, len(covE), 10):
        if sctSigma[i] == 0 or cerSigma[i] == 0:
            continue
        E.append(covE[i])
        sctNPE.append(sctPE[i])
        cerNPE.append(cerPE[i])
        sctSPE.append(sctSigma[i])
        cerSPE.append(cerSigma[i])
        sctSPEerr.append(sctSigmaErr[i])
        cerSPEerr.append(cerSigmaErr[i])
        covArr.append(cov[i])
        covErrArr.append(covErr[i])
        correlation.append(cov[i]/sctSigma[i]/cerSigma[i])
        correlation_err.append( np.sqrt( covErr[i]**2/(sctSigma[i]*cerSigma[i])**2  \
            + (cov[i]/sctSigma[i])**2 * cerSigmaErr[i]**2/cerSigma[i]**4 +          \
            + (cov[i]/cerSigma[i])**2 * sctSigmaErr[i]**2/sctSigma[i]**4 ) )

    covArr = np.array(covArr)
    cerNPE = np.array(cerNPE)
    sctNPE = np.array(sctNPE)
    cerSPE = np.array(cerSPE)
    cerSPEerr = np.array(cerSPEerr)
    sctSPE = np.array(sctSPE)
    sctSPEerr = np.array(sctSPEerr)
    correlation = np.array(correlation)
    correlation_err = np.array(correlation_err)

    ntot1 = sctNPE + cerNPE

    ke0, ntot0, stot0, nsct0, ssct0, ncer0, scer0, cov0 = readSPEfile("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/gamma/singles/gammaSPE.txt")
    ke2, ntot2, stot2, nsct2, ssct2, ncer2, scer2, cov2 = readSPEfile("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/positron/positronSPE.txt")

    popt, pcov = curve_fit(func, ntot2/es, cov2)
    print(popt)
    plt.plot(ntot2/es, cov2)
    plt.plot(ntot2/es, func(ntot2/es, *popt))

    plt.show()






