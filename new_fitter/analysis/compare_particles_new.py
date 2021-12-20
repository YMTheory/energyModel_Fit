import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import gridspec
import elecLoader as eloader

def func(x, a, b, c):
    return a+b*x+c*x**2

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

    # gamma and positrons
    ke0, ntot0, stot0, nsct0, ssct0, ncer0, scer0, cov0 = readSPEfile("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/gamma/singles/gammaSPE.txt")
    ke2, ntot2, stot2, nsct2, ssct2, ncer2, scer2, cov2 = readSPEfile("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/positron/positronSPE.txt")


    print(len(ke0), len(ntot1), len(ntot2))
    num = [5, 15, 25, 35, 45, 55, 65, 74]

    subntot0, subntot1, subntot2 = [], [], []
    subssct0, subssct1, subssct2 = [], [], []
    subscer0, subscer1, subscer2 = [], [], []
    subcov0, subcov1, subcov2 = [], [], []
    for i in num:
        subntot0.append(ntot0[i])
        subssct0.append(ssct0[i])
        subscer0.append(scer0[i])
        subcov0.append(cov0[i])
        subntot1.append(ntot1[i])
        subssct1.append(sctSPE[i])
        subscer1.append(cerSPE[i])
        subcov1.append(covArr[i])
        subntot2.append(ntot2[i])
        subssct2.append(ssct2[i])
        subscer2.append(scer2[i])
        subcov2.append(cov2[i])

    subntot0 = np.array(subntot0)
    subssct0 = np.array(subssct0)
    subscer0 = np.array(subscer0)
    subcov0  = np.array(subcov0)
    subntot1 = np.array(subntot1)
    subssct1 = np.array(subssct1)
    subscer1 = np.array(subscer1)
    subcov1  = np.array(subcov1)
    subntot2 = np.array(subntot2)
    subssct2 = np.array(subssct2)
    subscer2 = np.array(subscer2)
    subcov2  = np.array(subcov2)


    fig = plt.figure(constrained_layout=True, figsize=(6, 10))
    gs = fig.add_gridspec(3, 1)

    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[2, 0])

    es = 3134.078/2.223


    ax0.plot(subntot0/es, subssct0**2/subntot0**2,  "o-",  color="seagreen", lw=2, ms=7, label="Gamma")
    ax0.plot(subntot1/es, subssct1**2/subntot1**2, "X-",  color="coral",    lw=2, ms=7, label="Electron") 
    ax0.plot(subntot2/es, subssct2**2/subntot2**2,  "d-",  color="blue",     lw=2, ms=7, label="Positron")
    ax0.set_xlabel(r"$E_{vis}$[MeV]", fontsize=13)
    ax0.set_ylabel(r"$(\sigma_s/E_{vis})^2$", fontsize=13)
    ax0.legend(prop={'size':13})
    ax0.grid(True)
    ax0.tick_params(axis='both', which='major', labelsize=11)
    ax0.set_title("(a)", fontsize=12)
    #ax0.semilogy()
    ax0.ticklabel_format(style='sci', scilimits=(-1,2), axis='y')

    ax1.plot(subntot1/es, subscer1**2/subntot1**2, "X--", color="coral", lw=2, ms=7)
    ax1.plot(subntot0/es, subscer0**2/subntot0**2, "o--", color="seagreen", lw=2, ms=7)
    ax1.plot(subntot2/es, subscer2**2/subntot2**2, "d--", color="blue", lw=2, ms=7, label="Cherenkov")
    #ax1.plot(ntot1/es, 2*covArr/ntot1**2, "-.", color="coral", lw=2)
    ax1.set_xlabel(r"$E_{vis}$[MeV]", fontsize=13)
    ax1.set_ylabel(r"$(\sigma_C/E_{vis})^2$", fontsize=13)
    #ax1.legend(prop={'size':13})
    ax1.tick_params(axis='both', which='major', labelsize=11)
    ax1.grid(True)
    ax1.set_title("(b)", fontsize=12)
    ax1.ticklabel_format(style='sci', scilimits=(-1,2), axis='y')
    #ax1.semilogy()

    #ax2.plot(ntot2/es, 2*cov2/ntot2**2, "-.", color="blue", lw=2, label="Covariance")
    ax2.plot(subntot0/es, func(subntot0/es, -97.21230398, 193.87142506, 1.37056115)/subntot0**2, "o-.", color="seagreen", lw=2, ms=7, label="Gamma")
    ax2.plot(subntot2/es, func(subntot2/es, -65.41740609,  31.58792588,  13.31949795)/(subntot2)**2, "d-.", color="blue", lw=2, ms=7, label="Positron")
    ax2.plot(subntot1/es, func(subntot1/es, -78.97563212, 97.33498604, 7.603084)/(subntot1)**2, "X-.", color="coral", lw=2, ms=7, label="Electron")
    #par_name = [str(u) for u in range(1, 11, 1)]
    #ax2.set_xticks(np.arange(1, 11, 1))
    #ax2.set_xticklabels(par_name, fontsize=11)
    ax2.grid(True)
    ax2.set_xlabel(r"$E_{vis}$[MeV]", fontsize=13)
    ax2.set_ylabel(r"cov/$(E_{vis})^2$", fontsize=13)
    #ax2.legend(prop={'size':13}, ncol=2)
    ax2.tick_params(axis='both', which='major', labelsize=11)
    ax2.set_title("(c)", fontsize=12)
    #ax2.semilogy()


    plt.tight_layout()
    plt.savefig("decompResolThreePar.pdf")
    plt.show()




