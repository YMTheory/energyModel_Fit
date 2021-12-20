import numpy as np
from matplotlib import gridspec
import matplotlib.pyplot as plt
import ROOT

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


def param(npe):
    p1, p2 = 2.34299e+00, 4.95424e-03
    return p1*npe + p2*npe**2


def main():
    covE, covPE, cov, covErr = readFile("../data/electron/elecPECov1.txt")
    sctE, sctPE, sctSigma, sctSigmaErr = readFile("../data/electron/elecSctPEResol1.txt")
    cerE, cerPE, cerSigma, cerSigmaErr = readFile("../data/electron/elecCerPEResol1.txt")

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

    #plt.errorbar(E, correlation*E/cerNPE, yerr=correlation_err*E/cerNPE, fmt="-", color="chocolate")
    #plt.plot(E, correlation/np.sqrt(cerNPE)*50, "-", color="chocolate")

    #plt.plot(E, covArr/np.sqrt(sctNPE*cerNPE) )

    g1 = ROOT.TGraphErrors()
    for i in range(len(covArr)):
        g1.SetPoint(i, cerNPE[i]+sctNPE[i], covArr[i])
        g1.SetPointError(i, 0, covErrArr[i])

    f1 = ROOT.TF1("f1", "[0]+[1]*x+[2]*x*x", 0, 15000)
    f1.FixParameter(0, 0)
    g1.Fit(f1, "RE")
    q0 = f1.GetParameter(0)
    q1 = f1.GetParameter(1)
    q2 = f1.GetParameter(2)

    npe1 = np.arange(0, 15000, 100)
    fit1 = q0 + q1*npe1 + q2*npe1**2
    

    plt.errorbar(sctNPE+cerNPE, covArr, yerr=covErrArr, fmt="o")
    plt.plot(npe1, fit1, "-")
    plt.show()
    
    """
    pois = []
    for i in np.arange(sctNPE.min()-10, sctNPE.max()+10, 100):
        pois.append(i)
    poisC = []
    for i in np.arange(cerNPE.min()-10, cerNPE.max()+10, 100):
        poisC.append(i)
    poisC = np.array(poisC)

    fig = plt.figure(figsize=(12, 5))
    spec = gridspec.GridSpec(ncols=2, nrows=1)

    ax1 = fig.add_subplot(spec[0])
    ax3 = fig.add_subplot(spec[1])
    
    es = 3134.078/2.223
    sctNPE = np.array(sctNPE)
    sctSPE = np.array(sctSPE)
    cerNPE = np.array(cerNPE)
    cerSPE = np.array(cerSPE)
    pois = np.array(pois)
    poisC = np.array(poisC)

    ax1.errorbar(sctNPE/es, sctSPE**2/es**2, yerr=sctSPE*2*sctSPEerr/es/es, fmt="o", ms=4, color='blue')
    ax1.plot(pois/es, pois/es/es, "--", color="darkviolet", lw=2, label="Poisson")
    ax1.errorbar(cerNPE/es, cerSPE**2/es**2, yerr=cerSPE*2*cerSPEerr/es/es, fmt="o", ms=4, color='seagreen', zorder=1)
    #ax1.plot(poisC/es, poisC/es/es, "--", color="darkviolet", label="Poisson")
    ax1.plot(poisC/es, param(poisC)/es/es, "-", lw=2, color="orange", label="Parameterisation", zorder=2)
    ax1.text(4, 0.0022, r"$\sigma_s$", color="blue", fontsize=16)
    ax1.text(1, 0.0025, r"$\sigma_C$", color="seagreen", fontsize=16)
    #ax1.semilogx()
    #ax2 = ax1.twiny()
    #ax2.errorbar(cerNPE/es, cerSPE**2/es**2, yerr=cerSPE*2*cerSPEerr, fmt="o", ms=4, color='seagreen', label="Simulation", zorder=1)
    #ax2.plot(poisC/es, poisC/es, "--", color="darkviolet", label="Poisson")
    #ax2.plot(poisC, param(poisC), "-", lw=2, color="chocolate", label="Parameterisation", zorder=2)
    ax1.set_xlabel("NPE/A [MeV]", fontsize=17, color='black')
    ax1.set_ylabel(r"$\sigma^2$[MeV$^2$]", fontsize=17, labelpad=0.2)
    #ax2.set_xlabel(r"$N_{Cer}$", fontsize=14, color='seagreen')
    ax1.tick_params(axis='x', which='major', labelsize=16, labelcolor='black')
    ax1.tick_params(axis='y', which='major', labelsize=16, labelcolor='black')
    ax1.grid(axis='both')
    ax1.legend(loc="upper left", prop={'size' : 15})
    #ax2.tick_params(axis='x', which='major', labelsize=12, labelcolor='seagreen')
    #ax2.legend(loc="upper left", prop={'size' : 14})


    ax3.errorbar(E, correlation, yerr=correlation_err, fmt="o", ms=5, color="peru")

    ax3.grid(True)
    ax3.set_ylim(-0.1, 0.5)
    ax3.set_xlabel(r"$E_{dep}$/MeV", fontsize=17)
    ax3.set_ylabel("correlation coefficient", fontsize=17)
    ax3.tick_params(axis='both', which='major', labelsize=16)

    plt.tight_layout()
    plt.savefig("stddev+cov.pdf")

    plt.show()
    


    fig, ax = plt.subplots()
    ax.errorbar(E, correlation, yerr=correlation_err, fmt="o", color="#DA4510")
    #plt.errorbar(sctNPE+cerNPE, correlation, yerr=correlation_err, fmt="o", color="#FF7447")

    ax.grid(True)
    #plt.xlabel("NPE", fontsize=14)
    ax.set_xlabel(r"$E_{dep}$/MeV", fontsize=14)
    ax.set_ylabel("correlation coefficient", fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)

    plt.savefig("Corr2Edep.pdf")

    plt.show()

    """



if __name__ == "__main__":
    main()
