import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import ROOT
import elecLoader as eloader
from scipy.optimize import curve_fit

def func(x, p0, p1, p2):
    y = p0**2 + p1**2*x + p2**2*x**2
    return y
    #if y<=0:
    #    return 0
    #else:
    #    return y

def rfunc(x, a, b, c):
    return np.sqrt(a**2/x+b**2+c**2/x**2)

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


es = 3134.078/2.223
#c0, c1, c2 = -1.75e2, 4.46, 3.12e-8
#d0, d1, d2 = 0, 1.43288e-02, 9.34006e-06
c0, c1, c2 = 0, 1.774, 0.0577  #0, 5.318e-5, 9.52e-2
d0, d1, d2 = 0, 0.207, 2.376e-3

# No michel e-
#c0, c1, c2 = -76, 4.84, 4.74e-3
#d0, d1, d2 = 0, 4.5e-12, 3.6e-11
# Michel e-
#c0, c1, c2 = 0, 1.40598e+00, 3.06888e-03
#d0, d1, d2 = 0, 0, 1.94790e-05

def main():
    
    ### Load Truth
    covE, covPE, cov, covErr = readFile("../data/electron/elecPECov1.txt")
    sctE, sctPE, sctSigma, sctSigmaErr = readFile("../data/electron/elecSctPEResol1.txt")
    cerE, cerPE, cerSigma, cerSigmaErr = readFile("../data/electron/elecCerPEResol1.txt")
    resolE, totpe, resolData, resolerr = eloader.getResolArray()
    totpe = np.array(totpe)

    E, correlation, correlation_err = [], [], []

    sctNPE, cerNPE, covArr, covErrArr = [], [], [], []
    sctSPE, cerSPE, sctSPEerr, cerSPEerr = [], [], [], []

    for i in range(55, len(covE), 20):
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

    ## Fitting Results :
    #fig = plt.figure(figsize=(12, 12))
    #spec = gridspec.GridSpec(ncols=2, nrows=2)
    #fig = plt.figure(constrained_layout=True, figsize=(12, 10))
    #gs = fig.add_gridspec(3, 2)
    fig = plt.figure(constrained_layout=True, figsize=(6, 8))
    gs = fig.add_gridspec(3, 1)

    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[2, 0])
    #ax3 = fig.add_subplot(gs[:, 1])
    
    ax0.errorbar(sctNPE/es, sctSPE/sctNPE, yerr=sctSPEerr/sctNPE, fmt="o", \
                color="crimson", ms=7, lw=2, fillstyle='none', zorder=2, label="Simulation truth")
    ax0.plot(sctNPE/es, 1./np.sqrt(sctNPE), color="black", lw=2, zorder=1, label="Poisson")
    ax0.grid(True)
    ax0.set_xlabel(r"Electron $E_{sct}$[MeV]", fontsize=13)
    ax0.set_ylabel(r"$\sigma_{sct}/E_{sct}$", fontsize=13)
    ax0.legend(prop={'size':13})
    ax0.tick_params(axis='both', which='major', labelsize=11)
    ax0.set_title("(a)", fontsize=12)


    Ncer = []
    cerSigma2 = []
    for i in cerNPE:
        cc = func(i, c0, c1, c2)
        if i/es > 0.01:
            Ncer.append(i)
            cerSigma2.append(cc)
    Ncer = np.array(Ncer)
    cerSigma2 = np.array(cerSigma2)
    popt, pcov = curve_fit(rfunc, Ncer, np.sqrt(cerSigma2)/Ncer)
    print(popt)
   
    ax1.errorbar(cerNPE/es, cerSPE/cerNPE, yerr=cerSPEerr/cerNPE, fmt="o", \
                color="crimson", ms=7, lw=2, mfc="w", zorder=1, label="Simulation truth")
    ax1.plot(cerNPE/es, 1./np.sqrt(cerNPE), color="black", lw=2, zorder=2, label="Poisson")
    #ax1.plot(Ncer/es, np.sqrt(cerSigma2)/Ncer, color="blue", lw=2, zorder=3, label="Best fit")
    #ax1.plot(cerNPE/es, rfunc(cerNPE, *popt), color="blue", lw=2, zorder=3, label="Best fit")
    ax1.grid(True)
    ax1.set_xlabel(r"Electron $E_{Cer}$[MeV]", fontsize=13)
    ax1.set_ylabel(r"$\sigma_{Cer}/E_{Cer}$", fontsize=13)
    ax1.legend(prop={'size':13})
    ax1.tick_params(axis='both', which='major', labelsize=11)
    ax1.set_title("(b)", fontsize=12)
   


    Ntot = []
    covSigma2 = []
    for i, j in zip(sctNPE, cerNPE):
        cc = func(i+j, d0, d1, d2)
        if i/es > 0.1:
            Ntot.append(i)
            covSigma2.append(cc)
    Ntot = np.array(Ntot)
    covSigma2 = np.array(covSigma2)


    covArr = np.array(covArr)
    covErrArr = np.array(covErrArr)
    ax2.errorbar((sctNPE+cerNPE)/es, covArr/(sctSPE*cerSPE), yerr=np.sqrt((covErrArr/sctSPE/cerSPE)**2+(covArr/cerSPE/sctSPE/sctSPE*sctSPEerr)**2 + (covArr/sctSPE/cerSPE/cerSPE/cerSPEerr)**2), fmt="o", \
                color="crimson", ms=7, mfc="w", zorder=1, label="Simulation truth")
    #ax2.plot((sctNPE+cerNPE)/es, 1./np.sqrt(cerNPE+sctNPE), color="black", lw=2, zorder=2, label="Poisson")
    #ax2.plot(Ntot/es, covSigma2, color="blue", lw=2, zorder=3, label="Best fit")
    ax2.grid(True)
    ax2.set_xlabel(r"Electron $E_{vis}$[MeV]", fontsize=13)
    ax2.set_ylabel(r"cov$/\sigma_{Cer}/\sigma_{sct}$", fontsize=13)
    ax2.legend(loc="lower right", prop={'size':13})
    ax2.tick_params(axis='both', which='major', labelsize=11)
    ax2.set_ylim(-0.1, 0.5)
    ax2.set_title("(c)", fontsize=12)

    
    Emin, Emax = 0.5, 8
    num = 0
    totpe1, resolData1, resolErr1 = [], [], []
    for i, j, k in zip(totpe, resolData, resolerr):
        num += 1
        if num%5 != 0:
            continue
        if i/es > Emin and i/es<Emax:
            totpe1.append(i)
            resolData1.append(j)
            resolErr1.append(k)
    totpe1 = np.array(totpe1)
    resolData1 = np.array(resolData1)
    resolErr1 = np.array(resolErr1)

    #ax3.errorbar(totpe1/es, resolData1**2/totpe1**2, yerr=resolErr1**2/totpe1**2, fmt="o", color="red", ms=7, mfc="w", zorder=1, label="Simulation truth")

    #ax3.plot((sctNPE+cerNPE)/es, (sctNPE+rfunc(cerNPE, *popt)*cerNPE + 2*func(cerNPE+sctNPE, d0, d1, d2))/(sctNPE+cerNPE)**2, "-", lw=2, color="black", zorder=4, label="Best Fit" )
    #ax3.fill_between((sctNPE+cerNPE)/es, sctNPE/(sctNPE+cerNPE)**2+rfunc(cerNPE, *popt)*cerNPE/(sctNPE+cerNPE)**2+2*func(cerNPE+sctNPE, d0, d1, d2)/(sctNPE+cerNPE)**2, rfunc(cerNPE, *popt)*cerNPE/(sctNPE+cerNPE)**2+2*func(cerNPE+sctNPE, d0, d1, d2)/(sctNPE+cerNPE)**2, hatch="/", edgecolor="dimgray", color="dimgray", lw=1, fc="w", zorder=1, label="Scintillation")
    #ax3.fill_between((sctNPE+cerNPE)/es, rfunc(cerNPE, *popt)*cerNPE/(sctNPE+cerNPE)**2+2*func(cerNPE+sctNPE, d0, d1, d2)/(cerNPE+sctNPE)**2, 2*func(cerNPE+sctNPE, d0, d1, d2)/(sctNPE+cerNPE)**2, color="blue", fc="w", hatch="|", edgecolor="blue", lw=1, zorder=2, label="Cherenkov")
    #ax3.fill_between((sctNPE+cerNPE)/es, 2*func(cerNPE+sctNPE, d0, d1, d2)/(sctNPE+cerNPE)**2, color="seagreen", edgecolor="seagreen", fc="w", lw=1, hatch="-", zorder=3, label="Covariance")
    ##ax3.grid(True)
    #ax3.set_xlabel(r"Electron $E_{vis}$[MeV]", fontsize=13)
    #ax3.set_ylabel(r"$(\sigma_E/E_{vis})^2$", fontsize=13)
    #ax3.legend(prop={'size':13})
    #ax3.tick_params(axis='both', which='major', labelsize=11)
    #ax3.semilogy()
    #ax3.set_title("(d)", fontsize=12)

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=0.01, hspace=0.02)



    #fig1, ax4 = plt.subplots()

    #totNPE = sctNPE + cerNPE

    ##ax4.errorbar(totpe1/es, (resolData1/totpe1)**2, yerr=(resolErr1/totpe1)**2, fmt="o", color="red", ms=7, mfc="w", zorder=1, label="Simulation truth")
    ##ax4.plot((sctNPE+cerNPE)/es, (sctNPE+rfunc(cerNPE, *popt)*cerNPE + 2*func(cerNPE+sctNPE, d0, d1, d2))/(sctNPE+cerNPE)**2, "-", lw=2, color="black", zorder=2, label="Best Fit" )

    ##ax4.errorbar(totNPE/es, (cerSPE/totNPE)**2, yerr=(cerSPEerr/totNPE)**2, fmt="^", \
    ##            color="crimson", ms=7, lw=2, mfc="w", zorder=1, label="Simulation truth")
    ##ax4.plot(totNPE/es, rfunc(cerNPE, *popt)**2*cerNPE**2/totNPE**2, color="black", lw=2, zorder=3, label="Best fit")

    ##ax4.errorbar(totNPE/es, (sctSPE/totNPE)**2, yerr=(sctSPEerr/totNPE)**2, fmt="s", \
    ##            color="crimson", ms=7, lw=2, mfc="w", zorder=1, label="Simulation truth")
    ##ax4.plot(totNPE/es, sctNPE/totNPE**2, color="black", lw=2, zorder=3, label="Best fit")



    #ax4.errorbar(totpe1/es, (resolData1)**2, yerr=(resolErr1)**2, fmt="o", color="red", ms=7, mfc="w", zorder=1, label="Simulation truth")
    #ax4.plot((sctNPE+cerNPE)/es, (sctNPE+rfunc(cerNPE, *popt)*cerNPE + 2*func(cerNPE+sctNPE, d0, d1, d2)), "-", lw=2, color="black", zorder=2, label="Best Fit" )

    #ax4.errorbar(totNPE/es, (cerSPE)**2, yerr=(cerSPEerr)**2, fmt="^", \
    #            color="crimson", ms=7, lw=2, mfc="w", zorder=1, label="Simulation truth")
    #ax4.plot(totNPE/es, rfunc(cerNPE, *popt)**2*cerNPE**2, color="black", lw=2, zorder=3, label="Best fit")

    #ax4.errorbar(totNPE/es, (sctSPE)**2, yerr=(sctSPEerr)**2, fmt="s", \
    #            color="crimson", ms=7, lw=2, mfc="w", zorder=1, label="Simulation truth")
    #ax4.plot(totNPE/es, sctNPE, color="black", lw=2, zorder=3, label="Best fResol
    #ax4.errorbar((sctNPE+cerNPE)/es, covArr, yerr=covErrArr, fmt="*", \
    #            color="crimson", ms=7, zorder=1, label="Simulation truth")


    #ax4.semilogy()    

    plt.tight_layout()

    plt.savefig("ResolDecomp.pdf")
    plt.show()


if __name__ == "__main__":
    main()



