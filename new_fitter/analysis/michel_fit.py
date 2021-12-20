import numpy as np
import matplotlib.pyplot as plt
import ROOT
from scipy import special
from matplotlib import gridspec

def loadMichelPE():
    #ff = ROOT.TFile("../michel_spectrum.root", "read")  # water
    print("preparing loading michel electron spectrum")
    ff = ROOT.TFile("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/michel/michel_spectrum_v2.root", "read")   # LS
    hh = ff.Get("hData")
    print("Total entries in histogram = %d" %hh.GetEntries())
    npe, entries, error = [], [], []
    for i in range(hh.GetNbinsX()):
        npe.append(hh.GetBinCenter(i+1))
        entries.append(hh.GetBinContent(i+1))
        error.append(hh.GetBinError(i+1))
    npe = np.array(npe)
    entries = np.array(entries)
    error = np.array(error)
    return npe, entries, error



def resFunc(a, b, c, N):
    return np.sqrt(a+b*N+c*N**2)/N


def add_subplot_axes(ax, rect, axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    #subax = fig.add_axes([x,y,width,height],facecolor=facecolor)  # matplotlib 2.0+
    subax = fig.add_axes([x,y,width,height])
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax



def michelFit():
    npe, entries, error = loadMichelPE()
    #c0, c1, Emid, sigma = 8.53697e+02, 8.23519e+02, 8.18234e+02, 5.57008e+01     # water
    #npe1 = np.arange(750, 900, 1)
    #c0, c1, Emid, sigma = 8.48007e+01, 9.08933e+01, 7.76149e+04, 8.37823e+02     #LS
    c0, c1, Emid, sigma = 4.95464e+02, 5.75788e+02, 7.76588e+04, 8.19608e+02
    npe1 = np.arange(76800, 78500, 1)
    fit = c0 - c1*special.erf((npe1 - Emid)/sigma)



    global_es = 3134.078/2.223
    q01, q11, q21 = 0,  9.79299e-01, 6.12574e-05
    
    #npe, spe = [], []
    micN, micS, micNerr, micSerr = 7.76588e+04, 8.19608e+02, 5.52924e+01, 2.13277e+02

    nn = np.arange(1000, micN+1000, 100)

    fig = plt.figure(figsize=(6, 4))
    ax1 = fig.add_subplot(111)
    rect = [0.35,0.35,0.6,0.6]
    ax0 = add_subplot_axes(ax1,rect)
    ax0.errorbar(npe/global_es, entries, yerr=error, fmt="o", ms=4, color="coral", label="Simulation")
    ax0.plot(npe1/global_es, fit, "-", color="blue", label="Fitting")
    ax0.set_xlim(54, 56)

    chi2, n = 4.2617, 6
    c0, c0err = 4.95464e+02, 3.61427e+01
    c1, c1err =  5.75788e+02, 1.16838e+02
    Esigma, Esigmaerr = 8.19608e+02/global_es, 2.13277e+02/global_es
    Emid, Emiderr = 7.76588e+04/global_es, 5.52924e+01/global_es 
    textstr = '\n'.join((
    r'$chi^2/NDF = %.1f / %d$' %(chi2, n),
    #r'$c_0 = %.2f \pm %.2f$' %(c0, c0err),
    #r'$c_1 = %.2f \pm %.2f$' %(c1, c1err),
    r'$\sigma_E = %.1f \pm %.1f$' %(Esigma, Esigmaerr),
    r'$E_{M} = %.1f \pm %.1f$' %(Emid, Emiderr)
    ))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax0.text(54.1, 20, textstr, fontsize=11, bbox=props)

    ax0.set_xlabel(r"$E_{vis}$[MeV]", fontsize=10)
    ax0.set_ylabel("Entries", fontsize=10)
    ax0.tick_params(axis='both', which='major', labelsize=10)
    #ax0.legend(prop={'size' : 14})
    ax1.plot(nn/global_es, resFunc(q01, q11, q21, nn), "--", color="black", lw=2, label="Model")
    ax1.errorbar(micN/global_es, micS/micN, yerr=np.sqrt(micSerr**2/micN**2 + micNerr**2*micS**2/micN**4), fmt="o", color="crimson", label=r"Michel $e^-$")
    ax1.set_ylabel(r"$\sigma_E/E_{vis}$", fontsize=14)
    ax1.set_xlabel(r"Electron $E_{vis}$[MeV]", fontsize=14)
    ax1.tick_params(axis='both', which='major', labelsize=13)
    ax1.legend(prop={'size' : 14}, ncol=2)

    plt.tight_layout()
    plt.savefig("michelElectron.pdf")
    plt.show()

    """
    fig = plt.figure(figsize=(12, 5))
    spec = gridspec.GridSpec(ncols=2, nrows=1 )

    ax0 = fig.add_subplot(spec[0])
    ax1 = fig.add_subplot(spec[1])


    ax0.errorbar(npe, entries, yerr=error, fmt="o", ms=4, color="coral", label="Simulation")
    ax0.plot(npe1, fit, "-", color="blue", label="Fitting")

    #plt.xlim(0, 1000)
    ax0.set_xlim(76000, 79200)
    ax0.set_xlabel("NPE", fontsize=14)
    ax0.set_ylabel("Entries", fontsize=14)
    ax0.tick_params(axis='both', which='major', labelsize=13)
    ax0.legend(prop={'size' : 14})


    ax1.plot(nn/global_es, resFunc(q01, q11, q21, nn), "--", color="coral", label="")
    ax1.errorbar(micN/global_es, micS/micN, yerr=np.sqrt(micSerr**2/micN**2 + micNerr**2*micS**2/micN**4), fmt="o")
    ax1.set_xlabel(r"$N_{tot}$", fontsize=14)
    ax1.set_ylabel(r"$E_{vis}$/MeV", fontsize=14)
    ax1.tick_params(axis='both', which='major', labelsize=13)
    ax1.legend(prop={'size' : 14})
    
    plt.savefig("michelFitLS.pdf")
    plt.show()
    """

def loadCerSigma():
    Etrue, npe, sigma = [], [], []
    with open("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/water/electron/cerPESigma_water.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Etrue.append(float(data[0]))
            npe.append(float(data[1]))
            sigma.append(float(data[2]))

    Etrue = np.array(Etrue)
    npe = np.array(npe)
    sigma = np.array(sigma)

    return Etrue, npe, sigma

def michelEdgeResol():
    Esim, npesim, sigmasim = loadCerSigma()
    fitEmid = 8.18234e2
    fitSigma = 5.57008e1
    fitSigmaErr = 4.94324e0
    fitEmid = 8.07571e+02
    fitSigma = 6.43617e+01   
    fitSigmaErr = 1.62765e+01
    
    plt.plot(npesim, sigmasim, "o", ms=4, label="simulation")
    plt.errorbar(fitEmid, fitSigma, yerr=fitSigmaErr, fmt="v", label="fitting")

    plt.xlabel("NPE")
    plt.ylabel(r"$\sigma_{NPE}$")
    plt.grid(True)
    plt.legend()
    plt.savefig("resolWater.pdf")
    plt.show()
    




if __name__ == "__main__":
    michelFit()
    #michelEdgeResol()

