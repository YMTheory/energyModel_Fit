import numpy as np
from matplotlib import gridspec
import matplotlib.pyplot as plt
import ROOT
import prmBetaLoader as bloader


def main():

    es = 3134.078/2.223

    # Load Livermore Truth :
    Edep, totpe, totpeerr, sigma, sigmaerr = [], [], [], [], []
    with open("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/gamma/singles/fitNPE.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if float(data[0]) < 0.2:
                continue
            if float(data[0]) > 7.0:
                break
            Edep.append(float(data[0]))
            totpe.append(float(data[1]))
            totpeerr.append(float(data[2]))
            sigma.append(float(data[3]))
            sigmaerr.append(float(data[4]))

    Edep = np.array(Edep)
    totpe = np.array(totpe)
    totpeerr = np.array(totpeerr)
    sigma = np.array(sigma)
    sigmaerr = np.array(sigmaerr)
    Evis = totpe/es
    Eviserr = totpeerr/es
    sigmaE = sigma/es
    sigmaEerr = sigmaerr/es



    # Load Livermore Calculation
    EdepL, totpeL, sigmaL = [], [], []
    with open("Livermore_check.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if float(data[0]) < 0.2:
                continue
            if float(data[0]) > 7.0:
                break
            EdepL.append(float(data[0]))
            totpeL.append(float(data[1]))
            sigmaL.append(float(data[2]))

    EdepL = np.array(EdepL)
    totpeL = np.array(totpeL)
    sigmaL = np.array(sigmaL)
    
    

    # Load Penelope Calculation
    EdepP, totpeP, sigmaP = [], [], []
    with open("Penelope_check.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if float(data[0]) < 0.2:
                continue
            if float(data[0]) > 7.0:
                break
            EdepP.append(float(data[0]))
            totpeP.append(float(data[1]))
            sigmaP.append(float(data[2]))
    
    EdepP = np.array(EdepP)
    totpeP = np.array(totpeP)
    sigmaP = np.array(sigmaP)

    res = sigmaE / Evis
    reserr = np.sqrt(sigmaEerr**2/Evis**2 + sigmaE**2*Eviserr**2/Evis**4)
    resP = sigmaP / totpeP
    resL = sigmaL / totpeL


    fig = plt.figure(figsize=(8, 6))
    spec = gridspec.GridSpec(ncols=1, nrows=2)

    ax0 = fig.add_subplot(spec[0])
    ax1 = fig.add_subplot(spec[1])

    ax0.errorbar(Edep, (totpeP-Evis)/(Evis), yerr=totpeP**2*Eviserr**2/Evis**4, fmt="o", color="seagreen", ms=7, label="Penelope")
    ax0.errorbar(Edep, (totpeL-Evis)/(Evis), yerr=totpeL**2*Eviserr**2/Evis**4, fmt="d", color="blue", ms=7, label="Livermore")
    ax0.set_title("(a) Nonlinearity", fontsize=18)
    ax0.legend(prop={"size":14})
    ax0.set_xlabel(r"$E_{dep}$/MeV", fontsize=18)
    ax0.set_ylabel("Relative bias", fontsize=18)
    ax0.grid(True)
    ax0.set_ylim(-0.005, 0.005)
    ax0.fill_between([0, 7.1], [-0.003, -0.003], [0.003, 0.003], color="royalblue", alpha=0.3)
    ax0.hlines(0, 0, 7.1, linestyle="--", lw=2, color="coral")
    ax0.text(1.0, 0.004, "0.3% uncertainties band", color="darkviolet", fontsize=15)
    ax0.tick_params(axis='both', which='major', labelsize=16)


    ax1.errorbar(Edep, (resP-res)/res, yerr=np.sqrt(resP**2*reserr**2/res**4), fmt="o", color="seagreen", ms=7, label="Penelope")
    ax1.errorbar(Edep, (resL-res)/res, yerr=np.sqrt(resL**2*reserr**2/res**4), fmt="d", color="blue", ms=7, label="Livermore")
    ax1.set_title("(b) Resolution", fontsize=18)
    ax1.set_xlabel(r"$E_{dep}$/MeV", fontsize=18)
    ax1.set_ylabel("Relative bias", fontsize=18)
    ax1.grid(True)
    ax1.set_ylim(-0.08, 0.08)
    ax1.fill_between([0, 7.1], [-0.05, -0.05], [0.05, 0.05], color="royalblue", alpha=0.3)
    ax1.hlines(0, 0, 7.1, linestyle="--", lw=2, color="coral")
    ax1.text(1.0, 0.060, "5% uncertainties band", color="darkviolet", fontsize=15)
    ax1.tick_params(axis='both', which='major', labelsize=16)

    plt.tight_layout()
    plt.savefig("modelDiffGam.pdf")
    plt.show()

if __name__ == "__main__" :
    main()

