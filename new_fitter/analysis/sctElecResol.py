import numpy as np
import matplotlib.pyplot as plt


def readTruth(filename):
    Etrue, totpe, totpeErr, sigma, sigmaErr, resol, resolErr = [], [], [], [], [], [], []
    with open("/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/data/electron/"+filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if float(data[0]) <= 150:
                Etrue.append(float(data[0]))
                totpe.append(float(data[1]))
                totpeErr.append(float(data[2]))
                sigma.append(float(data[3]))
                sigmaErr.append(float(data[4]))
                resol.append(float(data[5]))
                resolErr.append(float(data[6]))
    
    Etrue = np.array(Etrue)
    totpe = np.array(totpe)
    totpeErr = np.array(totpeErr)
    sigma = np.array(sigma)
    sigmaErr = np.array(sigmaErr)
    resol = np.array(resol)
    resolErr = np.array(resolErr)

    return Etrue, totpe, totpeErr, sigma, sigmaErr, resol, resolErr

def readMC(energy):
    path = "/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/electron/"
    filename = path + energy + "/electron.root"
    ff = ROOT.TFile(filename, "read")
    h1 = ff.Get("gamRatio")
    mean = h1.GetMean()
    return mean


import ROOT


if __name__ == "__main__" :

    Y1 = 3134.078 / 2.223
    Y2 = 3021.21 / 2.223

    Etrue1, totpe1, totpeErr1, sigma1, sigmaErr1, resol1, resolErr1 = readTruth("elecResol4.txt")
    Etrue2, totpe2, totpeErr2, sigma2, sigmaErr2, resol2, resolErr2 = readTruth("sctResol.txt")

    g1 = ROOT.TGraphErrors()
    g2 = ROOT.TGraphErrors()

    f1 = ROOT.TF1("f1", "sqrt([0]*[0]*x+[1]*[1]*x*x)", 0, 100000)
    f2 = ROOT.TF1("f2", "sqrt([0]*[0]*x+[1]*[1]*x*x)", 0, 100000)
    
    for i in range(len(totpe1)):
        g1.SetPoint(i, totpe1[i], sigma1[i])
        g1.SetPointError(i, 0, sigmaErr1[i])
    for i in range(len(totpe2)):
        g2.SetPoint(i, totpe2[i], sigma2[i])
        g2.SetPointError(i, 0, sigmaErr2[i])


    g1.Fit(f1, "RE")
    g2.Fit(f2, "RE")


    dx = np.arange(100, 100000, 100)
    dy1, dy2 = [], []
    pois = []
    for i in dx:
        dy1.append(f1.Eval(i)/Y1)
        dy2.append(f2.Eval(i)/Y2)
        pois.append(np.sqrt(i)/Y2)

    dx  = np.array(dx)
    dy1 = np.array(dy1)
    dy2 = np.array(dy2)
    pois = np.array(pois)

    Edep = [1, 5, 10, 20, 30, 40, 52]
    arr = ["1MeV", "5MeV", "10MeV", "20MeV", "30MeV", "40MeV", "52MeV"]
    mean = []

    for i in arr:
        mean.append(readMC(i))

    fig, ax1 = plt.subplots()
    ax1.plot(dx/Y2, (dy2-pois)/pois, "-.", lw=2, color="slategrey", label="Diff")
    ax1.plot(Edep, mean, "-", lw=2, color="green", label="eBrem Ratio")
    ax1.set_xlabel(r"$E^{vis}$ [MeV]", fontsize=15)
    ax1.set_ylabel("Ratio", fontsize=15)
    ax1.tick_params(axis='both', which='major', labelsize=14, labelcolor="black")

    ax1.set_xlim(0.09, 70)
    ax1.grid(True)

    fig, ax = plt.subplots()

    ax.errorbar(totpe1/Y1, sigma1/Y1, yerr=sigmaErr1/Y1, fmt="o", ms=6, color="blue", label="Total PE")
    ax.plot(dx/Y1, dy1, "-", lw=2, color="black")
    ax.errorbar(totpe2/Y2, sigma2/Y2, yerr=sigmaErr2/Y2, fmt="o", ms=6, color="orange", label="Scintillation PE")
    ax.plot(dx/Y2, dy2, "--", lw=2, color="black")

    ax.plot(dx/Y2, pois, "-.", lw=2, color="slategrey", label="Poisson")

    ax.legend(prop={"size":15})
    ax.set_xlabel(r"$E^{vis}$ [MeV]", fontsize=15)
    ax.set_ylabel(r"$\sigma$ [MeV]", fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=14, labelcolor="black")
    ax.semilogx()
    ax.set_xlim(0.09, 70)
    ax.grid(True)

    plt.tight_layout()

    plt.savefig("totpe+sctpe_resol.pdf")
    plt.show()




