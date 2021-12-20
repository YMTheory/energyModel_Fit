import numpy as np
import matplotlib.pyplot as plt
import ROOT
from matplotlib import gridspec

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




if __name__ == "__main__":

    Y = 3134.078 / 2.223

    Etrue, totpe, totpeErr, sigma, sigmaErr, resol, resolErr = readTruth("elecResol4.txt")
    Etrue1, totpe1, totpeErr1, sigma1, sigmaErr1, resol1, resolErr1 = readTruth("elecResol1.txt")


    f1 = ROOT.TF1("f1", "sqrt(([1]*[1]/(3134.078/2.223)*x + [2]*[2]*pow(x, 1.2)))/x", 0, 66)
    f2 = ROOT.TF1("f2", "([0]*[0]*x + [1]*[1]*pow(x, [2]))", 0, 100000)
    f2.SetParameters(0.90, 0.11579, 1.42189)
    #f2 = ROOT.TF1("f2", "([0]*[0] + [1]*[1]/(3134.078/2.223)*x + [2]*[2]*pow(x, 1.18))", 0, 66);
    #f2 = ROOT.TF1("f2", "([0]*[0] + [1]*[1]/(3134.078/2.223)*x + [2]*[2]*pow(x, [3]))", 0, 66);
    #f2.SetParameter(2, 1.5)
    f3 = ROOT.TF1("f3", "[0]*pow(x, 0.8) + [1]*sqrt(x)", 0, 66)

    g1 = ROOT.TGraphErrors()
    g2 = ROOT.TGraphErrors()
    g3 = ROOT.TGraphErrors()

    for i in range(len(Etrue)):
        g1.SetPoint(i, totpe[i]/Y, resol[i])
        g1.SetPointError(i, 0, resolErr[i])

        g2.SetPoint(i, totpe[i], sigma[i]**2)
        g2.SetPointError(i, 0, 2*sigma[i]*sigmaErr[i])

        g3.SetPoint(i, totpe[i]/Y, sigma[i]/Y)
        g3.SetPointError(i, 0, sigmaErr[i]/Y)


    g1.Fit(f1, "RE")
    #g2.Fit(f2, "RE")
    g3.Fit(f3, "RE")


    dx = np.arange(0.1, 62, 0.1)
    dx2 = np.arange(100, 20000, 100)
    dy1, dy2, dy3 = [], [], []
    for i in dx:
        dy1.append(f1.Eval(i))
        dy3.append(f3.Eval(i))

    for i in dx2:
        dy2.append(f2.Eval(i))



    sigma2Diff, sigma2DiffErr = [], []
    sigmaDiff, sigmaDiffErr = [], []
    resolDiff, resolDiffErr = [], []
    for i in range(len(totpe)) :
        resolDiff.append((f1.Eval(totpe[i]/Y) - resol[i])/resol[i])
        resolDiffErr.append(f1.Eval(totpe[i]/Y)*resolErr[i]/resol[i]**2)

        sigma2Diff.append((f2.Eval(totpe[i]/Y)-sigma[i]**2/Y**2)/sigma[i]**2*Y**2)
        sigma2DiffErr.append(f2.Eval(totpe[i]/Y) * 2 /sigma[i]**3 * sigmaErr[i] * Y**2)

        sigmaDiff.append((f3.Eval(totpe[i]/Y)-sigma[i]/Y)/sigma[i]*Y)
        sigmaDiffErr.append(f3.Eval(totpe[i]/Y)*sigmaErr[i]/sigma[i]**2*Y)


    #fig1 = plt.figure(figsize=(10, 6))
    #spec = gridspec.GridSpec(ncols=2, nrows=1)

    #ax4 = fig1.add_subplot(spec[0])
    #ax5 = fig1.add_subplot(spec[1])

    #ax4.errorbar(Esub1, sigmasub1, yerr=sigmaerrsub1, fmt="o", color="peru")
    #ax4.plot(dx3, dy3, "-", color="blue")
    #ax5.errorbar(Esub2, sigmasub2, yerr=sigmaerrsub2, fmt="o", color="crimson")
    #ax5.plot(dx4, dy4, "-", color="blue")


    
    fig = plt.figure(figsize=(10, 6))
    spec = gridspec.GridSpec(ncols=2, nrows=2, height_ratios=[1, 2])

    ax2 = fig.add_subplot(spec[0])
    ax  = fig.add_subplot(spec[1])
    ax0 = fig.add_subplot(spec[2])
    ax1 = fig.add_subplot(spec[3])


    ax2.errorbar(totpe/Y, sigma2Diff, yerr=sigma2DiffErr, fmt="o", color="peru")
    ax2.set_xlabel(r"$E^{vis}$ [MeV]", fontsize=15)
    ax2.set_ylabel("Bias", fontsize=15)
    ax2.tick_params(axis='both', which='major', labelsize=14, labelcolor="black")
    ax2.grid(True)

    

    ax.errorbar(totpe/Y, resolDiff, yerr=resolDiffErr, fmt="o", color="crimson")
    ax.set_xlabel(r"$E^{vis}$ [MeV]", fontsize=15)
    ax.set_ylabel("Bias", fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=14, labelcolor="black")
    ax.grid(True)


    ax0.errorbar(totpe1, sigma1**2, yerr=2*sigma1*sigmaErr1, fmt="o", color="peru", label="Simulation Truth", zorder=1)
    #ax0.errorbar(totpe/Y, sigma/Y, yerr=sigmaErr/Y, fmt="o", color="peru", label="Simulation Truth")
    ax0.plot(dx2, dy2, "-", color="blue", lw=2, label="Parameterisation", zorder=2)
    ax0.legend(prop={"size":15})
    ax0.set_xlabel(r"$E^{vis}$ [MeV]", fontsize=15)
    ax0.set_ylabel(r"$\sigma^2$", fontsize=15)
    ax0.tick_params(axis='both', which='major', labelsize=14, labelcolor="black")
    ax0.grid(True)

    ax1.errorbar(totpe/Y, resol, yerr=resolErr, fmt="o", ms=5, color="crimson", label="Simulation Truth")
    ax1.plot(dx, dy1, "-",  color="blue", lw=2, label="Parameterisation")
    ax1.legend(prop={"size":15})
    ax1.set_xlabel(r"$E^{vis}$ [MeV]", fontsize=15)
    ax1.set_ylabel(r"$\sigma/E^{vis}$", fontsize=15)
    ax1.tick_params(axis='both', which='major', labelsize=14, labelcolor="black")
    ax1.grid(True)


    plt.tight_layout()
    plt.show()

















