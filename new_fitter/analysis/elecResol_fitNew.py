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


    f1 = ROOT.TF1("f1", "sqrt([0]*[0]*x + +[1]*[1]*x*x) ", 0, 100000)
    f1.SetParameter(0, 0.9386)
    #f1.SetParameter(1, 1e-3)
    f1.SetParameter(1, 0)
    f1.SetParameter(2, 0.1034)
    f1.SetParameter(3, 1.439)
    f2 = ROOT.TF1("f2", "sqrt([0]*[0]*x + [1]*[1]*x*x + [2]*[2]*pow(x, [3]))", 0, 100000)
    f2.SetParameter(0, 0.840)
    #f2.SetParameter(1, 1e-3)
    f2.SetParameter(1, 0)
    f2.SetParameter(2, 0.286)
    f2.SetParameter(3, 1.237)
    f3 = ROOT.TF1("f3", "sqrt([0]*[0]*x + [1]*[1]*pow(x, [2]))", 0, 100000)


    g1 = ROOT.TGraphErrors()
    g2 = ROOT.TGraphErrors()
    for i in range(len(totpe)):
        g1.SetPoint(i, totpe[i], sigma[i])
        g1.SetPointError(i, 0, sigmaErr[i])
        if totpe[i] / Y < 10:
            g2.SetPoint(i, totpe[i], sigma[i])
            g2.SetPointError(i, 0, sigmaErr[i])


    #g1.Fit(f1, "RE")
    g1.Fit(f3, "RE")
    g2.Fit(f3, "RE")
    #g1.Fit(f3, "RE")

    draw = np.arange(100, 92000, 100)
    y1, y2, y3 = [], [], []
    for i in draw:
        # sigma NPE
        y1.append(f1.Eval(i))
        y2.append(f2.Eval(i))
        y3.append(f3.Eval(i))
    y1 = np.array(y1)
    y2 = np.array(y2)
    y3 = np.array(y3)


    ext1, ext2, ext3 = [], [], []
    cal1, cal2, cal3 = [], [], []

    f2.SetParameters(0.939, 0, 0.103, 1.439)
    for i in draw:
        ext1.append(f2.Eval(i))
    for i in totpe:
        cal1.append(f2.Eval(i)/i)

    f2.SetParameters(0.840, 0, 0.286, 1.234)
    for i in draw:
        ext2.append(f2.Eval(i))
    for i in totpe:
        cal2.append(f2.Eval(i)/i)

    f2.SetParameters(0.855, 0, 0.262, 1.253)
    for i in draw:
        ext3.append(f2.Eval(i))
    for i in totpe:
        cal3.append(f2.Eval(i)/i)




    diff1, diff2, diff3 = [], [], []
    for i in range(len(totpe)):
        tru = sigma[i] / totpe[i]
        #cal1 = f1.Eval(totpe[i]) / totpe[i]
        #cal2 = f2.Eval(totpe[i]) / totpe[i]
        #cal3 = f3.Eval(totpe[i]) / totpe[i]
        diff1.append( (cal1[i] - tru) / (tru))
        diff2.append( (cal2[i] - tru) / (tru))
        diff3.append( (cal3[i] - tru) / (tru))

    diff1 = np.array(diff1)
    diff2 = np.array(diff2)
    diff3 = np.array(diff3)





    fig, ax = plt.subplots()

    ax.plot(totpe/Y, diff1, "o-", ms=6, color="crimson",   label="only gamma")
    ax.plot(totpe/Y, diff2, "o-", ms=6, color="slategray", label="gamma + B12 + Michel")
    ax.plot(totpe/Y, diff3, "o-", ms=6, color="royalblue", label="gamma + B12 + Michel (edge)")

    #ax.plot(totpe/Y, sigma/totpe, "o", ms=6, color="black", mfc="w", label="Simulation truth", zorder=4)
    #ax.plot(draw/Y, ext1/draw, "-", lw=2, color="crimson", label="only gamma")
    #ax.plot(draw/Y, ext2/draw, "--", lw=2, color="slategray", label="gamma + B12 + Michel")
    #ax.plot(draw/Y, ext3/draw, "-.", lw=2, color="royalblue", label="gamma + B12 + Michel (edge)")
    #ax.plot(draw/Y, y1/draw, "-", lw=2, color="crimson", label="Form 1")
    #ax.plot(draw/Y, y2/draw, "--", lw=2, color="slategray", label="Form 2")
    #ax.plot(draw/Y, y3/draw, "-.", lw=2, color="royalblue", label="Form 3")
    ax.legend(prop={"size":15})
    ax.set_xlabel(r"$E^{vis}$ [MeV]", fontsize=15)
    #ax.set_ylabel(r"$\sigma/E^{vis}$", fontsize=15)
    ax.set_ylabel("Diff", fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=14, labelcolor="black")
    ax.grid(True)
    #ax.set_ylim(-0.15, 0.1)
    ax.set_ylim(-0.2, 0.2)
    #ax.semilogx()
    #ax.semilogy()

    """
    
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


    """
    plt.tight_layout()
    plt.show()

















