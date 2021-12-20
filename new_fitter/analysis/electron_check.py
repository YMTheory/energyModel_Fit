import numpy as np
import matplotlib.pyplot as plt
import elecLoader as el
import ROOT

if __name__ == "__main__":
    
    arr = [146, 180, 210, 246, 300, 346, 400, 446, 500, 546, 600, 646, 700, 746, 846, 856, 866, 876, 886, 896, 916, 917, 918, 919] 
    totpe, Earr, resol, resolerr = [], [], [], []
    nn = 1
    filename = "../data/electron/elecResol1.txt"
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if nn in arr:
                Earr.append(float(data[0]))
                totpe.append(float(data[1]))
                resol.append(float(data[5]))
                resolerr.append(float(data[6]))
            nn += 1
    Earr  = np.array(Earr)
    resol = np.array(resol)
    resolerr = np.array(resolerr)


    g1 = ROOT.TGraphErrors()
    for i in range(len(Earr)):
        g1.SetPoint(i, Earr[i], resol[i])
        g1.SetPointError(i, 0, resolerr[i])

    A = 3134.078 / 2.223
    besta, bestb = 0.988, 7.89e-3

    f1 = ROOT.TF1("f1", "sqrt([0]*[0]/x+[1]*[1])", 0, 60, 2)

    g1.Fit(f1, "RE")
    print(f1.GetParameter(0), f1.GetParameter(1))


    dx = np.arange(1, 50, 1)
    dy = []; fit = []
    for i in dx:
        dy.append(f1.Eval(i))
        fit.append(np.sqrt(besta**2/A/i + bestb**2))

    plt.plot(Earr, resol, "o", label="simulation")
    plt.plot(dx, dy, "--", label="")
    plt.plot(dx, fit, "-", color="red")

    plt.show()
