import numpy as np
import matplotlib.pyplot as plt

def read(filename):
    Etrue, npe = [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Etrue.append(float(data[0]))
            npe.append(float(data[1]))
            
    Etrue = np.array(Etrue)
    npe = np.array(npe)

    return Etrue, npe

import ROOT

def main():
    E1, N1 = read("../data/electron/cerPE3.txt")
    E2, N2 = read("../data/electron/cerPESigma_water.txt")
    
    g1 = ROOT.TGraph()
    g2 = ROOT.TGraph()

    for i in range(len(E1)):
        g1.SetPoint(i, E1[i], N1[i])
    for i in range(len(E2)):
        g2.SetPoint(i, E2[i], N2[i])

    Ee = np.arange(1, 50, 1)
    ratio = []
    for i in Ee:
        ratio.append(g1.Eval(i)/g2.Eval(i))

    plt.plot(Ee, ratio)
    plt.show()


if __name__ == "__main__":
    main()
