import ROOT
import numpy  as np
import matplotlib.pyplot as plt

def loadQuench(filename):
    ff = ROOT.TFile(filename, "read")
    hh = ff.Get("kB65")
    Etrue, nonl = [], []
    for i in range(hh.GetNbinsX()):
        Etrue.append(hh.GetBinCenter(i))
        nonl.append(hh.GetBinContent(i))

    return Etrue, nonl

import elecLoader as eloader


def main():
    Etrue1, nonl1 = loadQuench("../data/electron/Quench5.root")
    
    sc = 1500
    nonl2 = []
    for i in Etrue1:
        nonl2.append(eloader.getNPE(i) / sc / i)



    plt.plot(Etrue1, nonl1, "o-", ms=2)
    plt.plot(Etrue1, nonl2, "o-", ms=2)
    plt.semilogx()
    plt.xlim(0, 15)
    plt.show()



if __name__ == "__main__":
    main()
