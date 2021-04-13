#!/usr/bin/env python
# coding=utf-8

import numpy as np
import uproot as up
import matplotlib.pyplot as plt
from ROOT import TFile, TGraph

def loadDiffGraph(filename):
    ff = TFile(filename, "read")
    g = ff.Get("PEDiff")
    return g


def Compare():
    g1 = loadDiffGraph("../GamPECheck_3.root")
    g2 = loadDiffGraph("../GamPECheck_7.root")
    g3 = loadDiffGraph("../GamPECheck_8.root")
    xarr, yarr1, yarr2, yarr3 =[], [], [], []
    for i in range(g1.GetN()):
        xarr.append(g1.GetPointX(i))
        yarr1.append(g1.GetPointY(i))
        yarr2.append(g2.GetPointY(i))
        yarr3.append(g3.GetPointY(i))

    #plt.plot(xarr, yarr1, "o-", color='royalblue',  label="w/o e+")
    #plt.plot(xarr, yarr2, "o-", color="darkviolet", label="w/ e+ KE")
    plt.plot(xarr, yarr3, "o-", color="hotpink", label="w/ e+ KE + annil gamma")

    plt.xlabel("Etrue/MeV")
    plt.ylabel("(calcPE-dataPE)/dataPE")

    plt.tight_layout()
    plt.legend()
    plt.grid(True)
    plt.savefig("gammaPEDiff_8.pdf")
    plt.show()

def singleSource(name, ratio):
    binning = [5, 8, 10, 20]
    plt.plot(binning, ratio, "o-", label=name)
    
    plt.legend()
    plt.xlabel("keV/bin")
    plt.ylabel("(calcpe-datape)/datape")
    plt.grid(True)



def main():
    #Compare()

    datape = 917.119
    ratio = [905.935/datape, 914.946/datape, 922.832/datape, 975.971/datape]
    singleSource("Cs137", ratio)

    datape = 1171.94
    ratio = [1158.8/datape, 1167.78/datape, 1175.68/datape,  1228.85/datape]
    singleSource("Mn54", ratio)

    datape = 11741.9
    ratio = [11712.1/datape, 11719.7/datape, 11726.2/datape, 11771.1/datape]
    singleSource("nFe56", ratio)

    plt.show()




if __name__ == '__main__':
    main()