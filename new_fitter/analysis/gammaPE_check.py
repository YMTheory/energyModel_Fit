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
    g2 = loadDiffGraph("../GamPECheck_4.root")
    g3 = loadDiffGraph("../GamPECheck_5.root")
    xarr, yarr1, yarr2, yarr3 =[], [], [], []
    for i in range(g1.GetN()):
        xarr.append(g1.GetPointX(i))
        yarr1.append(g1.GetPointY(i))
        yarr2.append(g2.GetPointY(i))
        yarr3.append(g3.GetPointY(i))

    plt.plot(xarr, yarr1, "o-", color='royalblue',  label="800 bins in [0, 8]MeV")
    plt.plot(xarr, yarr2, "o-", color="darkviolet", label="600 bins in [0, 12]MeV")
    plt.plot(xarr, yarr3, "o-", color="hotpink", label="1600 bins in [0, 8]MeV")

    plt.xlabel("Etrue/MeV")
    plt.ylabel("(calcPE-dataPE)/dataPE")

    plt.tight_layout()
    plt.legend()
    plt.grid(True)
    plt.savefig("gammaPEDiff.pdf")
    plt.show()

def main():
    Compare()


if __name__ == '__main__':
    main()