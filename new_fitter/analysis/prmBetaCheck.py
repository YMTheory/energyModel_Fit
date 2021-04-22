import matplotlib.pyplot as plt

import ROOT
def prmBetaDist(name, Etrue, secBeta, secAntiBeta):
    sumN = []
    sumE = []
    for a in range(len(secBeta)):
        tmpE = 0
        tmpN = 0
        for i in secBeta[a]:
            tmpE += i
            tmpN += 1
        for j in secAntiBeta[a]:
            tmpE += j + 0.511*2
            tmpN += 1

        sumE.append(tmpE)
        sumN.append(tmpN)

    plt.figure(0)
    plt.hist(sumE, bins=100, histtype='step')
    plt.vlines(Etrue, 0, 100, linestyle="--")
    plt.xlabel("Sum of primary e-/e+ energy/MeV")
    plt.savefig(name+"_sumE.pdf")

    plt.figure(1)
    plt.hist(sumN, bins=30, range=(0, 30), histtype='step')
    plt.xlabel("# Primary e-/e+")
    plt.savefig(name+"_sumN.pdf")

    plt.figure(2)
    plt.hist2d(sumN, sumE, bins=[30, 100], range=([0, 30], [4.425, 4.430]))
    plt.xlabel("# Primary e-/e+")
    plt.ylabel("Sum of primary e-/e+ energy/MeV")
    plt.savefig(name+"hist2d.pdf")
    

    print("Pdf has been saved !")
