import numpy as np
import uproot as up
import elecLoader as eloader
import matplotlib.pyplot as plt
import random

def readSim(filename):
    totpe = up.open(filename)["evt"]["totalPE"].array()
    edep  = up.open(filename)["evt"]["edep"].array()

    return totpe, edep


def main():
    totpe, edep = [], []
    for i in range(1000, 1300, 1):
        filename = "/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/B12/user-" +str(i) + ".root"
        print(filename)
        tmp_totpe, tmp_edep = readSim(filename)
        totpe.extend(tmp_totpe)
    for i in range(2000, 2100, 1):
        filename = "/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/B12/user-" +str(i) + ".root"
        print(filename)
        tmp_totpe, tmp_edep = readSim(filename)
        edep.extend(tmp_edep)

    pe_pred, pe_pred1 = [], []
    for i in edep:
        sample = eloader.getNPE(i)
        sample1 = 1.126*eloader.getSctNPE(i) - 0.812*eloader.getSPE(i)
        #sample = random.gauss(eloader.getNPE(i), eloader.getSPE(i))
        #sample1 = random.gauss((1.126* eloader.getSctNPE(i) - 0.812*eloader.getCerNPE(i)), eloader.getSPE(i))
        pe_pred.append(sample)
        pe_pred1.append(sample1)

    plt.hist(totpe,    bins=300, range=(0, 20000), density=True, histtype="step", label="simulation")
    plt.hist(pe_pred,  bins=300, range=(0, 20000), density=True, histtype="step", label="kA=1, kC=1 (nominal)")
    plt.hist(pe_pred1, bins=300, range=(0, 20000), density=True, histtype="step", label="kA=1.126, kC=-0.812")


    plt.legend()
    plt.xlabel("# P.E.")
    plt.tight_layout()
    #plt.savefig("B12_totpe.pdf")
    plt.show()


if __name__ == "__main__":
    main()
    

