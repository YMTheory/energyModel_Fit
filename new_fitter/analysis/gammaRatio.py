import uproot as up
import numpy as np
import matplotlib.pyplot as plt
import ROOT

def readTruth(energy):
    path = "/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/electron/"
    filename = path + energy + "/electron.root"
    ff = ROOT.TFile(filename, "read")
    h1 = ff.Get("gamRatio")
    mean = h1.GetMean()
    return mean



if __name__ == "__main__":
    Edep = [1, 5, 10, 20, 30, 40, 52]
    arr = ["1MeV", "5MeV", "10MeV", "20MeV", "30MeV", "40MeV", "52MeV"]
    mean = []

    for i in arr:
        mean.append(readTruth(i))


    fig, ax = plt.subplots()
    ax.plot(Edep, mean, "o-", ms=6, lw=2, color="crimson")
    ax.set_xlabel(r"$E^{dep}$ [MeV]", fontsize=15)
    ax.set_ylabel(r"$E^{eBrem}/E^{dep}$", fontsize=15)
    ax.tick_params(axis="both", which="major", labelsize=14)
    ax.grid(True)

    plt.tight_layout()
    plt.show()
    
