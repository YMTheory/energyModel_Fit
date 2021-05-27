import numpy as np
import matplotlib.pyplot as plt

def loadStopPow():
    filename = "../data/electron/ESTAR_GdLS.txt"
    Etrue, dEdx = [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Etrue.append(float(data[0]))
            dEdx.append(float(data[1]))

    return Etrue, dEdx

# load quenching from Geant4
def loadGeant4(name):
    ff = ROOT.TFile("../data/electron/Quench5.root", "read")
    hist2 = ff.Get(name)

    E1, n1, E2, n2, E3, n3 = [], [], [], [], [], []

    N = hist2.GetNbinsX()
    for i in range(N):
        E2.append(hist2.GetBinCenter(i))
        n2.append(hist2.GetBinContent(i))
    return E2, n2



from scipy import integrate
from scipy import interpolate

import ROOT
def loadBinCenter():
    ff = ROOT.TFile("../data/electron/Quench5.root", "read")
    hist = ff.Get("kB65")
    binCenter = []
    for i in range(hist.GetNbinsX()):
        binCenter.append(hist.GetBinCenter(i))
    return binCenter



#Etrue = loadBinCenter()
Etrue = np.arange(0.001, 15, 0.001)
#Etrue = [1, 2, 3]

def integral(kB):
    KE, dEdx = loadStopPow()
    fdEdx = interpolate.interp1d(KE, dEdx)
    nes = np.concatenate((np.array([0]), Etrue), axis=0)
    f = lambda x : 1./(1+kB*fdEdx(x))
    accum = 0.
    value = []
    for i in range(len(Etrue)):
        accum += integrate.quad(f, nes[i], nes[i+1])[0]
        value.append(accum/nes[i+1])
    return value
        


def main():
    #E_sim, nonl_sim = loadGeant4()
    #nonl_int = integral(0.0065)

    #plt.plot(E_sim, nonl_sim, "-",  label="Geant4")
    #plt.plot(Etrue, nonl_int, "--", label="integral")

    #plt.show()

    #ff = ROOT.TFile("Quench_NumInt.root", "recreate")

    #kB_arr = np.arange(50, 80, 1)
    #for i in kB_arr:
    #    print("kB = ", i)
    #    kB = i/1000.
    #    value = integral(kB)
    #    name = "kB" + str(i)
    #    hist = ROOT.TH1D(name, "", len(value), Etrue[0], Etrue[-1])
    #    for i in range( len(Etrue) ):
    #        hist.SetBinContent(i, value[i])
    #    hist.Write()

    #ff.Close()

    value0 = integral(0.0065)
    plt.plot(Etrue, value0, "-", color="blue", label=r"numerical integral, kB=0.0065g/cm$^2$/MeV")
    E2, n2 = loadGeant4("kB65")
    plt.plot(E2, n2, "-.", color="blue", label=r"Geant4, kB=0.0065g/cm$^2$/MeV")

    value1 = integral(0.0051)
    plt.plot(Etrue, value1, "-", color="orange", label=r"numerical integral, kB=0.0055g/cm$^2$/MeV")
    E3, n3 = loadGeant4("kB51")
    plt.plot(E3, n3, "-.", color="orange", label=r"Geant4, kB=0.0055g/cm$^2$/MeV")

    plt.grid(True)
    plt.xlabel(r"$E_{dep}$/MeV")
    plt.ylabel(r"$E_{vis}/{E_{dep}$")
    plt.legend()
    plt.semilogx()
    plt.ylim(0.5, 1)
    plt.show()

    



if __name__=="__main__":
    main()
