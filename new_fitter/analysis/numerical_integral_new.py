import numpy as np
import matplotlib.pyplot as plt

def loadStopPow():
    filename = "JUNOLS_StoppingPower.txt"
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
        E2.append(hist2.GetBinCenter(i+1))
        n2.append(hist2.GetBinContent(i+1))
    return E2, n2



from scipy import integrate
from scipy import interpolate

import ROOT
def loadBinCenter():
    ff = ROOT.TFile("../data/electron/Quench5.root", "read")
    hist = ff.Get("kB65")
    binCenter = []
    for i in range(hist.GetNbinsX()):
        binCenter.append(hist.GetBinCenter(i+1))
    return binCenter



def draw_stoppow(E, sp):
    plt.plot(E, sp, "-")
    plt.semilogy()
    plt.xlabel("electron kE/MeV")
    plt.ylabel("stopping power")



def write_data(arr, filename):
    with open(filename, "w") as f:
        for i in arr:
            f.write(i)
            f.write(" ")


#Etrue = loadBinCenter()
Etrue = np.arange(0.001, 15, 0.001)
#Etrue = [1, 2, 3]
#Etrue, t#mp = loadGeant4("kB65")
print(Etrue)

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
        

import elecLoader as el

def main():

    es = 3134.078 / 2.223
    
    fig, ax = plt.subplots()
    value0 = integral(0.0065)
    ax.plot(Etrue, value0, "-", lw=2, color="royalblue", label=r"ESTAR, $kB=6.50\times10^{-3}$g/cm$^2$/MeV")
    #for i, j in zip(Etrue, value0):
    #    print(i, j)
    E2, n2 = loadGeant4("kB65")
    ax.plot(E2, n2, "--", lw=2, color="royalblue", label=r"Geant4, $kB=6.50\times10^{-3}$g/cm$^2$/MeV")

    #value2 = integral(0.0063)
    #ax.plot(Etrue, value2, "-", lw=2, color="red", label=r"ESTAR, $kB=6.3\times10^{-3}$g/cm$^2$/MeV")
    #E4, n4 = loadGeant4("kB63")
    #ax.plot(E4, n4, "-.", lw=2, color="red", label=r"Geant4, $kB=6.3\times10^{-3}$g/cm$^2$/MeV")

    value1 = integral(0.0051)
    ax.plot(Etrue, value1, "-", lw=2, color="slategray", label=r"ESTAR, $kB=5.50\times10^{-3}$g/cm$^2$/MeV")
    E3, n3 = loadGeant4("kB51")
    ax.plot(E3, n3, "--", lw=2, color="slategray", label=r"Geant4, $kB=5.50\times10^{-3}$g/cm$^2$/MeV")

    As = 1416.1
    bestFit = 5.77e-3
    fitError = 0.18e-4

    n4, n4min, n4max = [], [], []
    for i in Etrue:
        n4.append(el.getFitNsct(i, bestFit, As, "Sim"))
        n4min.append(el.getFitNsct(i, bestFit-fitError, As, "Sim"))
        n4max.append(el.getFitNsct(i, bestFit+fitError, As, "Sim"))
    
    #n4 = integral(bestFit)
    #n4min = integral(bestFit-fitError)
    #n4max = integral(bestFit+fitError)

    n4 = np.array(n4)
    n4min = np.array(n4min)
    n4max = np.array(n4max)
    delmin = n4 - n4min
    delmax = n4max - n4
    n4min = n4 - delmin*50
    n4max = n4 + delmax*50

    ax.plot(Etrue, n4/As/Etrue, "--", lw=2, color="crimson", label=r"Geant4, $kB=5.77\times10^{-3}$g/cm$^2$/MeV")
    ax.fill_between(Etrue, n4min/As/Etrue, n4max/As/Etrue, color="crimson", alpha=0.3)
    #ax.plot(Etrue, n4, "--", lw=2, color="crimson", label=r"ESTAR, $kB=5.45\times10^{-3}$g/cm$^2$/MeV")
    #ax.fill_between(Etrue, n4min, n4max, color="crimson", alpha=0.3)



    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlabel("Electron $E$ [MeV]", fontsize=15)
    ax.set_ylabel("$E_{s}$/$E$", fontsize=15)
    ax.legend(prop={'size': 13})
    ax.semilogx()
    ax.set_xlim(1e-3, 15)
    ax.set_ylim(0.3, 1.1)
    plt.tight_layout()
    plt.savefig("Num+G4_Birk.pdf")

    plt.show()
    
    

    



if __name__=="__main__":
    main()
