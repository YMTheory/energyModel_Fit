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
        


def main():


    #E_sim, nonl_sim = loadGeant4()
    #nonl_int = integral(0.0065)

    #plt.plot(E_sim, nonl_sim, "-",  label="Geant4")
    #plt.plot(Etrue, nonl_int, "--", label="integral")

    #plt.show()

    """
    ff = ROOT.TFile("Quench_NumInt.root", "recreate")
    

    kB_arr = np.arange(51, 76, 1)
    for i in kB_arr:
        print("kB = ", i)
        kB = i/10000.
        value = integral(kB)
        name = "kB" + str(i)
        hist = ROOT.TH1D(name, "", len(value), 0, 15)
        for i in range( len(Etrue) ):
            hist.SetBinContent(i, value[i])
        hist.Write()
    #    graph = ROOT.TGraph()
    #    graph.SetName(name)
    #    for i in range( len(Etrue)):
    #        graph.SetPoint(i, Etrue[i], value[i])
    #    graph.Write()
    ff.Close()


    """
    
    fig, ax = plt.subplots()
    value0 = integral(0.0065)
    ax.plot(Etrue, value0, "-", lw=2, color="royalblue", label=r"ESTAR, $kB=6.5\times10^{-3}$g/cm$^2$/MeV")
    #for i, j in zip(Etrue, value0):
    #    print(i, j)
    E2, n2 = loadGeant4("kB65")
    ax.plot(E2, n2, "--", lw=2, color="royalblue", label=r"Geant4, $kB=6.5\times10^{-3}$g/cm$^2$/MeV")

    #value2 = integral(0.0063)
    #ax.plot(Etrue, value2, "-", lw=2, color="red", label=r"ESTAR, $kB=6.3\times10^{-3}$g/cm$^2$/MeV")
    E4, n4 = loadGeant4("kB63")
    ax.plot(E4, n4, "-.", lw=2, color="red", label=r"Geant4, $kB=6.3\times10^{-3}$g/cm$^2$/MeV")

    value1 = integral(0.0051)
    ax.plot(Etrue, value1, "-", lw=2, color="slategray", label=r"ESTAR, $kB=5.5\times10^{-3}$g/cm$^2$/MeV")
    E3, n3 = loadGeant4("kB51")
    ax.plot(E3, n3, "--", lw=2, color="slategray", label=r"Geant4, $kB=5.5\times10^{-3}$g/cm$^2$/MeV")

    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlabel("Electron $E_{dep}$ [MeV]", fontsize=15)
    ax.set_ylabel("$E_{sct}$/$E_{dep}$", fontsize=15)
    ax.legend(prop={'size': 13})
    ax.semilogx()
    ax.set_xlim(1e-3, 15)
    ax.set_ylim(0.3, 1.1)
    plt.tight_layout()
    plt.savefig("Num+G4_Birk.pdf")

    plt.show()
    
    

    



if __name__=="__main__":
    main()
