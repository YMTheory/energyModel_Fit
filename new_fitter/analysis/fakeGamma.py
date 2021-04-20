import numpy as np
import random
import matplotlib.pyplot as plt
import uproot as up

def loadPEFile(filename):
    Earr, PEarr = [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Earr.append(float(data[0]))
            PEarr.append(float(data[1]))
    return Earr, PEarr

def loadResFile(filename):
    Earr, resol = [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Earr.append(float(data[0]))
            resol.append(float(data[3]))
    return Earr, resol


sctE, sctPE = loadPEFile("../data/electron/scintillationPE.txt")
cerE, cerPE = loadPEFile("../data/electron/cerenkovPE.txt") 
resolE, resol = loadResFile("../data/electron/elecResol.txt")

def getPE(E):
    deltaE = 0.05
    lowbin = int(E/deltaE)
    higbin = lowbin+1
    bias = (E - lowbin*deltaE)/deltaE
    totpe = (1-bias)*( sctPE[lowbin] + cerPE[lowbin] ) + bias* (cerPE[higbin] + sctPE[higbin] )
    sctpe = (1-bias)*( sctPE[lowbin]  ) + bias* ( sctPE[higbin] )
    cerpe = (1-bias)*(  cerPE[lowbin] ) + bias* (cerPE[higbin]  )
    return cerpe, sctpe, totpe


def getResol(E):
    deltaE = 0.05
    lowbin = int(E/deltaE)
    higbin = lowbin+1
    bias = (E-lowbin*deltaE )/deltaE
    print(E, lowbin, higbin, bias, resol[lowbin])
    return (1-bias)*(resol[lowbin]) + bias* (resol[higbin])

def fakeGam():
    
    mom = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2]

    NPE = 0
    NPE_sigma = 0
    for i in mom:
        NPE += getPE(i)
        NPE_sigma += getResol(i)**2
        print(i, NPE, NPE_sigma)
    
    print(NPE, NPE_sigma)
    NPE_sigma = np.sqrt(NPE_sigma)

    sample_pe = []
    for i in range(5000):
        tmp_npe = random.gauss(NPE, NPE_sigma)
        sample_pe.append(tmp_npe)

    plt.hist(sample_pe, bins=100, range=(1200, 1500), histtype="step", label="calc")

    sim = []
    for i in range(10):
        filename = "/junofs/users/miaoyu/energy_model/production/electron/user-"+str(i)+".root"
        arr = up.open(filename)["evt"]["totalPE"].array()
        sim.extend(arr)

    sample_pe = np.array(sample_pe)
    sim = np.array(sim)
    
    print("Simulation Event # %d" %len(sim))
    print(sim.mean(), sample_pe.mean())
    plt.hist(sim, bins=100, range=(1200,1500), histtype="step", label="sim")

    plt.legend()
    plt.xlabel("total p.e.")

    plt.show()

secBetaArr = []
def loadPrmBeta(filename):
    with open(filename) as f:
        for lines in f.readlines():
            tmpE = 0
            oneEvt = []
            line = lines.strip("\n")
            data = line.split(" ")
            counta = 0
            countb = 0
            for i in data:
                #if "a" in i:
                #    counta+=1
                #    tmp = list(i)
                #    tmp.pop()
                #    j = ''.join(tmp)
                #    hh2.SetBinContent(evtid+1, counta, float(j))
                if "b" in i:
                    countb+=1
                    tmp = list(i)
                    tmp.pop()
                    j = ''.join(tmp)
                    oneEvt.append(float(j))
                    tmpE += float(j)
            secBetaArr.append(tmpE)
            #secBetaArr.append(oneEvt)
    return secBetaArr, tmpE

secBetaArr_old = []
def loadPrmBeta_old(filename):
    with open(filename) as f:
        for lines in f.readlines():
            tmpE = 0
            oneEvt = []
            line = lines.strip("\n")
            data = line.split(" ")
            for j in data:
                if j == "":
                    continue
                oneEvt.append(float(j))
                tmpE += float(j)
            #secBetaArr_old.append(oneEvt)
            secBetaArr_old.append(tmpE)
    return secBetaArr_old



def Cs137_source():
    secBetaArr = loadPrmBeta("../data/gamma/Cs137-nooptical.txt")
    data_arr = up.open("/junofs/users/miaoyu/energy_model/production/gamma/Cs137/user-sim-7000.root")["evt"]["totalPE"].array()
    calc_arr = []
    prmE_arr = []
    for i in range(500):
        tmpE = 0
        NPE, NPE_sigma = 0, 0
        for j in secBetaArr[i]:
            tmpE += j
            NPE += getPE(j)
            NPE_sigma += getResol(j)**2
        NPE_sigma = np.sqrt(NPE_sigma)
        calc_arr.append(random.gauss(NPE, NPE_sigma))
        prmE_arr.append(tmpE)

    data_arr = np.array(data_arr)
    calc_arr = np.array(calc_arr)

    print("Simulation NPE mean = %d" %data_arr.mean())
    print("Calculation NPE mean = %d" %calc_arr.mean())
    plt.plot(data_arr, "o", ms=0.5, color="blue", label="sim")
    plt.plot(calc_arr, "o", ms=0.5, color="red", label="calc")

    plt.legend()

    plt.show()


def Draw():
    Earr, sctarr, cerarr, pearr= [], [], [], []
    for i in range(100):
        E = 16./100 *i
        cerpe, sctpe, totpe  = getPE(E)
        Earr.append(E)
        sctarr.append(sctpe)
        cerarr.append(cerpe)
        pearr.append(totpe)
    plt.plot(Earr, cerarr, "-", label="cerenkov")
    plt.plot(Earr, sctarr, "--",  label="scintillation")
    plt.plot(Earr, pearr, "-.",  label="total")
    plt.legend()
    plt.grid(True)
    plt.xlabel("Etrue/MeV")
    plt.ylabel("# p.e.")
    plt.semilogy()
    plt.show()



def main():
    #Cs137_source()
    #fakeGam()
    #Draw()
    E1 = loadPrmBeta("/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/data/gamma/Cs137-nooptical.txt")
    E2 = loadPrmBeta_old("/junofs/users/miaoyu/energy_model/energyModel_Fit/Miao/data/naked_gamma/Cs137_all.txt")
    E3 = loadPrmBeta("/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/data/gamma/Cs137_newversion.txt")

    plt.hist(E1, bins=100, histtype="step", label="J20")
    #plt.hist(E2, bins=100, histtype="step", label="J19")
    plt.hist(E3, bins=100, histtype="step", label="J20 new")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
