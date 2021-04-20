import numpy as np
import matplotlib.pyplot as plt
import random
import uproot3

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

sctE, sctPE = loadPEFile("../data/electron/scintillationPE1.txt")
cerE, cerPE = loadPEFile("../data/electron/cerenkovPE1.txt") 
resolE, resol = loadResFile("../data/electron/elecResol.txt")

from ROOT import TGraph
g1 = TGraph()
g2 = TGraph()
for i in range(len(sctE)):
    g1.SetPoint(i, sctE[i], sctPE[i]+cerPE[i])
for i in range(len(resolE)):
    g2.SetPoint(i, resolE[i], resol[i])

import uproot as up
def loadSimTruth(filename):
    evt = up.open(filename)["evt"]
    totpe = evt["totalPE"].array()
    prmBetaPE = evt["prmBetaPE"].array()
    prmBetaE = evt["prmBetaEnergy"].array()
    return totpe, prmBetaPE, prmBetaE


def loadFitterCurve():
    # curves in fitter:
    Eetrue, sctpe, cerpe, simtotpe = [], [], [], []
    with open("/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/data/electron/scintillationPE.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Eetrue.append(float(data[0]))
            sctpe.append(float(data[1]))


    with open("/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/data/electron/cerenkovPE.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            cerpe.append(float(data[1]))

    for i in range(len(sctpe)):
        simtotpe.append(sctpe[i]+cerpe[i])

    plt.plot(Eetrue, simtotpe, "--", label="Fitter Curve")

from ROOT import TH1D, TF1

def gammaSource(name):
    path = "/junofs/users/miaoyu/energy_model/production/gamma/"
    totpe, prmBetaPE, prmBetaE = [], [], []
    for i in range(900, 910, 1):
        filename = path + name + "/user-sim-"+str(i) + ".root"
        print("------> Reading " + filename)
        tmptotpe, tmpprmBetaPE, tmpprmBetaE = loadSimTruth(filename)
        totpe.extend(tmptotpe)
        prmBetaPE.extend(tmpprmBetaPE)
        prmBetaE.extend(tmpprmBetaE)

    """ prmBeta-single-chekc """
    
    simPrmBetaE, simPrmBetaPE, calcPrmBetaPE = [], [], []
    binningSim, binningCalc = [[] for i in range(10)], [[] for i in range(10)]
    meanSim, meanCalc = [], []
    for i in range(len(totpe)):
        if len(prmBetaPE[i]) != len(prmBetaE[i]):
            continue
        for j in range(len(prmBetaPE[i])):
            simPrmBetaE.append(prmBetaE[i][j])
            simPrmBetaPE.append(prmBetaPE[i][j])
            calcPrmBetaPE.append(random.gauss(g1.Eval(prmBetaE[i][j], 0, "S"), g2.Eval(prmBetaE[i][j])))

    diff = []
    for i, j in zip(simPrmBetaPE, calcPrmBetaPE) :
        diff.append(j-i)
    diff = np.array(diff)
    plt.hist(diff, bins=100, histtype="step", label=name+" mean %.3f"%diff.mean() )



    ''' binning '''
    '''
    for a, b, e in zip(simPrmBetaPE, calcPrmBetaPE, simPrmBetaE):
        if e >= 0.01:
            continue
        binid = int(e/0.001)
        binningSim[binid].append(a)
        binningCalc[binid].append(b)

    meanFitter = []
    for i in range(10):
        tmpSim = np.array(binningSim[i])
        tmpCalc = np.array(binningCalc[i])
        meanSim.append(tmpSim.mean())
        meanCalc.append(tmpCalc.mean())
        #meanSim.append(np.std(tmpSim))
        #meanCalc.append(np.std(tmpCalc))

    Eeee = np.arange(0, 0.01, 0.001)
    for i in Eeee:
        meanFitter.append(g1.Eval(i+0.0005, 0, "S"))

    meanSim = np.array(meanSim)
    meanCalc = np.array(meanCalc)
    #plt.plot(Eeee+0.005, meanSim,  "o-", ms=0.3, label="Sim")
    #plt.plot(Eeee+0.005, meanCalc, "s-", ms=0.3, label="Calc")
    #plt.plot(Eeee, meanFitter, "--", label="Fitter Curve")
    #plt.plot(resolE, resol, "-")
    #plt.xlim(0, 0.5)
    #plt.ylim(0, 40)
    plt.plot(Eeee+0.0005, (meanCalc-meanSim)/meanSim, "-", label=name+"relative difference")
    plt.hlines(0, Eeee[0], Eeee[-1], colors="red", linestyle='--')
    plt.xlabel("Etrue/MeV")
    plt.ylabel("# P.E.")
    plt.grid(True)
    

    #plt.plot(simPrmBetaE, simPrmBetaPE, "o",  ms=0.5, alpha=0.5, label="Simulation")
    #plt.plot(simPrmBetaE, calcPrmBetaPE, "o", ms=0.5, alpha=0.5, label="Calculation")
    #plt.xlabel("Etrue/MeV")
    #plt.ylabel("# P.E.")
    '''

    #sumPreBetaPE = []
    #for i in range(len(totpe)):
    #    tmpsum = 0
    #    for j in range(len(prmBetaPE[i])):
    #        #tmpsum += g1.Eval(prmBetaE[i][j])
    #        tmpsum += random.gauss(g1.Eval(prmBetaE[i][j]), g2.Eval(prmBetaE[i][j]))    # consider resolution
    #    sumPreBetaPE.append(tmpsum)
    #
    #totpe = np.array(totpe)
    #sumPreBetaPE = np.array(sumPreBetaPE)
    #f1 = TF1("f1", "gaus", 750, 1050)
    #f2 = TF1("f2", "gaus", 750, 1050)
    #h1 = TH1D("h1", "", 100, 750, 1050)
    #h2 = TH1D("h2", "", 100, 750, 1050)
    #for i, j in zip(totpe, sumPreBetaPE):
    #    h1.Fill(i)
    #    h2.Fill(j)
    #h1.Fit(f1, "RE")
    #h2.Fit(f2, "RE")
    #
    #plotx1, ploty1, plotx2, ploty2 = [], [], [], []
    #for i in range(750, 1050, 1):
    #    plotx1.append(i)
    #    plotx2.append(i)
    #    ploty1.append(f1.Eval(i))
    #    ploty2.append(f2.Eval(i))
    
    #print("Simulation NPE: %.2f , Calculation NPE : %.2f" % (totpe.mean(), sumPreBetaPE.mean()) ) 
    ##plt.plot(totpe, "o", ms=0.5, label="Simulation")
    ##plt.plot(sumPreBetaPE, "o", ms=0.5, label="Calculation")
    #plt.hist(totpe,        bins=100, range=(750, 1050), alpha=0.5, color="royalblue", label="Simulation %.2f+-%.2f" %(totpe.mean(), np.std(totpe)) )
    #plt.hist(sumPreBetaPE, bins=100, range=(750, 1050), alpha=0.5, color="peru",      label="Calculation %.2f +- %.2f" %(sumPreBetaPE.mean(), np.std(sumPreBetaPE)) )
    #plt.plot(plotx1, ploty1, "-", color="red")
    #plt.plot(plotx2, ploty2, "-", color="red")
    #plt.xlabel("# P.E.")


def main():
    gammaSource("Cs137")
    gammaSource("Mn54")
    #loadFitterCurve()

    #plt.xlim(0, 0.5)
    #plt.ylim(0, 1000)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
