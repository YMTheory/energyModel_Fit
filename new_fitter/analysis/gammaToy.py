import numpy as np
import matplotlib.pyplot as plt
import random
import uproot3

# configureation :
saveHist = False
import ROOT
#ff = ROOT.TFile("totpe.root", "recreate")

import prmBetaCheck as pbe
import fileLoader as loader

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

sctE, sctPE = loadPEFile("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/electron/sctPE.txt")
cerE, cerPE = loadPEFile("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/electron/cerPE.txt") 
resolE, resol = loadResFile("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/electron/elecResol.txt")

def draw():
    plt.style.use("seaborn-muted")
    #totPE = np.array(sctPE) + np.array(cerPE)
    #plt.plot(sctE, sctPE, "-",  label="scintillation")
    #plt.plot(cerE, cerPE, "-.", label="cerenkov")
    #plt.plot(sctE, totPE, "--", label="total")
    #plt.ylabel("# P.E.")
    plt.plot(resolE, resol, "-")
    plt.xlabel("Etrue/MeV")
    plt.ylabel("PE sigma")
    plt.semilogy()
    plt.grid(True)

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
    #prmBetaPE = evt["prmBetaPE"].array()
    #prmBetaE = evt["prmBetaEnergy"].array()
    #return totpe, prmBetaPE, prmBetaE
    return totpe

def loadPrmBeta(filename):
    secBetaArr, secAntiBetaArr = [], []
    with open(filename) as f:
        for lines in f.readlines():
            oneEvtBeta, oneEvtAntiBeta = [], []
            line = lines.strip("\n")
            data = line.split(" ")
            counta = 0
            countb = 0
            for i in data:
                if "a" in i:
                    counta+=1
                    tmp = list(i)
                    tmp.pop()
                    j = ''.join(tmp)
                    oneEvtAntiBeta.append(float(j))
                #    hh2.SetBinContent(evtid+1, counta, float(j))
                if "b" in i:
                    countb+=1
                    tmp = list(i)
                    tmp.pop()
                    j = ''.join(tmp)
                    oneEvtBeta.append(float(j))
            #secBetaArr.append(tmpE)
            secAntiBetaArr.append(oneEvtAntiBeta)
            secBetaArr.append(oneEvtBeta)
    return secBetaArr, secAntiBetaArr


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

def GammaPePred(secBeta, secAntiBeta):
    sc = 1.00
    calcpe = []
    sample_id = np.random.randint(0, 5000, size=10000)
    for i in sample_id:
        tmppe = 0
        for j in secBeta[i]:
            tmppe += random.gauss(g1.Eval(j, 0, "S"), g2.Eval(j, 0, "S")*sc)
        for j in secAntiBeta[i]:
            tmppe += random.gauss(g1.Eval(j, 0, "S"), g2.Eval(j, 0, "S")*sc) + 2*random.gauss(660.8, 27.07*sc)
        calcpe.append(tmppe)
    return calcpe

import time
def gammaSource(name, Etrue):
    start = time.time()
    path = "./data/"
    totpe, prmBetaPE, prmBetaE = [], [], []
    filename = path + name + "_totpe.root"
    print("------> Reading " + filename)
    tmptotpe = loadSimTruth(filename)
    totpe.extend(tmptotpe)
    
    secBetaArr, secAntiBetaArr = loadPrmBeta("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/gamma/"+name+"/log-30.txt")

    calcpe = GammaPePred(secBetaArr, secAntiBetaArr)
    #pbe.prmBetaDist(name, Etrue, secBetaArr, secAntiBetaArr)

    totpe = np.array(totpe)
    calcpe = np.array(calcpe)
    low = totpe.min() - 50
    high = totpe.max() + 50

    f1 = TF1("f1", "gaus", low, high)
    f2 = TF1("f2", "gaus", low, high)
    h1 = TH1D(name+"h1", "", 100, low, high)
    h2 = TH1D(name+"h2", "", 100, low, high)
    for i, j in zip(totpe, calcpe):
        h1.Fill(i)
        h2.Fill(j)
    h1.Fit(f1, "RE")
    h2.Fit(f2, "RE")

    if saveHist:
        ff.cd()
        h1.Write()
        h2.Write()

    sim_mean, sim_sigma   = f1.GetParameter(1), f1.GetParameter(2)
    calc_mean, calc_sigma = f2.GetParameter(1), f2.GetParameter(2)
    sim_mean_err, sim_sigma_err   = f1.GetParError(1), f1.GetParError(2)
    calc_mean_err, calc_sigma_err = f2.GetParError(1), f2.GetParError(2)

    print(name + " -> Sim Mean: %.3f, Calc Mean: %.3f" %(totpe.mean(), calcpe.mean()) )

    #plt.hist(totpe,  bins=100, range=(low, high), density=True, histtype='step', color="royalblue", label="J19 Sim "+name)
    #plt.hist(calcpe, bins=100, range=(low, high), density=True, histtype='step', color="seagreen", label="J19 based calc "+name)
    #plt.hist(totpe,  bins=100, density=True, histtype='step', color="royalblue", label="J19 Sim "+name)
    plt.hist(calcpe, bins=100, density=True, histtype='step', color="seagreen", label="J19 based calc "+name)

    plt.xlabel("# P.E.")

    end = time.time()
    print("Total time consumed for " + name + " %s" % (end-start))
    return sim_mean, sim_mean_err, sim_sigma, sim_sigma_err, calc_mean, calc_mean_err, calc_sigma, calc_sigma_err

    """ prmBeta-single-chekc """
    
    #simPrmBetaE, simPrmBetaPE, calcPrmBetaPE = [], [], []
    #binningSim, binningCalc = [[] for i in range(50)], [[] for i in range(50)]
    #meanSim, meanCalc = [], []
    #for i in range(len(totpe)):
    #    if len(prmBetaPE[i]) != len(prmBetaE[i]):
    #        continue
    #    for j in range(len(prmBetaPE[i])):
    #        simPrmBetaE.append(prmBetaE[i][j])
    #        simPrmBetaPE.append(prmBetaPE[i][j])
    #        calcPrmBetaPE.append(random.gauss(g1.Eval(prmBetaE[i][j], 0, "S"), g2.Eval(prmBetaE[i][j])))
            #if 0.053 < simPrmBetaE[-1] < 0.054:
            #    print(simPrmBetaE[-1], simPrmBetaPE[-1], calcPrmBetaPE[-1])
    ''' 
    diff = []
    for i, j in zip(simPrmBetaPE, calcPrmBetaPE) :
        diff.append(j-i)
    diff = np.array(diff)
    plt.hist(diff, bins=100, histtype="step", label=name+" mean %.3f"%diff.mean() )

    with uproot3.recreate(name+"_bias.root") as f:
        f["t"] = uproot3.newtree({"diff" : "float64",
                                  "energy" : "float64"})
        f["t"].extend({"diff": diff,
                       "energy": simPrmBetaE})
    '''

    ''' binning '''
    '''
    for a, b, e in zip(simPrmBetaPE, calcPrmBetaPE, simPrmBetaE):
        if e > 0.5:
            continue
        binid = int(e/0.01)
        binningSim[binid].append(a)
        binningCalc[binid].append(b)

    meanFitter = []
    for i in range(50):
        tmpSim = np.array(binningSim[i])
        tmpCalc = np.array(binningCalc[i])
        meanSim.append(tmpSim.mean())
        meanCalc.append(tmpCalc.mean())
        #meanSim.append(np.std(tmpSim))
        #meanCalc.append(np.std(tmpCalc))

    Eeee = np.arange(0, 0.5, 0.01)
    for i in Eeee:
        meanFitter.append(g1.Eval(i+0.005, 0, "S"))

    meanSim = np.array(meanSim)
    meanCalc = np.array(meanCalc)
    #plt.plot(Eeee+0.005, meanSim,  "o-", ms=0.3, label="Sim")
    #plt.plot(Eeee+0.005, meanCalc, "s-", ms=0.3, label="Calc")
    #plt.plot(Eeee, meanFitter, "--", label="Fitter Curve")
    #plt.plot(resolE, resol, "-")
    #plt.xlim(0, 0.5)
    #plt.ylim(0, 40)
    plt.plot(Eeee+0.005, (meanCalc-meanSim)/meanSim, "-", label=name+"relative difference")
    plt.hlines(0, Eeee[0], Eeee[-1], colors="red", linestyle='--')
    plt.xlabel("Etrue/MeV")
    plt.ylabel("# P.E.")
    plt.grid(True)
    ''' 

    #plt.plot(simPrmBetaE, simPrmBetaPE, "o",  ms=0.5, alpha=0.5, label="Simulation")
    #plt.plot(simPrmBetaE, calcPrmBetaPE, "o", ms=0.5, alpha=0.5, label="Calculation")
    #plt.xlabel("Etrue/MeV")
    #plt.ylabel("# P.E.")


    ''' final totpe distribution '''
    ''''
    sumPreBetaPE = []
    sumPreBetaPE1 = []
    sumPreBetaPE2 = []
    for i in range(len(totpe)):
        tmpsum = 0
        tmpsum1 = 0
        tmpsum2 = 0
        for j in range(len(prmBetaPE[i])):
            #tmpsum += g1.Eval(prmBetaE[i][j])
            tmpsum  += random.gauss(g1.Eval(prmBetaE[i][j]), g2.Eval(prmBetaE[i][j]))    # consider resolution
            tmpsum1 += random.gauss(g1.Eval(prmBetaE[i][j]), g2.Eval(prmBetaE[i][j])*0.95)    # consider resolution
            tmpsum2 += random.gauss(g1.Eval(prmBetaE[i][j]), g2.Eval(prmBetaE[i][j])*1.05)    # consider resolution
        sumPreBetaPE.append(tmpsum)
        sumPreBetaPE1.append(tmpsum1)
        sumPreBetaPE2.append(tmpsum2)
    
    totpe = np.array(totpe)
    sumPreBetaPE = np.array(sumPreBetaPE)
    sumPreBetaPE1 = np.array(sumPreBetaPE1)
    sumPreBetaPE2 = np.array(sumPreBetaPE2)

    f1 = TF1("f1", "gaus", 750, 1050)
    f2 = TF1("f2", "gaus", 750, 1050)
    f3 = TF1("f4", "gaus", 750, 1050)
    f4 = TF1("f4", "gaus", 750, 1050)
    h1 = TH1D("h1", "", 100, 750, 1050)
    h2 = TH1D("h2", "", 100, 750, 1050)
    h3 = TH1D("h3", "", 100, 750, 1050)
    h4 = TH1D("h4", "", 100, 750, 1050)
    for i, j, k, p in zip(totpe, sumPreBetaPE, sumPreBetaPE1, sumPreBetaPE2):
        h1.Fill(i)
        h2.Fill(j)
        h3.Fill(k)
        h4.Fill(p)
    h1.Fit(f1, "RE")
    h2.Fit(f2, "RE")
    h3.Fit(f3, "RE")
    h4.Fit(f4, "RE")
    
    plotx1, ploty1, plotx2, ploty2 = [], [], [], []
    plotx3, ploty3, plotx4, ploty4 = [], [], [], []
    for i in range(750, 1050, 1):
        plotx1.append(i)
        plotx2.append(i)
        ploty1.append(f1.Eval(i))
        ploty2.append(f2.Eval(i))
        plotx3.append(i)
        plotx4.append(i)
        ploty3.append(f3.Eval(i))
        ploty4.append(f4.Eval(i))
    
    print("Simulation NPE: %.2f , Calculation NPE : %.2f" % (totpe.mean(), sumPreBetaPE.mean()) ) 
    #plt.plot(totpe, "o", ms=0.5, label="Simulation")
    #plt.plot(sumPreBetaPE, "o", ms=0.5, label="Calculation")
    plt.plot(plotx1, ploty1, "--", color="royalblue")
    plt.plot(plotx2, ploty2, "--", color="peru")
    plt.plot(plotx3, ploty3, "--", color="darkviolet")
    plt.plot(plotx4, ploty4, "--", color="lightseagreen")
    plt.hist(totpe,        bins=100, range=(750, 1050), histtype='step', color="royalblue", label="Simulation %.2f+-%.2f" %(f1.GetParameter(1), f1.GetParameter(2)) )
    plt.hist(sumPreBetaPE, bins=100, range=(750, 1050), histtype='step', color="peru",      label="Calculation %.2f +- %.2f" %(f2.GetParameter(1), f2.GetParameter(2)) )
    plt.hist(sumPreBetaPE1, bins=100, range=(750, 1050), histtype='step', color="darkviolet",      label="Calc->0.05 down e- resol %.2f +- %.2f" %(f3.GetParameter(1), f3.GetParameter(2)) )
    plt.hist(sumPreBetaPE2, bins=100, range=(750, 1050), histtype='step', color="lightseagreen",      label="Calc->0.05 up e- resol %.2f +- %.2f" %(f4.GetParameter(1), f4.GetParameter(2)) )
    plt.xlabel("# P.E.")
    '''

#name = ["Ge68", "Cs137", "Mn54", "Co60", "K40", "nH", "AmBe", "nC12", "AmC"]
name = ["Cs137", "Mn54", "Ge68", "K40", "nH", "Co60", "AmBe", "nC12", "AmC"]
Etrue = [0.662, 0.835, 1.022, 1.461, 2.223, 2.506, 4.43, 4.94, 6.13]
Etrue_nonl = [0.662, 0.835, 0.511, 1.461, 2.223, 1.253, 4.43, 4.94, 6.13]

def GammaCollection():
    start = time.time()

    simMeanArr      = [0 for i in range(9)]
    simMeanErrArr   = [0 for i in range(9)]
    simSigmaArr     = [0 for i in range(9)]
    simSigmaErrArr  = [0 for i in range(9)]
    calcMeanArr     = [0 for i in range(9)]
    calcMeanErrArr  = [0 for i in range(9)]
    calcSigmaArr    = [0 for i in range(9)]
    calcSigmaErrArr = [0 for i in range(9)]

    for i in range(9):
        simMeanArr[i], simMeanErrArr[i], simSigmaArr[i], simSigmaErrArr[i], \
        calcMeanArr[i], calcMeanErrArr[i], calcSigmaArr[i], calcSigmaErrArr[i] = gammaSource(name[i], Etrue[i])

    scale = simMeanArr[5] / Etrue[5]
    simnonl, simnonl_err, simres, simres_err = [], [], [], []
    calcnonl, calcnonl_err, calcres, calcres_err = [], [], [], []
    for i in range(9):
        if name[i] == "Ge68" or name[i]=="Co60":
            simnonl.append(simMeanArr[i]/2/scale/Etrue[i])
            simnonl_err.append(simMeanErrArr[i]/2/scale/Etrue[i])
            calcnonl.append(calcMeanArr[i]/2/scale/Etrue[i])
            calcnonl_err.append(calcMeanErrArr[i]/2/scale/Etrue[i])
        else:
            simnonl.append(simMeanArr[i]/scale/Etrue[i])
            simnonl_err.append(simMeanErrArr[i]/scale/Etrue[i])
            calcnonl.append(calcMeanArr[i]/scale/Etrue[i])
            calcnonl_err.append(calcMeanErrArr[i]/scale/Etrue[i])

        simres.append(simSigmaArr[i]/simMeanArr[i])
        simres_err.append(np.sqrt(simSigmaErrArr[i]**2/simMeanArr[i]**2 + simMeanErrArr[i]**2*simSigmaArr[i]**2/simMeanArr[i]**4))
        calcres.append(calcSigmaArr[i]/calcMeanArr[i])
        calcres_err.append(np.sqrt(calcSigmaErrArr[i]**2/calcMeanArr[i]**2 + calcMeanErrArr[i]**2*calcSigmaArr[i]**2/calcMeanArr[i]**4))

    with open("test.txt", "w") as f:
        for i in range(9):
            f.write("%.5f %.6f %.5f %.6f %.5f %.6f %.5f %.6f" %(simnonl[i], simnonl_err[i], calcnonl[i], calcnonl_err[i], simres[i], simres_err[i], calcres[i], calcres_err[i]))
            f.write("\n")
    
    end = time.time()
    print("All Gamma Soueces Prediction Finished with %s" %(end-start))
    


def readPredResults(filename):
    simnonl, simnonl_err, simres, simres_err = [], [], [], []
    calcnonl, calcnonl_err, calcres, calcres_err = [], [], [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            simnonl.append(float(data[0]))
            simnonl_err.append(float(data[1]))
            #simnonl_err.append(0.001)
            calcnonl.append(float(data[2]))
            calcnonl_err.append(float(data[3]))
            #calcnonl_err.append(0.001)
            simres.append(float(data[4]))
            simres_err.append(float(data[5]))
            calcres.append(float(data[6]))
            calcres_err.append(float(data[7]))


    # change nonlinearity data order :
    new_id = [2, 0, 1, 5, 3, 4, 6, 7, 8]
    simnonl_new, simnonl_new_err, calcnonl_new, calcnonl_new_err, Etrue_nonl_new = [], [], [], [], []
    for i in new_id:
        Etrue_nonl_new.append(Etrue_nonl[i])
        simnonl_new.append(simnonl[i]*Etrue[i]/Etrue_nonl[i])
        simnonl_new_err.append(simnonl_err[i]*Etrue[i]/Etrue_nonl[i])
        calcnonl_new.append(calcnonl[i]*Etrue[i]/Etrue_nonl[i])
        calcnonl_new_err.append(calcnonl_err[i]*Etrue[i]/Etrue_nonl[i])


    nonl_diff = (np.array(calcnonl_new) - np.array(simnonl_new)) / np.array(simnonl_new)
    nonl_diff_err = np.sqrt(np.array(calcnonl_new_err)**2/np.array(simnonl_new) + np.array(simnonl_new_err)**2 * np.array(calcnonl_new)**2 / np.array(simnonl_new)**4)
    res_diff = (np.array(calcres) - np.array(simres)) / np.array(simres)
    res_diff_err = np.sqrt(np.array(calcres_err)**2/np.array(simres) + np.array(simres_err)**2 * np.array(calcres)**2 / np.array(simres)**4)



    plt.figure(0, figsize=(6, 4))
    plt.errorbar(Etrue_nonl_new, simnonl_new, yerr=simnonl_new_err, fmt="o-", color="royalblue", label="Sim")
    plt.errorbar(Etrue_nonl_new, calcnonl_new, yerr=calcnonl_new_err, fmt="o-", color="lightseagreen", label="Calc")
    plt.grid(True)
    plt.xlabel("Etrue/MeV")
    plt.ylabel("nonlinearity")
    plt.legend()
    plt.savefig("nonl.pdf")


    plt.figure(1, figsize=(6, 4))
    plt.errorbar(Etrue, simres, yerr=simres_err, fmt="o-", color="royalblue", label="Sim")
    plt.errorbar(Etrue, calcres, yerr=calcres_err, fmt="o-", color="lightseagreen", label="Calc")
    plt.grid(True)
    plt.xlabel("Etrue/MeV")
    plt.ylabel("resolution")
    plt.legend()
    plt.savefig("res.pdf")


    plt.figure(2, figsize=(6, 2))
    plt.errorbar(Etrue_nonl_new, nonl_diff, yerr=nonl_diff_err, fmt="o", color="peru")
    plt.fill_between([0, 6.2], [-0.002, -0.002], [0.002, 0.002], alpha=0.3)
    plt.hlines(0, 0, 6.3, linestyle="--", color="red")
    plt.ylim(-0.005, 0.005)
    plt.grid(True)
    plt.xlabel("Etrue/MeV")
    plt.ylabel("relative different")
    plt.savefig("nonlDiff.pdf")


    plt.figure(3, figsize=(6, 2))
    plt.errorbar(Etrue, res_diff, yerr=res_diff_err, fmt="o", color="peru")
    plt.ylim(-0.1, 0.1)
    plt.fill_between([0, 6.2], [-0.03, -0.03], [0.03, 0.03], alpha=0.3)
    plt.hlines(0, 0, 6.3, linestyle="--", color="red")
    plt.grid(True)
    plt.xlabel("Etrue/MeV")
    plt.ylabel("relative different")
    plt.savefig("resDiff.pdf")



def resultsCompare():
    Etrue_nonl_new0, simnonl_new0, simnonl_new_err0, calcnonl_new0, calcnonl_new_err0, \
    Etrue0, simres0, simres_err0, calcres0, calcres_err0, resdiff0, resdiff_err0 = loader.ReadModelPrediction("./resdown10.txt")

    Etrue_nonl_new1, simnonl_new1, simnonl_new_err1, calcnonl_new1, calcnonl_new_err1, \
    Etrue1, simres1, simres_err1, calcres1, calcres_err1, resdiff1, resdiff_err1 = loader.ReadModelPrediction("./resdown5.txt")

    Etrue_nonl_new2, simnonl_new2, simnonl_new_err2, calcnonl_new2, calcnonl_new_err2, \
    Etrue2, simres2, simres_err2, calcres2, calcres_err2,resdiff2, resdiff_err2 = loader.ReadModelPrediction("./normal.txt")

    Etrue_nonl_new3, simnonl_new3, simnonl_new_err3, calcnonl_new3, calcnonl_new_err3, \
    Etrue3, simres3, simres_err3, calcres3, calcres_err3, resdiff3, resdiff_err3 = loader.ReadModelPrediction("./resup5.txt")
    
    Etrue_nonl_new4, simnonl_new4, simnonl_new_err4, calcnonl_new4, calcnonl_new_err4, \
    Etrue4, simres4, simres_err4, calcres4, calcres_err4, resdiff4, resdiff_err4 = loader.ReadModelPrediction("./resup10.txt")

    plt.figure(0)
    plt.plot(Etrue_nonl_new0, simnonl_new0,  "o-", label="Sim")
    plt.plot(Etrue_nonl_new0, calcnonl_new0, "o-", label="e- resol 0.10 down")
    plt.plot(Etrue_nonl_new1, calcnonl_new1, "o-", label="e- resol 0.05 down")
    plt.plot(Etrue_nonl_new2, calcnonl_new2, "o-", label="e- resol unchanged")
    plt.plot(Etrue_nonl_new3, calcnonl_new3, "o-", label="e- resol 0.05 up")
    plt.plot(Etrue_nonl_new4, calcnonl_new4, "o-", label="e- resol 0.10 up")
    plt.xlabel("Etrue/MeV")
    plt.ylabel("nonlinearity")
    plt.grid(True)

    plt.figure(1)
    plt.plot(Etrue0, simres0,  "o-", label="Sim")
    plt.plot(Etrue0, calcres0, "o-", label="e- resol 0.10 down")
    plt.plot(Etrue1, calcres1, "o-", label="e- resol 0.05 down")
    plt.plot(Etrue2, calcres2, "o-", label="e- resol unchanged")
    plt.plot(Etrue3, calcres3, "o-", label="e- resol 0.05 up")
    plt.plot(Etrue4, calcres4, "o-", label="e- resol 0.10 up")
    plt.xlabel("Etrue/MeV")
    plt.ylabel("resolution")
    plt.grid(True)

    plt.figure(2)
    plt.errorbar(Etrue0, resdiff0, yerr=resdiff_err0, fmt="o-", label="e- resol 0.10 down")
    plt.errorbar(Etrue1, resdiff1, yerr=resdiff_err1, fmt="o-", label="e- resol 0.05 down")
    plt.errorbar(Etrue2, resdiff2, yerr=resdiff_err2, fmt="o-", label="e- resol unchanged")
    plt.errorbar(Etrue3, resdiff3, yerr=resdiff_err3, fmt="o-", label="e- resol 0.05 up")
    plt.errorbar(Etrue4, resdiff4, yerr=resdiff_err4, fmt="o-", label="e- resol 0.10 up")
    plt.fill_between([0, 6.2], [-0.03, -0.03], [0.03, 0.03], alpha=0.3, color="green")
    plt.hlines(0, 0, 6.3, linestyle="--", color="red")
    plt.ylim(-0.2, 0.2)
    plt.xlabel("Etrue/MeV")
    plt.ylabel("resolution")
    plt.grid(True)




def main():
    GammaCollection()
    #draw()
    #readPredResults("./resdown5.txt")
    #gammaSource(name[0], Etrue[0])
    #resultsCompare()

    #ff.Close()
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
