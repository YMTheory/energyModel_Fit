import numpy as np
import time

def loadPEFile(filename):
    Earr, PEarr = [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Earr.append(float(data[0]))
            PEarr.append(float(data[1]))
    Earr = np.array(Earr)
    PEarr = np.array(PEarr)
    return Earr, PEarr

def loadResFile(filename):
    totpe = []
    Earr, resol, resolerr = [], [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Earr.append(float(data[0]))
            totpe.append(float(data[1]))
            resol.append(float(data[3]))
            resolerr.append(float(data[4]))
    Earr  = np.array(Earr)
    resol = np.array(resol)
    resolerr = np.array(resolerr)
    return Earr, totpe, resol, resolerr


import uproot as up
import ROOT

def loadQuench(filename):
    ff = ROOT.TFile(filename, "read")
    binCenter, graph_arr = [], []
    for i in range(51, 76, 1):
        histname = "kB" + str(i)
        tmphist = ff.Get(histname)
        tmp_arr = []
        for j in range(tmphist.GetNbinsX()):
            if i == 51:
                binCenter.append(tmphist.GetBinCenter(j+1))
            tmp_arr.append(tmphist.GetBinContent(j+1))

        graph_arr.append(tmp_arr)

    return binCenter, graph_arr



def loadStopPow(filename):
    print(" >>> Loading Stoppint Power Data <<< ")
    Etrue, sp = [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Etrue.append(float(data[0]))
            sp.append(float(data[1]))

    return Etrue, sp



st = time.time()

sctE, sctPE = loadPEFile("../data/electron/sctPE1.txt")
print(">>>>>>>>>> Load scintillation PE file <<<<<<<<<<")
cerE, cerPE = loadPEFile("../data/electron/cerPE1.txt") 
print(">>>>>>>>>> Load cerenkov PE file <<<<<<<<<<")
resolE, totpe, resol, resolerr = loadResFile("../data/electron/elecResol1.txt")
print(">>>>>>>>>> Load resolution PE file <<<<<<<<<<")
binCenter, graphArr = loadQuench("../data/electron/Quench5.root")

etr, sp = loadStopPow("../data/electron/StopPower.txt")

et = time.time()

print("elecLoader time: %.3f s" %(et-st))



from ROOT import TGraph
g1 = TGraph()
g2 = TGraph()
g3 = TGraph()
gsct = TGraph()
gcer = TGraph()
for i in range(len(sctE)):
    gsct.SetPoint(i, sctE[i], sctPE[i])
    gcer.SetPoint(i, cerE[i], cerPE[i])
    g1.SetPoint(i, sctE[i], sctPE[i]+cerPE[i])
for i in range(len(resolE)):
    g2.SetPoint(i, resolE[i], resol[i])
    g3.SetPoint(i, totpe[i], resol[i])

def getSctNPE(Etrue):
    return gsct.Eval(Etrue)

def getCerNPE(Etrue):
    return gcer.Eval(Etrue)

def getNPE(Etrue):
    return g1.Eval(Etrue, 0, "S")

def getSPE(Etrue):
    return g2.Eval(Etrue, 0, "S")


birkLow  = 0.0051
birkHigh = 0.0075
def getQPE(Etrue, kB, es):
    if kB < birkLow:
        kB = birkLow
    if kB > birkHigh:
        kB = birkHigh

    kBIdx = int(kB*1e4) - 51
    kBResid = kBIdx+52-kB*1e4;
    qnl_low  = graphArr[kBIdx]
    qnl_high = graphArr[kBIdx+1]


    if Etrue < 0.1:
        idx = int(Etrue/0.001)
    else:
        idx = int((Etrue-0.1)/0.01) + 100

    quenchNL_low  = qnl_low[idx-1] + (qnl_low[idx]-qnl_low[idx-1]) * (Etrue - binCenter[idx-1]) / (binCenter[idx] - binCenter[idx-1])
    quenchNL_high = qnl_high[idx-1] + (qnl_high[idx]-qnl_high[idx-1]) * (Etrue - binCenter[idx-1]) / (binCenter[idx] - binCenter[idx-1])

    quenchNL = kBResid * quenchNL_low + (1-kBResid) *quenchNL_high

    return quenchNL * es * Etrue


m_kA = 1
m_kB = 0.0065
def Integral_BirkLaw(E):
    num = len(etr)
    m_sum = 0
    integ_part1, integ_part2 = 0, 0
    for i in range(1, num, 1):
        integ_part1 = 1./(1+m_kB*sp[i-1])
        integ_part2 = 1./(1+m_kB*sp[i])
        if etr[i] <= E:
            m_sum += (integ_part1+integ_part2) * (etr[i] - etr[i-1]) /2.
        else:
            break

    return m_sum * m_kA/E



def getArray():
    sctE = np.array(sctE)
    sctPE = np.array(sctPE)
    cerPE = np.array(cerPE)
    resol = np.array(resol)
    return sctE, sctPE, cerPE, resol


def getResolArray():
    return resolE, totpe, resol, resolerr
