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


import ROOT
def loadQuenchIntFile(filename):
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





st = time.time()

sctE, sctPE = loadPEFile("../data/electron/sctPE2.txt")
print(">>>>>>>>>> Load scintillation PE file <<<<<<<<<<")
cerE, cerPE = loadPEFile("../data/electron/cerPE3.txt") 
print(">>>>>>>>>> Load cerenkov PE file <<<<<<<<<<")
resolE, totpe, resol, resolerr = loadResFile("../data/electron/elecResol1.txt")
print(">>>>>>>>>> Load resolution PE file <<<<<<<<<<")
binCenter, graphArr = loadQuench("../data/electron/Quench5.root")

intCenter, intContent = loadQuenchIntFile("../data/electron/Quench_NumInt.root")

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
def getQPE(Etrue, kB, es, quenchMode):
    if kB <= birkLow:
        kB = birkLow + 1e-4
    if kB >= birkHigh:
        kB = birkHigh -1e-4

    kBIdx = int(kB*1e4) - 51
    kBResid = kBIdx+52-kB*1e4;

    qnl_low  = graphArr[kBIdx]
    qnl_high = graphArr[kBIdx+1]

    if quenchMode == "Int":
        qnl_low  = intContent[kBIdx]
        qnl_high = intContent[kBIdx+1]



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




m_a, m_b, m_c = 0.990, 7.78e-3, 0
es = 3134.078 / 2.223
def getFitResol(E):
    return np.sqrt(m_a**2/es/E + m_b**2)


def getFitConstResol(E):
    return np.sqrt(m_b**2)


def readFile(filename):
    Etrue, totpe, sigma, sigmaErr = [], [], [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Etrue.append(float(data[0]))
            totpe.append(float(data[1]))
            sigma.append(float(data[2]))
            sigmaErr.append(float(data[3]))

    Etrue = np.array(Etrue)
    totpe = np.array(totpe)
    sigma = np.array(sigma)
    sigmaErr = np.array(sigmaErr)

    return Etrue, totpe, sigma, sigmaErr


covE1, covPE1, cov1, covErr1 = readFile("../data/electron/elecPECov1.txt")
sctE1, sctPE1, sctSigma1, sctSigmaErr1 = readFile("../data/electron/elecSctPEResol1.txt")
cerE1, cerPE1, cerSigma1, cerSigmaErr1 = readFile("../data/electron/elecCerPEResol1.txt")


gCov = ROOT.TGraph()
gSctSigma = ROOT.TGraph()
gCerSigma = ROOT.TGraph()

kB, Asct, kC = 6.26e-3, 1408, 0.996
aa, bb, cc = 0.990, 7.78e-3, 0
c0, c1, c2 = 0, 1.774, 0.0577  #0, 5.318e-5, 9.52e-2
d0, d1, d2 = 0, 0.207, 2.376e-3


def func(x, p0, p1, p2):
    y = p0**2 + p1**2*x + p2**2*x**2
    return y

def funcConst(x, p0, p1, p2):
    y = p2**2*x**2
    return y


for i in range(len(covE1)):
    gCov.SetPoint(i, covE1[i], cov1[i])

for i in range(len(sctE1)):
    gSctSigma.SetPoint(i, sctE1[i], sctSigma1[i])

for i in range(len(cerE1)):
    gCerSigma.SetPoint(i, cerE1[i], cerSigma1[i])



def getTruthCov(E):
    return gCov.Eval(E)/es/es



def getTruthSctSigma(E):
    return gSctSigma.Eval(E)/es


def getTruthCerSigma(E):
    return gCerSigma.Eval(E)/es



def getFitCov(E):
    NPE = getNPE(E)
    return func(NPE, d0, d1, d2)


def getFitCerSigma(E):
    NPE = getCerNPE(E) * kC
    return func(NPE, c0, c1, c2)



def getFitSctSigma(E, kB, Asct, mode):
    NPE = getQPE(E, kB, Asct, mode)
    return (NPE)




def getFitNsct(E, kB, Asct, mode):
    return getQPE(E, kB, Asct, mode)


def getFitNcer(E, kC):
    return kC * getCerNPE(E)


def getFitCerSigmaConst(E):
    NPE = getCerNPE(E) * kC
    return funcConst(NPE, c0, c1, c2)


def getFitCovConst(E):
    NPE = getNPE(E)
    return funcConst(NPE, d0, d1, d2)








