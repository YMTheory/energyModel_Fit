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
    Earr, resol = [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Earr.append(float(data[0]))
            resol.append(float(data[3]))
    Earr  = np.array(Earr)
    resol = np.array(resol)
    return Earr, resol

st = time.time()
sctE, sctPE = loadPEFile("../data/electron/sctPE1.txt")
print(">>>>>>>>>> Load scintillation PE file <<<<<<<<<<")
cerE, cerPE = loadPEFile("../data/electron/cerPE1.txt") 
print(">>>>>>>>>> Load cerenkov PE file <<<<<<<<<<")
resolE, resol = loadResFile("../data/electron/elecResol1.txt")
print(">>>>>>>>>> Load resolution PE file <<<<<<<<<<")
et = time.time()
print("elecLoader time: %.3f s" %(et-st))

from ROOT import TGraph
g1 = TGraph()
g2 = TGraph()
for i in range(len(sctE)):
    g1.SetPoint(i, sctE[i], sctPE[i]+cerPE[i])
for i in range(len(resolE)):
    g2.SetPoint(i, resolE[i], resol[i])


def getNPE(Etrue):
    return g1.Eval(Etrue, 0, "S")

def getSPE(Etrue):
    return g2.Eval(Etrue, 0, "S")


def getArray():
    sctE = np.array(sctE)
    sctPE = np.array(sctPE)
    cerPE = np.array(cerPE)
    resol = np.array(resol)
    return sctE, sctPE, cerPE, resol
