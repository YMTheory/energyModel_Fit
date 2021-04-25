from singleGamma import singleGamma
import matplotlib.pyplot as plt
import numpy as np

names = ["Cs137", "Mn54", "Ge68", "K40", "nH", "Co60", "AmBe", "nC12", "AmC"]
gEtrue = [0.662, 0.835, 1.022, 1.461, 2.223, 2.506, 4.43, 4.94, 6.13]
gEtrue_nonl = [ 0.662, 0.835, 0.511, 1.461, 2.223, 1.253, 4.43, 4.94, 6.13]
new_id = [2, 0, 1, 5, 3, 4, 6, 7, 8]

sGamArr = []
Cs137 = singleGamma("Cs137"); sGamArr.append(Cs137)
Mn54  = singleGamma("Mn54") ; sGamArr.append(Mn54)
Ge68  = singleGamma("Ge68") ; sGamArr.append(Ge68)
K40   = singleGamma("K40")  ; sGamArr.append(K40)
nH    = singleGamma("nH")   ; sGamArr.append(nH)
Co60  = singleGamma("Co60") ; sGamArr.append(Co60)
AmBe  = singleGamma("AmBe") ; sGamArr.append(AmBe)
nC12  = singleGamma("nC12") ; sGamArr.append(nC12)
AmC   = singleGamma("AmC")  ; sGamArr.append(AmC)

def DrawHisto():
    for i in sGamArr:
        i.Compare()

def DrawCurves():
    Etrue, Etrue_nonl = [], []
    nonlData, nonlData_err, nonlCalc = [], [], []
    resData, resData_err, resCalc = [], [], []

    for i in sGamArr:
        i.calcTruth()
        i.ModelPrediction()

    scale = nH.getNPESim() / gEtrue[4]

    for i in new_id:
        Etrue_nonl.append(gEtrue_nonl[i])
        nonlData.append(sGamArr[i].getNPESim()/scale/gEtrue[i])
        nonlData_err.append(sGamArr[i].getNPEErrSim()/scale/gEtrue[i])
        nonlCalc.append(sGamArr[i].getNPE()/scale/gEtrue[i])

    for i in range(len(sGamArr)):
        Etrue.append(gEtrue[i])
        resCalc.append(sGamArr[i].getSPE() / sGamArr[i].getNPE())
        resData.append(sGamArr[i].getSPESim()/sGamArr[i].getNPESim())
        resData_err.append(np.sqrt(sGamArr[i].getSPEErrSim()**2/sGamArr[i].getNPESim()**2 + \
                                   sGamArr[i].getNPEErrSim()**2*sGamArr[i].getSPE()**2 / sGamArr[i].getNPESim()**4))

    plt.style.use("seaborn-deep")

    plt.figure(0)
    plt.plot(Etrue_nonl, nonlCalc, "o-", label="calculation")
    plt.errorbar(Etrue_nonl, nonlData, yerr=nonlData_err, fmt="o-", label="simulation")
    plt.legend()
    plt.grid(True)

    plt.figure(1)
    plt.plot(Etrue, resCalc, "o-", label="calculation")
    plt.errorbar(Etrue, resData, yerr=resData_err, fmt="o-", label="simulation")
    plt.legend()
    plt.grid(True)


    plt.show()
