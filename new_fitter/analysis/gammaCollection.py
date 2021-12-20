from singleGamma import singleGamma
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

name = ["Cs137", "Mn54", "Ge68", "K40", "nH", "Co60", "AmBe", "nC12", "AmC"]
name_nonl = ["Ge68", "Cs137", "Mn54", "Co60", "K40", "nH", "AmBe", "nC12", "AmC"]
Etrue = [0.662, 0.835, 1.022, 1.461, 2.223, 2.506, 4.43, 4.94, 6.13]
Etrue_nonl = [0.511, 0.662, 0.835, 1.253, 1.461, 2.223, 4.43, 4.94, 6.13]
new_id = [2, 0, 1, 5, 3, 4, 6, 7, 8]
nonlid1 = [1, 2, 4, 5, 6, 8]
nonlid2 = [0, 3, 7]
name1 = ["Cs137", "Mn54", "K40", "nH",  "AmBe", "AmC"]
Etrue1 = [0.662, 0.835, 1.461, 2.223, 4.43, 6.13]
Etrue_nonl1 = [ 0.662, 0.835, 1.461, 2.223, 4.43, 6.13]
ingle
ame2 = ["Ge68", "Co60", "nC12"]
Etrue2 = [1.022, 2.506, 4.94]
Etrue_nonl2 = [0.511, 1.253, 4.94]
resid1 = [0, 1, 3, 4, 6, 8]
resid2 = [2, 5, 7]
    

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
    nonlCalc, resCalc = [], []
    nonlData1, nonlData_err1, nonlCalc1 = [], [], []
    nonlData2, nonlData_err2, nonlCalc2 = [], [], []
    resData1, resData_err1, resCalc1 = [], [], []
    resData2, resData_err2, resCalc2 = [], [], []

    for i in sGamArr:
        i.calcTruth()
        i.ModelPrediction()

    scale = nH.getNPESim() / Etrue[4]
    
    for i in new_id:
        nonlCalc.append(sGamArr[i].getNPE()/scale/Etrue[i])

    for i in range(9):
        resCalc.append(sGamArr[i].getSPE() / sGamArr[i].getNPE())

    for i in resid1:
        nonlData1.append(sGamArr[i].getNPESim()/scale/Etrue[i])
        nonlData_err1.append(sGamArr[i].getNPEErrSim()/scale/Etrue[i])
        nonlCalc1.append(sGamArr[i].getNPE()/scale/Etrue[i])

    for i in resid2:
        nonlData2.append(sGamArr[i].getNPESim()/scale/Etrue[i])
        nonlData_err2.append(sGamArr[i].getNPEErrSim()/scale/Etrue[i])
        nonlCalc2.append(sGamArr[i].getNPE()/scale/Etrue[i])

    nonlDiff1 = (np.array(nonlCalc1) - np.array(nonlData1)) / np.array(nonlData1)
    nonlDiff1_err = (np.array(nonlData_err1) * np.array(nonlCalc1) ) / (np.array(nonlData1))**2
    nonlDiff2 = (np.array(nonlCalc2) - np.array(nonlData2)) / np.array(nonlData2)
    nonlDiff2_err = (np.array(nonlData_err2) * np.array(nonlCalc2) ) / (np.array(nonlData2))**2

    for i in resid1:
        resCalc1.append(sGamArr[i].getSPE() / sGamArr[i].getNPE())
        resData1.append(sGamArr[i].getSPESim()/sGamArr[i].getNPESim())
        resData_err1.append(np.sqrt(sGamArr[i].getSPEErrSim()**2/sGamArr[i].getNPESim()**2 + \
                                   sGamArr[i].getNPEErrSim()**2*sGamArr[i].getSPE()**2 / sGamArr[i].getNPESim()**4))
    resDiff1 = (np.array(resCalc1) - np.array(resData1)) / np.array(resData1)
    resDiff1_err = (np.array(resData_err1) * np.array(resCalc1) ) / (np.array(resData1))**2

    for i in resid2:
        resCalc2.append(sGamArr[i].getSPE() / sGamArr[i].getNPE())
        resData2.append(sGamArr[i].getSPESim()/sGamArr[i].getNPESim())
        resData_err2.append(np.sqrt(sGamArr[i].getSPEErrSim()**2/sGamArr[i].getNPESim()**2 + \
                                   sGamArr[i].getNPEErrSim()**2*sGamArr[i].getSPE()**2 / sGamArr[i].getNPESim()**4))
    resDiff2 = (np.array(resCalc2) - np.array(resData2)) / np.array(resData2)
    resDiff2_err = (np.array(resData_err2) * np.array(resCalc2) ) / (np.array(resData2))**2


    print(nonlCalc1)
    print(nonlData1)
    print(nonlData_err1)
    print(nonlDiff1)
    print(nonlDiff1_err)
    print(resCalc1)
    print(resData1)
    print(resData_err1)
    print(resDiff1)
    print(resDiff1_err)


    print(nonlCalc2)
    print(nonlData2)
    print(nonlData_err2)
    print(nonlDiff2)
    print(nonlDiff2_err)
    print(resCalc2)
    print(resData2)
    print(resData_err2)
    print(resDiff2)
    print(resDiff2_err)


    """

    fig = plt.figure(figsize=(12, 6))
    spec = gridspec.GridSpec(ncols=2, nrows=2,
                         height_ratios=[1, 2])

    ax0 = fig.add_subplot(spec[0])
    ax1 = fig.add_subplot(spec[1])
    ax2 = fig.add_subplot(spec[2])
    ax3 = fig.add_subplot(spec[3])
    
    ax2.plot(Etrue_nonl, nonlCalc, color="blue", label="Simulation")
    ax2.errorbar(Etrue_nonl1, nonlData1, yerr=nonlData_err1, fmt="o", color="magenta", label="Fitting: single gamma")
    ax2.errorbar(Etrue_nonl2, nonlData2, yerr=nonlData_err2, fmt="d", color="magenta", label="Fitting: multiple gamma")
    for i in range(len(name_nonl)):
        if i == 6:
            ax2.text(Etrue_nonl[i]-0.3, nonlCalc[i]-0.01, name_nonl[i], color="blue")
        else:
            ax2.text(Etrue_nonl[i], nonlCalc[i]-0.01, name_nonl[i], color="blue")
    ax2.legend()
    ax2.set_xlabel(r"$E_{dep}$/MeV", fontsize=14)
    ax2.set_ylabel(r"$E_{vis}/E_{dep}$", fontsize=14)
    ax2.set_ylim(0.905, 1.04)
    ax2.set_xlim(0, 6.6)
    #ax1.tight_layout()
    ax2.grid(True)
    ax2.tick_params(axis='both', which='major', labelsize=13)
    #ax1.savefig("nonlFit.pdf")
    ax2.legend(prop={'size': 13})

    #ax0.figure(1, figsize=(6, 3))
    ax0.errorbar(Etrue_nonl1, nonlDiff1, yerr=nonlDiff1_err, fmt="o", color="peru", label="single gamma")
    ax0.errorbar(Etrue_nonl2, nonlDiff2, yerr=nonlDiff2_err, fmt="d", color="peru", label="multiple gamma")
    ax0.fill_between([0, 6.3], [-0.001, -0.001], [0.001, 0.001], color="royalblue", alpha=0.3)
    ax0.text(2.5, 0.002, "0.1% uncertainty band", color="darkviolet")
    ax0.hlines(0, 0, 6.3, linestyle="--", color="red")
    #ax0.set_xlabel("Etrue/MeV")
    ax0.set_ylabel("relative bias", fontsize=14)
    ax0.set_ylim(-0.003, 0.003)
    ax0.set_xlim(0, 6.6)
    #ax0.tight_layout()
    ax0.grid(True)
    ax0.tick_params(axis='both', which='major', labelsize=13)


    ax3.plot(Etrue, resCalc, color="blue", label="Simulation")
    ax3.errorbar(Etrue1, resData1, yerr=resData_err1, fmt="o", color="magenta", label="Fitting: single gamma")
    ax3.errorbar(Etrue2, resData2, yerr=resData_err2, fmt="d", color="magenta", label="Fitting: multiple gamma")
    for i in range(len(name)):
        if i != 6:
            ax3.text(Etrue[i]-0.05, resCalc[i]+0.001, name[i], color="blue")
        else:
            ax3.text(Etrue[i]-0.30, resCalc[i]+0.001, name[i], color="blue")
    ax3.legend()
    ax3.set_xlim(0, 6.6)
    ax3.set_xlabel(r"$E_{dep}$/MeV", fontsize=14)
    ax3.set_ylabel(r"$\sigma/N_{tot}$", fontsize=14)
    ax3.set_ylim(0.012, 0.040)
    ax3.grid(True)
    ax3.tick_params(axis='both', which='major', labelsize=13)
    ax3.legend(prop={'size': 13})

    ax1.errorbar(Etrue1, resDiff1, yerr=resDiff1_err, fmt="o", color="peru", label="single gamma")
    ax1.errorbar(Etrue2, resData2, yerr=resDiff2_err, fmt="d", color="peru", label="multiple gamma")
    ax1.set_ylabel("relative bias", fontsize=14)
    ax1.fill_between([0, 6.3], [-0.03, -0.03], [0.03, 0.03], color="royalblue", alpha=0.3)
    ax1.text(2.5, 0.035, "3% uncertainty band", color="darkviolet")
    ax1.hlines(0, 0, 6.3, linestyle="--", color="red")
    ax1.set_ylim(-0.04, 0.05)
    ax1.set_xlim(0, 6.6)
    ax1.grid(True)
    ax1.tick_params(axis='both', which='major', labelsize=13)



    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=0.02, hspace=0.02)
    plt.tight_layout()

    plt.savefig("modelVal_Penelope.pdf")
    plt.show()
    """
