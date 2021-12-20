import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

def readFile(filename):
    nonlCalc1, nonlData1, nonlData_err1, nonlDiff1, nonlDiff1_err= [], [], [], [], []
    nonlCalc2, nonlData2, nonlData_err2, nonlDiff2, nonlDiff2_err= [], [], [], [], []
    resCalc1, resData1, resData_err1, resDiff1, resDiff1_err= [], [], [], [], []
    resCalc2, resData2, resData_err2, resDiff2, resDiff2_err= [], [], [], [], []

    arr = [nonlCalc1, nonlData1, nonlData_err1, nonlDiff1, nonlDiff1_err,  \
           resCalc1, resData1, resData_err1, resDiff1, resDiff1_err,       \
           nonlCalc2, nonlData2, nonlData_err2, nonlDiff2, nonlDiff2_err,  \
           resCalc2, resData2, resData_err2, resDiff2, resDiff2_err        \
          ]

    num = 0 
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            for i in data:
                arr[num].append(float(i))
            num += 1

    return nonlCalc1, nonlData1, nonlData_err1, nonlDiff1, nonlDiff1_err, \
           nonlCalc2, nonlData2, nonlData_err2, nonlDiff2, nonlDiff2_err, \
           resCalc1, resData1, resData_err1, resDiff1, resDiff1_err, \
           resCalc2, resData2, resData_err2, resDiff2, resDiff2_err



if __name__ == "__main__" :

    nonlCalc1, nonlData1, nonlData_err1, nonlDiff1, nonlDiff1_err, \
    nonlCalc2, nonlData2, nonlData_err2, nonlDiff2, nonlDiff2_err, \
    resCalc1, resData1, resData_err1, resDiff1, resDiff1_err, \
    resCalc2, resData2, resData_err2, resDiff2, resDiff2_err  = readFile( "modelVal_Livermore.txt" )

    nonlCalc, resCalc = [], []
    nonlCalc.append(nonlCalc2[0])
    nonlCalc.append(nonlCalc1[0])
    nonlCalc.append(nonlCalc1[1])
    nonlCalc.append(nonlCalc2[1])
    nonlCalc.append(nonlCalc1[2])
    nonlCalc.append(nonlCalc1[3])
    nonlCalc.append(nonlCalc2[2])
    nonlCalc.append(nonlCalc1[4])
    nonlCalc.append(nonlCalc1[5])
    print(resCalc1)
    print(resCalc2)
    resCalc.append(resCalc1[0])
    resCalc.append(resCalc1[1])
    resCalc.append(resCalc2[0])
    resCalc.append(resCalc1[2])
    resCalc.append(resCalc1[3])
    resCalc.append(resCalc2[1])
    resCalc.append(resCalc1[4])
    resCalc.append(resCalc2[2])
    resCalc.append(resCalc1[5])


    nonlCalc3, nonlData3, nonlData_err3, nonlDiff3, nonlDiff3_err, \
    nonlCalc4, nonlData4, nonlData_err4, nonlDiff4, nonlDiff4_err, \
    resCalc3, resData3, resData_err3, resDiff3, resDiff3_err, \
    resCalc4, resData4, resData_err4, resDiff4, resDiff4_err  = readFile( "modelVal_Penelope.txt" )


    nonlCalcP, resCalcP = [], []
    nonlCalcP.append(nonlCalc4[0])
    nonlCalcP.append(nonlCalc3[0])
    nonlCalcP.append(nonlCalc3[1])
    nonlCalcP.append(nonlCalc4[1])
    nonlCalcP.append(nonlCalc3[2])
    nonlCalcP.append(nonlCalc3[3])
    nonlCalcP.append(nonlCalc4[2])
    nonlCalcP.append(nonlCalc3[4])
    nonlCalcP.append(nonlCalc3[5])
    print(resCalc1)
    print(resCalc2)
    resCalcP.append(resCalc3[0])
    resCalcP.append(resCalc3[1])
    resCalcP.append(resCalc4[0])
    resCalcP.append(resCalc3[2])
    resCalcP.append(resCalc3[3])
    resCalcP.append(resCalc4[1])
    resCalcP.append(resCalc3[4])
    resCalcP.append(resCalc4[2])
    resCalcP.append(resCalc3[5])





    name = ["Cs137", "Mn54", "Ge68", "K40", "nH", "Co60", "AmBe", "nC12", "AmC"]
    name_nonl = ["Ge68", "Cs137", "Mn54", "Co60", "K40", "nH", "nC12", "AmBe", "AmC"]
    Etrue = [0.662, 0.835, 1.022, 1.461, 2.223, 2.506,  4.43, 4.94, 6.13]
    Etrue_nonl = [0.511, 0.662, 0.835, 1.253, 1.461, 2.223, 4.15, 4.43,  6.13]
    Etrue_nonl1 = [ 0.662, 0.835, 1.461, 2.223, 4.43, 6.13]
    Etrue_nonl2 = [0.511, 1.253, 4.15]
    Etrue1 = [0.662, 0.835, 1.461, 2.223, 4.43, 6.13]
    Etrue2 = [1.022, 2.506, 4.94]
    name1 = ["Cs137", "Mn54", "K40", "nH", "AmBe", "AmC"]
    name2 = ["Ge68", "Co60", "nC12"]

    fig = plt.figure(figsize=(8, 6))
    spec = gridspec.GridSpec(ncols=1, nrows=2)

    ax0 = fig.add_subplot(spec[0])
    ax1 = fig.add_subplot(spec[1])


    ax0.errorbar(Etrue1, nonlDiff1, yerr=nonlDiff1_err, fmt="o", ms=8, color="royalblue", label="single gamma")
    ax0.errorbar(Etrue2, nonlDiff2, yerr=nonlDiff2_err, fmt="o", ms=8, color="royalblue", label="multiple gamma")
    ax0.errorbar(Etrue1, nonlDiff3, yerr=nonlDiff3_err, fmt="d", ms=8, color="seagreen")
    ax0.errorbar(Etrue2, nonlDiff4, yerr=nonlDiff4_err, fmt="d", ms=8, color="seagreen")
    ax0.fill_between([0, 6.3], [-0.001, -0.001], [0.001, 0.001], color="royalblue", alpha=0.3)
    ax0.text(2.5, 0.002, "0.1% uncertainty band", color="darkviolet", fontsize=13)
    ax0.hlines(0, 0, 6.3, linestyle="--", color="red")
    ax0.set_xlabel(r"$E_{dep}$/MeV", fontsize=14)
    ax0.set_ylabel("Relative bias", fontsize=14)
    ax0.set_ylim(-0.003, 0.003)
    ax0.set_xlim(0, 6.6)
    #ax0.tight_layout()
    ax0.grid(True)
    ax0.tick_params(axis='both', which='major', labelsize=13)
    for i in range(len(nonlDiff1)):
        if i == 1:
            ax0.text(Etrue1[i]-0.02, nonlDiff1[i]-0.001, name1[i], color="blue", fontsize=12)
            continue
        ax0.text(Etrue1[i]-0.05, nonlDiff1[i]+0.001, name1[i], color="blue", fontsize=12)
    for i in range(len(nonlDiff2)):
        ax0.text(Etrue2[i]-0.05, nonlDiff2[i]+0.001, name2[i], color="blue", fontsize=12)
    ax0.set_title("(a) Nonlinearity", fontsize=14)
  


    ax1.errorbar(Etrue1, resDiff1, yerr=resDiff1_err, fmt="o", ms=8, color="royalblue", label="Livermore")
    ax1.errorbar(Etrue2, resData2, yerr=resDiff2_err, fmt="o", ms=8, color="royalblue")
    ax1.errorbar(Etrue1, resDiff3, yerr=resDiff3_err, fmt="d", ms=8, color="seagreen", label="Penelope")
    ax1.errorbar(Etrue2, resData4, yerr=resDiff4_err, fmt="d", ms=8, color="seagreen")
    ax1.set_ylabel("Relative bias", fontsize=14)
    ax1.fill_between([0, 6.3], [-0.03, -0.03], [0.03, 0.03], color="royalblue", alpha=0.3)
    ax1.text(2.5, 0.035, "3% uncertainty band", color="darkviolet", fontsize=13)
    ax1.hlines(0, 0, 6.3, linestyle="--", color="red")
    ax1.set_ylim(-0.04, 0.07)
    ax1.set_xlim(0, 6.6)
    ax1.grid(True)
    ax1.tick_params(axis='both', which='major', labelsize=13)
    ax1.set_xlabel(r"$E_{dep}$/MeV", fontsize=14)
    ax1.legend(prop={"size":14})
    ax1.set_title("(b) Resolution", fontsize=14)

    plt.tight_layout()
    plt.savefig("modelDiff.pdf")
    plt.show()


    """
    fig = plt.figure(figsize=(12, 6))
    spec = gridspec.GridSpec(ncols=2, nrows=2,
                         height_ratios=[1, 2])

    ax0 = fig.add_subplot(spec[0])
    ax1 = fig.add_subplot(spec[1])
    ax2 = fig.add_subplot(spec[2])
    ax3 = fig.add_subplot(spec[3])
    
    ax2.plot(Etrue_nonl, nonlCalc, lw=2, color="blue")
    ax2.plot(Etrue_nonl, nonlCalcP, "--", lw=2, color="seagreen")
    ax2.errorbar(Etrue_nonl1, nonlData1, yerr=nonlData_err1, fmt="o", color="magenta", label="Fitting: single gamma")
    ax2.errorbar(Etrue_nonl2, nonlData2, yerr=nonlData_err2, fmt="d", color="magenta", label="Fitting: multiple gamma")
    for i in range(len(name_nonl)):
        if i == 6:
            ax2.text(Etrue_nonl[i]-0.3, nonlCalc[i]-0.01, name_nonl[i], color="blue", fontsize=11)
        else:
            ax2.text(Etrue_nonl[i], nonlCalc[i]-0.01, name_nonl[i], color="blue", fontsize=11)
    #ax2.legend()
    ax2.set_xlabel(r"$E_{dep}$/MeV", fontsize=15)
    ax2.set_ylabel(r"$E_{vis}/E_{dep}$", fontsize=15)
    ax2.set_ylim(0.905, 1.04)
    ax2.set_xlim(0, 6.6)
    #ax1.tight_layout()
    ax2.grid(True)
    ax2.tick_params(axis='both', which='major', labelsize=13)
    #ax1.savefig("nonlFit.pdf")
    #ax2.legend(prop={'size': 13})

    #ax0.figure(1, figsize=(6, 3))
    ax0.errorbar(Etrue_nonl1, nonlDiff1, yerr=nonlDiff1_err, fmt="o", color="royalblue", label="single gamma")
    ax0.errorbar(Etrue_nonl2, nonlDiff2, yerr=nonlDiff2_err, fmt="d", color="royalblue", label="multiple gamma")
    ax0.errorbar(Etrue_nonl1, nonlDiff3, yerr=nonlDiff3_err, fmt="o", color="seagreen")
    ax0.errorbar(Etrue_nonl2, nonlDiff4, yerr=nonlDiff4_err, fmt="d", color="seagreen")
    ax0.fill_between([0, 6.3], [-0.001, -0.001], [0.001, 0.001], color="royalblue", alpha=0.3)
    ax0.text(2.5, 0.002, "0.1% uncertainty band", color="darkviolet", fontsize=13)
    ax0.hlines(0, 0, 6.3, linestyle="--", color="red")
    #ax0.set_xlabel("Etrue/MeV")
    ax0.set_ylabel("Relative bias", fontsize=15)
    ax0.set_ylim(-0.003, 0.003)
    ax0.set_xlim(0, 6.6)
    #ax0.tight_layout()
    ax0.grid(True)
    ax0.tick_params(axis='both', which='major', labelsize=13)


    l1, = ax3.plot(Etrue, resCalc, "-", lw=2, color="blue")
    l2, = ax3.plot(Etrue, resCalcP, "--", lw=2, color="seagreen")
    ax3.errorbar(Etrue1, resData1, yerr=resData_err1, fmt="o", color="magenta", label="Simulation: single gamma")
    ax3.errorbar(Etrue2, resData2, yerr=resData_err2, fmt="d", color="magenta", label="Simulation: multiple gamma")
    for i in range(len(name)):
        if i != 6:
            ax3.text(Etrue[i]-0.05, resCalc[i]+0.001, name[i], color="blue", fontsize=11)
        else:
            ax3.text(Etrue[i]-0.30, resCalc[i]+0.001, name[i], color="blue", fontsize=11)
  
    first_legend = plt.legend(handles=[l1, l2], labels=["Livermore", "Penelope"], loc="center right", prop={'size': 13})
    plt.gca().add_artist(first_legend)
    
    ax3.legend()
    ax3.set_xlim(0, 6.6)
    ax3.set_xlabel(r"$E_{dep}$/MeV", fontsize=15)
    ax3.set_ylabel(r"$\sigma/N_{tot}$", fontsize=15)
    ax3.set_ylim(0.012, 0.040)
    ax3.grid(True)
    ax3.tick_params(axis='both', which='major', labelsize=13)
    ax3.legend(prop={'size': 13})

    ax1.errorbar(Etrue1, resDiff1, yerr=resDiff1_err, fmt="o", color="royalblue", label="single gamma")
    ax1.errorbar(Etrue2, resData2, yerr=resDiff2_err, fmt="d", color="royalblue", label="multiple gamma")
    ax1.errorbar(Etrue1, resDiff3, yerr=resDiff3_err, fmt="o", color="seagreen")
    ax1.errorbar(Etrue2, resData4, yerr=resDiff4_err, fmt="d", color="seagreen")
    ax1.set_ylabel("Relative bias", fontsize=15)
    ax1.fill_between([0, 6.3], [-0.03, -0.03], [0.03, 0.03], color="royalblue", alpha=0.3)
    ax1.text(2.5, 0.035, "3% uncertainty band", color="darkviolet", fontsize=13)
    ax1.hlines(0, 0, 6.3, linestyle="--", color="red")
    ax1.set_ylim(-0.04, 0.05)
    ax1.set_xlim(0, 6.6)
    ax1.grid(True)
    ax1.tick_params(axis='both', which='major', labelsize=13)



    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=0.02, hspace=0.02)
    plt.tight_layout()

    plt.savefig("modelVal_Livermore+Penelope.pdf")
    plt.show()

    """
