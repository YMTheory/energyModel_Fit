#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
import uproot as up

file = up.open("./gammaAllFit.root")
#print(dict(file.classes()))

def plot_nonl():
    ## Nonlineaity
    tge_nonlData = file["gNonlData"]
    tg_nonlCalc = file["gNonlCalc"]

    tg_nonlCalcX = tg_nonlCalc.xvalues
    tg_nonlCalcY = tg_nonlCalc.yvalues
    tge_nonlDataX = tge_nonlData.xvalues
    tge_nonlDataY = tge_nonlData.yvalues
    tge_nonlDataXerr = tge_nonlData._fEX
    tge_nonlDataYerr = tge_nonlData._fEY

    plt.errorbar(tge_nonlDataX, tge_nonlDataY, xerr=tge_nonlDataXerr, yerr=tge_nonlDataYerr, fmt="o--", color="blue", label="Nonlinearity Data")
    plt.plot(tg_nonlCalcX, tg_nonlCalcY, "s--", color="red",label="Nonliearity Calculation")
    plt.xlabel("gamma Etrue/MeV")
    plt.ylabel("Nonlinearity")
    plt.title("Nonlinearity Model")
    plt.legend()
    plt.show()
    #plt.savefig("gammaNonl.svg")

def plot_resol():
    # resolution
    tge_resData = file["gResData"]
    tge_resDataX = tge_resData.xvalues
    tge_resDataY = tge_resData.yvalues
    tge_resDataXerr = tge_resData._fEX
    tge_resDataYerr = tge_resData._fEY

    tg_resCalc = file["gResCalc"]
    tg_resCalcX = tg_resCalc.xvalues
    tg_resCalcY = tg_resCalc.yvalues

    plt.errorbar(tge_resDataX, tge_resDataY, xerr=tge_resDataXerr, yerr=tge_resDataYerr, fmt="o--", color="blue", label="Resolution Data")
    plt.plot(tg_resCalcX, tg_resCalcY, "s--", color="red",label="Resolution Calculation")
    plt.xlabel("gamma Etrue/MeV")
    plt.ylabel("Resolution")
    plt.title("Resolution Model")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    plot_nonl();



#tge_nonlData.matplotlib(showtitle="Nonliearity Data", show=True, fmt="o--", color="blue")
#tge_nonlCalc.matplotlib(showtitle="Nonliearity Pred", show=True, fmt='s--', color='red')