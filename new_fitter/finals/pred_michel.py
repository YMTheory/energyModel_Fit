import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from elecResponseLoader import elecResponseLoader
from outputHandler import outputHandler
import ROOT
import uproot as up
import random
plt.style.use("science")

if __name__ == "__main__" :

    ### Load simulation data
    nBinData, Emin, Emax = 100, 10, 57
    Y = 3111.44 / 2.223

    sim_file = "/junofs/users/miaoyu/simulation/LS_Sim/jobs/Det1/michel/michel_Det1.root"
    f = up.open(sim_file)
    totpe_sim = f["michel"]['totpe'].array()

    hsim = ROOT.TH1D("hsim", "", nBinData, Emin, Emax)    
    for i in totpe_sim:
        hsim.Fill(i / Y)

    ## generate calculation data
    Edep_file =  "/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/michel/michel_totpe_LS_v7.root"   
    f = up.open(Edep_file)
    edep = f["michel"]["edep"].array()

    qmode, cmode, rmode = "Int", "Ana1", "New"
    el = elecResponseLoader(qmode, cmode, rmode, "Det1")
    dirs = "Det1_kCerAnaNew1_kQInt_p2limits_B12"
    nonlfile = "/junofs/users/miaoyu/energy_model/fitter/energyModel_Fit/new_fitter/outputs/%s/nonlcov.txt"%dirs
    oh1 = outputHandler(nonlfile, 6)
    nonlpar, nonlparerr = oh1.loadBestFit() 
    el.setYs(nonlpar[0])
    el.setkB(nonlpar[1]) 
    if el.getCerenkovMode() == "Sim":
        el.setkC(nonlpar[2])
    if el.getCerenkovMode() == "Ana1":
        el.setp0(nonlpar[2])
        el.setp1(nonlpar[3])
        el.setp2(nonlpar[4])
        el.setE0(nonlpar[5])
    resfile = "/junofs/users/miaoyu/energy_model/fitter/energyModel_Fit/new_fitter/outputs/%s/rescov.txt"%dirs
    oh2 = outputHandler(resfile, 3)
    respar, resparerr = oh2.loadBestFit() 
    if rmode == "New":
        el.seta(respar[0])
        el.setb(respar[1])
        el.setn(respar[2])
    elif rmode == "New1":
        el.seta1(respar[0])
        el.setn1(respar[1])
    hcal = ROOT.TH1D("hcal", "", nBinData, Emin, Emax)
    for i in edep:
        evis = ( el.getScintillationNumber(i) + el.getCerenkovNumber(i) ) / Y
        esigma = el.getEvisSigma(evis)
        evis_new = random.gauss(evis, esigma)
        hcal.Fill(evis_new)



    x, y1, y2 = [], [], [] 
    for i in range(nBinData):
        x.append(hsim.GetBinCenter(i+1))
        y1.append(hsim.GetBinContent(i+1))
        y2.append(hcal.GetBinContent(i+1))
    x = np.array(x)
    y1 = np.array(y1); y1err = np.sqrt(y1)
    y2 = np.array(y2); y2err = np.sqrt(y2)
    dy = (y2 - y1)/y1
    dyerr = np.sqrt(y2err**2/y1**2 + y1err**2*y2**2/y1**4)

    
    fig = plt.figure(figsize=(10, 6))
    spec = gridspec.GridSpec(ncols=1, nrows=2,
                         height_ratios=[1, 2])

    ax0 = fig.add_subplot(spec[0])
    ax1 = fig.add_subplot(spec[1])

    ax0.errorbar(x, 100*dy, yerr=100*dyerr, fmt="o", ms=5, color="gray")
    ax0.hlines(0, 10, 57, color="red", lw=1.5)
    ax0.set_ylabel("Bias (\%)", fontsize=15)

    ax1.errorbar(x, y1, yerr=y1err, fmt="o", ms=5, label="Simulation")
    ax1.plot(x, y2, lw=2, label="Best fit")
    ax1.set_xlabel(r"$E_\mathrm{vis}$ [MeV]", fontsize=15)
    ax1.set_ylabel("counts per bin", fontsize=15)
    ax1.legend(loc="upper left", ncol=2, prop={"size":15})
    plt.tight_layout()
    plt.savefig("pred_michel.pdf")
    plt.show()



