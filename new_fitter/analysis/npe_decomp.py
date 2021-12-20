import numpy as np
import matplotlib.pyplot as plt
import ROOT
from singleGamma import singleGamma
import uproot as up
import elecLoader as el

def fileDealing(filename):
    tt = up.open(filename)["evt"]
    totpe = tt['totalPE'].array()
    cerpe = tt['cerPE'].array()
    sctpe = totpe - cerpe
    edep = tt['edep'].array()

    meanSct, meanCer = np.mean(sctpe), np.mean(cerpe)
    stdSct, stdCer = np.std(sctpe), np.std(cerpe)
    stdTot = np.std(totpe)

    return stdSct, stdCer, stdTot, np.mean(edep), meanSct+meanCer, meanCer



def fileLoading(filename):
    tt = up.open(filename)["evt"]
    totpe = tt['totalPE'].array()
    cerpe = tt['cerPE'].array()
    sctpe = totpe - cerpe
    edep = tt['edep'].array()

    return sctpe, cerpe, totpe, edep


def main():
    #########################################
    ######  Truth 
    # Electrons:
    Etrue_elec, tot_elec, sct_elec, cer_elec, npe_elec, ncer_elec = [], [], [], [], [], []
    for i in [100, 200, 320, 400, 501, 583, 702, 789]:
        filename = "/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/electron/user-"+str(i)+".root"
        stdSct, stdCer, stdTot, edep, totpe, cerpe = fileDealing(filename)
        Etrue_elec.append(edep)
        tot_elec.append(stdTot)
        sct_elec.append(stdSct)
        cer_elec.append(stdCer)
        npe_elec.append(totpe)
        ncer_elec.append(cerpe)

    Etrue_elec = np.array(Etrue_elec)
    tot_elec   = np.array(tot_elec)
    sct_elec   = np.array(sct_elec)
    cer_elec   = np.array(cer_elec)
    npe_elec   = np.array(npe_elec)
    ncer_elec  = np.array(ncer_elec)


    # Positron :
    Etrue_posi, tot_posi, sct_posi, cer_posi, npe_posi, ncer_posi = [], [], [], [], [], []
    for i in range(1, 80, 10):
        filename = "/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/positron/user-"+str(i)+".root"
        stdSct, stdCer, stdTot, edep, totpe, cerpe = fileDealing(filename)
        Etrue_posi.append(edep)
        tot_posi.append(stdTot)
        sct_posi.append(stdSct)
        cer_posi.append(stdCer)
        npe_posi.append(totpe)
        ncer_posi.append(cerpe)

    Etrue_posi = np.array(Etrue_posi)
    tot_posi   = np.array(tot_posi)
    sct_posi   = np.array(sct_posi)
    cer_posi   = np.array(cer_posi)
    npe_posi   = np.array(npe_posi)
    ncer_posi  = np.array(ncer_posi)


    # Gamma :
    # From gamma_decomp.py outputs:
    sigma1_gam = [135.75754411290345, 214.08054787476925, 701.0206474514603, 1465.5137067971593, 3233.109245421706, 4032.931805211035]
    sigma2_gam = [856.3847286529768, 1100.8440713142763, 2069.100274912067, 3385.8515058172984, 7797.8901612275795, 11827.917014913259]
    Etrue_gam, tot_gam, sct_gam, cer_gam, npe_gam, ncer_gam = [], [], [], [], [], []
    names = ['Cs137', "Mn54", "K40", 'nH', 'AmBe', 'AmC']
    for j in range(6):
        Sct, Cer, Tot, edep, npe = [], [], [], [], []
        for i in range(10):
            filename = "/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/gamma/"+names[j]+"/user-"+str(i)+".root"
            tmpSct, tmpCer, tmpTot, tmpedep = fileLoading(filename)
            for i in range(len(tmpSct)):
                Sct.append(tmpSct[i])
                Cer.append(tmpCer[i])
                Tot.append(tmpTot[i])
                edep.append( tmpedep[i])

        Etrue_gam.append(np.mean(edep))
        tot_gam.append(np.std(Tot))
        sct_gam.append(np.std(Sct))
        cer_gam.append(np.std(Cer))
        npe_gam.append(np.mean(Tot))
        ncer_gam.append(np.mean(Cer))

    Etrue_gam = np.array(Etrue_gam)
    tot_gam   = np.array(tot_gam)
    sct_gam   = np.array(sct_gam)
    cer_gam   = np.array(cer_gam)
    npe_gam   = np.array(npe_gam)
    ncer_gam  = np.array(ncer_gam)


    global_es = 3134.078/2.223

    fig, ax0 = plt.subplots()

    ax0.plot(Etrue_elec, ncer_elec/npe_elec, "^-", lw=2, ms=8, color="coral", label=r"$e^-$")
    ax0.plot(Etrue_posi, ncer_posi/npe_posi, "^-.", lw=2, ms=8, color="blue", label=r"$e^+$")
    ax0.plot(Etrue_gam, ncer_gam/npe_gam, "^--", lw=2, ms=8, color="seagreen", label=r"Single $\gamma$")


    ######################## Fitting
    kB, Asct, kC = 6.26e-3, 1408, 0.996

    # electron
    elecFitNPE, elecFitCNPE = [], []
    for i in npe_elec:
        E = i/global_es
        elecFitNPE.append(el.getCerNPE(E)*kC + el.getQPE(E, kB, Asct))
        print(i, elecFitNPE[-1])
        elecFitCNPE.append(el.getCerNPE(E)*kC)

    elecFitNPE = np.array(elecFitNPE)
    elecFitCNPE = np.array(elecFitCNPE)
    ax0.plot(Etrue_elec, elecFitCNPE/elecFitNPE, "o", lw=2, mfc="w", ms=6, color="coral")



    # gamma
    name = ["Cs137", "Mn54", "K40", "nH",  "AmBe", "AmC"]
    Etrue = [0.662, 0.835, 1.461, 2.223, 4.43, 6.13]
    
    gamFitNPE, gamFitCNPE = [], []

    sGamArr = []
    Cs137 = singleGamma("Cs137"); sGamArr.append(Cs137)
    Mn54  = singleGamma("Mn54") ; sGamArr.append(Mn54)
    K40   = singleGamma("K40")  ; sGamArr.append(K40)
    nH    = singleGamma("nH")   ; sGamArr.append(nH)
    AmBe  = singleGamma("AmBe") ; sGamArr.append(AmBe)
    AmC   = singleGamma("AmC")  ; sGamArr.append(AmC)
    
    sigmaTot2, sigmaPart12, sigmaPart22 = [], [], []
    
    for i in sGamArr:
        i.ModelPrediction()
        #sigmaTot2.append(i.getSPE()**2)
        #sigmaPart12.append(i.getSigmaPart1())
        #sigmaPart22.append(i.getSigmaPart2())
    
        #print(sigmaTot2[-1], sigmaPart12[-1], sigmaPart22[-1])
        gamFitNPE.append(i.getNPE())
        gamFitCNPE.append(i.getCNPE())
        print(i.getName(), gamFitCNPE[-1], gamFitNPE[-1])
   
    gamFitNPE = np.array(gamFitNPE) 
    gamFitCNPE = np.array(gamFitCNPE)

    ax0.plot(Etrue_gam, gamFitCNPE/gamFitNPE, "o", mfc="w", lw=2, ms=6, color="seagreen")


    # positron
    posiFitNPE, posiFitCNPE = [], []
    for E in Etrue_posi:
        E -= 2*0.511
        posiFitNPE.append(el.getCerNPE(E)*kC + el.getQPE(E, kB, Asct) + 2*660.8)
        print(i, posiFitNPE[-1])
        posiFitCNPE.append(el.getCerNPE(E)*kC + 2.49)

    posiFitNPE = np.array(posiFitNPE)
    posiFitCNPE = np.array(posiFitCNPE)
    ax0.plot(Etrue_posi, posiFitCNPE/posiFitNPE, "o", lw=2, mfc="w", ms=6, color="blue")



    ax0.tick_params(axis='both', which='major', labelsize=13)
    ax0.set_xlabel(r"$E_{dep}$[MeV]", fontsize=16)
    ax0.set_ylabel(r"$E_{Cer}/E_{dep}$", fontsize=16)
    #axmin = npe_gam.min()-100
    #axmax = npe_posi.max()+100
    #ax0.set_xlim(axmin, axmax)
    ax0.grid(True)
    ax0.set_xlim(0, 9)
    ax0.set_ylim(0., 0.07)

    p1, = ax0.plot(12, 1, "o", mfc="w", ms=6, color="dimgray")
    p2, = ax0.plot(12, 1, "^", ms=8, color="dimgray")
    
    legend1 = ax0.legend([p1, p2], ["Model", "Simulation truth"], loc="center right", prop={"size":14},  ncol=1)
    ax0.add_artist(legend1)
    ax0.legend(loc="lower right", prop={"size" : 14})

    plt.tight_layout()

    plt.savefig("CerRatioThreePar.pdf")
    #plt.show()

if __name__ == "__main__" :
    main()




