import numpy as np
import matplotlib.pyplot as plt
import uproot as up
from matplotlib import gridspec

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




if __name__ == "__main__" :

    #fig = plt.figure(figsize=(12, 5))
    #spec = gridspec.GridSpec(ncols=2, nrows=1 )

    #ax0 = fig.add_subplot(spec[0])
    #ax1 = fig.add_subplot(spec[1])

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

    ax0.plot(npe_elec, ncer_elec/npe_elec, "X-", lw=2, ms=8, color="coral", label=r"$e^-$")
    ax0.plot(npe_posi, ncer_posi/npe_posi, "d-.", lw=2, ms=8, color="blue", label=r"$e^+$")
    ax0.plot(npe_gam, ncer_gam/npe_gam, "o--", lw=2, ms=8, color="seagreen", label=r"Single $\gamma$")
    ax0.tick_params(axis='both', which='major', labelsize=13)
    ax0.set_xlabel(r"$N_{tot}$", fontsize=16)
    ax0.set_ylabel(r"$N_{Cer}/N_{tot}$", fontsize=16)
    axmin = npe_gam.min()-100
    axmax = npe_posi.max()+100
    ax0.set_xlim(axmin, axmax)
    ax0.grid(True)
    ax0.legend(loc="lower right", prop={"size" : 14})
    ax2 = ax0.twiny()
    ax2.plot(npe_elec/global_es, ncer_elec/npe_elec, "o", color="coral", ms=0.01, zorder=1)
    ax2.set_xlabel(r"$E_{vis}$/MeV", fontsize=16)
    ax2.set_xlim(axmin/global_es, axmax/global_es)
    ax0.legend(prop={'size' : 14})


    #ax1.bar(npe_elec, height=tot_elec**2, width=300, edgecolor="darkviolet", color="royalblue", label="cov")
    #ax1.bar(npe_elec, height=sct_elec**2, width=300, edgecolor="darkviolet", color="orange", label=r"$\sigma^2_{sct}$")
    #ax1.bar(npe_elec, height=cer_elec**2, width=300, edgecolor="darkviolet", bottom=sct_elec**2, color="lightseagreen", label=r"$\sigma^2_{Cer}$")
    #line_elec, = ax1.plot(npe_elec, tot_elec**2, "X-", ms=6, color="coral")


    #ax1.bar(npe_posi, height=tot_posi**2, width=300, edgecolor="darkviolet", color="royalblue")
    #ax1.bar(npe_posi, height=sct_posi**2, width=300, edgecolor="darkviolet", color="orange")
    #ax1.bar(npe_posi, height=cer_posi**2, width=300, edgecolor="darkviolet", bottom=sct_posi**2, color="lightseagreen")
    #line_posi, = ax1.plot(npe_posi, tot_posi**2, "d-", ms=6, color="blue")
    #
    ##ax1.bar(npe_gam, height=tot_gam**2, width=300, edgecolor="darkviolet", color="royalblue")
    ##ax1.bar(npe_gam, height=sct_gam**2, width=300, edgecolor="darkviolet", color="orange")
    ##ax1.bar(npe_gam, height=cer_gam**2, width=300, edgecolor="darkviolet", bottom=sct_gam**2, color="lightseagreen")
    #bar1 = ax1.bar(npe_gam, height=tot_gam**2, width=300, edgecolor="darkviolet", color="peru")
    #bar2 = ax1.bar(npe_gam, height=sigma2_gam, width=300, edgecolor="darkviolet", color="slateblue")
    #line_gam, = ax1.plot(npe_gam, tot_gam**2, "v-", ms=6, color="seagreen")
    #

    #ax1.tick_params(axis='both', which='major', labelsize=13)
    #ax1.set_xlabel(r"$N_{tot}$", fontsize=14)
    #ax1.set_ylabel(r"$\sigma^2$", fontsize=14)

    #legend1 = ax1.legend([bar1, bar2], [r"$\sigma_{N_{single}}$", r"$\overline{(\sigma^\gamma)^2}$"], loc="center left", prop={"size":14})
    #ax1.add_artist(legend1)
    #ax1.legend(loc="upper left", prop={"size" : 14})
    
    plt.tight_layout()

    plt.savefig("CerRatio.pdf")

    plt.show()









