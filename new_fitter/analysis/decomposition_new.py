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

    fig = plt.figure(figsize=(8, 12))
    spec = gridspec.GridSpec(ncols=1, nrows=3 )

    ax2 = fig.add_subplot(spec[2])
    ax0 = fig.add_subplot(spec[0], sharex = ax2)
    ax1 = fig.add_subplot(spec[1], sharex = ax2)

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

    #ax0.plot(npe_elec, ncer_elec/npe_elec, "X-", color="coral", label="Electron")
    #ax0.plot(npe_posi, ncer_posi/npe_posi, "d-.", color="blue", label="Positron")
    #ax0.plot(npe_gam, ncer_gam/npe_gam, "o--", color="seagreen", label="Gamma")
    #ax0.tick_params(axis='both', which='major', labelsize=13)
    #ax0.set_xlabel(r"$N_{tot}$", fontsize=14)
    #ax0.set_ylabel(r"$N_{Cer}/N_{tot}$", fontsize=14)
    #axmin = npe_gam.min()-100
    #axmax = npe_posi.max()+100
    #ax0.set_xlim(axmin, axmax)
    #ax0.grid(True)
    #ax0.legend(loc="lower right", prop={"size" : 14})
    #ax2 = ax0.twiny()
    #ax2.plot(npe_elec/global_es, ncer_elec/npe_elec, "o", color="coral", ms=0.01, zorder=1)
    #ax2.set_xlabel(r"$E_{vis}$[MeV]", fontsize=14)
    #ax2.set_xlim(axmin/global_es, axmax/global_es)
    #ax0.legend(prop={'size' : 14})


    ax1.bar(npe_elec/global_es, height=tot_elec**2/npe_elec**2, width=500/global_es, edgecolor="darkviolet", color="royalblue", label="cov")
    ax1.bar(npe_elec/global_es, height=sct_elec**2/npe_elec**2, width=500/global_es, edgecolor="darkviolet", color="orange", label=r"$\sigma^2_{s}$")
    ax1.bar(npe_elec/global_es, height=cer_elec**2/npe_elec**2, width=500/global_es, edgecolor="darkviolet", bottom=sct_elec**2/npe_elec**2, color="lightseagreen", label=r"$\sigma^2_{C}$")
    #ax1.bar(npe_elec/global_es, height=sct_elec**2/global_es**2, width=500/global_es, facecolor="white", lw=2, edgecolor="orange", label=r"$\sigma^2_{s}$")
    #ax1.bar(npe_elec/global_es, height=cer_elec**2/global_es**2, width=500/global_es, facecolor="white", lw=2, edgecolor="lightseagreen", bottom=sct_elec**2/global_es**2, label=r"$\sigma^2_{C}$")
    line_elec, = ax1.plot(npe_elec/global_es, tot_elec**2/npe_elec**2, "X-", ms=6, color="coral")
    ax1.legend(loc="upper right", prop={"size" : 14}, ncol=3)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #ax1.set_ylabel(r"$\sigma^2_E$", fontsize=14)
    ax1.set_ylabel(r"$(\sigma_E/E_{vis})^{2}$", fontsize=17)
    ax1.tick_params(axis='x', which='major', labelsize=0)
    ax1.tick_params(axis='y', which='major', labelsize=16)
    ax1.text(0.4, 0.7, "Electron", fontsize=17, horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
    ax1.set_title("(b)", fontsize=17)
    ax1.grid(True)


    ax2.bar(npe_posi/global_es, height=tot_posi**2/npe_posi**2, width=500/global_es, edgecolor="darkviolet", color="royalblue", label="cov")
    ax2.bar(npe_posi/global_es, height=sct_posi**2/npe_posi**2, width=500/global_es, edgecolor="darkviolet", color="orange", label=r"$\sigma^2_s$")
    ax2.bar(npe_posi/global_es, height=cer_posi**2/npe_posi**2, width=500/global_es, edgecolor="darkviolet", bottom=sct_posi**2/npe_posi**2, color="lightseagreen", label=r"$\sigma^2_C$")
    line_posi, = ax2.plot(npe_posi/global_es, tot_posi**2/npe_posi**2, "d-", ms=6, color="blue")
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax2.legend(loc="upper right", prop={"size" : 14}, ncol=3)
    ax2.grid(True)
    ax2.text(0.4, 0.7, "Positron", fontsize=17, horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
    
    #ax1.bar(npe_gam, height=tot_gam**2, width=300/global_es, edgecolor="darkviolet", color="royalblue")
    #ax1.bar(npe_gam, height=sct_gam**2, width=300/global_es, edgecolor="darkviolet", color="orange")
    #ax1.bar(npe_gam, height=cer_gam**2, width=300/global_es, edgecolor="darkviolet", bottom=sct_gam**2, color="lightseagreen")
    npe_gam = np.array(npe_gam)
    tot_gam = np.array(tot_gam)
    sigma2_gam = np.array(sigma2_gam)

    wid_gam = [200/global_es, 200/global_es, 500/global_es, 500/global_es, 500/global_es, 500/global_es]
    bar1 = ax0.bar(npe_gam/global_es, height=tot_gam**2/ (npe_gam)**2, width=wid_gam, edgecolor="darkviolet", color="peru")
    bar2 = ax0.bar(npe_gam/global_es, height=sigma2_gam/npe_gam**2, width=wid_gam, edgecolor="darkviolet", color="slateblue")
    #bar1 = ax0.bar(npe_gam/global_es, height=tot_gam**2/global_es**2, width=500/global_es, facecolor="white", linewidth=2, edgecolor="peru", hatch="|")
    #bar2 = ax0.bar(npe_gam/global_es, height=sigma2_gam/global_es**2, width=500/global_es, facecolor="white", linewidth=2, edgecolor="slateblue", hatch="/")
    line_gam, = ax0.plot(npe_gam/global_es, tot_gam**2/npe_gam**2, "v-", ms=6, color="seagreen")
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax0.set_ylabel(r"$(\sigma_E/E_{vis})^{2}$", fontsize=17)
    ax0.tick_params(axis='x', which='major', labelsize=0)
    ax0.tick_params(axis='y', which='major', labelsize=16)
    ax0.set_title("(a)", fontsize=17)
    ax0.grid(True)
    ax0.text(0.4, 0.7, "Gamma", fontsize=17, horizontalalignment='center', verticalalignment='center', transform=ax0.transAxes)

    

    ax2.tick_params(axis='both', which='major', labelsize=16)
    ax2.set_xlabel(r"$E_{vis}$[MeV]", fontsize=17)
    ax2.set_ylabel(r"$(\sigma_E/E_{vis})^{2}$", fontsize=17)
    ax2.set_title("(c)", fontsize=17)

    legend1 = ax0.legend([bar1, bar2], [r"$(\sigma^\gamma)^2_{mean}$", r"$(\sigma^\gamma)^2_{ave}$"], loc="upper right", prop={"size":14}, ncol=2)
    ax0.add_artist(legend1)
    
    #plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
    #            wspace=0.05, hspace=0.05)
    plt.tight_layout()

    plt.savefig("ResolDecomp_new1.pdf")

    plt.show()









