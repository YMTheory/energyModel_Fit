import numpy as np
import ROOT
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.stats import norm
import uproot as up

def loadOutput(name):
    print("Loading source " + name )
    ff = ROOT.TFile("./"+name+"hist.root", "read")
    hh = ff.Get(name+"_data")
    f1 = ff.Get("func")

    totpe, totentry = [], []
    for i in range(hh.GetNbinsX()):
        totpe.append(hh.GetBinCenter(i))
        totentry.append(hh.GetBinContent(i))

    A  = f1.GetParameter(0)
    mu = f1.GetParameter(1)
    sg = f1.GetParameter(2)

    totpe = np.array(totpe)
    totentry = np.array(totentry)

    return totpe, totentry, A, mu, sg


def fitOutput(name):
    print("Loading source " + name )
    ff = ROOT.TFile("./"+name+"hist.root", "read")
    hh = ff.Get(name+"_data")
    hh.Fit("gaus", "Q")
    f2 = hh.GetFunction("gaus")
    
    return f2.GetParameter(1), f2.GetParError(1), f2.GetParameter(2), f2.GetParError(2)




name = ["Cs137", "Mn54", "Ge68", "K40", "nH", "Co60", "AmBe", "nC12", "AmC"]
name_nonl = ["Ge68", "Cs137", "Mn54", "Co60", "K40", "nH", "AmBe", "nC12", "AmC"]
Etrue = [0.662, 0.835, 1.022, 1.461, 2.223, 2.506, 4.43, 4.94, 6.13]
Etrue_nonl = [0.511, 0.662, 0.835, 1.253, 1.461, 2.223, 4.43, 4.94, 6.13]
nonlid1 = [1, 2, 4, 5, 6, 8]
nonlid2 = [0, 3, 7]
name1 = ["Cs137", "Mn54", "K40", "nH",  "AmBe", "AmC"]
Etrue1 = [0.662, 0.835, 1.461, 2.223, 4.43, 6.13]
Etrue_nonl1 = [ 0.662, 0.835, 1.461, 2.223, 4.43, 6.13]
name2 = ["Ge68", "Co60", "nC12"]
Etrue2 = [1.022, 2.506, 4.94]
Etrue_nonl2 = [0.511, 1.253, 4.94]
resid1 = [0, 1, 3, 4, 6, 8]
resid2 = [2, 5, 7]

def gammaSource():

    Cs137pe, Cs137entry, Cs137A, Cs137mu, Cs137sg = loadOutput("Cs137")
    Mn54pe, Mn54entry, Mn54A, Mn54mu, Mn54sg = loadOutput("Mn54")
    Ge68pe, Ge68entry, Ge68A, Ge68mu, Ge68sg = loadOutput("Ge68")
    K40pe, K40entry, K40A, K40mu, K40sg = loadOutput("K40")
    nHpe, nHentry, nHA, nHmu, nHsg = loadOutput("nH")
    Co60pe, Co60entry, Co60A, Co60mu, Co60sg = loadOutput("Co60")
    AmBepe, AmBeentry, AmBeA, AmBemu, AmBesg = loadOutput("AmBe")
    nC12pe, nC12entry, nC12A, nC12mu, nC12sg = loadOutput("nC12")
    AmCpe, AmCentry, AmCA, AmCmu, AmCsg = loadOutput("AmC")



    Cs137mean, Cs137meanerr, Cs137sigma, Cs137sigmaerr = fitOutput("Cs137")
    Mn54mean, Mn54meanerr, Mn54sigma, Mn54sigmaerr = fitOutput("Mn54")
    Ge68mean, Ge68meanerr, Ge68sigma, Ge68sigmaerr = fitOutput("Ge68")
    K40mean, K40meanerr, K40sigma, K40sigmaerr = fitOutput("K40")
    nHmean, nHmeanerr, nHsigma, nHsigmaerr = fitOutput("nH")
    Co60mean, Co60meanerr, Co60sigma, Co60sigmaerr = fitOutput("Co60")
    AmBemean, AmBemeanerr, AmBesigma, AmBesigmaerr = fitOutput("AmBe")
    nC12mean, nC12meanerr, nC12sigma, nC12sigmaerr = fitOutput("nC12")
    AmCmean, AmCmeanerr, AmCsigma, AmCsigmaerr = fitOutput("AmC")


    # nonlinearity:
    scale = nHmean/Etrue[4]
    nonlData, nonlCalc, nonlDataErr = [], [], []
    nonlData1, nonlCalc1, nonlDataErr1 = [], [], []
    nonlData2, nonlCalc2, nonlDataErr2 = [], [], []
    nonlData.append(Ge68mean/scale/Etrue[2])
    nonlCalc.append(Ge68mu/scale/Etrue[2])
    nonlDataErr.append(Ge68meanerr/scale/Etrue[2])
    nonlData.append(Cs137mean/scale/Etrue[0])
    nonlCalc.append(Cs137mu/scale/Etrue[0])
    nonlDataErr.append(Cs137meanerr/scale/Etrue[0])
    nonlData.append(Mn54mean/scale/Etrue[1])
    nonlCalc.append(Mn54mu/scale/Etrue[1])
    nonlDataErr.append(Mn54meanerr/scale/Etrue[1])
    nonlData.append(Co60mean/scale/Etrue[5])
    nonlCalc.append(Co60mu/scale/Etrue[5])
    nonlDataErr.append(Co60meanerr/scale/Etrue[5])
    nonlData.append(K40mean/scale/Etrue[3])
    nonlCalc.append(K40mu/scale/Etrue[3])
    nonlDataErr.append(K40meanerr/scale/Etrue[3])
    nonlData.append(nHmean/scale/Etrue[4])
    nonlCalc.append(nHmu/scale/Etrue[4])
    nonlDataErr.append(nHmeanerr/scale/Etrue[4])
    nonlData.append(AmBemean/scale/Etrue[6])
    nonlCalc.append(AmBemu/scale/Etrue[6])
    nonlDataErr.append(AmBemeanerr/scale/Etrue[6])
    nonlData.append(nC12mean/scale/Etrue[7])
    nonlCalc.append(nC12mu/scale/Etrue[7])
    nonlDataErr.append(nC12meanerr/scale/Etrue[7])
    nonlData.append(AmCmean/scale/Etrue[8])
    nonlCalc.append(AmCmu/scale/Etrue[8])
    nonlDataErr.append(AmCmeanerr/scale/Etrue[8])
    for i in nonlid1:
        nonlData1.append(nonlData[i])
        nonlCalc1.append(nonlCalc[i])
        nonlDataErr1.append(nonlDataErr[i])
    for i in nonlid2:
        nonlData2.append(nonlData[i])
        nonlCalc2.append(nonlCalc[i])
        nonlDataErr2.append(nonlDataErr[i])

    nonlData = np.array(nonlData)
    nonlDataErr = np.array(nonlDataErr)
    nonlCalc = np.array(nonlCalc)
    nonlData1 = np.array(nonlData1)
    nonlDataErr1 = np.array(nonlDataErr1)
    nonlCalc1 = np.array(nonlCalc1)
    nonlData2 = np.array(nonlData2)
    nonlDataErr2 = np.array(nonlDataErr2)
    nonlCalc2 = np.array(nonlCalc2)

    relnonloDiff = (nonlCalc - nonlData) / nonlData
    relnonlDiffErr = nonlDataErr * nonlCalc / nonlData/ nonlData
    relnonloDiff1 = (nonlCalc1 - nonlData1) / nonlData1
    relnonlDiffErr1 = nonlDataErr1 * nonlCalc1 / nonlData1 / nonlData1
    relnonloDiff2 = (nonlCalc2 - nonlData2) / nonlData2
    relnonlDiffErr2 = nonlDataErr2 * nonlCalc2 / nonlData2/ nonlData2

    
    #plt.figure(0, figsize=(6, 4))
    #plt.plot(Etrue_nonl, nonlCalc, "o-", label="Fitting")
    #plt.errorbar(Etrue_nonl, nonlData, yerr=nonlDataErr, label="Simulation")
    #plt.legend()
    #plt.xlabel("Etrue/MeV")
    #plt.ylabel(r"$E_{vis}/E_{true}$")
    #plt.tight_layout()
    #plt.grid(True)
    #plt.savefig("nonlFit.pdf")

    #plt.figure(1, figsize=(6, 3))
    #plt.errorbar(Etrue_nonl, relnonloDiff, yerr=relnonlDiffErr, fmt="o", color="peru")
    #plt.fill_between([0, 6.3], [-0.001, -0.001], [0.001, 0.001], color="royalblue", alpha=0.3)
    #plt.hlines(0, 0, 6.3, linestyle="--", color="red")
    #plt.xlabel("Etrue/MeV")
    #plt.ylabel("relative bias")
    #plt.ylim(-0.003, 0.003)
    #plt.tight_layout()
    #plt.grid(True)
    #plt.savefig("nonlBias.pdf")



    # resolution:
    resData, resCalc, resDataErr = [], [], []
    resData1, resCalc1, resDataErr1 = [], [], []
    resData2, resCalc2, resDataErr2 = [], [], []
    resData.append(Cs137sigma/Cs137mean)
    resCalc.append(Cs137sg/Cs137mu)
    resDataErr.append(np.sqrt(Cs137sigmaerr**2/Cs137mean**2 + Cs137meanerr**2*Cs137sigma**2/Cs137mean**4))
    resData.append(Mn54sigma/Mn54mean)
    resCalc.append(Mn54sg/Mn54mu)
    resDataErr.append(np.sqrt(Mn54sigmaerr**2/Mn54mean**2 + Mn54meanerr**2*Mn54sigma**2/Mn54mean**4))
    resData.append(Ge68sigma/Ge68mean)
    resCalc.append(Ge68sg/Ge68mu)
    resDataErr.append(np.sqrt(Ge68sigmaerr**2/Ge68mean**2 + Ge68meanerr**2*Ge68sigma**2/Ge68mean**4))
    resData.append(K40sigma/K40mean)
    resCalc.append(K40sg/K40mu)
    resDataErr.append(np.sqrt(K40sigmaerr**2/K40mean**2 + K40meanerr**2*K40sigma**2/K40mean**4))
    resData.append(nHsigma/nHmean)
    resCalc.append(nHsg/nHmu)
    resDataErr.append(np.sqrt(nHsigmaerr**2/nHmean**2 + nHmeanerr**2*nHsigma**2/nHmean**4))
    resData.append(Co60sigma/Co60mean)
    resCalc.append(Co60sg/Co60mu)
    resDataErr.append(np.sqrt(Co60sigmaerr**2/Co60mean**2 + Co60meanerr**2*Co60sigma**2/Co60mean**4))
    resData.append(AmBesigma/AmBemean)
    resCalc.append(AmBesg/AmBemu)
    resDataErr.append(np.sqrt(AmBesigmaerr**2/AmBemean**2 + AmBemeanerr**2*AmBesigma**2/AmBemean**4))
    resData.append(nC12sigma/nC12mean)
    resCalc.append(nC12sg/nC12mu)
    resDataErr.append(np.sqrt(nC12sigmaerr**2/nC12mean**2 + nC12meanerr**2*nC12sigma**2/nC12mean**4))
    resData.append(AmCsigma/AmCmean)
    resCalc.append(AmCsg/AmCmu)
    resDataErr.append(np.sqrt(AmCsigmaerr**2/AmCmean**2 + AmCmeanerr**2*AmCsigma**2/AmCmean**4))
    for i in resid1:
        resData1.append(resData[i])
        resCalc1.append(resCalc[i])
        resDataErr1.append(resDataErr[i])
    for i in resid2:
        resData2.append(resData[i])
        resCalc2.append(resCalc[i])
        resDataErr2.append(resDataErr[i])
    
    resData = np.array(resData)
    resDataErr = np.array(resDataErr)
    resCalc = np.array(resCalc)
    resData1 = np.array(resData1)
    resDataErr1 = np.array(resDataErr1)
    resCalc1 = np.array(resCalc1)
    resData2 = np.array(resData2)
    resDataErr2 = np.array(resDataErr2)
    resCalc2 = np.array(resCalc2)

    relresoDiff = (resCalc - resData) / resData
    relresDiffErr = resDataErr * resCalc / resData/ resData
    relresoDiff1 = (resCalc1 - resData1) / resData1
    relresDiffErr1 = resDataErr1 * resCalc1 / resData1 / resData1
    relresoDiff2 = (resCalc2 - resData2) / resData2
    relresDiffErr2 = resDataErr2 * resCalc2 / resData2/ resData2


    print("=====> Nonlinearity <=====")
    print(Etrue_nonl)
    print(nonlData)
    print(nonlDataErr)
    print(nonlCalc)
    print("=====> Resolution <=====")
    print(Etrue)
    print(resData)
    print(resDataErr)
    print(resCalc)

    #plt.figure(2, figsize=(6, 4))
    #plt.plot(Etrue, resCalc, "o-", label="Fitting")
    #plt.errorbar(Etrue, resData, yerr=resDataErr, label="Simulation")
    #plt.legend()
    #plt.xlabel("Etrue/MeV")
    #plt.ylabel("p.e. resolution")
    #plt.tight_layout()
    #plt.grid(True)
    #plt.savefig("resFit.pdf")

    #plt.figure(3, figsize=(6, 3))
    #plt.errorbar(Etrue, relresoDiff, yerr=relresDiffErr, fmt="o", color="peru")
    #plt.xlabel("Etrue/MeV")
    #plt.ylabel("relative bias")
    #plt.fill_between([0, 6.3], [-0.02, -0.02], [0.02, 0.02], color="royalblue", alpha=0.3)
    #plt.hlines(0, 0, 6.3, linestyle="--", color="red")
    #plt.tight_layout()
    #plt.ylim(-0.04, 0.04)
    #plt.grid(True)
    #plt.savefig("resBias.pdf")

    #plt.show()


    fig = plt.figure(figsize=(12, 6))
    spec = gridspec.GridSpec(ncols=2, nrows=2,
                         height_ratios=[1, 2])

    ax0 = fig.add_subplot(spec[0])
    ax1 = fig.add_subplot(spec[1])
    ax2 = fig.add_subplot(spec[2])
    ax3 = fig.add_subplot(spec[3])
    
    ax2.errorbar(Etrue_nonl, nonlData, yerr=nonlDataErr, color="blue", label="Simulation")
    ax2.plot(Etrue_nonl1, nonlCalc1, "o", color="magenta", label="Fitting: single gamma")
    ax2.plot(Etrue_nonl2, nonlCalc2, "d", color="magenta", label="Fitting: multiple gamma")
    for i in range(len(name_nonl)):
        if i == 6:
            ax2.text(Etrue_nonl[i]-0.3, nonlCalc[i]-0.01, name_nonl[i], color="blue")
        else:
            ax2.text(Etrue_nonl[i], nonlCalc[i]-0.01, name_nonl[i], color="blue")
    ax2.legend()
    ax2.set_xlabel(r"$E_{dep}$/MeV")
    ax2.set_ylabel(r"$E_{vis}/E_{dep}$")
    ax2.set_ylim(0.905, 1.04)
    ax2.set_xlim(0, 6.6)
    #ax1.tight_layout()
    ax2.grid(True)
    #ax1.savefig("nonlFit.pdf")

    #ax0.figure(1, figsize=(6, 3))
    ax0.errorbar(Etrue_nonl1, relnonloDiff1, yerr=relnonlDiffErr1, fmt="o", color="peru", label="single gamma")
    ax0.errorbar(Etrue_nonl2, relnonloDiff2, yerr=relnonlDiffErr2, fmt="d", color="peru", label="multiple gamma")
    ax0.fill_between([0, 6.3], [-0.001, -0.001], [0.001, 0.001], color="royalblue", alpha=0.3)
    ax0.text(2.5, 0.002, "0.1% uncertainty band", color="darkviolet")
    ax0.hlines(0, 0, 6.3, linestyle="--", color="red")
    #ax0.set_xlabel("Etrue/MeV")
    ax0.set_ylabel("relative bias")
    ax0.set_ylim(-0.003, 0.003)
    ax0.set_xlim(0, 6.6)
    #ax0.tight_layout()
    ax0.grid(True)


    ax3.errorbar(Etrue, resData, yerr=resDataErr, color="blue", label="Simulation")
    ax3.plot(Etrue1, resCalc1, "o", color="magenta", label="Fitting: single gamma")
    ax3.plot(Etrue2, resCalc2, "d", color="magenta", label="Fitting: multiple gamma")
    for i in range(len(name)):
        if i != 6:
            ax3.text(Etrue[i]-0.05, resCalc[i]+0.001, name[i], color="blue")
        else:
            ax3.text(Etrue[i]-0.30, resCalc[i]+0.001, name[i], color="blue")
    ax3.legend()
    ax3.set_xlim(0, 6.6)
    ax3.set_xlabel(r"$E_{dep}$/MeV")
    ax3.set_ylabel(r"$\sigma/N_{tot}$")
    ax3.set_ylim(0.012, 0.040)
    ax3.grid(True)

    ax1.errorbar(Etrue1, relresoDiff1, yerr=relresDiffErr1, fmt="o", color="peru", label="single gamma")
    ax1.errorbar(Etrue2, relresoDiff2, yerr=relresDiffErr2, fmt="d", color="peru", label="multiple gamma")
    ax1.set_ylabel("relative bias")
    ax1.fill_between([0, 6.3], [-0.02, -0.02], [0.02, 0.02], color="royalblue", alpha=0.3)
    ax1.text(2.5, 0.025, "2% uncertainty band", color="darkviolet")
    ax1.hlines(0, 0, 6.3, linestyle="--", color="red")
    ax1.set_ylim(-0.04, 0.04)
    ax1.set_xlim(0, 6.6)
    ax1.grid(True)



    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=0.02, hspace=0.02)
    plt.tight_layout()
    plt.show()









    """

    fig, ax = plt.subplots(2, 5, figsize=(14, 6))
    ax0 = ax[0, 0]
    ax1 = ax[0, 1]
    ax2 = ax[0, 2]
    ax3 = ax[0, 3]
    ax4 = ax[0, 4]
    ax5 = ax[1, 0]
    ax6 = ax[1, 1]
    ax7 = ax[1, 2]
    ax8 = ax[1, 3]
    ax9 = ax[1, 4]


    #plt.figure(0)
    ax0.errorbar(Cs137pe, Cs137entry, yerr=np.sqrt(Cs137entry), fmt="o",  \
                 ms=2.0, lw=1, color="indianred", label="Simulation")
    Cs137X = np.arange(np.min(Cs137pe), np.max(Cs137pe), 1)
    Cs137pdf = norm.pdf(Cs137X, loc=Cs137mu, scale=Cs137sg)
    ax0.plot(Cs137X, Cs137pdf*Cs137A*np.sqrt(2*np.pi)*Cs137sg, "-", label="Fitting")
    ax0.legend()
    ax0.set_xlabel("# P.E.")
    ax0.set_title("Cs137")
    #ax0.tight_layout()
    #plt.savefig('Cs137_spec.pdf')
    
    #plt.figure(1)
    ax1.errorbar(Mn54pe, Mn54entry, yerr=np.sqrt(Mn54entry), fmt="o",  \
                 ms=2.0, lw=1, color="indianred", label="Simulation")
    Mn54X = np.arange(np.min(Mn54pe), np.max(Mn54pe), 1)
    Mn54pdf = norm.pdf(Mn54X, loc=Mn54mu, scale=Mn54sg)
    ax1.plot(Mn54X, Mn54pdf*Mn54A*np.sqrt(2*np.pi)*Mn54sg, "-", label="Fitting")
    ax1.legend()
    ax1.set_xlabel("# P.E.")
    ax1.set_title("Mn54")
    #ax1.tight_layout()
    #plt.savefig('Mn54_spec.pdf')

    #plt.figure(2)
    ax2.errorbar(Ge68pe, Ge68entry, yerr=np.sqrt(Ge68entry), fmt="o",  \
                 ms=2.0, lw=1, color="indianred", label="Simulation")
    Ge68X = np.arange(np.min(Ge68pe), np.max(Ge68pe), 1)
    Ge68pdf = norm.pdf(Ge68X, loc=Ge68mu, scale=Ge68sg)
    ax2.plot(Ge68X, Ge68pdf*Ge68A*np.sqrt(2*np.pi)*Ge68sg, "-", label="Fitting")
    ax2.legend()
    ax2.set_xlabel("# P.E.")
    ax2.set_title("Ge68")
    #ax2.tight_layout()
    #ax2.savefig('Ge68_spec.pdf')

    #plt.figure(3)
    ax3.errorbar(K40pe, K40entry, yerr=np.sqrt(K40entry), fmt="o",  \
                 ms=2.0, lw=1, color="indianred", label="Simulation")
    K40X = np.arange(np.min(K40pe), np.max(K40pe), 1)
    K40pdf = norm.pdf(K40X, loc=K40mu, scale=K40sg)
    ax3.plot(K40X, K40pdf*K40A*np.sqrt(2*np.pi)*K40sg, "-", label="Fitting")
    ax3.legend()
    ax3.set_xlabel("# P.E.")
    ax3.set_title("K40")
    #ax3.tight_layout()
    #plt.savefig('K40_spec.pdf')

    #plt.figure(4)
    ax4.errorbar(nHpe, nHentry, yerr=np.sqrt(nHentry), fmt="o",  \
                 ms=2.0, lw=1, color="indianred", label="Simulation")
    nHX = np.arange(np.min(nHpe), np.max(nHpe), 1)
    nHpdf = norm.pdf(nHX, loc=nHmu, scale=nHsg)
    ax4.plot(nHX, nHpdf*nHA*np.sqrt(2*np.pi)*nHsg, "-", label="Fitting")
    ax4.legend()
    ax4.set_xlabel("# P.E.")
    ax4.set_title("nH")
    #ax4.tight_layout()
    #plt.savefig('nH_spec.pdf')

    #plt.figure(5)
    ax5.errorbar(Co60pe, Co60entry, yerr=np.sqrt(Co60entry), fmt="o",  \
                 ms=2.0, lw=1, color="indianred", label="Simulation")
    Co60X = np.arange(np.min(Co60pe), np.max(Co60pe), 1)
    Co60pdf = norm.pdf(Co60X, loc=Co60mu, scale=Co60sg)
    ax5.plot(Co60X, Co60pdf*Co60A*np.sqrt(2*np.pi)*Co60sg, "-", label="Fitting")
    ax5.legend()
    ax5.set_xlabel("# P.E.")
    ax5.set_title("Co60")
    #ax5.tight_layout()
    #plt.savefig('Co60_spec.pdf')

    #plt.figure(6)
    ax6.errorbar(AmBepe, AmBeentry, yerr=np.sqrt(AmBeentry), fmt="o",  \
                 ms=2.0, lw=1, color="indianred", label="Simulation")
    AmBeX = np.arange(np.min(AmBepe), np.max(AmBepe), 1)
    AmBepdf = norm.pdf(AmBeX, loc=AmBemu, scale=AmBesg)
    ax6.plot(AmBeX, AmBepdf*AmBeA*np.sqrt(2*np.pi)*AmBesg, "-", label="Fitting")
    ax6.legend()
    ax6.set_xlabel("# P.E.")
    ax6.set_title("AmBe")
    #ax6.tight_layout()
    #plt.savefig('AmBe_spec.pdf')

    #plt.figure(7)
    ax7.errorbar(nC12pe, nC12entry, yerr=np.sqrt(nC12entry), fmt="o",  \
                 ms=2.0, lw=1, color="indianred", label="Simulation")
    nC12X = np.arange(np.min(nC12pe), np.max(nC12pe), 1)
    nC12pdf = norm.pdf(nC12X, loc=nC12mu, scale=nC12sg)
    ax7.plot(nC12X, nC12pdf*nC12A*np.sqrt(2*np.pi)*nC12sg, "-", label="Fitting")
    ax7.legend()
    ax7.set_xlabel("# P.E.")
    ax7.set_title("nC12")
    #ax7.tight_layout()
    #plt.savefig('nC12_spec.pdf')

    #ax8.figure(8)
    ax8.errorbar(AmCpe, AmCentry, yerr=np.sqrt(AmCentry), fmt="o",  \
                 ms=2.0, lw=1, color="indianred", label="Simulation")
    AmCX = np.arange(np.min(AmCpe), np.max(AmCpe), 1)
    AmCpdf = norm.pdf(AmCX, loc=AmCmu, scale=AmCsg)
    ax8.plot(AmCX, AmCpdf*AmCA*np.sqrt(2*np.pi)*AmCsg, "-", label="Fitting")
    ax8.legend()
    ax8.set_xlabel("# P.E.")
    ax8.set_title("AmC")
    #ax8.tight_layout()
    #ax8.savefig('AmC_spec.pdf')
    
    ax9.get_xaxis().set_visible(False)
    ax9.get_yaxis().set_visible(False)
    ax9.spines['top'].set_visible(False)
    ax9.spines['right'].set_visible(False)
    ax9.spines['bottom'].set_visible(False)
    ax9.spines['left'].set_visible(False)


    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=0.1, hspace=0.1)

    plt.tight_layout()
    plt.savefig("AllGamFit.pdf")
    plt.show()



def readB12():
    scale = 1386.48

    f1 = ROOT.TFile("spectrum.root", "read")
    hData = f1.Get("hData")
    hTheo = f1.Get("hTheo")
    hRela = f1.Get("hRela")
    binCenter, binData, binTheo, binTheoErr, binRela, binDataErr = [], [], [], [], [], []
    for i in range(hData.GetNbinsX()):
        binCenter.append(hData.GetBinCenter(i))
        binData.append(hData.GetBinContent(i))
        binDataErr.append(hData.GetBinError(i))
        binTheo.append(hTheo.GetBinContent(i))
        binTheoErr.append(hTheo.GetBinError(i))
        binRela.append(hRela.GetBinContent(i))

    binCenter = np.array(binCenter) / 1386.48
    binData = np.array(binData)
    binDataErr = np.array(binDataErr)
    binTheo = np.array(binTheo)
    binTheoErr = np.array(binTheoErr)
    binRela = np.array(binRela)

    return binCenter, binData, binDataErr, binTheo, binTheoErr, binRela
    """



def B12():
    binCenter, binData, binDataErr, binTheo, binTheoErr, binRela = readB12()


    fig = plt.figure()
    spec = gridspec.GridSpec(ncols=1, nrows=2,
                         height_ratios=[1, 2])

    ax1 = fig.add_subplot(spec[0])
    ax2 = fig.add_subplot(spec[1])


    ax2.errorbar(binCenter, binData, yerr=binDataErr, fmt="o", ms=2, color="royalblue", label='simulation')
    ax2.errorbar(binCenter, binTheo, yerr=binTheoErr, fmt="o", ms=2, color="chocolate", label="fitting")
    ax2.legend()
    ax2.set_xlabel(r"$E_{vis}$")
    ax2.set_ylabel("a.u.")
    #plt.savefig("B12_spec.pdf")

    #plt.figure(1, figsize=(6, 3))
    binRelaErr = np.sqrt(binTheoErr**2/binData**2 + binDataErr**2*binTheo**2/binData**4)
    ax1.errorbar(binCenter, binRela, yerr=binRelaErr, fmt="o", ms=3, color="peru")
    ax1.grid(True)
    #ax1.set_xlabel(r"$E_{vis}$")
    ax1.set_ylabel("relative ratio")
    ax1.hlines(1, 0, 15, linestyle="--", color='red')
    ax1.set_xlim(3, 12)
    ax1.set_ylim(0.7, 1.3)

    plt.savefig("B12_spec_fit.pdf")
    plt.show()




def main():
    gammaSource()
    #B12()



if __name__ == "__main__":
    main()
