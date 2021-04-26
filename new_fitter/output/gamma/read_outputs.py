import numpy as np
import ROOT
import matplotlib.pyplot as plt
from scipy.stats import norm

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
Etrue = [0.662, 0.835, 1.022, 1.461, 2.223, 2.506, 4.43, 4.94, 6.13]
Etrue_nonl = [0.511, 0.662, 0.835, 1.253, 1.461, 2.223, 4.43, 4.94, 6.13]

def main():

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

    nonlData = np.array(nonlData)
    nonlDataErr = np.array(nonlDataErr)
    nonlCalc = np.array(nonlCalc)

    relnonloDiff = (nonlCalc - nonlData) / nonlData
    relnonlDiffErr = nonlDataErr * nonlCalc / nonlData/ nonlData
    
    plt.figure(0, figsize=(6, 4))
    plt.plot(Etrue_nonl, nonlCalc, "o-", label="Fitting")
    plt.errorbar(Etrue_nonl, nonlData, yerr=nonlDataErr, label="Simulation")
    plt.legend()
    plt.xlabel("Etrue/MeV")
    plt.ylabel(r"$E_{vis}/E_{true}$")
    plt.tight_layout()
    plt.grid(True)
    plt.savefig("nonlFit.pdf")

    plt.figure(1, figsize=(6, 3))
    plt.errorbar(Etrue_nonl, relnonloDiff, yerr=relnonlDiffErr, fmt="o", color="peru")
    plt.fill_between([0, 6.3], [-0.001, -0.001], [0.001, 0.001], color="royalblue", alpha=0.3)
    plt.hlines(0, 0, 6.3, linestyle="--", color="red")
    plt.xlabel("Etrue/MeV")
    plt.ylabel("relative bias")
    plt.ylim(-0.003, 0.003)
    plt.tight_layout()
    plt.grid(True)
    plt.savefig("nonlBias.pdf")



    # resolution:
    resData, resCalc, resDataErr = [], [], []
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
    
    resData = np.array(resData)
    resDataErr = np.array(resDataErr)
    resCalc = np.array(resCalc)

    relresoDiff = (resCalc - resData) / resData
    relresDiffErr = resDataErr * resCalc / resData/ resData


    plt.figure(2, figsize=(6, 4))
    plt.plot(Etrue, resCalc, "o-", label="Fitting")
    plt.errorbar(Etrue, resData, yerr=resDataErr, label="Simulation")
    plt.legend()
    plt.xlabel("Etrue/MeV")
    plt.ylabel("p.e. resolution")
    plt.tight_layout()
    plt.grid(True)
    plt.savefig("resFit.pdf")

    plt.figure(3, figsize=(6, 3))
    plt.errorbar(Etrue, relresoDiff, yerr=relresDiffErr, fmt="o", color="peru")
    plt.xlabel("Etrue/MeV")
    plt.ylabel("relative bias")
    plt.fill_between([0, 6.3], [-0.02, -0.02], [0.02, 0.02], color="royalblue", alpha=0.3)
    plt.hlines(0, 0, 6.3, linestyle="--", color="red")
    plt.tight_layout()
    plt.ylim(-0.04, 0.04)
    plt.grid(True)
    plt.savefig("resBias.pdf")

    plt.show()

    '''
    plt.figure(0)
    plt.errorbar(Cs137pe, Cs137entry, yerr=np.sqrt(Cs137entry), fmt="o",  \
                 ms=3.0, lw=1.5, color="indianred", label="Simulation")
    Cs137X = np.arange(np.min(Cs137pe), np.max(Cs137pe), 1)
    Cs137pdf = norm.pdf(Cs137X, loc=Cs137mu, scale=Cs137sigma)
    plt.plot(Cs137X, Cs137pdf*Cs137A*np.sqrt(2*np.pi)*Cs137sigma, "-", label="Fitting")
    plt.legend()
    plt.xlabel("# P.E.")
    plt.tight_layout()
    plt.savefig('Cs137_spec.pdf')
    
    plt.figure(1)
    plt.errorbar(Mn54pe, Mn54entry, yerr=np.sqrt(Mn54entry), fmt="o",  \
                 ms=3.0, lw=1.5, color="indianred", label="Simulation")
    Mn54X = np.arange(np.min(Mn54pe), np.max(Mn54pe), 1)
    Mn54pdf = norm.pdf(Mn54X, loc=Mn54mu, scale=Mn54sigma)
    plt.plot(Mn54X, Mn54pdf*Mn54A*np.sqrt(2*np.pi)*Mn54sigma, "-", label="Fitting")
    plt.legend()
    plt.xlabel("# P.E.")
    plt.tight_layout()
    plt.savefig('Mn54_spec.pdf')

    plt.figure(2)
    plt.errorbar(Ge68pe, Ge68entry, yerr=np.sqrt(Ge68entry), fmt="o",  \
                 ms=3.0, lw=1.5, color="indianred", label="Simulation")
    Ge68X = np.arange(np.min(Ge68pe), np.max(Ge68pe), 1)
    Ge68pdf = norm.pdf(Ge68X, loc=Ge68mu, scale=Ge68sigma)
    plt.plot(Ge68X, Ge68pdf*Ge68A*np.sqrt(2*np.pi)*Ge68sigma, "-", label="Fitting")
    plt.legend()
    plt.xlabel("# P.E.")
    plt.tight_layout()
    plt.savefig('Ge68_spec.pdf')

    plt.figure(3)
    plt.errorbar(K40pe, K40entry, yerr=np.sqrt(K40entry), fmt="o",  \
                 ms=3.0, lw=1.5, color="indianred", label="Simulation")
    K40X = np.arange(np.min(K40pe), np.max(K40pe), 1)
    K40pdf = norm.pdf(K40X, loc=K40mu, scale=K40sigma)
    plt.plot(K40X, K40pdf*K40A*np.sqrt(2*np.pi)*K40sigma, "-", label="Fitting")
    plt.legend()
    plt.xlabel("# P.E.")
    plt.tight_layout()
    plt.savefig('K40_spec.pdf')

    plt.figure(4)
    plt.errorbar(nHpe, nHentry, yerr=np.sqrt(nHentry), fmt="o",  \
                 ms=3.0, lw=1.5, color="indianred", label="Simulation")
    nHX = np.arange(np.min(nHpe), np.max(nHpe), 1)
    nHpdf = norm.pdf(nHX, loc=nHmu, scale=nHsigma)
    plt.plot(nHX, nHpdf*nHA*np.sqrt(2*np.pi)*nHsigma, "-", label="Fitting")
    plt.legend()
    plt.xlabel("# P.E.")
    plt.tight_layout()
    plt.savefig('nH_spec.pdf')

    plt.figure(5)
    plt.errorbar(Co60pe, Co60entry, yerr=np.sqrt(Co60entry), fmt="o",  \
                 ms=3.0, lw=1.5, color="indianred", label="Simulation")
    Co60X = np.arange(np.min(Co60pe), np.max(Co60pe), 1)
    Co60pdf = norm.pdf(Co60X, loc=Co60mu, scale=Co60sigma)
    plt.plot(Co60X, Co60pdf*Co60A*np.sqrt(2*np.pi)*Co60sigma, "-", label="Fitting")
    plt.legend()
    plt.xlabel("# P.E.")
    plt.tight_layout()
    plt.savefig('Co60_spec.pdf')

    plt.figure(6)
    plt.errorbar(AmBepe, AmBeentry, yerr=np.sqrt(AmBeentry), fmt="o",  \
                 ms=3.0, lw=1.5, color="indianred", label="Simulation")
    AmBeX = np.arange(np.min(AmBepe), np.max(AmBepe), 1)
    AmBepdf = norm.pdf(AmBeX, loc=AmBemu, scale=AmBesigma)
    plt.plot(AmBeX, AmBepdf*AmBeA*np.sqrt(2*np.pi)*AmBesigma, "-", label="Fitting")
    plt.legend()
    plt.xlabel("# P.E.")
    plt.tight_layout()
    plt.savefig('AmBe_spec.pdf')

    plt.figure(7)
    plt.errorbar(nC12pe, nC12entry, yerr=np.sqrt(nC12entry), fmt="o",  \
                 ms=3.0, lw=1.5, color="indianred", label="Simulation")
    nC12X = np.arange(np.min(nC12pe), np.max(nC12pe), 1)
    nC12pdf = norm.pdf(nC12X, loc=nC12mu, scale=nC12sigma)
    plt.plot(nC12X, nC12pdf*nC12A*np.sqrt(2*np.pi)*nC12sigma, "-", label="Fitting")
    plt.legend()
    plt.xlabel("# P.E.")
    plt.tight_layout()
    plt.savefig('nC12_spec.pdf')

    plt.figure(8)
    plt.errorbar(AmCpe, AmCentry, yerr=np.sqrt(AmCentry), fmt="o",  \
                 ms=3.0, lw=1.5, color="indianred", label="Simulation")
    AmCX = np.arange(np.min(AmCpe), np.max(AmCpe), 1)
    AmCpdf = norm.pdf(AmCX, loc=AmCmu, scale=AmCsigma)
    plt.plot(AmCX, AmCpdf*AmCA*np.sqrt(2*np.pi)*AmCsigma, "-", label="Fitting")
    plt.legend()
    plt.xlabel("# P.E.")
    plt.tight_layout()
    plt.savefig('AmC_spec.pdf')

    '''


if __name__ == "__main__":
    main()
