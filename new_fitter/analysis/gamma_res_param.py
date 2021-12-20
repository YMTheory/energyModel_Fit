import numpy as np
import matplotlib.pyplot as plt
import uproot as up
import ROOT


def loadOutput(name):
    print("Loading source " + name )
    ff = ROOT.TFile("../output/gamB12NewMicwhole/"+name+"hist.root", "read")
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
    ff = ROOT.TFile("../output/gamB12NewMicwhole/"+name+"hist.root", "read")
    hh = ff.Get(name+"_data")
    hh.Fit("gaus", "Q")
    f2 = hh.GetFunction("gaus")
    
    return f2.GetParameter(1), f2.GetParError(1), f2.GetParameter(2), f2.GetParError(2)




if __name__ == "__main__" :

    Y = 3134.078 / 2.223

    name = ["Cs137", "Mn54", "Ge68", "K40", "nH", "Co60", "AmBe", "nC12", "AmC"]

    name_nonl = ["Ge68", "Cs137", "Mn54", "Co60", "K40", "nH", "nC12", "AmBe",  "AmC"]
    Etrue = [0.662, 0.835, 1.022, 1.461, 2.223, 2.506, 4.43, 4.94, 6.13]
    Etrue_nonl = [0.511, 0.662, 0.835, 1.253, 1.461, 2.223, 4.15, 4.43, 6.13]
    nonlid1 = [1, 2, 4, 5, 7, 8]
    nonlid2 = [0, 3, 6]
    name1 = ["Cs137", "Mn54", "K40", "nH",  "AmBe", "AmC"]
    Etrue1 = [0.662, 0.835, 1.461, 2.223, 4.43, 6.13]
    Etrue_nonl1 = [ 0.662, 0.835, 1.461, 2.223, 4.43, 6.13]
    name2 = ["Ge68", "Co60", "nC12"]
    Etrue2 = [1.022, 2.506, 4.94]
    Etrue_nonl2 = [0.511, 1.253, 4.15]
    resid1 = [0, 1, 3, 4, 6, 8]
    resid2 = [2, 5, 7]


    #### Loading & Fitting ####
    #####################################################
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


    # resolution:
    Evis, Evis1, Evis2 = [], [], []
    resData, resCalc, resDataErr = [], [], []
    resData1, resCalc1, resDataErr1 = [], [], []
    resData2, resCalc2, resDataErr2 = [], [], []
    Evis.append(Cs137mean/Y)
    resData.append(Cs137sigma/Cs137mean)
    resCalc.append(Cs137sg/Cs137mu)
    resDataErr.append(np.sqrt(Cs137sigmaerr**2/Cs137mean**2 + Cs137meanerr**2*Cs137sigma**2/Cs137mean**4))
    Evis.append(Mn54mean/Y)
    resData.append(Mn54sigma/Mn54mean)
    resCalc.append(Mn54sg/Mn54mu)
    resDataErr.append(np.sqrt(Mn54sigmaerr**2/Mn54mean**2 + Mn54meanerr**2*Mn54sigma**2/Mn54mean**4))
    Evis.append(Ge68mean/Y)
    resData.append(Ge68sigma/Ge68mean)
    resCalc.append(Ge68sg/Ge68mu)
    resDataErr.append(np.sqrt(Ge68sigmaerr**2/Ge68mean**2 + Ge68meanerr**2*Ge68sigma**2/Ge68mean**4))
    Evis.append(K40mean/Y)
    resData.append(K40sigma/K40mean)
    resCalc.append(K40sg/K40mu)
    resDataErr.append(np.sqrt(K40sigmaerr**2/K40mean**2 + K40meanerr**2*K40sigma**2/K40mean**4))
    Evis.append(nHmean / Y)
    resData.append(nHsigma/nHmean)
    resCalc.append(nHsg/nHmu)
    resDataErr.append(np.sqrt(nHsigmaerr**2/nHmean**2 + nHmeanerr**2*nHsigma**2/nHmean**4))
    Evis.append(Co60mean / Y)
    resData.append(Co60sigma/Co60mean)
    resCalc.append(Co60sg/Co60mu)
    resDataErr.append(np.sqrt(Co60sigmaerr**2/Co60mean**2 + Co60meanerr**2*Co60sigma**2/Co60mean**4))
    Evis.append(AmBemean / Y)
    resData.append(AmBesigma/AmBemean)
    resCalc.append(AmBesg/AmBemu)
    resDataErr.append(np.sqrt(AmBesigmaerr**2/AmBemean**2 + AmBemeanerr**2*AmBesigma**2/AmBemean**4))
    Evis.append(nC12mean / Y)
    resData.append(nC12sigma/nC12mean)
    resCalc.append(nC12sg/nC12mu)
    resDataErr.append(np.sqrt(nC12sigmaerr**2/nC12mean**2 + nC12meanerr**2*nC12sigma**2/nC12mean**4))
    Evis.append(AmCmean / Y)
    resData.append(AmCsigma/AmCmean)
    resCalc.append(AmCsg/AmCmu)
    resDataErr.append(np.sqrt(AmCsigmaerr**2/AmCmean**2 + AmCmeanerr**2*AmCsigma**2/AmCmean**4))
    for i in resid1:
        Evis1.append(Evis[i])
        resData1.append(resData[i])
        resCalc1.append(resCalc[i])
        resDataErr1.append(resDataErr[i])
    for i in resid2:
        Evis2.append(Evis[i])
        resData2.append(resData[i])
        resCalc2.append(resCalc[i])
        resDataErr2.append(resDataErr[i])


    Etrue1.append(10)
    Etrue1.append(20)
    Etrue1.append(30)
    Evis1.append(14660.276247614132/Y)
    resData1.append(0.011165031552469784)
    resDataErr1.append(0.00011838695308845511)
    Evis1.append(29524.22795852505/Y)
    resData1.append(0.00837763619662347)
    resDataErr1.append(8.109382814377594e-05)
    Evis1.append(44370.491814067784/Y)
    resData1.append( 0.007270069440771564)
    resDataErr1.append(7.601395521851264e-05)


    
    Evis = np.array(Evis)
    Evis1 = np.array(Evis1)
    Evis2 = np.array(Evis2)
    resData = np.array(resData)
    resDataErr = np.array(resDataErr)
    resCalc = np.array(resCalc)
    resData1 = np.array(resData1)
    resDataErr1 = np.array(resDataErr1)
    resCalc1 = np.array(resCalc1)
    resData2 = np.array(resData2)
    resDataErr2 = np.array(resDataErr2)
    resCalc2 = np.array(resCalc2)

    #relresoDiff = (resCalc - resData) / resData
    #relresDiffErr = resDataErr * resCalc / resData/ resData
    #relresoDiff1 = (resCalc1 - resData1) / resData1
    #relresDiffErr1 = resDataErr1 * resCalc1 / resData1 / resData1
    #relresoDiff2 = (resCalc2 - resData2) / resData2
    #relresDiffErr2 = resDataErr2 * resCalc2 / resData2/ resData2





    g1 = ROOT.TGraphErrors()
    for i in range(len(Evis1)):
        g1.SetPoint(i, Evis1[i], resData1[i])
        g1.SetPointError(i, 0, resDataErr1[i])

    f1 = ROOT.TF1("f1", "sqrt([0]**2/x + [1]**2 + [2]**2/x**2)", 0, 8)

    g1.Fit(f1, "RE")


    dx = np.arange(0.1, 35, 0.1)
    dy = []
    for i in dx:
        dy.append(f1.Eval(i))



    ## New Production 
    Edep_new, Eave_new, Evis_new, resol_new, resolErr_new = [], [], [], [], []
    with open("../data/gamma/gammaResol.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Edep_new.append(float(data[0]))
            Eave_new.append(float(data[1]))
            Evis_new.append(float(data[2]) / Y)
            resol_new.append(float(data[6]))
            resolErr_new.append(float(data[7]))





    fig, ax = plt.subplots()
    #ax.errorbar(Etrue1, resData1, yerr=resDataErr1, fmt="o", color="blue", ms=5, label=r"single $\gamma$s Sim")
    ax.errorbar(Evis1, resData1, yerr=resDataErr1, fmt="o", color="blue", ms=5, label=r"single $\gamma$s Sim")

    ax.plot(dx, dy, "--", color="coral", label="ABC model")

    ax.errorbar(Evis_new, resol_new, yerr=resolErr_new, fmt="v", color="black", ms=5, label=r"single $\gamma$s Sim New")

    ax.legend(prop={"size":14})

    ax.set_xlabel(r"$\gamma$ $E^{vis}$ [MeV]", fontsize=15)
    ax.set_ylabel(r"$\sigma/E^{vis}$", fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=14)

    plt.tight_layout()
    plt.show()




















