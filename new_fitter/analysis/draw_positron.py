import uproot as up
import matplotlib.pyplot as plt
import numpy as np
import ROOT
import elecLoader as el

def loadTotalPE(name):
    print("Loading source " + name )
    totpe = up.open(name)["evt"]["totalPE"]
    totpe = np.array(totpe)

    return totpe



def gauss(x, mu, sigma):
    return 1/(np.sqrt(2*np.pi) * sigma) * np.exp(-(x-mu)**2/2/sigma**2)


def PESigma(x, ra, rb, rc):
    sigma2 = ra+rb*x+rc*x**2
    if sigma2 < 0:
        return 0
    else:
        return np.sqrt(sigma2)

import ROOT
def gaussFit(arr, low, high):
    hist = ROOT.TH1D("hist", "", 200, low, high)
    for i in arr:
        hist.Fill(i)

    hist.Fit("gaus", "Q0")
    ff = hist.GetFunction("gaus")

    return ff.GetParameter(1),  ff.GetParError(1), ff.GetParameter(2), ff.GetParError(2)

    


kA, kB, kC, es, ra, rb, rc = 0.99879, 6.34823e-3, 0.984026, 1409.61, -11.0053, 1484.17, 121.885

gau_mu, gam_sigma = 1321.197, 38.484

def main():

    path = []
    path.append("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/positron/e+1MeV/")
    path.append("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/positron/e+2MeV/")
    path.append("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/positron/e+3MeV/")
    path.append("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/positron/e+5MeV/")
    path.append("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/positron/e+8MeV/")

    KE = [1, 2, 3, 5, 8]
    totE = [2.022, 3.022, 4.022, 6.022, 9.022]
    low  = [2000, 3500, 5000, 8000, 12000]
    high = [4000, 5500, 7000, 9000, 14000]

    gam_mu, gam_sigma = 1321.197, 38.484

    nonl_sim, nonlerr_sim, nonl_calc, nonlerr_calc, nonlmaxerr_calc = [], [], [], [[], []], []
    res_sim, reserr_sim, res_calc, reserr_calc, resmaxerr_calc     = [], [], [], [[], []], []
    nonldiff, nonldiff_err = [], []
    resdiff, resdiff_err = [], []

    nonlmin_calc = [1431.8381389757549+gam_mu, 2913.0388007609768+gam_mu, 4401.233268426745+gam_mu, 7371.814214576139+gam_mu, 11823.149593869492+gam_mu]
    nonlmax_calc = [1437.8121664943255+gam_mu, 2925.5246537477137+gam_mu, 4420.31365383517+gam_mu, 7404.033416866311+gam_mu, 11875.027720320764+gam_mu]
    nonlmin_calc = np.array(nonlmin_calc)
    nonlmax_calc = np.array(nonlmax_calc)

    resmin_calc = [1498.913928982478, 3235.138681521096, 5155.421562764596, 9495.534628055557, 16740.93191173247]
    resmax_calc = [1681.2479261177136, 3604.857079158652, 5855.433378370916, 11364.529454522733, 22447.7253381445]
    resmin_calc = np.array(resmin_calc)
    resmax_calc = np.array(resmax_calc)
    resmin_calc = np.sqrt(resmin_calc+gam_sigma**2)
    resmax_calc = np.sqrt(resmax_calc+gam_sigma**2)

    global_es = 3134.078/2.223

    for i in range(5):
    #for i in range(len(path)):

        filename = path[i] + "user.root"
        totpe = loadTotalPE(filename)

        npe, npe_err, sigma, sigma_err = gaussFit(totpe, low[i], high[i])
        #print(npe, npe_err, sigma ,sigma_err)
        nonl_sim.append(npe/(KE[i]+1.022)/global_es)
        nonlerr_sim.append(npe_err/(KE[i]+1.022)/global_es)
        res_sim.append(sigma/npe)
        reserr_sim.append(np.sqrt(sigma_err**2/npe**2 + npe_err**2*sigma**2/npe**4))
        #plt.hist(totpe, bins=100, density=True, label=str(KE[i])+" MeV")

        sctpe_elec = kA*el.getQPE(KE[i], kB, es)
        cerpe_elec = kC * el.getCerNPE(KE[i])
        pesigma_elec = PESigma(KE[i], ra, rb, rc)

        ntotpe = sctpe_elec + cerpe_elec + gam_mu
        totsigma = np.sqrt(pesigma_elec**2 + gam_sigma**2)
        nonl_calc.append(ntotpe/(KE[i]+1.022)/global_es)
        res_calc.append(totsigma/ntotpe)

        nonlerr_calc[0].append(nonl_calc[-1] - nonlmin_calc[i]/(KE[i]+1.022)/global_es)
        nonlerr_calc[1].append(nonlmax_calc[i]/(KE[i]+1.022)/global_es - nonl_calc[-1])
        if nonlerr_calc[0][-1] < nonlerr_calc[1][-1]:
            nonlmaxerr_calc.append(nonlerr_calc[1][-1])
        else:
            nonlmaxerr_calc.append(nonlerr_calc[0][-1])


        nonldiff.append((nonl_calc[-1] - nonl_sim[-1])/nonl_sim[-1])
        nonldiff_err.append( np.sqrt(nonlerr_sim[-1]**2 * nonl_calc[-1]**2 / nonl_sim[-1]**4 + nonlmaxerr_calc[-1]**2/nonl_sim[-1]**2 ) )

        reserr_calc[0].append( res_calc[-1] - resmin_calc[i]/ntotpe )
        reserr_calc[1].append( resmax_calc[i]/ntotpe - res_calc[-1] )
        if reserr_calc[0][-1] < reserr_calc[1][-1]:
            resmaxerr_calc.append(reserr_calc[1][-1])
        else:
            resmaxerr_calc.append(reserr_calc[0][-1])

        resdiff.append((res_calc[-1] - res_sim[-1])/res_sim[-1])
        resdiff_err.append( np.sqrt(reserr_sim[-1]**2 * res_calc[-1]**2 / res_sim[-1]**4 + resmaxerr_calc[-1]**2/res_sim[-1]**2 ) )

        #plotx = np.arange(totpe-10*totsigma, totpe+10*totsigma, 1)
        #ploty = []
        #for k in plotx:
        #    ploty.append(gauss(k, totpe, totsigma))

        #plt.plot(plotx, ploty, "-")

    plt.figure(0, figsize=(6, 4))
    plt.errorbar(totE, nonl_sim, yerr=nonlerr_sim, fmt="-", color='royalblue', label="Simulation")
    plt.errorbar(totE, nonl_calc, yerr=nonlerr_calc, fmt="o", ms=5, color='magenta', label="Calculation")
    plt.grid(True)
    plt.xlabel("Etrue/MeV")
    plt.ylabel("Nonlinearity")
    plt.legend()
    plt.savefig("Positron_nonl.pdf")

    plt.figure(1, figsize=(6, 3))
    plt.errorbar(totE, nonldiff, yerr=nonldiff_err, fmt="o", color="peru")
    plt.xlabel("Etrue/MeV")
    plt.ylabel("relative bias")
    plt.fill_between([0, 10], [-0.003, -0.003], [0.003, 0.003], color="royalblue", alpha=0.3)
    plt.hlines(0, 0, 10, linestyle="--", color="red")
    plt.tight_layout()
    plt.ylim(-0.01, 0.01)
    plt.grid(True)
    plt.savefig("Positron_nonlBias.pdf")

    plt.figure(2, figsize=(6, 4))
    plt.errorbar(totE, res_sim, yerr=reserr_sim, fmt="-", color="royalblue", label="Simulation")
    plt.errorbar(totE, res_calc, yerr=reserr_calc, fmt="o", ms=5, color='magenta', label="Calculation")
    plt.grid(True)
    plt.xlabel("Etrue/MeV")
    plt.ylabel("NPE Resolution")
    plt.legend()
    plt.savefig("Positron_res.pdf")

    plt.figure(3, figsize=(6, 3))
    plt.errorbar(totE, resdiff, yerr=resdiff_err, fmt="o", color="peru")
    plt.xlabel("Etrue/MeV")
    plt.ylabel("relative bias")
    plt.fill_between([0, 10], [-0.03, -0.03], [0.03, 0.03], color="royalblue", alpha=0.3)
    plt.hlines(0, 0, 10, linestyle="--", color="red")
    plt.tight_layout()
    plt.ylim(-0.1, 0.1)
    plt.grid(True)
    plt.savefig("Positron_resBias.pdf")

    #plt.xlabel("N.P.E")
    #plt.legend()
    #plt.title("Positron NPE distribution")
    #plt.savefig("PositronNPE.pdf")
    

    plt.show()


if __name__ == "__main__":
    main()


