import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerLine2D
plt.style.use("science")

from elecResponseLoader import elecResponseLoader
from outputHandler import outputHandler
from truthLoader import truthLoader


def fit_res(dirs, det, qmode, cmode, rmode, Nnonl, Nres, KE, par, kBflag=False, method="prop"):

    el = elecResponseLoader(qmode, cmode, rmode, det)
    el.IskBHigh(kBflag)
    err_method = method
    
    #### Nonlinearity #### 
    nonlfile = "/junofs/users/miaoyu/energy_model/fitter/energyModel_Fit/new_fitter/outputs/%s/nonlcov.txt"%dirs
    oh1 = outputHandler(nonlfile, Nnonl)
    nonlpar, nonlparerr = oh1.loadBestFit() 
    el.setYs(nonlpar[0])
    el.setkB(nonlpar[1]) 
    if el.getCerenkovMode() == "Sim":
        el.setkC(nonlpar[2])
    if el.getCerenkovMode() == "Ana1":
        el.setp0(nonlpar[2])
        #el.setp1(nonlpar[3])
        el.setp1(0)
        el.setp2(nonlpar[4])
        el.setE0(nonlpar[5])
    
    totpe_sim = []
    for i in KE:
        if par == "e-":
            totpe_sim.append(el.getScintillationNumber(i) + el.getCerenkovNumber(i))
        if par == "e+":
            totpe_sim.append(el.getScintillationNumber(i) + el.getCerenkovNumber(i) + el.getNPE_Ge68())
    totpe_sim = np.array(totpe_sim)
    
    ### calculate errorbar
    if err_method == "sampling" :
        dtotpe_sim = [-10000 for i in range(len(KE))]
        nonlpar_sigma = oh1.sample_corelation()
        for i in range(oh1.getSampleSize()):
            el.setYs(nonlpar[0] + nonlpar_sigma[0, i])
            el.setkB(nonlpar[1] + nonlpar_sigma[1, i]) 
            if el.getCerenkovMode() == "Sim":
                el.setkC(nonlpar[2] + nonlpar_sigma[2, i])
            if el.getCerenkovMode() == "Ana1":
                el.setp0(nonlpar[2] + nonlpar_sigma[2, i])
                el.setp1(0)
                #el.setp1(nonlpar[3] + nonlpar_sigma[3, i])
                el.setp2(nonlpar[4] + nonlpar_sigma[4, i])
                el.setE0(nonlpar[5] + nonlpar_sigma[5, i])
            tmp_totpe = []
            for j in KE:
                if par == "e-":
                    tmp_totpe.append(el.getScintillationNumber(j) + el.getCerenkovNumber(j))
                if par == "e+":
                    tmp_totpe.append(el.getScintillationNumber(j) + el.getCerenkovNumber(j) + el.getNPE_Ge68())

            for k in range(len(KE)):
                if dtotpe_sim[k] < abs(tmp_totpe[k] - totpe_sim[k]):
                    dtotpe_sim[k] = abs(tmp_totpe[k] - totpe_sim[k])
                    if dtotpe_sim[k] / totpe_sim[k] > 0.05:
                        dtotpe_sim[k] = 0
        dtotpe_sim = np.array(dtotpe_sim)

    if err_method == "prop":
        cov = oh1.loadCov()
        dtotpe_sim = []
        Ysct, kB, p0, p1, p2, E0 = nonlpar
        Ysct_err, kB_err, p0_err, p1_err, p2_err, E0_err = nonlparerr
        for i in KE:
            deltakB = 0.05e-3
            kBmin = kB - deltakB
            kBmax = kB + deltakB
            el.setkB(kBmin)
            Nsmin = el.getScintillationNumber(i)
            el.setkB(kBmax)
            Nsmax = el.getScintillationNumber(i)
            dfdkB   = (Nsmax - Nsmin) / (2*deltakB)
            el.setkB(kB)
            dfdYsct = el.getScintillationNumber(i) / Ysct
            E = i - E0
            dfdp0 = E**2/(-E0+np.exp(-p2*E)*p1+i)**2
            dfdp1 = -(np.exp(-p2*E)*p0 * E**2) / (-E0+np.exp(-p2*E)*p1+i)**2
            dfdp2 = -(np.exp(-p2*E)*p0*p1*(-E) * E**2) / (-E0+np.exp(-p2*E)*p1+i)**2
            dfdE0 = -(p0*(-1+np.exp(-p2*E)*p1*p2)*E**2)/(-E0+np.exp(-p2*E)*p1+i)**2 - (2*p0*E)/(-E0+np.exp(-p2*E)*p1+i)

            sigma_nonl = np.sqrt(Ysct_err**2*dfdYsct**2 + kB_err**2*dfdkB**2 + p0_err**2*dfdp0**2 + p1_err**2*dfdp1**2 + p2_err**2*dfdp2**2+E0_err**2*dfdE0**2
                       + 2* cov[0, 1] * dfdkB*dfdYsct + 2 * cov[0, 2] * dfdYsct*dfdp0 + 2 * cov[0, 3]*dfdYsct*dfdp1 + 2*cov[0,4]*dfdYsct*dfdp2 + 2*cov[0,5]*dfdYsct*dfdE0
                       + 2*cov[1, 2]*dfdkB*dfdp0 + 2*cov[1, 3]*dfdkB*dfdp1 + 2*cov[1,4]*dfdkB*dfdp2 + 2*cov[1, 5]*dfdkB*dfdE0
                       + 2*cov[2, 3]*dfdp0*dfdp1 + 2*cov[2,4]*dfdp0*dfdp2 + 2*cov[2, 5]*dfdp0*dfdE0
                       + 2*cov[3,4]*dfdp1*dfdp2 + 2*cov[3, 5]*dfdp1*dfdE0
                       +2*cov[4, 5]*dfdp2*dfdE0)
            dtotpe_sim.append(sigma_nonl)
        dtotpe_sim = np.array(dtotpe_sim)


    #### Resolution ####
    resfile = "/junofs/users/miaoyu/energy_model/fitter/energyModel_Fit/new_fitter/outputs/%s/rescov.txt"%dirs
    oh2 = outputHandler(resfile, Nres) 
    respar, resparerr = oh2.loadBestFit() 
    if rmode == "New":
        el.seta(respar[0])
        el.setb(respar[1])
        el.setn(respar[2])
    elif rmode == "New1":
        el.seta1(respar[0])
        el.setn1(respar[1])
    Evis_sigma = []
    for i in KE:
        Evis = (el.getScintillationNumber(i) + el.getCerenkovNumber(i) ) / el.getY()
        Esigma = el.getEvisSigma(Evis)
        if par == "e-":
            Evis_sigma.append(Esigma)
        if par == "e+":
            Evis_sigma.append(np.sqrt(Esigma**2 + (el.getSigma_Ge68()/el.getY())**2))
    Evis_sigma = np.array(Evis_sigma)

    ### calculate errorbar
    if err_method == "sampling":
        dEvissigma_sim = [-10000 for i in range(len(KE))]
        respar_sigma = oh2.sample_corelation()
        for i in range(oh2.getSampleSize()):
            if rmode == "New":
                el.seta(respar[0] + respar_sigma[0, i])
                el.setb(respar[1] + respar_sigma[1, i])
                el.setn(respar[2] + respar_sigma[2, i])
            elif rmode == "New1":
                el.seta1(respar[0] + respar_sigma[0, i])
                el.setn1(respar[1] + respar_sigma[1, i])
            tmp_Evissigma = []
            for m in KE:
                j = (el.getScintillationNumber(m)+el.getCerenkovNumber(m) ) / el.getY()
                Esigma = el.getEvisSigma(j)
                if par == "e-":
                    tmp_Evissigma.append(Esigma)
                if par == "e+":
                    tmp_Evissigma.append(np.sqrt(Esigma**2 + (el.getSigma_Ge68()/el.getY())**2))
            for k in range(len(KE)):
                if dEvissigma_sim[k] < abs(tmp_Evissigma[k] - Evis_sigma[k]):
                    dEvissigma_sim[k] = abs(tmp_Evissigma[k] - Evis_sigma[k])

        dEvissigma_sim = np.array(dEvissigma_sim)

    if err_method == "prop":
        cov = oh2.loadCov()
        dEvissigma_sim = []
        Y = el.getY()
        m_npeGe68 = el.getNPE_Ge68()
        m_sigmaGe68 = el.getSigma_Ge68()
        for j in KE:
            i = (el.getScintillationNumber(i) + el.getCerenkovNumber(i) ) / Y
            a, b, n = respar
            aerr, berr, nerr = resparerr

            dfda = (a*i*Y) / ( (m_npeGe68 + i * Y) * np.sqrt(m_sigmaGe68**2 + a**2 * i * Y + b**2 * np.power(i*Y, n)) )
            dfdb = (b*np.power(i*Y, n)) / ( (m_npeGe68 + i * Y) * np.sqrt(m_sigmaGe68**2 + a**2 * i * Y + b**2 * np.power(i*Y, n)) )
            dfdn = (b**2*np.power(i*Y, n)*np.log(i*Y)) / ( (m_npeGe68 + i * Y) * np.sqrt(m_sigmaGe68**2 + a**2 * i * Y + b**2 * np.power(i*Y, n)) )
    
            sigma_res = np.sqrt( dfda**2 * aerr**2 + dfdb**2*berr**2 + dfdn**2*nerr**2 + 2*dfda*dfdb*cov[0, 1] + 2*dfda*dfdn*cov[0, 2] * 2*dfdb*dfdn*cov[1, 2] )
            dEvissigma_sim.append(sigma_res)
        dEvissigma_sim = np.array(dEvissigma_sim)

    return el.getY(), el.getNPE_Ge68(), el.getSigma_Ge68(), totpe_sim, dtotpe_sim, Evis_sigma, dEvissigma_sim



def get_one_fit(dirs, det, particle, KE, KE1, kBflag, qmode, cmode, rmode, Nnonl, Nres):

    Y, m_npeGe68, m_sigmaGe68, totpe_sim, dtotpe_sim, Evis_sigma, dEvissigma_sim = fit_res(dirs,
                                                                                           det,
                                                                                           qmode,
                                                                                           cmode,
                                                                                           rmode,
                                                                                           Nnonl,
                                                                                           Nres,
                                                                                           KE, 
                                                                                           particle,
                                                                                           kBflag,
                                                                                           "sampling"
                                                                                           )

    Y1, m_npeGe68, m_sigmaGe68, totpe_sim1, dtotpe_sim1, Evis_sigma1, dEvissigma_sim1 = fit_res(dirs,
                                                                                                det,
                                                                                                qmode,
                                                                                                cmode,
                                                                                                rmode,
                                                                                                Nnonl,
                                                                                                Nres,
                                                                                                KE1, 
                                                                                                particle,
                                                                                                kBflag,
                                                                                                "sampling"
                                                                                                )




    #el = truthLoader("kB15.2e-3kC0.1")
    if det == "DYBnonl":
        det = "kB15.2e-3kC0.5"
    el = truthLoader(det)
    if particle == "e-":
        el.fitElecNPESimTruth(KE1)
    elif particle == "e+":
        el.fitPosiNPESimTruth(KE1)
    NPE_truth, NPEerr_truth, sigma_truth, sigmaerr_truth = [], [], [], []
    for ke in KE1:
        ke = str(ke) + "MeV"
        NPE_truth.append(el.getNPE(ke))
        NPEerr_truth.append(el.getNPEerr(ke))
        sigma_truth.append(el.getsigma(ke))
        sigmaerr_truth.append(el.getsigmaerr(ke))
    NPE_truth = np.array(NPE_truth)
    NPEerr_truth = np.array(NPEerr_truth)
    sigma_truth = np.array(sigma_truth)
    sigmaerr_truth = np.array(sigmaerr_truth)


    ### for plotting :
    Edep = KE
    Edep1 = KE1
    if particle == "e+":
        Edep = KE + 1.022
        Edep1 = KE1 + 1.022

    r_nonl = 10
    Evis = totpe_sim / Y
    nonl_sim = NPE_truth / Y / Edep1
    nonlerr_sim = NPEerr_truth / Y / Edep1
    nonl_fit = totpe_sim1 / Y / Edep1
    dnonl_fit_min = (totpe_sim1 - dtotpe_sim1) / Y / Edep1
    dnonl_fit_max = (totpe_sim1 + dtotpe_sim1) / Y / Edep1
    nonl_draw = totpe_sim / Y / Edep
    dnonl_draw_min = (totpe_sim - r_nonl*dtotpe_sim) / Y / Edep
    dnonl_draw_max = (totpe_sim + r_nonl*dtotpe_sim) / Y / Edep

    r_res = 10
    Evis1 = totpe_sim1 / Y
    res_sim = sigma_truth / NPE_truth
    reserr_sim = np.sqrt(sigmaerr_truth**2/NPE_truth**2 + sigma_truth**2*NPEerr_truth**2/NPE_truth**4)
    res_fit = Evis_sigma1 * Y / totpe_sim1
    dres_fit_min = (Evis_sigma1 - dEvissigma_sim1) * Y / totpe_sim1
    dres_fit_max = (Evis_sigma1 + dEvissigma_sim1) * Y / totpe_sim1
    res_draw = Evis_sigma * Y / totpe_sim
    dres_draw_min = (Evis_sigma - r_res*dEvissigma_sim) * Y / totpe_sim
    dres_draw_max = (Evis_sigma + r_res*dEvissigma_sim) * Y / totpe_sim

    return Edep, nonl_fit, dnonl_fit_min, dnonl_fit_max, Edep1, nonl_draw, dnonl_draw_min, dnonl_draw_max, nonl_sim, nonlerr_sim, Evis, res_fit, dres_fit_min, dres_fit_max, Evis1, res_draw, dres_draw_min, dres_draw_max, res_sim, reserr_sim



def compare_betas():
    #dirs        = "DYBnonl_kCer1p2limit_kQInt_kNew_onlygam"
    dirs        = "Det5_kCerAnaNew1_kQInt_p2limits"
    KE          = np.arange(0.1, 13, 0.1)
    KE1         = np.arange(0.5, 12.5, 0.5)
    #KE1         = np.append(KE1, 20)
    #KE1         = np.append(KE1, 30)
    #KE1         = np.append(KE1, 40)
    #KE1         = np.append(KE1, 50)
    #KE1         = np.append(KE1, 60)
    particle    = "e-"
    det         = "Det5"
    kBflag      = False
    qmode       = "Int"
    cmode       = "Ana1"
    rmode       = "New"
    Nnonl       = 6
    Nres        = 3
    Edep, nonl_fit, dnonl_fit_min, dnonl_fit_max, Edep1, nonl_draw, dnonl_draw_min, dnonl_draw_max, nonl_sim, nonlerr_sim, Evis, res_fit, dres_fit_min, dres_fit_max, Evis1, res_draw, dres_draw_min, dres_draw_max, res_sim, reserr_sim = get_one_fit(dirs, det, particle, KE, KE1, kBflag, qmode, cmode, rmode, Nnonl, Nres)

    #dirs        = "Det1_kCerAnaNew1_kQInt_p2limits_B12"
    #Edepm, nonl_fitm, dnonl_fit_minm, dnonl_fit_maxm, Edep1m, nonl_drawm, dnonl_draw_minm, dnonl_draw_maxm, nonl_simm, nonlerr_simm, Evism, res_fitm, dres_fit_minm, dres_fit_maxm, Evis1m, res_drawm, dres_draw_minm, dres_draw_maxm, res_simm, reserr_simm = get_one_fit(dirs, det, particle, KE, KE1, kBflag, qmode, cmode, rmode, Nnonl, Nres)

    flag = True
    if flag :
        fig = plt.figure(figsize=(12, 6))
        spec = gridspec.GridSpec(ncols=2, nrows=2,
                             height_ratios=[1, 2])

        ax1 = fig.add_subplot(spec[2])
        ax0 = fig.add_subplot(spec[0])
        ax2 = fig.add_subplot(spec[1])
        ax3 = fig.add_subplot(spec[3])

        l0, = ax0.plot(Edep1, 100*(nonl_fit - nonl_sim) / nonl_sim, "o-", color="blue")
        ax0.fill_between(Edep1, 100*(dnonl_fit_min-nonl_sim) / nonl_sim, 100*(dnonl_fit_max-nonl_sim)/nonl_sim, color="blue", alpha=0.3)
        pc0 = mpatches.Patch(facecolor="blue", alpha=0.3, linewidth=0)
        ##l0m, = ax0.plot(Edep1m, 100*(nonl_fitm - nonl_simm) / nonl_simm, "o-", color="green")
        ##ax0.fill_between(Edep1m, 100*(dnonl_fit_minm-nonl_simm) / nonl_simm, 100*(dnonl_fit_maxm-nonl_simm)/nonl_simm, color="green", alpha=0.3)
        ##pc0m = mpatches.Patch(facecolor="green", alpha=0.3, linewidth=0)
        #ax0.legend([(pc0, l0), (pc0m, l0m)], [r"1$\sigma$ errorbar w/o Michel $e^-$", r"1$\sigma$ errorbar w/ Michel $e^-$"], handler_map = {l0: HandlerLine2D(marker_pad = 0), l0m:HandlerLine2D(marker_pad = 0)}, loc="upper right", prop={"size":16})
        ax0.tick_params(axis='x', which='major', labelsize=0, labelcolor="black")
        ax0.tick_params(axis='y', which='major', labelsize=16, labelcolor="black")
        ax0.set_ylabel("Bias (\%)", fontsize=16)

        p1 = ax1.errorbar(Edep1, nonl_sim, yerr=nonlerr_sim, fmt="o", color="red", fillstyle="none")
        l1, = ax1.plot(Edep, nonl_draw, lw=2, color="blue")
        ax1.fill_between(Edep, dnonl_draw_min, dnonl_draw_max, color="blue", alpha=0.3)
        pc1 = mpatches.Patch(facecolor="blue", alpha=0.3, linewidth=0)
        #l1m, = ax1.plot(Edepm, nonl_drawm, "--", lw=2, color="green")
        ##ax1.fill_between(Edepm, dnonl_draw_minm, dnonl_draw_maxm, color="green", alpha=0.3)
        ##pc1m = mpatches.Patch(facecolor="green", alpha=0.3, linewidth=0)
        ax1.set_xlabel(r"E [MeV]", fontsize=16)
        ax1.set_ylabel(r"$E_\mathrm{vis}/E$", fontsize=16)
        #ax1.legend([(pc1, l1), p1, l1m], [r"$10\times 1\sigma$ errorbar ", "simulation", "nominal"], handler_map = {l1: HandlerLine2D(marker_pad = 0), l1m:HandlerLine2D(marker_pad=0)}, loc="lower right", prop={"size":16})
        ax1.legend([(pc1, l1), p1], [r"$10\times 1\sigma$ errorbar ", "simulation"], handler_map = {l1: HandlerLine2D(marker_pad = 0)}, loc="lower right", prop={"size":16})
        #ax1.legend([(pc1, l1), (pc1m, l1m), p1], [r"$10\times 1\sigma$ errorbar w/o Michel $e^-$", r"$10\times 1\sigma$ errorbar w/ Michel $e^-$","simulation"], handler_map = {l1: HandlerLine2D(marker_pad = 0), l1m: HandlerLine2D(marker_pad=0)}, loc="lower right", prop={"size":16})
        ax1.tick_params(axis='both', which='major', labelsize=16, labelcolor="black")

        l2, = ax2.plot(Evis1, 100*(res_fit - res_sim)/res_sim, "o-", color="blue")
        ax2.fill_between(Evis1, 100*(dres_fit_min-res_sim) / res_sim, 100*(dres_fit_max-res_sim)/res_sim, color="blue", alpha=0.3)
        pc2 = mpatches.Patch(facecolor="blue", alpha=0.3, linewidth=0)
        ##l2m, = ax2.plot(Evis1m, 100*(res_fitm - res_simm)/res_simm, "o-", color="green")
        ##ax2.fill_between(Evis1m, 100*(dres_fit_minm-res_simm) / res_simm, 100*(dres_fit_maxm-res_simm)/res_simm, color="green", alpha=0.3)
        ##pc2m = mpatches.Patch(facecolor="blue", alpha=0.3, linewidth=0)
        ax2.legend([(pc2, l2)], [r"1$\sigma$ errorbar"], handler_map = {l0: HandlerLine2D(marker_pad = 0)}, loc="upper center", prop={"size":16})
        ##ax2.legend([(pc2, l2), (pc2m, l2m)], [r"1$\sigma$ errorbar w/o Michel $e^-$", r"1$\sigma$ errorbar w/ Michel $e^-$"], handler_map = {l0: HandlerLine2D(marker_pad = 0), l0m:HandlerLine2D(marker_pad = 0)}, loc="upper left", prop={"size":16})
        ax2.tick_params(axis='x', which='major', labelsize=0, labelcolor="black")
        ax2.tick_params(axis='y', which='major', labelsize=16, labelcolor="black")
        ax2.set_ylabel("Bias (\%)", fontsize=16)

        p3 = ax3.errorbar(Evis1, res_sim, yerr=reserr_sim, fmt="o", color="red", fillstyle="none")
        l3, = ax3.plot(Evis, res_draw, lw=2, color="blue")
        pc3 = mpatches.Patch(facecolor="blue", alpha=0.3, linewidth=0)
        ax3.fill_between(Evis, dres_draw_min, dres_draw_max, color="blue", alpha=0.3)
        #l3m, = ax3.plot(Evism, res_drawm, "--", lw=2, color="green")
        ##pc3m = mpatches.Patch(facecolor="blue", alpha=0.3, linewidth=0)
        ##ax3.fill_between(Evism, dres_draw_minm, dres_draw_maxm, color="green", alpha=0.3)
        ax3.set_xlabel(r"E$_\mathrm{vis}$ [MeV]", fontsize=16)
        ax3.set_ylabel(r"$\sigma/E_\mathrm{vis}$", fontsize=16)
        #ax3.legend([(pc3, l3), p3], [r"$10\times 1\sigma$ errorbar", "simulation"], handler_map = {l1: HandlerLine2D(marker_pad = 0)}, loc="upper right", prop={"size":16})
        ax3.tick_params(axis='both', which='major', labelsize=16, labelcolor="black")

        plt.tight_layout()
        plt.savefig("./results/Det5_kQInt_kCAna1_kNew_%s.pdf"%particle)
        #plt.savefig("./results/Det1_kQInt_kCAna1_kNew_gammaB12Mic_%s.pdf"%particle)
        #plt.show()

def compare_particles():
    KE = np.arange(0.5, 12.3, 0.1)
    Edep_elec = KE

    par = "e-"
    Y, m_npeGe68, m_sigmaGe68, totpe_sim_elec, dtotpe_sim_elec, Evis_sigma_elec, dEvissigma_sim_elec = fit_res("Det1_kCerAnaNew1_kQInt_p2limits"    ,   ## output directory
                                                                                                               "Det1",                                  ## detconfig
                                                                                                               "Int",                                   ## quench_mode 
                                                                                                               "Ana1",                                  ## cerenkov_mode
                                                                                                               "New1",                                  ## resolution_mode
                                                                                                               6,                                       ## number of nonlinearity parameters
                                                                                                               3,                                       ## number of resolution parameters
                                                                                                               KE,                                      ## KE array
                                                                                                               par)                                     ## particle definition

    KE = np.arange(0.1, 12.3, 0.1)
    Edep_posi = KE + 1.022
    par = "e+"
    Y, m_npeGe68, m_sigmaGe68, totpe_sim_posi, dtotpe_sim_posi, Evis_sigma_posi, dEvissigma_sim_posi = fit_res("Det1_kCerAnaNew1_kQInt_p2limits"    ,   ## output directory
                                                                                                               "Det1",                                  ## detconfig
                                                                                                               "Int",                                   ## quench_mode 
                                                                                                               "Ana1",                                  ## cerenkov_mode
                                                                                                               "New1",                                  ## resolution_mode
                                                                                                               6,                                       ## number of nonlinearity parameters
                                                                                                               3,                                       ## number of resolution parameters
                                                                                                               KE,                                      ## KE array
                                                                                                               par)                                     ## particle definition

    el = truthLoader()
    el.fitGammaNPESimTruth()
    NPE_gamma, NPEerr_gamma, sigma_gamma, sigmaerr_gamma = [], [], [], []
    gams = ["Cs137", "Mn54", "K40", "nH", "AmBe", "AmC"]
    Edep_gamma = [0.662, 0.834, 1.461, 2.223, 4.430, 6.130]
    for ke in gams:
        NPE_gamma.append(el.getNPE(ke))
        NPEerr_gamma.append(el.getNPEerr(ke))
        sigma_gamma.append(el.getsigma(ke))
        sigmaerr_gamma.append(el.getsigmaerr(ke))
    NPE_gamma = np.array(NPE_gamma)
    NPEerr_gamma = np.array(NPEerr_gamma)
    sigma_gamma = np.array(sigma_gamma)
    sigmaerr_gamma = np.array(sigmaerr_gamma)

    ## for plotting
    ### nonlinearity
    r = 5
    elec_nonl = totpe_sim_elec / Y / Edep_elec
    elec_nonl_min = (totpe_sim_elec - r * dtotpe_sim_elec) / Y /  Edep_elec
    elec_nonl_max = (totpe_sim_elec + r * dtotpe_sim_elec) / Y /  Edep_elec
    posi_nonl = totpe_sim_posi / Y / Edep_posi
    posi_nonl_min = (totpe_sim_posi - r * dtotpe_sim_posi) / Y /  Edep_posi
    posi_nonl_max = (totpe_sim_posi + r * dtotpe_sim_posi) / Y /  Edep_posi
    gam_nonl = NPE_gamma / Y / Edep_gamma
    gam_nonl_err = NPEerr_gamma / Y / Edep_gamma
    ### resolution
    Evis_elec = totpe_sim_elec / Y
    elec_res = Evis_sigma_elec * Y / totpe_sim_elec
    elec_res_min = (Evis_sigma_elec - r * dEvissigma_sim_elec) * Y / totpe_sim_elec
    elec_res_max = (Evis_sigma_elec + r * dEvissigma_sim_elec) * Y / totpe_sim_elec
    Evis_posi = totpe_sim_posi / Y
    posi_res = Evis_sigma_posi * Y / totpe_sim_posi
    posi_res_min = (Evis_sigma_posi - r * dEvissigma_sim_posi) * Y / totpe_sim_posi
    posi_res_max = (Evis_sigma_posi + r * dEvissigma_sim_posi) * Y / totpe_sim_posi
    Evis_gamma = NPE_gamma / Y
    gam_res = sigma_gamma / NPE_gamma
    gam_res_err = np.sqrt(sigmaerr_gamma**2/NPE_gamma**2 + sigma_gamma**2*NPEerr_gamma**2/NPE_gamma**4)


    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    lposi, = ax1.plot(Edep_posi, posi_nonl, "-.", lw=2.0, color="red", label=r"e$^+$")
    patchposi = mpatches.Patch(facecolor='crimson', alpha=0.5, linewidth=0)
    ax1.fill_between(Edep_posi, posi_nonl_min, posi_nonl_max, color="crimson", alpha=0.5, label=r"Fit $10\times$ error")
    lgam = ax1.errorbar(Edep_gamma, gam_nonl, yerr=gam_nonl_err, fmt="o--", fillstyle="none", color="slategray", lw=2.0, label=r"$\gamma$")
    lelec, = ax1.plot(Edep_elec, elec_nonl, "-", lw=2.0, color="blue", label=r"e$^-$")
    ax1.fill_between(Edep_elec, elec_nonl_min, elec_nonl_max, color="blue", alpha=0.3, label=r"Fit $10\times$ error")
    ax1.plot(1.022, m_npeGe68/Y/1.022, "*", ms=10, fillstyle="none", color="black")
    patchelec = mpatches.Patch(facecolor='royalblue', alpha=0.5, linewidth=0)
    #ax1.legend([(patchelec, lelec), (patchposi, lposi), lgam], [r"$e^-$", r"$e^+$", r"$\gamma$"], handler_map = {lelec: HandlerLine2D(marker_pad = 0), lposi: HandlerLine2D(marker_pad = 0)}, loc="lower right", prop={"size":15})
    ax1.text(1.022+0.20, m_npeGe68/Y/1.022-0.004, r"$^{68}$Ge", fontsize=14, color="black")
    ax1.grid(True)
    ax1.set_title("(a)", y =-0.3, fontsize=17)
    ax1.set_xlabel(r"$E$ [MeV]", fontsize=17)
    ax1.set_ylabel(r"$E_\mathrm{vis}/E$", fontsize=17)
    ax1.text(8.1, 0.96, r"$%d\times 1\sigma$ error band"%r, fontsize=15, color="black", bbox=dict(boxstyle="round", ec="black", fc="white"))
    
    gam_name = [r"$^{137}$Cs", r"$^{54}$Mn", r"$^{40}$K", "n-H", "AmBe", "AmC"]
    nonlx = [Edep_gamma[0]-0.95, Edep_gamma[1]-0.45, Edep_gamma[2]-0.45, Edep_gamma[3]-0.2, Edep_gamma[4]-0.1, Edep_gamma[5]-0.1 ]
    nonly = [gam_nonl[0]-0.01, gam_nonl[1]+0.01, gam_nonl[2]+0.01, gam_nonl[3]+0.006, gam_nonl[4]+0.006, gam_nonl[5]+0.005 ]
    resx  = [Evis_gamma[0], Evis_gamma[1], Evis_gamma[2]+0.1, Evis_gamma[3], Evis_gamma[4], Evis_gamma[5] ]
    resy  = [gam_res[0], gam_res[1], gam_res[2]-0.001, gam_res[3]+0.001, gam_res[4]+0.001, gam_res[5]+0.001 ]
    for i, j in enumerate(gam_name):
        ax1.text(nonlx[i], nonly[i], j, fontsize=14, color="dimgray")



    lposi, = ax2.plot(Edep_posi, posi_res, "-.", lw=2.0, color="red", label=r"e$^+$")
    lelec, = ax2.plot(Edep_elec, elec_res, "-", lw=2.0, color="blue", label=r"e$^-$")
    ax2.fill_between(Evis_elec, elec_res_min, elec_res_max, color="blue", alpha=0.5)
    ax2.fill_between(Evis_posi, posi_res_min, posi_res_max, color="crimson", alpha=0.5)
    lgam = ax2.errorbar(Evis_gamma, gam_res, yerr=gam_res_err, fmt="o--", lw=2, fillstyle="none", color="slategray")
    ax2.plot(m_npeGe68/Y, m_sigmaGe68/m_npeGe68, "*", ms=12, fillstyle="none", color="black")
    ax2.legend([(patchelec, lelec), (patchposi, lposi), lgam], [r"$e^-$", r"$e^+$", r"$\gamma$"], handler_map = {lelec: HandlerLine2D(marker_pad = 0), lposi: HandlerLine2D(marker_pad = 0)}, loc="upper right", prop={"size":15})
    ax2.grid(True)
    ax2.set_title("(b)", y =-0.3, fontsize=17)
    ax2.set_xlabel(r"$E_\mathrm{vis}$ [MeV]", fontsize=17)
    ax2.set_ylabel(r"$\sigma / E_\mathrm{vis}$", fontsize=17)
    ax2.text(1.7, m_sigmaGe68/m_npeGe68, r"$^{68}$Ge", fontsize=14, color="black")
    ax2.text(8.1, m_sigmaGe68/m_npeGe68, r"$%d\times 1\sigma$ error band"%r, fontsize=15, color="black", bbox=dict(boxstyle="round", ec="black", fc="white"))
    for i, j in enumerate(gam_name):
        ax2.text(resx[i], resy[i], j, fontsize=14, color="dimgray")
        


    plt.tight_layout()
    plt.savefig("Det1_kCerAnaNew1_kQInt_p2limits_allPar.pdf")
    plt.show()









if __name__ == "__main__" :
    compare_betas()

    #compare_particles()













