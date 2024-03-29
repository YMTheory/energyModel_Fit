import numpy as np
import matplotlib.pyplot as plt
import elecLoader as eloader
from matplotlib import gridspec
import random
import ROOT


# no correlation in fitter...

global_es = 3134.078/2.223

def cerpe_func(E, A2, A3, A4):
    E0 = 0.2
    A1 = 0
    x = np.log(1+E/E0)
    return (A1*x+A2*x**2+A3*x**3) * (1/E+A4) * E 


def cerpe_sim(E, kC):
    return kC * eloader.getCerNPE(E)


def sctpe_sim(E, kB, scale):
    return eloader.getFitNsct(E, kB, scale, "Sim")


def resFunc(x, a, b, c):
    s2 = a + b*x + c*x**2
    if s2<0:
        return 0
    else:
        return np.sqrt(s2)

def cerFunc(x, a, b, c):
    s2 = a + b*x + c*x**2
    if s2<0:
        return 0
    else:
        return np.sqrt(s2)


def sctFunc(x, a):
    s2 = a*x
    if s2 <= 0:
        return 0
    else:
        return np.sqrt(s2)


def bestfit(p0, p1, p2):
    m_p0 = p0
    m_p1 = p1
    m_p2 = p2

    m_nonl = []
    for i in Etrue:
        m_nonl.append(resFunc(i, m_p0, m_p1, m_p2) )

    m_nonl = np.array(m_nonl)
    return m_nonl

def nominal():
    m_p0 = -2.17203
    m_p1 = 1.31498e3
    #m_p2 = 3134.078/2.223
    m_p2 = 1.60508e2

    m_nonl = []
    for i in Etrue:
        m_nonl.append(resFunc(i, m_p0, m_p1, m_p2) )

    m_nonl = np.array(m_nonl)
    return m_nonl

def loadCov(filename, num):
    row, col = 0, 0
    cov_mat = np.ones((num, num))
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            col = 0
            for i in data:
               cov_mat[row, col] = float(i) 
               cov_mat[col, row] = float(i)
               col+=1
            row += 1
    

    return cov_mat


def loadSimTruth():
    filename = '/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/positron/positronNPE.txt'
    ke, mu, muerr, sigma, sigmaerr = [], [], [], [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")

            ke.append(float(data[0]))
            mu.append(float(data[1]))
            muerr.append(float(data[2]))
            sigma.append(float(data[3]))
            sigmaerr.append(float(data[4]))

    ke = np.array(ke)
    mu = np.array(mu)
    muerr = np.array(muerr)
    sigma = np.array(sigma)
    sigmaerr = np.array(sigmaerr)
    return ke, mu, muerr, sigma, sigmaerr




from scipy.linalg import eigh, cholesky
from scipy.stats import norm


def sample_corelation(filename, num, sampleSize):
    method = 'eigenvectors'
    
    num_sample = sampleSize

    cov = loadCov(filename, num)
    #print(cov)
    
    x = norm.rvs(size=(num, num_sample))

    if method == 'cholesky':
        # Compute the Cholesky decomposition.
        c = cholesky(cov, lower=True)
    else:
        # Compute the eigenvalues and eigenvectors.
        evals, evecs = eigh(cov)
        # Construct c, so c*c^T = r.
        c = np.dot(evecs, np.diag(np.sqrt(evals)))

    y = np.dot(c, x)

    return y    


def cerpe_sim(E, kC):
    return kC * eloader.getCerNPE(E)


def sctpe_sim(E, kB, scale):
    return eloader.getQPE(E, kB, scale, "Sim")


def wenNcer(E, p0, p1, p2, p3, p4):
    if E<=0.2:
        return 0
    if E>0.2:
        E = E - 0.2
        #return p3*E**2/(p0+p4*E+p1*np.exp(-p2*E))
        return p3*E**2/(p4+p0*E+p1*np.exp(-p2*E))




def main():

    fig = plt.figure(figsize=(12, 6))
    spec = gridspec.GridSpec(ncols=2, nrows=2,
                         height_ratios=[1, 2])

    ax0 = fig.add_subplot(spec[0])
    ax1 = fig.add_subplot(spec[1])
    ax2 = fig.add_subplot(spec[2])
    ax3 = fig.add_subplot(spec[3])


    ## nonlinearity part
    Etrue = np.arange(0.1, 8, 0.5)
    Etrue_posi = Etrue + 1.022
    dnum = len(Etrue)
    #filename4 = "../output/articles/gam_kSimulation_kSimulationCer_kNPE_fixedc_nonlcov.txt"
    #filename5 = "../tmp_results/gam+B12_kIntegralCalc_kAnalyticalCer_kNPE_fixedp0_nonlcov.txt"
    #filename4 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/gam_kSimulation_kNewAnaCer_kNPE_nonlcov.txt"
    #filename5 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/gamB12_kSimulation_kNewAnaCer_kNPE_nonlcov.txt"
    filename4 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/articles_gam_kSimQ_kNewAnaC1_kNPE/gam_kSimulation_kNewAnaC1_kNPE_nonlcov.txt"
    filename5 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/articles_gamB12_kSimQ_kNewAnaC1_kNPE/gamB12_kSimulation_kNewAnaC1_kNPE_nonlcov.txt"
    #filename5 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/gamB12_kIntegral_kNewAnaCer_kNPE_nonlcov.txt"

    sampleSize = 50000
    #sigma4 = sample_corelation(filename4, 3, sampleSize)
    sigma4 = sample_corelation(filename4, 7, sampleSize)
    sigma5 = sample_corelation(filename5, 7, sampleSize)

    #kA4, kB4, scale4, kC4 = 1, 6.29208e-03, 1.40845e+03, 9.90770e-01
    #kA5, kB5, scale5, A25, A35, A45 = 1, 5.20018e-03, 1.41741e+03, -1.34659e+01 ,  1.46295e+01 , 1.08040e-02
    
    #kA4, kB4, scale4, kC4 = 1, 6.31e-3, 1.409e3, 0.989
    #kA5, kB5, scale5, kC5 = 1, 6.26e-3, 1.408e3, 0.996
    
    #kA4, kB4, scale4, A24, A34, A44 = 1, 5.08e-3, 1403.89, -8.97, 14.18, 0.028
    #kA5, kB5, scale5, A25, A35, A45 = 1, 5.08e-3, 1403.56, -8.75, 14.11, 2.90e-2
    
    # kSimQ_kNewAnaCer
    #scale4, kB4, p04, p14, p24, p34, p44 = 1407.55, 6.18e-3, -0.409, 0.421, -0.0224, 3.02, 0.0189
    #scale5, kB5, p05, p15, p25, p35, p45 = 1409.64, 6.30e-3, -0.409, 0.422, -0.0327, 3.14, 0.0155
    # kIntegral_kNewAnaCer
    #scale4, kB4, p04, p14, p24, p34, p44 = 1410.19, 5.45e-3, -0.407, 0.423, -0.0340, 3.40, 0.0192
    #scale5, kB5, p05, p15, p25, p35, p45 = 1410.19, 5.45e-3, -0.409, 0.423, -0.0327, 3.14, 0.0155
    # kSimQ_kNewAnaCer1
    #scale4, kB4, p04, p14, p24, p34, p44 = 1408.07, 6.21e-3,  4.20 , 3.64 , 0.115  , 405., -2.09 
    #scale5, kB5, p05, p15, p25, p35, p45 = 1411.39, 6.40e-3,  4.33 , 3.76 , 0.192  , 394., -2.01 
    # kIntQ_kNewAnaCer1
    #scale4, kB4, p04, p14, p24, p34, p44 = 1415.58, 5.73e-3,  5.03 , 1.37 , 1.161  , 452.,  1.76 
    #scale5, kB5, p05, p15, p25, p35, p45 = 1413.73, 5.63e-3,  5.17 , 3.83 , 0.283  , 435., -1.28 
    # kSimQ_kNewAnaCer1
    scale4, kB4, p04, p14, p24, p34, p44 = 1408.07, 6.21e-3,  4.20 , 3.64 , 0.115  , 405.,  -2.1
    scale5, kB5, p05, p15, p25, p35, p45 = 1411.39, 6.40e-3,  4.33 , 3.76 , 0.192  , 394., -2.01 

    global_es = 3134.078/2.223

    ymin4, ymax4 = [], []
    ymin5, ymax5 = [], []
    ymin6, ymax6 = [], []
    for i in range(dnum):
        ymin4.append(100)
        ymax4.append(-100)
        ymin5.append(100)
        ymax5.append(-100)
    

    for i in range(sampleSize):

        pe_arr1, pe_arr2, pe_arr3, pe_arr4, pe_arr5, pe_arr6 = [], [], [], [], [], []

        m_scale4 = scale4 + sigma4[0, i]
        m_kB4 = kB4 + sigma4[1, i]
        m_p04 = p04 + sigma4[2, i]
        m_p14 = p14 + sigma4[3, i]
        m_p24 = p24 + sigma4[4, i]
        m_p34 = p34 + sigma4[5, i]
        m_p44 = p44 + sigma4[6, i]

        m_scale5 = scale5 + sigma5[0, i]
        m_kB5 = kB5 + sigma5[1, i]
        m_p05 = p05 + sigma5[2, i]
        m_p15 = p15 + sigma5[3, i]
        m_p25 = p25 + sigma5[4, i]
        m_p35 = p35 + sigma5[5, i]
        m_p45 = p45 + sigma5[6, i]

        #m_kA4 = kA4
        #m_scale4 = scale4 + sigma4[0, i]
        #m_kB4 = kB4 + sigma4[1, i]
        #m_A24 = A24 + sigma4[2, i]
        #m_A34 = A34 + sigma4[3, i]
        #m_A44 = A44 + sigma4[4, i]

        #m_kA5 = kA5
        #m_scale5 = scale5 + sigma5[0, i]
        #m_kB5 = kB5 + sigma5[1, i]
        #m_A25 = A25 + sigma5[2, i]
        #m_A35 = A35 + sigma5[3, i]
        #m_A45 = A45 + sigma5[4, i]

        for i, j in zip(Etrue,  Etrue_posi):
            #pe_arr4.append( (sctpe_sim(i, m_kB4, m_scale4) + cerpe_sim(i, m_kC4) + 2*660.8 ) / j /global_es)
            #pe_arr5.append( (sctpe_sim(i, m_kB5, m_scale5) + cerpe_sim(i, m_kC5) + 2*660.8 ) / j /global_es)
            pe_arr4.append( (sctpe_sim(i, m_kB4, m_scale4) + wenNcer(i, m_p04, m_p14, m_p24, m_p34, m_p44) + 2*660.8 ) / j /global_es)
            pe_arr5.append( (sctpe_sim(i, m_kB5, m_scale5) + wenNcer(i, m_p05, m_p15, m_p25, m_p35, m_p45) + 2*660.8 ) / j /global_es)
            
            #pe_arr4.append( (sctpe_sim(i, m_kB4, m_scale4) + cerpe_func(i, m_A24, m_A34, m_A44) + 2*660.8) / j / global_es)
            #pe_arr5.append( (sctpe_sim(i, m_kB5, m_scale5) + cerpe_func(i, m_A25, m_A35, m_A45) + 2*660.8) / j / global_es)

        for n in range(dnum):
            if pe_arr4[n] < ymin4[n] :
                ymin4[n] = pe_arr4[n]
            if pe_arr4[n] > ymax4[n]:
                ymax4[n] = pe_arr4[n]
            if pe_arr5[n] < ymin5[n] :
                ymin5[n] = pe_arr5[n]
            if pe_arr5[n] > ymax5[n]:
                ymax5[n] = pe_arr5[n]


    kesim, musim, muerrsim, sigmasim, sigmaerrsim = loadSimTruth()
    gSim = ROOT.TGraph()
    for i in range(len(musim)):
        gSim.SetPoint(i, kesim[i], musim[i])
    kesim1, musim1, muerrsim1, sigmasim1, sigmaerrsim1 = [], [], [], [], []
    for i in range(len(musim)):
        if kesim[i] > Etrue[0] and kesim[i] < Etrue[-1]-0.3:
            kesim1.append(kesim[i])
            musim1.append(musim[i])
            muerrsim1.append(muerrsim[i])
            sigmasim1.append(sigmasim[i])
            sigmaerrsim1.append(sigmaerrsim[i])
    kesim1 = np.array(kesim1)
    musim1 = np.array(musim1)
    muerrsim1 = np.array(muerrsim1)
    sigmasim1 = np.array(sigmasim1)
    sigmaerrsim1 = np.array(sigmaerrsim1)


    best1, best2, best3 = [], [], []
    best4, best5, best6 = [], [], []
    simData = []
    diff, diff_err = [], []
    diffSim, diffSimerr = [], []
    diffAna, diffAnaerr = [], []
    nominal = []
    nn = 0
    nonl_elec, nonl_gam = [], []
    for i, j in zip(Etrue, Etrue_posi):
        #nonl_elec.append( (sctpe_sim(i, kB4, scale4) + cerpe_sim(i, kC4) ) /i/global_es)
        #best4.append( (sctpe_sim(i, kB4, scale4) + cerpe_sim(i, kC4) + 2*660.8) /j/global_es)
        #best5.append( (sctpe_sim(i, kB5, scale5) + cerpe_sim(i, kC5) + 2*660.8) /j/global_es)
        #best4.append( (sctpe_sim(i, kB4, scale4) + cerpe_func(i, A24, A34, A44) + 2*660.8) / j / global_es)
        #best5.append( (sctpe_sim(i, kB5, scale5) + cerpe_func(i, A25, A35, A45) + 2*660.8) / j / global_es)
        best4.append( (sctpe_sim(i, kB4, scale4) + wenNcer(i, p04, p14, p24, p34, p44) + 2*660.8 ) / j /global_es)
        #print(i, j, sctpe_sim(i, kB4, scale4), wenNcer(i, p04, p14, p24, p34, p44))
        best5.append( (sctpe_sim(i, kB5, scale5) + wenNcer(i, p05, p15, p25, p35, p45) + 2*660.8 ) / j /global_es)
        simData.append(gSim.Eval(i, 0, "S")/global_es/j)
        diff.append( (best5[-1] - best4[-1]) / best4[-1] )
        err4 = ymax4[nn] - best4[nn]
        err5 = ymax5[nn] - best5[nn]
        diff_err.append( np.sqrt(err5**2/best4[nn]**2 + err4**2*best5[nn]**2/best4[nn]**2) )
        diffSim.append( (best4[-1] - gSim.Eval(i, 0, "S")/global_es/j) / (gSim.Eval(i, 0, "S")/global_es/j) )
        diffSimerr.append(err4/(gSim.Eval(i, 0, "S")/global_es/j))
        diffAna.append( (best5[-1] - gSim.Eval(i, 0, "S")/global_es/j) / (gSim.Eval(i, 0, "S")/global_es/j) )
        diffAnaerr.append(err5/(gSim.Eval(i, 0, "S")/global_es/j))
        nn += 1

    diff = np.array(diff)
    diff_err = np.array(diff_err)
    diffSim = np.array(diffSim)
    diffSimerr = np.array(diffSimerr)
    diffAna = np.array(diffAna)
    diffAnaerr = np.array(diffAnaerr)

    best4 = np.array(best4)
    ymin4 = np.array(ymin4)
    ymax4 = np.array(ymax4)
    best5 = np.array(best5)
    ymin5 = np.array(ymin5)
    ymax5 = np.array(ymax5)

    #ax0.fill_between(Etrue_posi, diff+diff_err, diff-diff_err, alpha=0.3, color="#E49D22")
    ax0.fill_between(Etrue_posi, diffSim+5*diffSimerr, diffSim-5*diffSimerr, alpha=0.3, color="dimgray")
    #ax0.fill_between(Etrue_posi, diffAna+5*diffAnaerr, diffAna-5*diffAnaerr, alpha=0.5, color="blue")
    #ax0.plot(Etrue_posi, diff, "-", ms=4, color="#E49D22")
    ax0.plot(Etrue_posi, diffSim, "-.", lw=2, ms=4, color="black")
    #ax0.plot(Etrue_posi, diffAna, "--", lw=2, ms=4, color="blue")
    #ax0.hlines(0.001, 1, 9, color='red')
    ax0.grid(True)
    ticknames = ["-0.015", "-0.01", "-0.005",  "0", "0.005", "0.01", "0.015"]
    ax0.set_yticks([ -0.015, -0.01, -0.005, 0, 0.005, 0.01, 0.015])
    ax0.set_yticklabels(ticknames)
    ax0.set_ylabel("Bias", fontsize=17)
    ax0.tick_params(axis='x', which='major', labelsize=15)
    ax0.tick_params(axis='y', which='major', labelsize=14)
    #ax0.set_ylim(-0.015, 0.015)

    subkesim1, submusim1, submuerrsim1 = [], [], []
    for i in [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]:
        subkesim1.append(kesim1[i])
        submusim1.append(musim1[i])
        submuerrsim1.append(muerrsim1[i])

    subkesim1 = np.array(subkesim1)
    submusim1 = np.array(submusim1)
    submuerrsim1 = np.array(submuerrsim1)
    print(len(subkesim1), len(submusim1), len(submuerrsim1))
    
    ax2.plot(Etrue_posi, best4, "-.", lw=2, color="black",    alpha=0.8,  label="Fitting: simulation-based", zorder=1)
    #ax2.plot(Etrue_posi, best5, "--", lw=2, color="black", alpha=0.8,    label="Fitting: analytical", zorder=2)
    ax2.errorbar(subkesim1+1.022, submusim1/global_es/(subkesim1+1.022), yerr=submuerrsim1/global_es/(subkesim1+1.022), fmt="o", ms=8, alpha=1.0, color="crimson", fillstyle="none", label="Simulation truth", zorder=3)
    ax2.fill_between(Etrue_posi, best4-5*(best4-ymin4), best4+5*(ymax4-best4), alpha=0.3, color="dimgray")
    #ax2.fill_between(Etrue_posi, best5-5*(best5-ymin5), best5+5*(ymax5-best5), alpha=0.5, color="dimgray")
    #ax2.fill_between(Etrue_posi, ymin5, ymax5, alpha=0.3, color="dimgray")
    ax2.set_xlabel(r"Positron $E_{dep}$[MeV]", fontsize=17)
    ax2.set_ylabel(r"$E_{vis}/E_{dep}$", fontsize=17)
    ax2.tick_params(axis='both', which='major', labelsize=16)
    #plt.plot(Etrue_posi, best6, "-",  color="royalblue",   label="only gamma")
    
    #ax2.legend(loc="lower right", prop={'size':14})
    ax2.grid(True)


    ## resolution part
    Etrue = np.arange(1.5, 8, 0.1)
    Etrue_elec = np.arange(0.1, 8, 0.1)
    Etrue_posi = Etrue_elec + 1.022
    #npe = np.arange(400, 10000, 100)
    npe = Etrue_elec * global_es
    npe_posi = npe + 2*660.8
    dnum = len(Etrue_posi)
    sampleSize = 50000
    #sigma1 = sample_corelation("../output/articles/gam_kSimulation_kSimulationCer_kNPE_fixedc_rescov.txt", 2, sampleSize)
    #sigma1 = sample_corelation("../output/articles/gam+B12_kIntegral_kAnalytical_kNPE_rescov.txt", 2, sampleSize)
    #sigma2 = sample_corelation("../gam_kIntegral_kAnalytical_kNPE_fixedc0_rescov.txt", 2, sampleSize)
    #sigma2 = sample_corelation("../output/articles/gam_kIntegral_kAnalytical_kNPE_fixedc0_rescov.txt", 2, sampleSize)
    #filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/gam_kSimulation_kNewAnaCer_kNPE_rescov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/gamB12_kSimulation_kNewAnaCer_kNPE_rescov.txt"
    #filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/gam_kIntegral_kNewAnaCer_kNPE_rescov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/gam_kIntegral_kNewAnaCer_kNPE_rescov.txt"
    filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/articles_gam_kSimQ_kNewAnaC1_kNPE/gam_kSimulation_kNewAnaC1_kNPE_rescov.txt"
    filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/articles_gamB12_kSimQ_kNewAnaC1_kNPE/gamB12_kSimulation_kNewAnaC1_kNPE_rescov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/gamB12_kIntegral_kNewAnaCer_kNPE_rescov.txt"
    sigma1 = sample_corelation(filename1, 2, sampleSize)
    sigma2 = sample_corelation(filename2, 2, sampleSize)

    #pa1, pb1, pc1 = 0.990, 7.92e-3, 0
    #pa2, pb2, pc2 = 0.990, 7.92e-3, 0
    q03, q13, q23 = 0, 1.00710, 3.98696e-05

    # kInt_kNewAnaCer
    #pa1, pb1, pc1 = 0.993, 7.74e-3, 0
    #pa2, pb2, pc2 = 0.988, 8.026e-3, 0
    
    #kSim_kNewAnaCer1
    pa1, pb1, pc1 = 0.988, 7.89e-3, 0
    pa2, pb2, pc2 = 0.988, 8.09e-3, 0




    ymin1, ymax1 = [], []
    ymin2, ymax2 = [], []
    for i in range(dnum):
        ymin1.append(1000000)
        ymax1.append(-100)
        ymin2.append(1000000)
        ymax2.append(-100)


    for i in range(sampleSize):

        #a1 = q01 
        #b1 = q11 + sigma1[0, i]
        #c1 = q21 + sigma1[1, i]

        #a2 = q02 
        #b2 = q12 + sigma2[0, i]
        #c2 = q22 + sigma2[1, i]

        a1 = pa1 + sigma1[0, i]
        b1 = pb1 + sigma1[1, i]
        c1 = pc1 

        a2 = pa2 + sigma2[0, i]
        b2 = pb2 + sigma2[1, i]
        c2 = pc2 


        totsigma1 = []
        totsigma2 = []

        for i in npe:

            totsigma1.append( np.sqrt(resFunc(i, c1**2, a1**2, b1**2)**2 + 2*27.07**2 ) )
            totsigma2.append( np.sqrt(resFunc(i, c2**2, a2**2, b2**2)**2 + 2*27.07**2 ) )
        

        
        for n in range(dnum):
            if totsigma1[n] < ymin1[n]:
                ymin1[n] = totsigma1[n]
            if totsigma1[n] > ymax1[n]:
                ymax1[n] = totsigma1[n]

            if totsigma2[n] < ymin2[n]:
                ymin2[n] = totsigma2[n]
            if totsigma2[n] > ymax2[n]:
                ymax2[n] = totsigma2[n]
    
    ymin1 = np.array(ymin1)
    ymax1 = np.array(ymax1)
    ymin2 = np.array(ymin2)
    ymax2 = np.array(ymax2)


    kesim, musim, muerrsim, sigmasim, sigmaerrsim = loadSimTruth()
    musim1, muerrsim1, sigmasim1, sigmaerrsim1 = [], [], [], []
    gSim1 = ROOT.TGraph()
    for i in range(len(musim)):
        gSim1.SetPoint(i, musim[i], sigmasim[i])
        if kesim[i] > Etrue_elec[0] and kesim[i] < Etrue_elec[-1]-0.3:
            musim1.append(musim[i])
            muerrsim1.append(muerrsim[i])
            sigmasim1.append(sigmasim[i])
            sigmaerrsim1.append(sigmaerrsim[i])
    musim1 = np.array(musim1)
    muerrsim1 = np.array(muerrsim1)
    sigmasim1 = np.array(sigmasim1)
    sigmaerrsim1 = np.array(sigmaerrsim1)


    simData = []
    best, best_posi1, best_posi2 = [] , [], []
    diff, diff_err = [], []
    diffSim, diffSimerr = [], []
    diffAna, diffAnaerr = [], []


    nn = 0
    for i, j in zip(npe, npe_posi):
        #best.append( resFunc(i, q01, q11, q21))
        #best_posi1.append( np.sqrt(resFunc(i, q01, q11, q21)**2 + 2*27.07**2 ) )
        #best_posi2.append( np.sqrt(resFunc(i, q02, q12, q22)**2 + 2*27.07**2 ) )
        best_posi1.append( np.sqrt(resFunc(i, pc1, pa1**2, pb1**2)**2 + 2*27.07**2 ) )
        best_posi2.append( np.sqrt(resFunc(i, pc2, pa2**2, pb2**2)**2 + 2*27.07**2 ) )
        simData.append(resFunc(j, q03, q13, q23))

        diff.append( (best_posi1[nn] - best_posi2[nn])/best_posi1[nn] )
        err1 = ymax1[nn] - best_posi1[nn]
        err2 = ymax2[nn] - best_posi2[nn]
        diff_err.append( np.sqrt(err2**2/best_posi1[nn]**2 + err1**2*best_posi2[nn]**2/best_posi1[nn]**4) )

        diffSim.append((best_posi1[nn] - simData[-1])/simData[-1])
        diffSimerr.append( err1/simData[-1])
        diffAna.append((best_posi2[nn] - simData[-1])/simData[-1])
        diffAnaerr.append( err2/simData[-1])

        
        nn += 1


    simData = np.array(simData)

    best = np.array(best)
    best_posi1 = np.array(best_posi1)
    best_posi2 = np.array(best_posi2)
    diff = np.array(diff)
    diff_err = np.array(diff_err)

    diffSim = np.array(diffSim)
    diffSimerr = np.array(diffSimerr)
    diffAna = np.array(diffAna)
    diffAnaerr = np.array(diffAnaerr)

    subsigmasim1, subsigmaerrsim1 = [], []
    for i in [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]:
        subsigmasim1.append(sigmasim1[i])
        subsigmaerrsim1.append(sigmaerrsim1[i])

    subsigmasim1 = np.array(subsigmasim1)
    subsigmaerrsim1 = np.array(subsigmaerrsim1)

    #ax0.fill_between(npe_posi/global_es, diff+diff_err, diff-diff_err, alpha=0.3, color="#E49D22")
    ax1.fill_between(Etrue_posi, diffSim+diffSimerr, diffSim-diffSimerr, alpha=0.3, color="dimgray")
    #ax1.fill_between(Etrue_posi, diffAna+diffAnaerr, diffAna-diffAnaerr, alpha=0.5, color="blue")
    #ax0.plot(npe_posi/global_es, diff, "-", ms=4, color="#E49D22")
    ax1.plot(Etrue_posi, diffSim, "-.", lw=2, ms=4, color="black")
    #ax1.plot(Etrue_posi, diffAna, "--", lw=2, ms=4, color="blue")
    ticknames = ["-0.1", "-0.05",  "0", "0.05", "0.1", "0.15", "0.2"]
    ax1.set_yticks([ -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2])
    ax1.set_yticklabels(ticknames)
    #ax1.hlines(0.10, 1, 9, color='red')
    ax1.grid(True)
    ax1.set_ylabel("Bais", fontsize=17)
    ax1.tick_params(axis='both', which='major', labelsize=15)
    ax1.set_ylim(-0.10, 0.20)
    ax1.tick_params(axis='x', which='major', labelsize=15)
    ax1.tick_params(axis='y', which='major', labelsize=14)

    #plt.plot(Etrue, nominal(), "--", color="darkviolet", label="nominal")
    ax3.plot(npe_posi/global_es, best_posi1/npe_posi, "-.", lw=2, color="black", alpha=1.0, label="Fitting: simulation-based", zorder = 1)
    ax3.fill_between(npe_posi/global_es, ymin1/npe_posi, ymax1/npe_posi, color="dimgray", alpha=0.3)
    #ax3.plot(npe_posi/global_es, best_posi2/npe_posi, "--", lw=2, color="blue", alpha=1.0, label="Fitting", zorder = 2)
    #ax3.fill_between(npe_posi/global_es, ymin2/npe_posi, ymax2/npe_posi, color="blue", alpha=0.5)
    ax3.errorbar(submusim1/global_es, subsigmasim1/submusim1, yerr=np.sqrt(subsigmaerrsim1**2/submusim1**2 + submuerrsim1**2*subsigmasim1**2/submusim1**4) ,fmt="o", ms=8, alpha=1.0, color='crimson', fillstyle="none", label="Simulation truth", zorder = 3)
    #ax1.plot(npe_posi/global_es, simData/npe_posi, "-", color="red")
    ax3.set_xlabel(r"Positron $E_{vis}$[MeV]", fontsize=17)
    ax3.set_ylabel(r"$\sigma_E/E_{vis}$", fontsize=17)
    ax3.legend(prop={"size":15})
    ax3.tick_params(axis='both', which='major', labelsize=16)
    ax3.grid(True)

    """
    plt.plot(npe/global_es, best/npe, "-", color="#FF8D70",  label="electron")
    plt.plot(npe_posi/global_es, best_posi1/npe_posi, "-", color="#FF60F1", label="positron")
    plt.fill_between(npe_posi/global_es, ymin1/npe_posi, ymax1/npe_posi, color="lightskyblue", label=r"positron: $1\sigma$ zone")
    plt.plot(gamEtrue/global_es, (gamResCalc), "--", color="#7054FF", label="single gamma: fitting")
    plt.errorbar(gamEtrue/global_es, (gamResData), yerr=gamResDataErr, fmt="o", color="#90F7C4", label="single gamma: simulation")



    #plt.title("Electron NPE Sigma")

    plt.xlabel(r"$E_{vis}/MeV$")
    plt.ylabel(r"$\sigma_{N_{tot}}/N_{tot}$")
    plt.legend()
    plt.grid(True)
    """

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=0.3, hspace=0.2)
    plt.tight_layout()
    plt.savefig("compare_positron2.pdf")
    #plt.show()



if __name__ == "__main__":
    main()
