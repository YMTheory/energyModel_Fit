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
    return eloader.getQPE(E, kB, scale)


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


def loadGamTruth():
    gamE, gammu, gammuerr, gamsigma, gamsigmamu = [], [], [], [], []
    with open("../data/gamma/gamma_J19.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if data[0] != "Ge68" and data[0]!="Co60" and data[0]!="nC12" and data[0]!="gamma3215" and data[0]!="gamma4440":
                gamE.append(float(data[1]))
                gammu.append(float(data[2]))
                gammuerr.append(float(data[3]))
                gamsigma.append(float(data[4]))
                gamsigmamu.append(float(data[5]))
                print(data[0], gamE[-1], gammu[-1]/global_es/gamE[-1])

    gamE = np.array(gamE)
    gammu = np.array(gammu)
    gammuerr = np.array(gammuerr)
    gamsigma = np.array(gamsigma)
    gamsigmamu = np.array(gamsigmamu)

    return gamE, gammu, gammuerr, gamsigma, gamsigmamu



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
    print(cov)
    
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
    return eloader.getQPE(E, kB, scale)



def main():

    fig = plt.figure(figsize=(12, 5))
    spec = gridspec.GridSpec(ncols=2, nrows=1 )

    ax0 = fig.add_subplot(spec[0])
    ax1 = fig.add_subplot(spec[1])


    ## nonlinearity part
    Etrue = np.arange(0.1, 8, 0.5)
    Etrue_posi = Etrue + 1.022
    dnum = len(Etrue)
    filename4 = "../tmp_results/gam+B12_kSimulation_kSimulationCer_kNPE_fixedq0_nonlcov1.txt"
    filename5 = "../tmp_results/gam+B12_kIntegralCalc_kAnalyticalCer_kNPE_fixedp0_nonlcov.txt"


    sampleSize = 5000
    sigma4 = sample_corelation(filename4, 3, sampleSize)
    sigma5 = sample_corelation(filename5, 5, sampleSize)

    kA4, kB4, scale4, kC4 = 1, 6.29208e-03, 1.40845e+03, 9.90770e-01
    kA5, kB5, scale5, A25, A35, A45 = 1, 5.20018e-03, 1.41741e+03, -1.34659e+01 ,  1.46295e+01 , 1.08040e-02
    
    global_es = 3134.078/2.223

    ymin4, ymax4 = [], []
    ymin5, ymax5 = [], []
    ymin6, ymax6 = [], []
    for i in range(dnum):
        ymin4.append(100)
        ymax4.append(-100)
        ymin5.append(100)
        ymax5.append(-100)
        ymin6.append(100)
        ymax6.append(-100)
    

    for i in range(sampleSize):

        pe_arr1, pe_arr2, pe_arr3, pe_arr4, pe_arr5, pe_arr6 = [], [], [], [], [], []

        m_kA4 = kA4
        m_kB4 = kB4 + sigma4[1, i]
        m_kC4 = kC4 + sigma4[2, i]
        m_scale4 = scale4 + sigma4[0, i]

        m_kA5 = kA5
        m_scale5 = scale5 + sigma5[0, i]
        m_kB5 = kB5 + sigma5[1, i]
        m_A25 = A25 + sigma5[2, i]
        m_A35 = A35 + sigma5[3, i]
        m_A45 = A45 + sigma5[4, i]

        for i, j in zip(Etrue,  Etrue_posi):
            pe_arr4.append( (sctpe_sim(i, m_kB4, m_scale4) + cerpe_sim(i, m_kC4) + 2*660.8 ) / j /global_es)
            pe_arr5.append( (sctpe_sim(i, m_kB5, m_scale5) + cerpe_func(i, m_A25, m_A35, m_A45) + 2*660.8) / j / global_es)
            pe_arr6.append( (sctpe_sim(i, m_kB4, m_scale4) + cerpe_sim(i, m_kC4) ) / i /global_es)

        for n in range(dnum):
            if pe_arr4[n] < ymin4[n] :
                ymin4[n] = pe_arr4[n]
            if pe_arr4[n] > ymax4[n]:
                ymax4[n] = pe_arr4[n]
            if pe_arr5[n] < ymin5[n] :
                ymin5[n] = pe_arr5[n]
            if pe_arr5[n] > ymax5[n]:
                ymax5[n] = pe_arr5[n]
            if pe_arr6[n] < ymin6[n] :
                ymin6[n] = pe_arr6[n]
            if pe_arr6[n] > ymax6[n]:
                ymax6[n] = pe_arr6[n]


    kesim, musim, muerrsim, sigmasim, sigmaerrsim = loadSimTruth()
    gamE, gammu, gammuerr, gamsigma, gamsigmaerr = loadGamTruth()

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
        nonl_elec.append( (sctpe_sim(i, kB4, scale4) + cerpe_sim(i, kC4) ) /i/global_es)
        best4.append( (sctpe_sim(i, kB4, scale4) + cerpe_sim(i, kC4) + 2*660.8) /j/global_es)
        best5.append( (sctpe_sim(i, kB5, scale5) + cerpe_func(i, A25, A35, A45) + 2*660.8) / j / global_es)
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


    ax0.plot(Etrue_posi, best4, "-.",  color="royalblue", lw=2, label="Positron", zorder=3)
    ax0.fill_between(Etrue_posi, ymin4, ymax4, alpha=1.0, color="lightskyblue")
    ax0.plot(Etrue, nonl_elec, "-", lw=2, color="coral", label="Electron")
    ax0.fill_between(Etrue, ymin6, ymax6, alpha=1.0, color="bisque")
    ax0.errorbar(gamE, gammu/global_es/gamE, yerr=gammuerr/global_es/gamE,fmt="o--", ms=5, color="seagreen", label="Gamma")
    ax0.set_xlabel(r"$E_{dep}$/MeV", fontsize=14)
    ax0.set_ylabel(r"$E_{vis}/E_{dep}$", fontsize=14)
    ax0.tick_params(axis='both', which='major', labelsize=13)
    #plt.plot(Etrue_posi, best6, "-",  color="royalblue",   label="only gamma")
    
    #ax2.legend(loc="lower right", prop={'size':14})
    ax0.grid(True)


    ## resolution part
    Etrue = np.arange(1.5, 8, 0.1)
    Etrue_elec = np.arange(0.1, 8, 0.1)
    Etrue_posi = Etrue_elec + 1.022
    #npe = np.arange(400, 10000, 100)
    npe = Etrue_elec * global_es
    npe_posi = npe + 2*660.8
    dnum = len(Etrue_posi)
    sampleSize = 5000
    #sigma = sample_corelation("../gam+B12_kSimulationQuench_kSimulationCer_kNPE_rescov.txt", 3, sampleSize)
    sigma1 = sample_corelation("../tmp_results/gam+B12_kSimulation_kSimulationCer_kNPE_fixedq0_rescov1.txt", 2, sampleSize)
    sigma2 = sample_corelation("../tmp_results/gam+B12_kIntegralCalc_kAnalyticalCer_kNPE_fixedp0_rescov.txt", 2, sampleSize)
    sigma4 = sample_corelation("../gam+B12+mic_kSimulationQuench_kSimulationCer_kNPE_resolcov.txt", 2, sampleSize)

    kA4, kB4, scale4, kC4 = 1, 6.29208e-03, 1.40845e+03, 9.90770e-01
    kA5, kB5, scale5, A25, A35, A45 = 1, 5.20018e-03, 1.41741e+03, -1.34659e+01 ,  1.46295e+01 , 1.08040e-02

    q01, q11, q21 = 0,  9.79299e-01, 6.12574e-05
    q02, q12, q22 = 0,  9.80275e-01, 6.68470e-05
    q03, q13, q23 = 0, 1.00710, 3.98696e-05
    q04, q14, q24 = 0, 0.9787, 6.2306e-5


    ymin1, ymax1 = [], []
    ymin2, ymax2 = [], []
    for i in range(dnum):
        ymin1.append(1000000)
        ymax1.append(-100)
        ymin2.append(1000000)
        ymax2.append(-100)


    for i in range(sampleSize):

        a1 = q01 
        b1 = q11 + sigma1[0, i]
        c1 = q21 + sigma1[1, i]

        a2 = q02 
        b2 = q12 + sigma2[0, i]
        c2 = q22 + sigma2[1, i]

        totsigma1 = []
        totsigma2 = []

        for i in npe:

            totsigma1.append( np.sqrt(resFunc(i, a1, b1, c1)**2 + 2*27.07**2 ) )
            totsigma2.append( np.sqrt(resFunc(i, a2, b2, c2)**2 + 2*27.07**2 ) )
        

        
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
        best.append( resFunc(i, q01, q11, q21))
        best_posi1.append( np.sqrt(resFunc(i, q01, q11, q21)**2 + 2*27.07**2 ) )
        best_posi2.append( np.sqrt(resFunc(i, q02, q12, q22)**2 + 2*27.07**2 ) )
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


    # Electron Resolution
    resE_elec, res_elec, res_elec2 = [], [], []
    for i in np.arange(0.5*global_es, 8*global_es, 0.1*global_es):
        resE_elec.append(i)
        res_elec.append(resFunc(i, q01, q11, q21))
        res_elec2.append(resFunc(i, q04, q14, q24))

    resE_elec = np.array(resE_elec)
    res_elec = np.array(res_elec)
    res_elec2 = np.array(res_elec2)



    ymin7, ymax7 = [], []
    for i in range(len(resE_elec)):
        ymin7.append(1000000)
        ymax7.append(-100)
    ymin9, ymax9 = [], []
    for i in range(len(resE_elec)):
        ymin9.append(1000000)
        ymax9.append(-100)


    for i in range(sampleSize):

        a1 = q01 
        b1 = q11 + sigma1[0, i]
        c1 = q21 + sigma1[1, i]

        a2 = q04 
        b2 = q14 + sigma4[0, i]
        c2 = q24 + sigma4[1, i]

        totsigmae = []
        totsigmae1 = []

        for i in resE_elec:

            totsigmae.append( np.sqrt(resFunc(i, a1, b1, c1)**2 ) )
            totsigmae1.append( np.sqrt(resFunc(i, a2, b2, c2)**2 ) )
        

        for n in range(len(resE_elec)):
            if totsigmae[n] < ymin7[n]:
                ymin7[n] = totsigmae[n]
            if totsigmae[n] > ymax7[n]:
                ymax7[n] = totsigmae[n]
        for n in range(len(resE_elec)):
            if totsigmae1[n] < ymin9[n]:
                ymin9[n] = totsigmae1[n]
            if totsigmae1[n] > ymax9[n]:
                ymax9[n] = totsigmae1[n]

    ymin7 = np.array(ymin7)
    ymax7 = np.array(ymax7)
    ymin9 = np.array(ymin9)
    ymax9 = np.array(ymax9)
    
    print(ymin7)
    print(ymin9)


    #ax1.plot(npe_posi/global_es, best_posi1/npe_posi, "-.", color="blue", label="Positron", zorder = 2)
    #ax1.fill_between(npe_posi/global_es, ymin1/npe_posi, ymax1/npe_posi, color="lightskyblue")
    #ax1.errorbar(gammu/global_es, gamsigma/gammu, yerr=np.sqrt(gammuerr**2*gamsigma**2/gammu**4 + gamsigmaerr**2/gammu**2), fmt="o--", ms=5, color="seagreen", label="Gamma")
    ax1.plot(resE_elec/global_es, res_elec/(resE_elec), color="coral", lw=2, label="Electron")
    ax1.fill_between(resE_elec/global_es, ymin7/resE_elec, ymax7/resE_elec, alpha=0.6, color="bisque")
    ax1.plot(resE_elec/global_es, res_elec2/(resE_elec), color="seagreen", lw=2, label="Electron")
    ax1.fill_between(resE_elec/global_es, ymin9/resE_elec, ymax9/resE_elec, alpha=0.2, color="green")
    ax1.set_xlabel(r"$E_{vis}$/MeV", fontsize=14)
    ax1.set_ylabel(r"$\sigma/N_{tot}$", fontsize=14)
    ax1.legend(prop={"size":13})
    ax1.tick_params(axis='both', which='major', labelsize=13)
    ax1.grid(True)

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=0.02, hspace=0.02)
    plt.tight_layout()
    #plt.savefig("compare_particles.pdf")
    plt.show()



if __name__ == "__main__":
    main()
