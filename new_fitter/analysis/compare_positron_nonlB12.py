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

    fig = plt.figure(figsize=(6, 6))
    spec = gridspec.GridSpec(ncols=1, nrows=1)


    ax2 = fig.add_subplot(spec[0])


    ## nonlinearity part
    Etrue = np.arange(0.1, 14, 0.5)
    Etrue_posi = Etrue + 1.022
    dnum = len(Etrue)
    #filename4 = "../B12/gam_kIntegral_kAnalytical_kNPE_fixedres_nonlcov.txt"
    #filename5 = "../B12/gam+B12_kIntegral_kAnalytical_kNPE_fixedres_nonlcov.txt"
    filename4 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/articles/gam_kIntegral_kAnalytical_kNPE_fixedc0_nonlcov.txt"
    filename5 = "../B12/gam+B12_kIntegral_kAnalytical_kNPE_fixedres_nonlcov.txt"

    sampleSize = 50000
    sigma4 = sample_corelation(filename4, 5, sampleSize)
    sigma5 = sample_corelation(filename5, 5, sampleSize)

    #kA4, kB4, scale4, A24, A34, A44 = 1, 5.20e-3, 1407.88, -9.79, 14.05, 0.026
    kA4, kB4, scale4, A24, A34, A44 = 1, 5.08e-3, 1403.56, -8.75, 14.11, 0.029
    kA5, kB5, scale5, A25, A35, A45 = 1, 5.08e-3, 1404.28, -9.23, 14.28, 2.64e-2
    

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

        m_kA4 = kA4
        m_scale4 = scale4 + sigma4[0, i]
        m_kB4 = kB4 + sigma4[1, i]
        m_A24 = A24 + sigma4[2, i]
        m_A34 = A34 + sigma4[3, i]
        m_A44 = A44 + sigma4[4, i]

        m_kA5 = kA5
        m_scale5 = scale5 + sigma5[0, i]
        m_kB5 = kB5 + sigma5[1, i]
        m_A25 = A25 + sigma5[2, i]
        m_A35 = A35 + sigma5[3, i]
        m_A45 = A45 + sigma5[4, i]

        for i, j in zip(Etrue,  Etrue_posi):
            pe_arr4.append( (sctpe_sim(i, m_kB4, m_scale4) + cerpe_func(i, m_A24, m_A34, m_A44) + 2*660.8) / j / global_es)
            pe_arr5.append( (sctpe_sim(i, m_kB5, m_scale5) + cerpe_func(i, m_A25, m_A35, m_A45) + 2*660.8) / j / global_es)

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
        best4.append( (sctpe_sim(i, kB4, scale4) + cerpe_func(i, A24, A34, A44) + 2*660.8) / j / global_es)
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

    best4 = np.array(best4)
    ymin4 = np.array(ymin4)
    ymax4 = np.array(ymax4)
    best5 = np.array(best5)
    ymin5 = np.array(ymin5)
    ymax5 = np.array(ymax5)

    #ax0.fill_between(Etrue_posi, diff+diff_err, diff-diff_err, alpha=0.3, color="#E49D22")
    #ax0.fill_between(Etrue_posi, diffSim+5*diffSimerr, diffSim-5*diffSimerr, alpha=0.3, color="blue")
    #ax0.fill_between(Etrue_posi, diffAna+5*diffAnaerr, diffAna-5*diffAnaerr, alpha=0.5, color="dimgray")
    ##ax0.plot(Etrue_posi, diff, "-", ms=4, color="#E49D22")
    #ax0.plot(Etrue_posi, diffSim, "-.", lw=2, ms=4, color="blue")
    #ax0.plot(Etrue_posi, diffAna, "--", lw=2, ms=4, color="black")
    ##ax0.hlines(0.001, 1, 9, color='red')
    #ax0.grid(True)
    #ticknames = ["-0.015", "-0.01", "-0.005",  "0", "0.005", "0.01", "0.015"]
    #ax0.set_yticks([ -0.015, -0.01, -0.005, 0, 0.005, 0.01, 0.015])
    #ax0.set_yticklabels(ticknames)
    #ax0.set_ylabel("Bias", fontsize=17)
    #ax0.tick_params(axis='x', which='major', labelsize=15)
    #ax0.tick_params(axis='y', which='major', labelsize=14)
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
    
    ax2.plot(Etrue_posi, best4, "-.", lw=2, color="blue",    alpha=0.8,  label=r"Fitting w/ $\gamma$", zorder=1)
    ax2.plot(Etrue_posi, best5, "--", lw=2, color="black", alpha=0.8,    label=r"Fitting w/ $\gamma+{\rm ^{12}B}$", zorder=2)
    ax2.errorbar(subkesim1+1.022, submusim1/global_es/(subkesim1+1.022), yerr=submuerrsim1/global_es/(subkesim1+1.022), fmt="o", ms=8, alpha=1.0, color="crimson", fillstyle="none", label="Simulation truth", zorder=3)
    ax2.fill_between(Etrue_posi, best4-5*(best4-ymin4), best4+5*(ymax4-best4), alpha=0.3, color="blue")
    ax2.fill_between(Etrue_posi, best5-5*(best5-ymin5), best5+5*(ymax5-best5), alpha=0.5, color="dimgray")
    #ax2.fill_between(Etrue_posi, ymin5, ymax5, alpha=0.3, color="dimgray")
    ax2.set_xlabel(r"Positron $E_{dep}$[MeV]", fontsize=17)
    ax2.set_ylabel(r"$E_{vis}/E_{dep}$", fontsize=17)
    ax2.tick_params(axis='both', which='major', labelsize=16)
    #plt.plot(Etrue_posi, best6, "-",  color="royalblue",   label="only gamma")
    
    ax2.legend(loc="lower right", prop={'size':14})
    ax2.grid(True)


    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=0.3, hspace=0.2)
    plt.tight_layout()
    plt.savefig("compare_positron_B12.pdf")
    plt.show()



if __name__ == "__main__":
    main()
