import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import elecLoader as eloader
import random
import ROOT

Etrue = np.arange(0.1, 8, 0.5)
orue_posi = Etrue + 1.022
dnum = len(Etrue)

def cerpe_func(E, A2, A3, A4):
    E0 = 0.2
    A1 = 0
    x = np.log(1+E/E0)
    return (A1*x+A2*x**2+A3*x**3) * (1/E+A4) * E 


def cerpe_sim(E, kC):
    return kC * eloader.getCerNPE(E)


def sctpe_sim(E, kB, scale):
    return eloader.getQPE(E, kB, scale)


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


from scipy.linalg import eigh, cholesky
from scipy.stats import norm


def sample_corelation(filename, num, sampleSize):
    method = 'eigenvectors'
    
    num_sample = sampleSize

    cov = loadCov(filename, num)
    
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


def sample_uncorelation(filename, num, sampleSize):

    num_sample = sampleSize

    cov = loadCov(filename, num)

    pars = [ np.zeros((1, sampleSize)) for i in range(num) ]
    for i in range(num):
        pars[i] = np.random.normal(loc=0, scale=np.sqrt(scale=cov[i ,i]), size=sampleSize )

    return pars
        


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



def main():

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


    best1, best2, best3 = [], [], []
    best4, best5, best6 = [], [], []
    simData = []
    diff, diff_err = [], []
    diffSim, diffSimerr = [], []
    diffAna, diffAnaerr = [], []
    nominal = []
    nn = 0
    for i, j in zip(Etrue, Etrue_posi):
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

    fig = plt.figure(figsize=(8, 6))
    spec = gridspec.GridSpec(ncols=1, nrows=2,
                         height_ratios=[1, 2])

    ax0 = fig.add_subplot(spec[0])
    ax1 = fig.add_subplot(spec[1])

    #ax0.fill_between(Etrue_posi, diff+diff_err, diff-diff_err, alpha=0.3, color="#E49D22")
    ax0.fill_between(Etrue_posi, diffSim+diffSimerr, diffSim-diffSimerr, alpha=0.3, color="blue")
    ax0.fill_between(Etrue_posi, diffAna+diffAnaerr, diffAna-diffAnaerr, alpha=0.3, color="#FF60F1")
    #ax0.plot(Etrue_posi, diff, "-", ms=4, color="#E49D22")
    ax0.plot(Etrue_posi, diffSim, "-", ms=4, color="blue")
    ax0.plot(Etrue_posi, diffAna, "-", ms=4, color="#FF60F1")
    ax0.grid(True)
    ax0.set_ylabel("Relaive bias", fontsize=15)
    ax0.tick_params(axis='both', which='major', labelsize=13)
    ax0.set_ylim(-0.02, 0.01)

    ax1.plot(Etrue_posi, best4, "-.",  color="blue",   label="Fitting: simulation-based")
    ax1.plot(Etrue_posi, best5, "--", color="#FF60F1",     label="Fitting: analytical")
    ax1.plot(Etrue_posi, simData, color="lightseagreen", label="Simulation truth")
    ax1.fill_between(Etrue_posi, ymin4, ymax4, alpha=1.0, color="lightskyblue")
    ax1.fill_between(Etrue_posi, ymin5, ymax5, alpha=0.5, color="pink")
    ax1.set_xlabel(r"$E_{dep}$/MeV", fontsize=15)
    ax1.set_ylabel(r"$E_{vis}/E_{dep}$", fontsize=15)
    ax1.tick_params(axis='both', which='major', labelsize=13)
    #plt.plot(Etrue_posi, best6, "-",  color="royalblue",   label="only gamma")
    
    ax1.legend(loc="lower right", prop={'size':14})
    ax1.grid(True)

    """
    plt.xlabel(r"$E_{dep}/MeV$")
    plt.ylabel(r"$E_{vis}/E_{dep}$")
    #plt.title("Electron NPE Nonlinearity")
    #plt.semilogy()
    plt.grid(True)
    plt.legend()
    """
    plt.savefig("Compare_positron_nonl.pdf")
    plt.show()


if __name__ == "__main__":
    main()
