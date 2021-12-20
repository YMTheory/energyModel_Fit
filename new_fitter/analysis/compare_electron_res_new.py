import numpy as np
import matplotlib.pyplot as plt
import elecLoader as el
from matplotlib import gridspec
import random
import ROOT


# no correlation in fitter...

global_es = 3134.078/2.223

def wenNcer(E, p0, p1, p2, p3, p4):
    if E<=0.2:
        return 0
    if E>0.2:
        E = E - 0.2
        return p3*E**2/(p4+p0*E+p1*np.exp(-p2*E))



def fitNsct(E, As, kB):
    if E > 15:
        return el.getFitNsct(15.9, kB, As, "Sim") / 15.9 * E
    else:
        return el.getFitNsct(E, kB, As, "Sim")


def sigma2Func(a, b, c, A, x):
    s2 = a**2/A*x + b**2*x**2 + c**2/A**2
    return s2 if s2>0  else 0


def sigma2FuncNew(a, b, n, A, x):
    s2 = a**2/A*x + b**2*np.pow(x, n)
    return s2 if s2>0  else 0


def loadBestFit(filename, num):
    bestFit, err = [], []
    row = 0
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if row == num:
                for kk in range(num):
                    bestFit.append(float(data[kk]))
            if row == num+1:
                for kk in range(num):
                    err.append(float(data[kk]))
            row += 1

    return bestFit, err


def loadCov(filename, num):
    row, col = 0, 0
    cov_mat = np.ones((num, num))
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if row == num:
                break
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



def loadSimTruth():
    totpe, totpeerr, Earr, resol, resolerr = [], [], [], [], []
    nn = 1
    filename = "../data/electron/elecResol4.txt"
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Earr.append(float(data[0]))
            totpe.append(float(data[1]))
            totpeerr.append(float(data[2]))
            resol.append(float(data[5]))
            resolerr.append(float(data[6]))
    Earr  = np.array(Earr)
    totpe = np.array(totpe)
    totpeerr = np.array(totpeerr)
    resol = np.array(resol)
    resolerr = np.array(resolerr)
    return Earr, totpe, totpeerr, resol, resolerr






def main():

    fig = plt.figure(figsize=(12, 6))
    spec = gridspec.GridSpec(ncols=2, nrows=2,
                         height_ratios=[1, 2])

    ax0 = fig.add_subplot(spec[0])
    ax1 = fig.add_subplot(spec[1])
    ax2 = fig.add_subplot(spec[2])
    ax3 = fig.add_subplot(spec[3])


    ## nonlinearity part
    Etrue = np.arange(0.1, 52.1, 0.5)
    dnum = len(Etrue)
    filename4 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/articles_gamB12Mic_kSimQ_kNewAnaC1_kNPE/gamB12Mic_kSimulation_kNewAnaCer_kNPE_nonlcov.txt"

    sampleSize = 50000
    sigma4 = sample_corelation(filename4, 7, sampleSize)

    kA4, scale4, kB4, p04, p14, p24, p34, p44 = 1, 1406.61, 6.06e-3, 0.91, 10.2, 0.02, 77.22, -9.9
    a4, b4, c4 = 1.018, 3.3e-3, 0

    global_es = 3134.078/2.223

    ymin4, ymax4 = [], []
    for i in range(dnum):
        ymin4.append(100)
        ymax4.append(-100)
    

    for j in range(sampleSize):

        pe_arr4 = []

        m_scale4 = scale4 + sigma4[0, j]
        m_kB4 = kB4 + sigma4[1, j]
        m_p04 = p04 + sigma4[2, j]
        m_p14 = p14 + sigma4[3, j]
        m_p24 = p24 + sigma4[4, j]
        m_p34 = p34 + sigma4[5, j]
        m_p44 = p44 + sigma4[6, j]


        for i in Etrue:
            pe_arr4.append( (fitNsct(i, scale4, kB4) + wenNcer(i, m_p04, m_p14, m_p24, m_p34, m_p44)) / i /global_es)
            

        for n in range(dnum):
            if pe_arr4[n] < ymin4[n] :
                ymin4[n] = pe_arr4[n]
            if pe_arr4[n] > ymax4[n]:
                ymax4[n] = pe_arr4[n]

    kesim, totpesim, totpeerrsim, resolsim, resolsimerr = loadSimTruth()
    gNonl = ROOT.TGraph()
    gRes  = ROOT.TGraph()
    for i in range(len(kesim)):
        gNonl.SetPoint(i, kesim[i], totpesim[i]/global_es/kesim[i])
        gRes.SetPoint(i, totpesim[i]/global_es, resolsim[i])

    best4 = []
    simData = []
    diff4, differr4 = [], []
    nn = 0
    for i in Etrue:
        best4.append( (fitNsct(i, scale4, kB4) + wenNcer(i, p04, p14, p24, p34, p44)) / i /global_es)
        #print(i, fitNsct(i, scale4, kB4), wenNcer(i, p04, p14, p24, p34, p44), best4[-1])
        simData.append(gNonl.Eval(i, 0, "S"))
        diff4.append( (best4[-1] - simData[-1]) / (simData[-1]) )
        err4 = ymax4[nn] - best4[nn]
        differr4.append(err4/(gNonl.Eval(i, 0, "S")))
        nn += 1

    best4    = np.array(best4)
    ymin4    = np.array(ymin4)
    ymax4    = np.array(ymax4)
    diff4    = np.array(diff4)
    differr4 = np.array(differr4)




    ax0.fill_between(Etrue, diff4+1*differr4, diff4-1*differr4, alpha=0.3, color="dimgray")
    ax0.plot(Etrue, diff4, "-.", lw=2, ms=4, color="black")
    ax0.grid(True)
    #ticknames = ["-0.015", "-0.01", "-0.005",  "0", "0.005", "0.01", "0.015"]
    #ax0.set_yticks([ -0.015, -0.01, -0.005, 0, 0.005, 0.01, 0.015])
    #ax0.set_yticklabels(ticknames)
    ax0.set_ylabel("Bias", fontsize=17)
    ax0.tick_params(axis='x', which='major', labelsize=15)
    ax0.tick_params(axis='y', which='major', labelsize=14)
    #ax0.set_ylim(-0.015, 0.015)

    
    ax2.plot(Etrue, best4, "-.", lw=2, color="black",    alpha=0.8,  label="Fitting", zorder=1)
    ax2.errorbar(kesim, totpesim/global_es/kesim, yerr=totpeerrsim/global_es/kesim , fmt="o", ms=5, color="crimson", label="Simulation truth", zorder=2)
    #ax2.errorbar(subkesim1+1.022, submusim1/global_es/(subkesim1+1.022), yerr=submuerrsim1/global_es/(subkesim1+1.022), fmt="o", ms=8, alpha=1.0, color="crimson", fillstyle="none", label="Simulation truth", zorder=3)
    ax2.fill_between(Etrue, best4-1*(best4-ymin4), best4+1*(ymax4-best4), alpha=0.3, color="dimgray")
    ax2.set_xlabel(r"Positron $E_{dep}$[MeV]", fontsize=17)
    ax2.set_ylabel(r"$E_{vis}/E_{dep}$", fontsize=17)
    ax2.tick_params(axis='both', which='major', labelsize=16)
    #plt.plot(Etrue_posi, best6, "-",  color="royalblue",   label="only gamma")
    
    #ax2.legend(loc="lower right", prop={'size':14})
    ax2.grid(True)




    ## resolution part
    Etrue = np.arange(0.1, 56.1, 0.5)
    dnum = len(Etrue)
    sampleSize = 50000
    filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/articles_gamB12Mic_kSimQ_kNewAnaC1_kNPE/gamB12Mic_kSimulation_kNewAnaCer_kNPE_rescov.txt"
    #filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12NewMicwhole_kSimQ_kNewAnaCer_kNew/NewgamB12NewMic_kSimQ_kNewCerAna_kNew_nonlcov.txt"
    #sigma1 = sample_corelation(filename1, 2, sampleSize)
    pa1, pb1, pc1 = 1.019, 3.297e-3, 0
    qa1, qb1, qc1 = 1.086, 4.321e-3, 0

    ymin1, ymax1 = [], []
    for i in range(dnum):
        ymin1.append(1000000)
        ymax1.append(-100)


    for j in range(sampleSize):

        a1 = pa1 + sigma1[0, j]
        b1 = pb1 + sigma1[1, j]
        c1 = pc1 

        totsigma1 = []
        for i in Etrue:
            totsigma1.append( np.sqrt(sigma2Func(a1, b1, c1, global_es, i))/i )
        
        
        for n in range(len(Etrue)):
            if totsigma1[n] < ymin1[n]:
                ymin1[n] = totsigma1[n]
            if totsigma1[n] > ymax1[n]:
                ymax1[n] = totsigma1[n]
    
    ymin1 = np.array(ymin1)
    ymax1 = np.array(ymax1)



    best1 = []
    diff, diff_err = [], []

    nn = 0
    for i in Etrue:
        best1.append( np.sqrt(sigma2Func(pa1, pb1, pc1, global_es, i))/i)
        #simData.append(np.sqrt(sigma2Func(qa1, qb1, qc1, global_es, i))/i)
        simData.append(gRes.Eval(i))

        diff.append((best1[-1] - simData[-1]) / simData[-1])
        err = ymax1[nn] - best1[nn]
        diff_err.append(np.sqrt(err**2/simData[-1]**2))

        nn += 1


    simData = np.array(simData)

    best1 = np.array(best1)
    diff = np.array(diff)
    diff_err = np.array(diff_err)


    ax1.fill_between(Etrue, diff+diff_err, diff-diff_err, alpha=0.3, color="dimgray")
    ax1.plot(Etrue, diff, "-.", lw=2, ms=4, color="black")
    #ax1.plot(Etrue_posi, diffAna, "--", lw=2, ms=4, color="blue")
    #ticknames = ["-0.1", "-0.05",  "0", "0.05", "0.1", "0.15", "0.2"]
    #ax1.set_yticks([ -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2])
    #ax1.set_yticklabels(ticknames)
    #ax1.hlines(0.10, 1, 9, color='red')
    ax1.grid(True)
    ax1.set_ylabel("Bais", fontsize=17)
    ax1.tick_params(axis='both', which='major', labelsize=15)
    #ax1.set_ylim(-0.10, 0.20)
    ax1.tick_params(axis='x', which='major', labelsize=15)
    ax1.tick_params(axis='y', which='major', labelsize=14)

    ax3.plot(Etrue, best1, "-.", lw=2, color="black", alpha=1.0, label="Fitting: simulation-based", zorder = 1)
    ax3.fill_between(Etrue, ymin1, ymax1, color="dimgray", alpha=0.3)
    ax3.errorbar(totpesim/global_es, resolsim, yerr=resolsimerr , fmt="o", ms=5, color="crimson", label="Simulation truth", zorder=2)
    ax3.set_xlabel(r"Positron $E_{vis}$[MeV]", fontsize=17)
    ax3.set_ylabel(r"$\sigma_E/E_{vis}$", fontsize=17)
    ax3.semilogy()
    ax3.legend(prop={"size":15})
    ax3.tick_params(axis='both', which='major', labelsize=16)
    ax3.grid(True)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
