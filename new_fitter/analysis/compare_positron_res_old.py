import numpy as np
import matplotlib.pyplot as plt
import elecLoader as eloader
from matplotlib import gridspec
import random


# no correlation in fitter...

global_es = 3134.078/2.223
Etrue = np.arange(0.5, 8, 0.1)
Etrue_elec = np.arange(0.1, 8, 0.1)
Etrue_posi = Etrue_elec + 1.022
#npe = np.arange(400, 10000, 100)
npe = Etrue_elec * global_es
npe_posi = npe + 2*660.8
dnum = len(Etrue_posi)


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

    sampleSize = 5000
    #sigma = sample_corelation("../gam+B12_kSimulationQuench_kSimulationCer_kNPE_rescov.txt", 3, sampleSize)
    sigma1 = sample_corelation("../tmp_results/gam+B12_kSimulation_kSimulationCer_kNPE_fixedq0_rescov1.txt", 2, sampleSize)
    sigma2 = sample_corelation("../tmp_results/gam+B12_kIntegralCalc_kAnalyticalCer_kNPE_fixedp0_rescov.txt", 2, sampleSize)

    kA, kB, scale, kC = 1, 6.32175e-3, 1408.910, 0.986594
    c0, c1, c2, s0 = -119.07, 3.35634, 0.00390735, 1
    ra, rb, rc = -1.10053e1, 1.48417e3, 1.21885e2

    #q01, q11, q21 = -4.74923, 1.0457, 4.83731e-5
    #q02, q12, q22 = -2.28558, 1.01722, 5.87549e-5

    q01, q11, q21 = 0,  9.79299e-01, 6.12574e-05
    q02, q12, q22 = 0,  9.80275e-01, 6.68470e-05



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

    #print(ymin[9]**2, ymin[19]**2, ymin[29]**2, ymin[49]**2, ymin[79]**2)
    #print(ymax[9]**2, ymax[19]**2, ymax[29]**2, ymax[49]**2, ymax[79]**2)

    


    gamEtrue = [870.26, 1115.197, 2017.245, 3134.078, 6392.879, 8916.649]
    gamResData = [0.03667723, 0.03270304, 0.02632371, 0.02269925, 0.01655166,  0.01399196]
    gamResDataErr = [0.00037148, 0.0003263,  0.00028073, 0.00023905, 0.00017784, 0.00014814]
    gamResCalc = [0.0363443,  0.03284537, 0.02644983, 0.02247591, 0.01648638, 0.01400519]
    gamEtrue = np.array(gamEtrue)
    gamResCalc = np.array(gamResCalc)
    gamResData = np.array(gamResData)
    gamResDataErr = np.array(gamResDataErr)


    best, best_posi1, best_posi2 = [] , [], []
    diff, diff_err = [], []
    #bc0, bc1, bc2 = -158.41, 4.333, 0.00181082
    ##for i, j in zip(npe, npe_posi):
    #for i, j in zip(Etrue, Etrue_posi):
    #    sctpe = sctpe_sim(i, kA, kB, scale)
    #    cerpe = cerpe_sim(i, kC)
    #    pp = sctpe / (sctpe + cerpe)
    #    cersigma2 = cerFunc(cerpe, c0, c1, c2)**2
    #    sctsigma2 = sctFunc(sctpe, s0)**2
    #    #totsigma2 = (sctsigma2+cersigma2) / (1-2*pp*(1-pp))
    #    totsigma2 = (2-pp)/pp*sctsigma2 + cersigma2
    #  
    #    best.append(np.sqrt(totsigma2))

        #sctpe = sctpe_sim(i, m_kA, m_kB, m_scale)
        #cerpe = cerpe_sim(i, m_kC)
        #pp = sctpe / (sctpe + cerpe)
        #cersigma2 = cerFunc(cerpe, m_c0, m_c1, m_c2)**2
        #sctsigma2 = sctFunc(sctpe, m_s0)**2
        #totsigma2 = (sctsigma2+cersigma2) / (1-2*pp*(1-pp))
        

        
            #totsigma.append(totsigma)
        #best.append(resFunc(i, q0, q1, q2)) 
        
        #best_posi.append( np.sqrt(best[-1]**2 + 2*27.07**2) )


    nn = 0
    for i in npe:
        best.append( resFunc(i, q01, q11, q21))
        best_posi1.append( np.sqrt(resFunc(i, q01, q11, q21)**2 + 2*27.07**2 ) )
        best_posi2.append( np.sqrt(resFunc(i, q02, q12, q22)**2 + 2*27.07**2 ) )

        diff.append( (best_posi1[nn] - best_posi2[nn])/best_posi1[nn] )
        err1 = ymax1[nn] - best_posi1[nn]
        err2 = ymax2[nn] - best_posi2[nn]
        diff_err.append( np.sqrt(err2**2/best_posi1[nn]**2 + err1**2*best_posi2[nn]**2/best_posi1[nn]**4) )
        
        nn += 1


    best = np.array(best)
    best_posi1 = np.array(best_posi1)
    best_posi2 = np.array(best_posi2)
    diff = np.array(diff)
    diff_err = np.array(diff_err)

    print(len(diff))
    print(len(diff_err))
    print(len(Etrue_posi))

    fig = plt.figure(figsize=(8, 6))
    spec = gridspec.GridSpec(ncols=1, nrows=2,
                         height_ratios=[1, 2])

    ax0 = fig.add_subplot(spec[0])
    ax1 = fig.add_subplot(spec[1])

    ax0.fill_between(npe_posi/global_es, diff+diff_err, diff-diff_err, alpha=0.3, color="#E49D22")
    ax0.plot(npe_posi/global_es, diff, "-", ms=4, color="#E49D22")
    ax0.grid(True)
    ax0.set_ylabel("relaive difference", fontsize=14)
    ax0.tick_params(axis='both', which='major', labelsize=13)
    ax0.grid(True)

    #plt.plot(Etrue, nominal(), "--", color="darkviolet", label="nominal")
    ax1.plot(npe_posi/global_es, best_posi1/npe_posi, "-", color="blue", label="simulation based")
    ax1.fill_between(npe_posi/global_es, ymin1/npe_posi, ymax1/npe_posi, color="lightskyblue")
    ax1.plot(npe_posi/global_es, best_posi2/npe_posi, "--", color="coral", label="analytical")
    ax1.fill_between(npe_posi/global_es, ymin2/npe_posi, ymax2/npe_posi, color="pink", alpha=0.5)
    ax1.set_xlabel(r"$E_{vis}$/MeV", fontsize=15)
    ax1.set_ylabel(r"$\sigma/N_{tot}$", fontsize=15)
    ax1.legend(prop={"size":13})
    ax1.tick_params(axis='both', which='major', labelsize=13)
    ax1.grid(True)

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

    plt.savefig("compare_positron_res.pdf")
    plt.show()



if __name__ == "__main__":
    main()
