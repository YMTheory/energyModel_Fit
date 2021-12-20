import numpy as np
import matplotlib.pyplot as plt
import ROOT
import elecLoader as el
from matplotlib import gridspec


def sigma2Func(a, b, c, A, x):
    s2 = a**2/A*x + b**2*x**2 + c**2/A**2
    return s2 if s2>0  else 0

def sigma2FuncNew(a, b, n, x):
    s2 = a**2*x + b**2*np.power(x, n)
    return s2 if s2>0  else 0


def loadTruth():
    resolE1, resolData1, resolerr1 = [], [], []
    with open("../data/electron/elecResol4.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if float(data[0]) < 0.2:
                continue
            resolE1.append(float(data[1]))
            resolData1.append(float(data[5]))
            resolerr1.append(float(data[6]))

    resolE1    = np.array(resolE1)
    resolData1 = np.array(resolData1)
    resolerr1  = np.array(resolerr1)

    return resolE1, resolData1, resolerr1

def resFunc(x, a, b, c):
    s2 = a + b*x + c*x**2
    if s2<0:
        return 0
    else:
        return np.sqrt(s2)

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


if __name__ == "__main__":

    Y = 3134.078 / 2.223


    filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/Newgam_kSimQ_kNewAnaCer_kNew/Newgam_kSimQ_kNewCerAna_kNew_rescov.txt"
    filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12_kSimQ_kNewAnaCer_kNew/NewgamB12_kSimQ_kNewCerAna_kNew_rescov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12NewMicedge_kSimQ_kNewAnaCer_kNew/NewgamB12NewMicedge_kSimQ_kNewCerAna_kNew_rescov.txt"
    filename3 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12NewMicwhole_kSimQ_kNewAnaCer_kNew/NewgamB12NewMic_kSimQ_kNewCerAna_kNew_rescov.txt"


    par1, parerr1 = loadBestFit(filename1, 3)
    a1, b1, n1 = par1
    print(par1)
    print(parerr1)

    par2, parerr2 = loadBestFit(filename2, 3)
    a2, b2, n2 = par2
    print(par2)
    print(parerr2)

    par3, parerr3 = loadBestFit(filename3, 3)
    a3, b3, n3 = par3
    print(par3)
    print(parerr3)

    Evis = np.arange(0.1, 9, 0.1)
    Evis_posi = []
    res1, res2, res3 = [], [], []


    kesim, musim, muerrsim, sigmasim, sigmaerrsim = loadSimTruth()

    q03, q13, q23 = 0, 1.00710, 3.98696e-05
    simData = []
    diff1, diff2, diff3 = [], [], []

    diff13 = []
    for i in Evis:
        ## Best Fit
        Evistot = i*Y 

        sigma_elec1 = sigma2FuncNew(a1, b1, n1, Y*i)
        sigma_elec2 = sigma2FuncNew(a2, b2, n2, Y*i)
        sigma_elec3 = sigma2FuncNew(a3, b3, n3, Y*i)

        Evistot = i*Y + 660.8 * 2 
        Evis_posi.append(Evistot/Y)

        simData.append(resFunc(Evistot, q03, q13, q23) / Evistot)

        res1.append(np.sqrt(sigma_elec1 + (2*27.07**2)) / Evistot )
        res2.append(np.sqrt(sigma_elec2 + (2*27.07**2)) / Evistot )
        res3.append(np.sqrt(sigma_elec3 + (2*27.07**2)) / Evistot )

        diff1.append((res1[-1] - simData[-1])/simData[-1])
        diff2.append((res2[-1] - simData[-1])/simData[-1])
        diff3.append((res3[-1] - simData[-1])/simData[-1])


        #res1.append(np.sqrt(sigma_elec1) / Evistot )
        #res2.append(np.sqrt(sigma_elec2) / Evistot )
        #res3.append(np.sqrt(sigma_elec3) / Evistot )

        diff13.append((res1[-1] - res3[-1])/res3[-1])


    simData = np.array(simData)
    diff13 = np.array(diff13)


    sampleSize = 5000
    sigma1 = sample_corelation(filename1, 3, sampleSize)
    sigma2 = sample_corelation(filename2, 3, sampleSize)
    sigma3 = sample_corelation(filename3, 3, sampleSize)

    ymin1, ymax1, ymin2, ymax2, ymin3, ymax3 = [], [], [], [], [], []
    for i in range(len(Evis)):
        ymin1.append(1000000)
        ymax1.append(-100)
        ymin2.append(1000000)
        ymax2.append(-100)
        ymin3.append(1000000)
        ymax3.append(-100)


    for i in range(sampleSize):

        m_a1 = a1 + sigma1[0, i]
        m_b1 = b1 + sigma1[1, i]
        m_n1 = n1 + sigma1[2, i]

        m_a2 = a2 + sigma2[0, i]
        m_b2 = b2 + sigma2[1, i]
        m_n2 = n2 + sigma2[2, i]

        m_a3 = a3 + sigma3[0, i]
        m_b3 = b3 + sigma3[1, i]
        m_n3 = n3 + sigma3[2, i]

        tmp_res1, tmp_res2, tmp_res3 = [], [], []
        for j in Evis:
            #Evistot = j*Y
            Evistot = j*Y + 660.8 * 2 

            sigma_elec1 = sigma2FuncNew(m_a1, m_b1, m_n1, Y*j)
            sigma_elec2 = sigma2FuncNew(m_a2, m_b2, m_n2, Y*j)
            sigma_elec3 = sigma2FuncNew(m_a3, m_b3, m_n3, Y*j)

            tmp_res1.append(np.sqrt(sigma_elec1 + (2*27.07**2)) / Evistot )
            tmp_res2.append(np.sqrt(sigma_elec2 + (2*27.07**2)) / Evistot )
            tmp_res3.append(np.sqrt(sigma_elec3 + (2*27.07**2)) / Evistot )
            #tmp_res1.append(np.sqrt(sigma_elec1) / Evistot )
            #tmp_res2.append(np.sqrt(sigma_elec2) / Evistot )
            #tmp_res3.append(np.sqrt(sigma_elec3) / Evistot )

        for j in range(len(Evis)):
            if tmp_res1[j] > ymax1[j]:
                ymax1[j] = tmp_res1[j]
            if tmp_res1[j] < ymin1[j]:
                ymin1[j] = tmp_res1[j]
            if tmp_res2[j] > ymax2[j]:
                ymax2[j] = tmp_res2[j]
            if tmp_res2[j] < ymin2[j]:
                ymin2[j] = tmp_res2[j]
            if tmp_res3[j] > ymax3[j]:
                ymax3[j] = tmp_res3[j]
            if tmp_res3[j] < ymin3[j]:
                ymin3[j] = tmp_res3[j]

    Evis = np.array(Evis)
    res1 = np.array(res1)
    ymin1 = np.array(ymin1)
    ymax1 = np.array(ymax1)
    res2 = np.array(res2)
    ymin2 = np.array(ymin2)
    ymax2 = np.array(ymax2)
    res3 = np.array(res3)
    ymin3 = np.array(ymin3)
    ymax3 = np.array(ymax3)

    del1min = res1 - ymin1
    del1max = ymax1 - res1
    del2min = res2 - ymin2
    del2max = ymax2 - res2
    del3min = res3 - ymin3
    del3max = ymax3 - res3



    #fig, ax = plt.subplots()
    fig = plt.figure(figsize=(6, 8))
    spec = gridspec.GridSpec(ncols=1, nrows=2,
                         height_ratios=[1, 2])

    ax = fig.add_subplot(spec[1])
    ax0 = fig.add_subplot(spec[0])

    #ax0.plot(elecE/Y, diffe1, "o-", ms=6, lw=2, color="crimson", label="only gamma")
    #ax0.plot(elecE/Y, diffe2, "o-", ms=6, lw=2, color="black"  , label="gamma + B12")
    #ax0.set_xlabel(r"Electron  $E^{vis}$ [MeV]", fontsize=15)
    ax0.plot(Evis_posi, diff1, "-", color="crimson", label="only gamma")
    ax0.plot(Evis_posi, diff2, "-", color="black", label="gamma + B12")
    ax0.fill_between(Evis_posi, (res1-1*del1min-simData)/simData, (res1+1*del1max-simData)/simData, color="crimson", alpha=0.3)
    ax0.fill_between(Evis_posi, (res2-1*del2min-simData)/simData, (res2+1*del2max-simData)/simData, color="slategray", alpha=0.3)
    ax0.set_ylabel("Bias", fontsize=15, color="black")
    ax0.tick_params(axis='both', which='major', labelsize=14, labelcolor="black")
    ax0.set_ylim(-0.05, 0.05)
    ax0.grid(True)

    ax.errorbar(musim/Y, sigmasim/musim, yerr=np.sqrt(sigmaerrsim**2/musim**2+muerrsim**2*sigmasim**2/musim**4), fmt="o", color="blue", ms=6, mfc="w", label="Simulation truth", zorder=1)
    #ax.errorbar(elecE/Y, elecRes, yerr=elecReserr, fmt="o", ms=6, color="blue", label="Simulation", zorder=2)
    ax.plot(Evis_posi, res1, "--", lw=2, color="crimson", label="only gamma", zorder=3)
    ax.plot(Evis_posi, res2, "-", lw=2, color="black", label="gamma + B12", zorder=2)
    #ax.plot(Evis, res3, "-", lw=2, color="royalblue", label="gamma + B12 + Michel", zorder=3)
    ax.fill_between(Evis_posi, res1-10*del1min, res1+10*del1max, color="crimson", alpha=0.3)
    ax.fill_between(Evis_posi, res2-10*del2min, res2+10*del2max, color="slategray", alpha=0.3)
    #ax.set_xlim(0.9, 2)
    #ax.set_ylim(0.02, 0.03)

    ax.legend(prop={"size" : 15})
    ax.set_xlabel(r"Positron $E^{vis}$ [MeV]", fontsize=15)
    ax.set_ylabel(r"$\sigma/E^{vis}$", fontsize=15, color="black")
    ax.tick_params(axis='both', which='major', labelsize=14, labelcolor="black")
    ax.grid(True)
    #ax.semilogy()

    plt.tight_layout()

    #plt.savefig("posiResol_GamB12.pdf")

    plt.show()










