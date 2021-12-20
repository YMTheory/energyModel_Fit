import numpy as np
import matplotlib.pyplot as plt
import ROOT
import elecLoader as el
from matplotlib import gridspec


def wenNcer(E, p0, p1, p2, p3, p4):
    if E<=0.2:
        return 0
    if E>0.2:
        E = E - 0.2
        return p3*E**2/(p4+p0*E+p1*np.exp(-p2*E))

def sigma2FuncNew(a, b, n, x):
    s2 = a**2*x + b**2*np.power(x, n)
    return s2 if s2>0  else 0



def fitNsct(E, As, kB):
    if E > 15.9:
        return el.getFitNsct(15.9, kB, As, "Sim") / 15.9 * E
    else:
        return el.getFitNsct(E, kB, As, "Sim")



def sigma2Func(a, b, c, A, x):
    s2 = a**2/A*x + b**2*x**2 + c**2/A**2
    return s2 if s2>0  else 0


def resFunc(x, a, b, c):
    s2 = a + b*x + c*x**2
    if s2<0:
        return 0
    else:
        return np.sqrt(s2)


def loadTruth():
    Y = 3134.078 / 2.223
    elecE, elecPE, elecPEerr = [], [], []
    with open("../data/electron/elecResol4.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if float(data[0]) > 0.1:
                elecE.append(float(data[0]))
                elecPE.append(float(data[1]) / Y)
                elecPEerr.append(float(data[2]) / Y)

    elecE = np.array(elecE)
    elecPE = np.array(elecPE)
    elecPEerr = np.array(elecPEerr)

    elecNonl = elecPE / elecE
    elecNonlerr = elecPEerr / elecE

    return elecE, elecNonl, elecNonlerr


def loadSimTruth():
    filename = '/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/positron/positronNPE1.txt'
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

    elecE, elecNonl, elecNonlerr = loadTruth()

    filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/Newgam_kSimQ_kNewAnaCer_kNew/Newgam_kSimQ_kNewCerAna_kNew_nonlcov.txt"
    filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamNewB12_kSimQ_kNewAnaCer_kNew/NewgamNewB12_kSimQ_kNewAnaCer_kNew_nonlcov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12NewMicwhole_kSimQ_kNewAnaCer_kNew/NewgamB12NewMic_kSimQ_kNewAnaCer_kNew_nonlcov.txt"
    filename3 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamNewB12NewMic_kSimQ_kNewAnaCer_kNew/NewgamNewB12NewMic_kSimQ_kNewAnaCer_kNew_nonlcov.txt"
    filename3 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamNewB12_kSimQ_kNewAnaCer_kNew/NewgamNewB12_kSimQ_kNewAnaCer_kNew_nonlcov.txt"

    par1, parerr1 = loadBestFit(filename1, 7)
    scale1, kB1, p01, p11, p21, p31, p41 = par1

    par2, parerr2 = loadBestFit(filename2, 7)
    scale2, kB2, p02, p12, p22, p32, p42 = par2

    par3, parerr3 = loadBestFit(filename3, 7)
    scale3, kB3, p03, p13, p23, p33, p43 = par3

    kesim, musim, muerrsim, sigmasim, sigmaerrsim = loadSimTruth()

    Edep = np.arange(0.05, 9, 0.1)
    Edep_posi = Edep + 1.022
    nonl1, nonl2, nonl3 = [], [], []

    for i, j in zip(Edep, Edep_posi):
        nonl1.append((fitNsct(i, scale1, kB1) + wenNcer(i, p01, p11, p21, p31, p41) + 2*660.8) / j / Y )
        nonl2.append((fitNsct(i, scale2, kB2) + wenNcer(i, p02, p12, p22, p32, p42) + 2*660.8) / j / Y )
        nonl3.append((fitNsct(i, scale3, kB3) + wenNcer(i, p03, p13, p23, p33, p43) + 2*660.8) / j / Y )

    sampleSize = 5000
    sigma1 = sample_corelation(filename1, 7, sampleSize)
    sigma2 = sample_corelation(filename2, 7, sampleSize)
    sigma3 = sample_corelation(filename3, 7, sampleSize)

    diffmin1, diffmax1, diffmin2, diffmax2 = [], [], [], []
    for i in range(len(kesim)):
        diffmin1.append(1000000)
        diffmax1.append(-100)
        diffmin2.append(1000000)
        diffmax2.append(-100)
    for i in range(sampleSize):

        m_scale1 = scale1 + sigma1[0, i]
        m_kB1    = kB1 + sigma1[1, i]
        m_p01    = p01 + sigma1[2, i]
        m_p11    = p11 + sigma1[3, i]
        m_p21    = p21 + sigma1[4, i]
        m_p31    = p31 + sigma1[5, i]
        m_p41    = p41 + sigma1[6, i]

        m_scale2 = scale2 + sigma2[0, i]
        m_kB2    = kB2 + sigma2[1, i]
        m_p02    = p02 + sigma2[2, i]
        m_p12    = p12 + sigma2[3, i]
        m_p22    = p22 + sigma2[4, i]
        m_p32    = p32 + sigma2[5, i]
        m_p42    = p42 + sigma2[6, i]
        tmp_nonl1, tmp_nonl2, tmp_nonl3 = [], [], []
        for k, j in zip(Edep_posi, Edep):
            tmp_nonl1.append((fitNsct(j, m_scale1, m_kB1) + wenNcer(j, m_p01, m_p11, m_p21, m_p31, m_p41) + 2*660.8) / k / Y )
            tmp_nonl2.append((fitNsct(j, m_scale2, m_kB2) + wenNcer(j, m_p02, m_p12, m_p22, m_p32, m_p42) + 2*660.8) / k / Y )
        for j in range(len(kesim)):
            if tmp_nonl1[j] > diffmax1[j]:
                diffmax1[j] = tmp_nonl1[j]
            if tmp_nonl1[j] < diffmin1[j]:
                diffmin1[j] = tmp_nonl1[j]
            if tmp_nonl2[j] > diffmax2[j]:
                diffmax2[j] = tmp_nonl2[j]
            if tmp_nonl2[j] < diffmin2[j]:
                diffmin2[j] = tmp_nonl2[j]

    diffmin1 = np.array(diffmin1)
    diffmax1 = np.array(diffmax1)
    diffmin2 = np.array(diffmin2)
    diffmax2 = np.array(diffmax2)

    err1, err2 = [], []
    for i in range(len(kesim)):
        del1min = nonl1[i] - diffmin1[i]
        del1max = diffmax1[i] - nonl1[i]
        del1 = del1min > del1max and del1min or del1max
        err1.append(del1)
        del2min = nonl2[i] - diffmin2[i]
        del2max = diffmax2[i] - nonl2[i]
        del2 = del2min > del2max and del2min or del2max
        err2.append(del2)
    err2 = np.array(err2)


    ### Comparing with MC truth :
    diff1, diff2, diff3 = [], [], []
    for i in range(len(musim)):
        ke = kesim[i]
        te = ke + 1.022
        tru = musim[i] / Y / te
        nl1 = ((fitNsct(ke, scale1, kB1) + wenNcer(ke, p01, p11, p21, p31, p41) + 2*660.8) / te / Y )
        nl2 = ((fitNsct(ke, scale2, kB2) + wenNcer(ke, p02, p12, p22, p32, p42) + 2*660.8) / te / Y )
        nl3 = ((fitNsct(ke, scale3, kB3) + wenNcer(ke, p03, p13, p23, p33, p43) + 2*660.8) / te / Y )
        
        diff1.append(nl1)
        diff2.append(nl2)
        diff3.append(nl3)



    ymin1, ymax1, ymin2, ymax2, ymin3, ymax3 = [], [], [], [], [], []
    for i in range(len(Edep)):
        ymin1.append(1000000)
        ymax1.append(-100)
        ymin2.append(1000000)
        ymax2.append(-100)
        ymin3.append(1000000)
        ymax3.append(-100)

    for i in range(sampleSize):

        m_scale1 = scale1 + sigma1[0, i]
        m_kB1    = kB1 + sigma1[1, i]
        m_p01    = p01 + sigma1[2, i]
        m_p11    = p11 + sigma1[3, i]
        m_p21    = p21 + sigma1[4, i]
        m_p31    = p31 + sigma1[5, i]
        m_p41    = p41 + sigma1[6, i]

        m_scale2 = scale2 + sigma2[0, i]
        m_kB2    = kB2 + sigma2[1, i]
        m_p02    = p02 + sigma2[2, i]
        m_p12    = p12 + sigma2[3, i]
        m_p22    = p22 + sigma2[4, i]
        m_p32    = p32 + sigma2[5, i]
        m_p42    = p42 + sigma2[6, i]

        m_scale3 = scale3 + sigma3[0, i]
        m_kB3    = kB3 + sigma3[1, i]
        m_p03    = p03 + sigma3[2, i]
        m_p13    = p13 + sigma3[3, i]
        m_p23    = p23 + sigma3[4, i]
        m_p33    = p33 + sigma3[5, i]
        m_p43    = p43 + sigma3[6, i]


        tmp_nonl1, tmp_nonl2, tmp_nonl3 = [], [], []
        for k, j in zip(Edep_posi, Edep):
            tmp_nonl1.append((fitNsct(j, m_scale1, m_kB1) + wenNcer(j, m_p01, m_p11, m_p21, m_p31, m_p41) + 2*660.8) / k / Y )
            tmp_nonl2.append((fitNsct(j, m_scale2, m_kB2) + wenNcer(j, m_p02, m_p12, m_p22, m_p32, m_p42) + 2*660.8) / k / Y )
            tmp_nonl3.append((fitNsct(j, m_scale3, m_kB3) + wenNcer(j, m_p03, m_p13, m_p23, m_p33, m_p43) + 2*660.8) / k / Y )
    
        for j in range(len(Edep)):
            if tmp_nonl1[j] > ymax1[j]:
                ymax1[j] = tmp_nonl1[j]
                #print("ymax1" , ymax1[j])
            if tmp_nonl1[j] < ymin1[j]:
                ymin1[j] = tmp_nonl1[j]
                #print("ymin1" , ymin1[j])
            if tmp_nonl2[j] > ymax2[j]:
                ymax2[j] = tmp_nonl2[j]
            if tmp_nonl2[j] < ymin2[j]:
                ymin2[j] = tmp_nonl2[j]
            if tmp_nonl3[j] > ymax3[j]:
                ymax3[j] = tmp_nonl3[j]
            if tmp_nonl3[j] < ymin3[j]:
                ymin3[j] = tmp_nonl3[j]

    Edep = np.array(Edep)
    nonl1 = np.array(nonl1)
    ymin1 = np.array(ymin1)
    ymax1 = np.array(ymax1)
    nonl2 = np.array(nonl2)
    ymin2 = np.array(ymin2)
    ymax2 = np.array(ymax2)
    nonl3 = np.array(nonl3)
    ymin3 = np.array(ymin3)
    ymax3 = np.array(ymax3)


    # Plotting

    #fig, ax = plt.subplots()
    fig = plt.figure(figsize=(12, 6))
    spec = gridspec.GridSpec(ncols=2, nrows=2,
                         height_ratios=[1, 2])

    ax1 = fig.add_subplot(spec[2])
    ax = fig.add_subplot(spec[0])
    ax2 = fig.add_subplot(spec[1])
    ax3 = fig.add_subplot(spec[3])

    nonlSim = musim / Y / (kesim+1.022)
    nonlSimErr = muerrsim / Y / (kesim+1.022)

    #ax.plot(kesim+1.022, (diff1-nonlSim)/nonlSim, "-", lw=2, color="crimson", label="only gamma")
    ax.plot(kesim+1.022, (diff2-nonlSim)/nonlSim*100, "--", lw=2, color="black")
    #ax.fill_between(kesim+1.022, (diff1-del1-nonlSim)/nonlSim, (diff1+del1-nonlSim)/nonlSim, color="crimson", alpha=0.3)
    ax.fill_between(kesim+1.022, (diff2-10*err2-nonlSim)/nonlSim*100, (diff2+10*err2-nonlSim)/nonlSim*100, color="slategrey", alpha=0.5, label=r"$10\times 1\sigma$ errorbar")
    #ax.errorbar(elecE, elecNonl, yerr=elecNonlerr, fmt="o", color="black", mfc="w", label="Simulation truth")
    #ax.errorbar(musim/Y, musim/Y/(kesim+1.022), yerr=muerrsim/Y/(kesim+1.022), fmt="o", color="black", mfc="w", label="Simulation")

    #ax.plot(Edep_posi, nonl1, "-", color="crimson",   label="gamma only")
    #ax.plot(Edep_posi, nonl2, "-", color="royalblue", label="Fitting+Michel")

    #ax.fill_between(Edep_posi, nonl1-10*(nonl1-ymin1), nonl1+10*(ymax1-nonl1), color="crimson", alpha=0.3)
    #ax.fill_between(Edep_posi, nonl2-10*(nonl2-ymin2), nonl2+10*(ymax2-nonl2), color="royalblue", alpha=0.3)

    ticknames = ["-0.6", "-0.4", "-0.2", "0", "0.2", "0.4", "0.6"]
    ax.set_yticks([-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6])
    ax.set_yticklabels(ticknames)
    ax.legend(loc="lower left", prop={"size":15})
    #ax.set_xlabel(r"$E^{dep}$ [MeV]", fontsize=16)
    ax.set_ylabel("Bias (%)", fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_ylim(-0.6, 0.6)
    ax.grid(True)


    #ax1.errorbar(elecE, elecNonl, yerr=elecNonlerr, fmt="o", color="black", mfc="w", label="Simulation truth")
    ax1.errorbar(kesim+1.022, musim/Y/(kesim+1.022), yerr=muerrsim/Y/(kesim+1.022), fmt="o", ms=6, color="crimson", mfc="w", label="Simulation", zorder=1)

    #ax1.plot(Edep_posi, nonl1, "-", color="crimson",   label="gamma only", zorder=2)
    ax1.plot(Edep_posi, nonl2, "-", lw=2, color="black", label=r"Fitting: $10\times1\sigma$ errorbar", zorder=3)

    #ax1.fill_between(Edep_posi, nonl1-20*(nonl1-ymin1), nonl1+20*(ymax1-nonl1), color="crimson", alpha=0.3)
    ax1.fill_between(Edep_posi, nonl2-10*(nonl2-ymin2), nonl3+10*(ymax2-nonl2), color="slategrey", alpha=0.5)

    ax1.legend(prop={"size":15})
    ax1.set_xlabel(r"Positron $E^{dep}$ [MeV]", fontsize=16)
    ax1.set_ylabel(r"$E^{vis}/E^{dep}$", fontsize=16)
    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax1.grid(True)
    ax1.set_title("(a)", y=-0.3, fontsize=15)



    filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamNewB12_kSimQ_kNewAnaCer_kNew/NewgamNewB12_kSimQ_kNewAnaCer_kNew_rescov.txt"
    par2, parerr2 = loadBestFit(filename2, 3)
    a2, b2, n2 = par2
    Evis = np.arange(0.1, 9, 0.1)
    Evis_posi = []
    res1, res2, res3 = [], [], []


    q03, q13, q23 = 0, 1.00710, 3.98696e-05
    simData = []

    diff2 = []
    for i in Evis:
        ## Best Fit
        Evistot = i*Y 

        sigma_elec2 = sigma2FuncNew(a2, b2, n2, Y*i)

        Evistot = i*Y + 660.8 * 2 
        Evis_posi.append(Evistot/Y)

        simData.append(resFunc(Evistot, q03, q13, q23) / Evistot)

        res2.append(np.sqrt(sigma_elec2 + (2*27.07**2)) / Evistot )

        diff2.append((res2[-1] - simData[-1])/simData[-1])

    simData = np.array(simData)
    diff2 = np.array(diff2)


    sampleSize = 5000
    sigma2 = sample_corelation(filename2, 3, sampleSize)

    ymin2, ymax2 = [], []
    for i in range(len(Evis)):
        ymin2.append(1000000)
        ymax2.append(-100)


    for i in range(sampleSize):

        m_a2 = a2 + sigma2[0, i]
        m_b2 = b2 + sigma2[1, i]
        m_n2 = n2 + sigma2[2, i]

        tmp_res2 = []
        for j in Evis:
            #Evistot = j*Y
            Evistot = j*Y + 660.8 * 2 

            sigma_elec2 = sigma2FuncNew(m_a2, m_b2, m_n2, Y*j)

            tmp_res2.append(np.sqrt(sigma_elec2 + (2*27.07**2)) / Evistot )

        for j in range(len(Evis)):
            if tmp_res2[j] > ymax2[j]:
                ymax2[j] = tmp_res2[j]
            if tmp_res2[j] < ymin2[j]:
                ymin2[j] = tmp_res2[j]

    Evis = np.array(Evis)
    res2 = np.array(res2)
    ymin2 = np.array(ymin2)
    ymax2 = np.array(ymax2)


    err2 = []
    for i in range(len(Evis)):
        del2min = res2[i] - ymin2[i]
        del2max = ymax2[i] - res2[i]
        del2 = del2min > del2max and del2min or del2max
        err2.append(del2)

    err2 = np.array(err2)


    ax2.plot(Evis_posi, diff2*100, "--",lw=2, color="black")
    ax2.fill_between(Evis_posi, (res2-err2-simData)/simData*100, (res2+err2-simData)/simData*100, color="slategray", alpha=0.5, label=r"$1\sigma$ errorbar")
    ax2.set_ylabel("Bias (%)", fontsize=15, color="black")
    ax2.tick_params(axis='both', which='major', labelsize=14, labelcolor="black")
    ax2.set_ylim(-6, 6)
    ax2.legend(loc="lower left", prop={"size" : 15})
    ticknames = ["-6", "-4", "-2", "0", "2", "4", "6"]
    ax2.set_yticks([-6, -4, -2, 0, 2, 4, 6])
    ax2.set_yticklabels(ticknames)
    ax2.grid(True)

    ax3.errorbar(musim/Y, sigmasim/musim, yerr=np.sqrt(sigmaerrsim**2/musim**2+muerrsim**2*sigmasim**2/musim**4), fmt="o", color="crimson", ms=6, mfc="w", label="Simulation", zorder=1)
    #ax.plot(Evis_posi, res1, "--", lw=2, color="crimson", label="only gamma", zorder=3)
    ax3.plot(Evis_posi, res2, "-", lw=2, color="black", label=r"Fitting:$10\times 1\sigma$ errorbar", zorder=2)
    #ax.plot(Evis, res3, "-", lw=2, color="royalblue", label="gamma + B12 + Michel", zorder=3)
    #ax.fill_between(Evis_posi, res1-10*del1min, res1+10*del1max, color="crimson", alpha=0.3)
    ax3.fill_between(Evis_posi, res2-err2*10, res2+10*err2, color="slategray", alpha=0.5)

    ax3.legend(prop={"size" : 15})
    ax3.set_xlabel(r"Positron $E^{vis}$ [MeV]", fontsize=15)
    ax3.set_ylabel(r"$\sigma/E^{vis}$", fontsize=15, color="black")
    ax3.tick_params(axis='both', which='major', labelsize=14, labelcolor="black")
    ax3.grid(True)
    #ax.semilogy()
    ax3.set_title("(b)", y=-0.3, fontsize=15)





    plt.subplots_adjust(left=None, bottom=None, right=None, top=None , wspace=None, hspace=None)

    plt.tight_layout()
    plt.savefig("compare_positron_nonlres.pdf")
    plt.show()
















