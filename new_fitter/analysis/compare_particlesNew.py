import numpy as np
import matplotlib.pyplot as plt
import elecLoader as el
from matplotlib import gridspec
import random
import ROOT


# no correlation in fitter...

Y = 3134.078/2.223

def sigma2FuncNew(a, b, n, x):
    s2 = a**2*x + b**2*np.power(x, n)
    return s2 if s2>0  else 0

def wenNcer(E, p0, p1, p2, p3, p4):
    if E<=0.2:
        return 0
    if E>0.2:
        E = E - 0.2
        return p3*E**2/(p4+p0*E+p1*np.exp(-p2*E))



def fitNsct(E, As, kB):
    if E > 15.9:
        return el.getFitNsct(15.9, kB, As, "Sim") / 15.9 * E
    else:
        return el.getFitNsct(E, kB, As, "Sim")

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
                print(data[0], gamE[-1], gammu[-1]/Y/gamE[-1])

    gamE = np.array(gamE)
    gammu = np.array(gammu)
    gammuerr = np.array(gammuerr)
    gamsigma = np.array(gamsigma)
    gamsigmamu = np.array(gamsigmamu)

    return gamE, gammu, gammuerr, gamsigma, gamsigmamu



def loadGe68():
    gamE, gammu, gammuerr, gamsigma, gamsigmamu = 0, 0, 0, 0, 0
    with open("../data/gamma/gamma_J19.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if data[0] == "Ge68":
                gamE = (float(data[1]))
                gammu = (float(data[2]))
                gammuerr = (float(data[3]))
                gamsigma = (float(data[4]))
                gamsigmamu = (float(data[5]))


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



def main():

    fig = plt.figure(figsize=(12, 5))
    spec = gridspec.GridSpec(ncols=2, nrows=1)

    ax0 = fig.add_subplot(spec[0])
    ax1 = fig.add_subplot(spec[1])


    ## nonlinearity part
    Edep = np.arange(0.05, 9, 0.1)
    Edep_posi = Edep + 1.022

    dnum = len(Edep)

    filename4 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamNewB12_kSimQ_kNewAnaCer_kNew/NewgamNewB12_kSimQ_kNewAnaCer_kNew_nonlcov.txt"

    par1, parerr1 = loadBestFit(filename4, 7)
    scale1, kB1, p01, p11, p21, p31, p41 = par1

    nonl_elec, nonl_posi = [], []
    for i, j in zip(Edep, Edep_posi):
        nonl_elec.append((fitNsct(i, scale1, kB1) + wenNcer(i, p01, p11, p21, p31, p41) ) / i / Y )
        nonl_posi.append((fitNsct(i, scale1, kB1) + wenNcer(i, p01, p11, p21, p31, p41) + 2*660.8) / j / Y )

    sampleSize = 5000
    sigma1 = sample_corelation(filename4, 7, sampleSize)
    ymin1, ymax1, ymin2, ymax2 = [], [] ,[], []
    for i in range(dnum):
        ymin1.append(1000000)
        ymax1.append(-100)
        ymin2.append(1000000)
        ymax2.append(-100)
    for i in range(sampleSize):
        m_scale1 = scale1 + sigma1[0, i]
        m_kB1    = kB1 + sigma1[1, i]
        m_p01    = p01 + sigma1[2, i]
        m_p11    = p11 + sigma1[3, i]
        m_p21    = p21 + sigma1[4, i]
        m_p31    = p31 + sigma1[5, i]
        m_p41    = p41 + sigma1[6, i]


        tmp1, tmp2 = [], []
        for k, j in zip(Edep, Edep_posi):
            tmp1.append((fitNsct(k, m_scale1, m_kB1) + wenNcer(k, m_p01, m_p11, m_p21, m_p31, m_p41) ) / k / Y )
            tmp2.append((fitNsct(k, m_scale1, m_kB1) + wenNcer(k, m_p01, m_p11, m_p21, m_p31, m_p41) + 2*660.8) / j / Y )

        for j in range(len(Edep)):
            if tmp1[j] > ymax1[j]:
                ymax1[j] = tmp1[j]
                #print("ymax1" , ymax1[j])
            if tmp1[j] < ymin1[j]:
                ymin1[j] = tmp1[j]
                #print("ymin1" , ymin1[j])
            if tmp2[j] > ymax2[j]:
                ymax2[j] = tmp2[j]
            if tmp2[j] < ymin2[j]:
                ymin2[j] = tmp2[j]


    ymin1 = np.array(ymin1)
    ymax1 = np.array(ymax1)
    ymin2 = np.array(ymin2)
    ymax2 = np.array(ymax2)

    kesim, musim, muerrsim, sigmasim, sigmaerrsim = loadSimTruth()
    gamE, gammu, gammuerr, gamsigma, gamsigmaerr = loadGamTruth()
    Ge68E, Ge68mu, Ge68muerr, Ge68sigma, Ge68sigmaerr = loadGe68()


    filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamNewB12_kSimQ_kNewAnaCer_kNew/NewgamNewB12_kSimQ_kNewAnaCer_kNew_rescov.txt"
    par1, parerr1 = loadBestFit(filename1, 3)
    a1, b1, n1 = par1
    
    Evis = np.arange(0.5, 9, 0.1)
    Evis_posi = []
    res_elec, res_posi = [], []
    for i in Evis:
        sigma_elec1 = sigma2FuncNew(a1, b1, n1, Y*i)
        res_elec.append(np.sqrt(sigma_elec1) / Y / i)
        Evistot = (i-0.5)*Y + 660.8 * 2 
        sigma_elec2 = sigma2FuncNew(a1, b1, n1, Y*(i-0.5))
        Evis_posi.append(Evistot/Y)
        res_posi.append(np.sqrt(sigma_elec2 + (2*27.07**2)) / Evistot )
    
    dnum = len(Evis)

    sigma1 = sample_corelation(filename1, 3, sampleSize)
    ymin3, ymax3, ymin4, ymax4 = [], [] ,[], []
    for k in range(dnum):
        ymin3.append(1000000)
        ymax3.append(-100)
        ymin4.append(1000000)
        ymax4.append(-100)
    for k in range(sampleSize):
        m_a1 = a1 + sigma1[0, k]
        m_b1 = b1 + sigma1[1, k]
        m_n1 = n1 + sigma1[2, k]

        tmp1, tmp2 = [], []
        for i in Evis:
            sigma_elec1 = sigma2FuncNew(m_a1, m_b1, m_n1, Y*i)
            tmp1.append(np.sqrt(sigma_elec1) / Y / i)
            Evistot = (i-0.5)*Y + 660.8 * 2 
            sigma_elec2 = sigma2FuncNew(m_a1, m_b1, m_n1, Y*(i-0.5))
            tmp2.append(np.sqrt(sigma_elec2 + (2*27.07**2)) / Evistot )

        for j in range(dnum):
            if tmp1[j] > ymax3[j]:
                ymax3[j] = tmp1[j]
                #print("ymax1" , ymax1[j])
            if tmp1[j] < ymin3[j]:
                ymin3[j] = tmp1[j]
                #print("ymin3" , ymin3[j])
            if tmp2[j] > ymax4[j]:
                ymax4[j] = tmp2[j]
            if tmp2[j] < ymin4[j]:
                ymin4[j] = tmp2[j]


    ymin3 = np.array(ymin3)
    ymax3 = np.array(ymax3)
    ymin4 = np.array(ymin4)
    ymax4 = np.array(ymax4)

    gam_name = [r"$^{137}$Cs", r"$^{54}$Mn", r"$^{40}$K", "n-H", "AmBe", "AmC"]
    nonlx = [gamE[0]-0.95, gamE[1]-0.45, gamE[2]-0.45, gamE[3]-0.1, gamE[4]-0.1, gamE[5]-0.1 ]
    nonly = [gammu[0]/Y/gamE[0]-0.02, gammu[1]/Y/gamE[1]+0.01, gammu[2]/Y/gamE[2]+0.01, gammu[3]/Y/gamE[3]+0.006, gammu[4]/Y/gamE[4]+0.006, gammu[5]/Y/gamE[5]+0.005 ]
    resx  = [-0.95, -0.8, gammu[2]/Y+0.1, gammu[3]/Y, gammu[4]/Y, gammu[5]/Y ]
    resy  = [gamsigma[0]/gammu[0], gamsigma[1]/gammu[1], gamsigma[2]/gammu[2]-0.001, gamsigma[3]/gammu[3]+0.001, gamsigma[4]/gammu[4]+0.001, gamsigma[5]/gammu[5]+0.001 ]
    

    ax0.plot(Edep, nonl_elec, color="blue", lw=2, label=r"$e^-$ with $1\sigma$ errorbar $\times 10$")
    ax0.fill_between(Edep, nonl_elec- 10*(nonl_elec - ymin1), nonl_elec + 10*(ymax1-nonl_elec), color="royalblue", alpha=0.5)
    ax0.plot(Edep_posi, nonl_posi, color="crimson", lw=2, label=r"$e^+$ with $1\sigma$ errorbar $\times 10$")
    ax0.fill_between(Edep_posi, nonl_posi- 10*(nonl_posi - ymin2), nonl_posi + 10*(ymax2-nonl_posi), color="crimson", alpha=0.5)
    ax0.errorbar(gamE, gammu/Y/gamE, yerr=10*gammuerr/Y/gamE,fmt="o--", lw=2, ms=6, color="dimgray", label=r"Single $\gamma$ with $1\sigma$ errorbars $\times 10$ ")
    for i, j in enumerate(gam_name):
        ax0.text(nonlx[i], nonly[i], j, fontsize=14, color="dimgray")
    ax0.errorbar(Ge68E, Ge68mu/Y/Ge68E, yerr=10*Ge68muerr/Y/Ge68E,fmt="*", lw=2, ms=7, color="black", label=r"$^{68}$Ge with $1\sigma$ errorbars $\times 10$ ")
    ax0.text(Ge68E+0.20, Ge68mu/Y/Ge68E-0.004, r"$^{68}$Ge", fontsize=14, color="black")
    ax0.set_xlabel(r"$E^{dep}$ [MeV]", fontsize=17)
    ax0.set_ylabel(r"$E^{vis}/E^{dep}$", fontsize=17)
    ax0.set_ylim(0.80, 1.06)
    ax0.tick_params(axis='both', which='major', labelsize=16)
    ax0.legend(loc="lower right", prop={"size":15})
    ax0.grid(True)
    ax0.set_title("(a)", y =-0.3, fontsize=17)

    ax1.plot(Evis, res_elec, color="blue", lw=2, label=r"$e^-$ with $1\sigma$ errorbar $\times 5$")
    ax1.plot(Evis_posi, res_posi, color="crimson", lw=2, label=r"$e^+$ with $1\sigma$ errorbar $\times 5$")
    ax1.fill_between(Evis, res_elec- 5*(res_elec - ymin3), res_elec + 5*(ymax3-res_elec), color="royalblue", alpha=0.5)
    ax1.fill_between(Evis_posi, res_posi- 5*(res_posi - ymin4), res_posi + 5*(ymax4-res_posi), color="crimson", alpha=0.5)
    ax1.errorbar(gammu/Y, gamsigma/gammu, yerr=5*np.sqrt(gammuerr**2*gamsigma**2/gammu**4 + gamsigmaerr**2/gammu**2), fmt="o--", lw=2, ms=6, color="dimgray", label=r"Single $\gamma$ with $1\sigma$ errorbars $\times$ 5")
    for i, j in enumerate(gam_name):
        plt.text(resx[i], resy[i], j, fontsize=14, color="dimgray")
    ax1.errorbar(Ge68mu/Y, Ge68sigma/Ge68mu, yerr=5*np.sqrt(Ge68muerr**2*Ge68sigma**2/Ge68mu**4 + Ge68sigmaerr**2/Ge68mu**2), fmt="*", lw=2, ms=7, color="black", label=r"$^{68}$Ge with $1\sigma$ errorbars $\times$ 5")
    ax1.text(-0.4, Ge68sigma/Ge68mu-0.003, r"$^{68}$Ge", fontsize=14, color="black")
    ax1.set_xlabel(r"$E^{vis}$ [MeV]", fontsize=17)
    ax1.set_ylabel(r"$\sigma / E^{vis}$", fontsize=17)
    ax1.tick_params(axis='both', which='major', labelsize=16)
    ax1.legend(prop={"size":15})
    ax1.grid(True)
    ax1.set_xlim(-1, 12)
    ax1.set_ylim(0.008, 0.050)
    ax1.set_title("(b)", y =-0.3, fontsize=17)

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=0.02, hspace=0.02)
    plt.tight_layout()
    plt.savefig("compare_particles.pdf")
    plt.show()



if __name__ == "__main__":
    main()
