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



def fitNsct(E, As, kB):
    if E > 15.9:
        return el.getFitNsct(15.9, kB, As, "Sim") / 15.9 * E
    else:
        return el.getFitNsct(E, kB, As, "Sim")



def sigma2Func(a, b, c, A, x):
    s2 = a**2/A*x + b**2*x**2 + c**2/A**2
    return s2 if s2>0  else 0


def sigma2FuncNew(a, b, n, x):
    s2 = a**2*x + b**2*np.power(x, n)
    return s2 if s2>0  else 0


def loadTruth():
    Y = 3134.078 / 2.223
    elecE, elecPE, elecPEerr = [], [], []
    elecRes, elecReserr = [], []
    with open("../data/electron/elecResol4.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if float(data[0]) > 0.1:
                elecE.append(float(data[0]))
                elecPE.append(float(data[1]) / Y)
                elecPEerr.append(float(data[2]) / Y)
                elecRes.append(float(data[5]))
                elecReserr.append(float(data[6]))

    elecE = np.array(elecE)
    elecPE = np.array(elecPE)
    elecPEerr = np.array(elecPEerr)
    elecRes = np.array(elecRes)
    elecReserr = np.array(elecReserr)

    elecNonl = elecPE / elecE
    elecNonlerr = elecPEerr / elecE

    return elecE, elecNonl, elecNonlerr, elecRes, elecReserr





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

    elecE, elecNonl, elecNonlerr, elecRes, elecReserr = loadTruth()

    #filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/articles_gam_kSimQ_kNewAnaC1_kNPE/gam_kSimulation_kNewAnaC1_kNPE_nonlcov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/articles_gamB12Mic_kSimQ_kNewAnaC1_kNPE/gamB12Mic_kSimulation_kNewAnaCer_kNPE_nonlcov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/articles_gamB12_kSimQ_kNewAnaC1_kNPE/gamB12_kSimulation_kNewAnaC1_kNPE_nonlcov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/articles_gamB12Mic_kSimQ_kNewAnaC1_kNPE_2/gamB12Mic_kSimulation_kNewAnaCer_kNPE_nonlcov2.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/articles_Mic_kSimQ_kNewAnaC1_kNPE/Mic_kSimulation_kNewAnaCer_kNPE_nonlcov3.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/gamB12NewMicwhole/gamB12NewMicwhole_kSimQ_kCerNewAna_kNPE_nonlcov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/gamB12NewMicedge/gamB12NewMicedge_kSimQ_kCerNewAna_kNPE_nonlcov.txt"
    #filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/Newgam/Newgam_kSimQ_kNewCerAna_kNPE_nonlcov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12/NewgamB12_kSimQ_kNewCerAna_kNPE_nonlcov.txt"
    #filename3 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12NewMic1whole/NewgamB12NewMicwhole1_nonlcov.txt"
    #filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/Newgam_kSimQ_kNewAnaCer_kNew/Newgam_kSimQ_kNewCerAna_kNew_nonlcov.txt"
    filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12_kSimQ_kNewAnaCer_kNew/NewgamB12_kSimQ_kNewCerAna_kNew_nonlcov.txt"
    filename3 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12_kSimQ_kNewAnaCer_kNew/NewgamB12_kSimQ_kNewCerAna_kNew_nonlcov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12NewMicwhole_kSimQ_kNewAnaCer_kNew/NewgamB12NewMic_kSimQ_kNewCerAna_kNew_nonlcov.txt"
    filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12NewMicwhole_kSimQ_kNewAnaCer_kNew2/NewgamB12NewMicwhole_kSimQ_kNewCerAna_kNew2_nonlcov.txt"

    par1, parerr1 = loadBestFit(filename1, 7)
    scale1, kB1, p01, p11, p21, p31, p41 = par1

    par2, parerr2 = loadBestFit(filename2, 7)
    scale2, kB2, p02, p12, p22, p32, p42 = par2

    par3, parerr3 = loadBestFit(filename3, 7)
    scale3, kB3, p03, p13, p23, p33, p43 = par3


    Edep = np.arange(0.1, 65, 0.1)
    nonl1, nonl2, nonl3 = [], [], []

    for i in Edep:
        nonl1.append((fitNsct(i, scale1, kB1) + wenNcer(i, p01, p11, p21, p31, p41)) / i / Y )
        nonl2.append((fitNsct(i, scale2, kB2) + wenNcer(i, p02, p12, p22, p32, p42)) / i / Y )
        nonl3.append((fitNsct(i, scale3, kB3) + wenNcer(i, p03, p13, p23, p33, p43)) / i / Y )

    sampleSize = 5000
    sigma1 = sample_corelation(filename1, 7, sampleSize)
    sigma2 = sample_corelation(filename2, 7, sampleSize)
    sigma3 = sample_corelation(filename3, 7, sampleSize)

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
        for j in Edep:
            tmp_nonl1.append((fitNsct(j, m_scale1, m_kB1) + wenNcer(j, m_p01, m_p11, m_p21, m_p31, m_p41)) / j / Y )
            tmp_nonl2.append((fitNsct(j, m_scale2, m_kB2) + wenNcer(j, m_p02, m_p12, m_p22, m_p32,  m_p42)) / j / Y )
            tmp_nonl3.append((fitNsct(j, m_scale3, m_kB3) + wenNcer(j, m_p03, m_p13, m_p23, m_p33,  m_p43)) / j / Y )
    
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
    #fig, (ax, ax1) = plt.subplots(1, 2, figsize=(10, 6))


    #ax.set_title("(a)", y=-0.3)
    ax.errorbar(elecE, elecNonl, yerr=elecNonlerr, fmt="o", ms=6, color="blue", mfc="w", label="Simulation")

    ax.plot(Edep, nonl1, "-", color="crimson",   label=r"$\gamma$ + $^{12}$B")
    ax.plot(Edep, nonl2, "-", color="black", label=r"$\gamma$ + $^{12}$B + Michel $e^-$")

    ax.fill_between(Edep, nonl1-5*(nonl1-ymin1), nonl1+5*(ymax1-nonl1), color="crimson", alpha=0.3)
    ax.fill_between(Edep, nonl2-5*(nonl2-ymin2), nonl2+5*(ymax2-nonl2), color="slategrey", alpha=0.5)

    ax.legend(prop={"size":15})
    ax.set_xlabel(r"Electron $E^{dep}$ [MeV]", fontsize=16)
    ax.set_ylabel(r"$E^{vis}/E^{dep}$", fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True)


    #ax1.set_title("(b)", y=-0.3)
    #ax1.errorbar(elecE, elecNonl, yerr=elecNonlerr, fmt="o", color="black", mfc="w", label="Simulation truth")

    #ax1.plot(Edep, nonl1, "-", color="crimson",   label="gamma only")
    #ax1.plot(Edep, nonl2, "-", color="black", label="gamma+B12+Michel")

    #ax1.fill_between(Edep, nonl1-5*(nonl1-ymin1), nonl1+5*(ymax1-nonl1), color="crimson", alpha=0.3)
    #ax1.fill_between(Edep, nonl2-5*(nonl2-ymin2), nonl3+5*(ymax2-nonl2), color="slategrey", alpha=0.5)

    #ax1.legend(prop={"size":15})
    #ax1.set_xlabel(r"$E^{dep}$ [MeV]", fontsize=16)
    #ax1.set_ylabel(r"$E^{vis}/E^{dep}$", fontsize=16)
    #ax1.tick_params(axis='both', which='major', labelsize=14)


    
    elecE, elecNonl, elecNonlerr, elecRes, elecReserr = loadTruth()
    filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/Newgam_kSimQ_kNewAnaCer_kNew/Newgam_kSimQ_kNewCerAna_kNew_rescov.txt"
    filename3 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12_kSimQ_kNewAnaCer_kNew/NewgamB12_kSimQ_kNewCerAna_kNew_rescov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12NewMicedge_kSimQ_kNewAnaCer_kNew/NewgamB12NewMicedge_kSimQ_kNewCerAna_kNew_rescov.txt"
    filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12NewMicwhole_kSimQ_kNewAnaCer_kNew2/NewgamB12NewMicwhole_kSimQ_kNewCerAna_kNew2_rescov.txt"

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

    Evis = np.arange(0.1, 65, 0.1)
    Evis_posi = []
    res1, res2, res3 = [], [], []


    diffe1, diffe2, diffe3 = [], [], []
    for i in range(len(elecE)):
        Ntot = elecE[i] 
        sigma_elec1 = sigma2FuncNew(a1, b1, n1, Ntot)
        sigma_elec2 = sigma2FuncNew(a2, b2, n2, Ntot)
        sigma_elec3 = sigma2FuncNew(a3, b3, n3, Ntot)

        diffe1.append( (np.sqrt(sigma_elec1)/Ntot - elecRes[i])/elecRes[i] )
        diffe2.append( (np.sqrt(sigma_elec2)/Ntot - elecRes[i])/elecRes[i] )
        diffe3.append( (np.sqrt(sigma_elec3)/Ntot - elecRes[i])/elecRes[i] )



    
    diff13 = []
    for i in Evis:
        ## Best Fit
        Evistot = i*Y 
        #Evis_posi.append(Evistot/Y)

        sigma_elec1 = sigma2FuncNew(a1, b1, n1, Y*i)
        sigma_elec2 = sigma2FuncNew(a2, b2, n2, Y*i)
        sigma_elec3 = sigma2FuncNew(a3, b3, n3, Y*i)

        #Evistot = i*Y + 660.8 * 2 
        #res1.append(np.sqrt(sigma_elec1 + (2*27.07**2)) / Evistot )
        #res2.append(np.sqrt(sigma_elec2 + (2*27.07**2)) / Evistot )
        #res3.append(np.sqrt(sigma_elec3 + (2*27.07**2)) / Evistot )
        res1.append(np.sqrt(sigma_elec1) / Evistot )
        res2.append(np.sqrt(sigma_elec2) / Evistot )
        res3.append(np.sqrt(sigma_elec3) / Evistot )

        diff13.append((res1[-1] - res3[-1])/res3[-1])


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
            Evistot = j*Y
            #Evistot = j*Y + 660.8 * 2 

            sigma_elec1 = sigma2FuncNew(m_a1, m_b1, m_n1, Y*j)
            sigma_elec2 = sigma2FuncNew(m_a2, m_b2, m_n2, Y*j)
            sigma_elec3 = sigma2FuncNew(m_a3, m_b3, m_n3, Y*j)

            #tmp_res1.append(np.sqrt(sigma_elec1 + (2*27.07**2)) / Evistot )
            #tmp_res2.append(np.sqrt(sigma_elec2 + (2*27.07**2)) / Evistot )
            #tmp_res3.append(np.sqrt(sigma_elec3 + (2*27.07**2)) / Evistot )
            tmp_res1.append(np.sqrt(sigma_elec1) / Evistot )
            tmp_res2.append(np.sqrt(sigma_elec2) / Evistot )
            tmp_res3.append(np.sqrt(sigma_elec3) / Evistot )

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

    print(elecE)

    ax1.errorbar(elecNonl*elecE, elecRes, yerr=elecReserr, fmt="o", color="blue", ms=6, mfc="w", label="Simulation", zorder=2)
    ax1.plot(Evis, res3, "-", lw=2, color="crimson", label=r"$\gamma$ + $^{12}$B", zorder=1)
    ax1.plot(Evis, res2, "--", lw=2, color="black", label=r"$\gamma$ + $^{12}$B + Michel $e^-$", zorder=2)
    ax1.fill_between(Evis, res1-10*del1min, res1+10*del1max, color="crimson", alpha=0.3)
    ax1.fill_between(Evis, res2-10*del2min, res2+10*del2max, color="slategray", alpha=0.5)
    ax1.legend(prop={"size" : 15})
    ax1.set_xlabel(r"Electron $E^{vis}$ [MeV]", fontsize=15)
    ax1.set_ylabel(r"$\sigma/E^{vis}$", fontsize=15, color="black")
    ax1.tick_params(axis='both', which='major', labelsize=14, labelcolor="black")
    ax1.grid(True)




    plt.subplots_adjust(left=None, bottom=None, right=None, top=None , wspace=None, hspace=None)

    plt.tight_layout()
    plt.savefig("compareElec_gammaB12Mic.pdf")
    plt.show()
















