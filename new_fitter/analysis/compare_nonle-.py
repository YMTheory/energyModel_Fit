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

    filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamNewB12_kSimQ_kNewAnaCer_kNew/NewgamNewB12_kSimQ_kNewAnaCer_kNew_nonlcov.txt"
    filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamNewB12NewMicedge_kSimQ_kNewAnaCer_kNew/NewgamNewB12NewMicedge_kSimQ_kNewAnaCer_kNew_nonlcov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamNewB12NewMic_kSimQ_kNewAnaCer_kNew/NewgamNewB12NewMic_kSimQ_kNewAnaCer_kNew_nonlcov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/Newgam_kSimQ_kNewAnaCer_kNew/Newgam_kSimQ_kNewAnaCer_kNew_nonlcov.txt"

    par1, parerr1 = loadBestFit(filename1, 7)
    scale1, kB1, p01, p11, p21, p31, p41 = par1

    par2, parerr2 = loadBestFit(filename2, 7)
    scale2, kB2, p02, p12, p22, p32, p42 = par2


    Edep = np.arange(0.1, 65, 0.1)
    nonl1, nonl2 = [], []

    for i in Edep:
        nonl1.append((fitNsct(i, scale1, kB1) + wenNcer(i, p01, p11, p21, p31, p41)) / i / Y )
        nonl2.append((fitNsct(i, scale2, kB2) + wenNcer(i, p02, p12, p22, p32, p42)) / i / Y )

    sampleSize = 5000
    sigma1 = sample_corelation(filename1, 7, sampleSize)
    sigma2 = sample_corelation(filename2, 7, sampleSize)

    diffmin1, diffmax1, diffmin2, diffmax2 = [], [], [], []
    for i in range(len(Edep)):
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

        tmp_nonl1, tmp_nonl2 = [], []
        for k in Edep:
            tmp_nonl1.append((fitNsct(k, m_scale1, m_kB1) + wenNcer(k, m_p01, m_p11, m_p21, m_p31, m_p41)) / k / Y )
            tmp_nonl2.append((fitNsct(k, m_scale2, m_kB2) + wenNcer(k, m_p02, m_p12, m_p22, m_p32, m_p42)) / k / Y )
        for j in range(len(Edep)):
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
    for i in range(len(Edep)):
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
    nonl_mc1, nonl_mc2 = [], []
    for i in elecE:
        ke = i
        nl1 = ((fitNsct(ke, scale1, kB1) + wenNcer(ke, p01, p11, p21, p31, p41) ) / ke / Y )
        nl2 = ((fitNsct(ke, scale2, kB2) + wenNcer(ke, p02, p12, p22, p32, p42) ) / ke / Y )
        
        nonl_mc1.append(nl1)
        nonl_mc2.append(nl2)

    nonl_mc1 = np.array(nonl_mc1)
    nonl_mc2 = np.array(nonl_mc2)


    ymin1, ymax1, ymin2, ymax2 = [], [], [], []
    for i in range(len(elecE)):
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

        m_scale2 = scale2 + sigma2[0, i]
        m_kB2    = kB2 + sigma2[1, i]
        m_p02    = p02 + sigma2[2, i]
        m_p12    = p12 + sigma2[3, i]
        m_p22    = p22 + sigma2[4, i]
        m_p32    = p32 + sigma2[5, i]
        m_p42    = p42 + sigma2[6, i]

        tmp_nonl1, tmp_nonl2 = [], []
        for j in elecE:
            tmp_nonl1.append((fitNsct(j, m_scale1, m_kB1) + wenNcer(j, m_p01, m_p11, m_p21, m_p31, m_p41)) / j / Y )
            tmp_nonl2.append((fitNsct(j, m_scale2, m_kB2) + wenNcer(j, m_p02, m_p12, m_p22, m_p32, m_p42)) / j / Y )
    
        for j in range(len(elecE)):
            if tmp_nonl1[j] > ymax1[j]:
                ymax1[j] = tmp_nonl1[j]
            if tmp_nonl1[j] < ymin1[j]:
                ymin1[j] = tmp_nonl1[j]
            if tmp_nonl2[j] > ymax2[j]:
                ymax2[j] = tmp_nonl2[j]
            if tmp_nonl2[j] < ymin2[j]:
                ymin2[j] = tmp_nonl2[j]

    Edep = np.array(Edep)
    ymin1 = np.array(ymin1)
    ymax1 = np.array(ymax1)
    ymin2 = np.array(ymin2)
    ymax2 = np.array(ymax2)

    derr1, derr2 = [], []
    for i in range(len(elecE)):
        del1min = nonl_mc1[i] - ymin1[i]
        del1max = ymax1[i] - nonl_mc1[i]
        del1 = del1min > del1max and del1min or del1max
        derr1.append(del1)
        del2min = nonl_mc2[i] - ymin2[i]
        del2max = ymax2[i] - nonl_mc2[i]
        del2 = del2min > del2max and del2min or del2max
        derr2.append(del2)

    derr1 = np.array(derr1)
    derr2 = np.array(derr2)

    # Plotting

    #fig, ax = plt.subplots()
    fig = plt.figure(figsize=(12, 6))
    spec = gridspec.GridSpec(ncols=2, nrows=2,
                         height_ratios=[1, 2])

    ax1 = fig.add_subplot(spec[2])
    ax = fig.add_subplot(spec[0])
    ax2 = fig.add_subplot(spec[1])
    ax3 = fig.add_subplot(spec[3])


    #ax.plot(kesim+1.022, (diff1-nonlSim)/nonlSim, "-", lw=2, color="crimson", label="only gamma")
    ax.plot(elecE, (nonl_mc1-elecNonl)/elecNonl*100, "-", lw=2, color="crimson")
    ax.plot(elecE, (nonl_mc2-elecNonl)/elecNonl*100, "--", lw=2, color="black")
    ax.fill_between(elecE, (nonl_mc1-1*derr1-elecNonl)/elecNonl*100, (nonl_mc1+1*derr1-elecNonl)/elecNonl*100, color="crimson", alpha=0.3)
    ax.fill_between(elecE, (nonl_mc2-1*derr2-elecNonl)/elecNonl*100, (nonl_mc2+1*derr2-elecNonl)/elecNonl*100, color="slategrey", alpha=0.5, label=r"$1\sigma$ errorbar")
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
    #ax.set_ylim(-0.01, 0.01)
    ax.grid(True)


    ax1.errorbar(elecE, elecNonl, yerr=elecNonlerr, fmt="o", ms=6, color="blue", mfc="w", zorder=3, label="Simulation")

    ax1.plot(Edep, nonl1, "-",  lw=2, color="crimson", label=r"$\gamma$+$^{12}$B: $10\times1\sigma$ errorbar", zorder=1)
    ax1.plot(Edep, nonl2, "--", lw=2, color="black",  label=r"$\gamma$+$^{12}$B+Michel $e^-$:" "\n"  r"$10\times1\sigma$ errorbar", zorder=2)

    ax1.fill_between(Edep, nonl1-10*(nonl1-diffmin1), nonl1+10*(diffmax1-nonl1), color="crimson", alpha=0.3)
    ax1.fill_between(Edep, nonl2-10*(nonl2-diffmin2), nonl2+10*(diffmax2-nonl2), color="slategrey", alpha=0.5)
    #ax1.set_xlim(-5, 70)

    ax1.legend(prop={"size":15})
    ax1.set_xlabel(r"Electron $E^{dep}$ [MeV]", fontsize=16)
    ax1.set_ylabel(r"$E^{vis}/E^{dep}$", fontsize=16)
    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax1.grid(True)
    ax1.set_title("(a)", y=-0.3, fontsize=15)



    ######################################################################################
    ################################### Resolution #######################################
    ######################################################################################
    
    filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamNewB12_kSimQ_kNewAnaCer_kNew/NewgamNewB12_kSimQ_kNewAnaCer_kNew_rescov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12NewMicwhole_kSimQ_kNewAnaCer_kNew2/NewgamB12NewMicwhole_kSimQ_kNewAnaCer_kNew2_rescov.txt"
    filename3 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12NewMic1whole/NewgamB12NewMicwhole1_rescov.txt"
    filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamNewB12NewMicedge_kSimQ_kNewAnaCer_kNew/NewgamNewB12NewMicedge_kSimQ_kNewAnaCer_kNew_rescov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamNewB12NewMic_kSimQ_kNewAnaCer_kNew/NewgamNewB12NewMic_kSimQ_kNewAnaCer_kNew_rescov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/Newgam_kSimQ_kNewAnaCer_kNew/Newgam_kSimQ_kNewAnaCer_kNew_rescov.txt"

    par1, parerr1 = loadBestFit(filename1, 3)
    a1, b1, n1 = par1

    par2, parerr2 = loadBestFit(filename2, 3)
    a2, b2, n2 = par2

    par3, parerr3 = loadBestFit(filename3, 2)
    a3, b3 = par3

    Evis = np.arange(0.1, 65, 0.1)
    res1, res2 = [], []

    res3 = []

    for i in Evis:
        ## Best Fit
        Evistot = i*Y 

        sigma_elec1 = sigma2FuncNew(a1, b1, n1, Y*i)
        sigma_elec2 = sigma2FuncNew(a2, b2, n2, Y*i)

        res1.append(np.sqrt(sigma_elec1) / Evistot )
        res2.append(np.sqrt(sigma_elec2) / Evistot )

        sigma_elec3 = sigma2FuncNew(a3, b3, 2, Y*i)
        res3.append(np.sqrt(sigma_elec3) / Evistot )


    sampleSize = 5000
    sigma1 = sample_corelation(filename1, 3, sampleSize)
    sigma2 = sample_corelation(filename2, 3, sampleSize)

    ymin1, ymax1 = [], []
    ymin2, ymax2 = [], []
    for i in range(len(Evis)):
        ymin1.append(1000000)
        ymax1.append(-100)
        ymin2.append(1000000)
        ymax2.append(-100)


    for i in range(sampleSize):

        m_a1 = a1 + sigma1[0, i]
        m_b1 = b1 + sigma1[1, i]
        m_n1 = n1 + sigma1[2, i]

        m_a2 = a2 + sigma2[0, i]
        m_b2 = b2 + sigma2[1, i]
        m_n2 = n2 + sigma2[2, i]

        tmp_res1 = []
        tmp_res2 = []
        for j in Evis:
            Evistot = j*Y

            sigma_elec1 = sigma2FuncNew(m_a1, m_b1, m_n1, Y*j)
            sigma_elec2 = sigma2FuncNew(m_a2, m_b2, m_n2, Y*j)

            tmp_res1.append(np.sqrt(sigma_elec1) / Evistot )
            tmp_res2.append(np.sqrt(sigma_elec2) / Evistot )

        for j in range(len(Evis)):
            if tmp_res1[j] > ymax1[j]:
                ymax1[j] = tmp_res1[j]
            if tmp_res1[j] < ymin1[j]:
                ymin1[j] = tmp_res1[j]
            if tmp_res2[j] > ymax2[j]:
                ymax2[j] = tmp_res2[j]
            if tmp_res2[j] < ymin2[j]:
                ymin2[j] = tmp_res2[j]


    res1 = np.array(res1)
    ymin1 = np.array(ymin1)
    ymax1 = np.array(ymax1)
    res2 = np.array(res2)
    ymin2 = np.array(ymin2)
    ymax2 = np.array(ymax2)


    # compare with MC Truth :
    res_mc1, res_mc2 = [], []
    res_mc3 = []
    for i in range(len(elecE)):
        evis = elecNonl[i] * elecE[i]
        sigma_elec1 = sigma2FuncNew(a1, b1, n1, evis*Y)
        sigma_elec2 = sigma2FuncNew(a2, b2, n2, evis*Y)
        sigma_elec3 = sigma2FuncNew(a3, b3, 2, evis*Y)

        res_mc1.append(np.sqrt(sigma_elec1) / evis / Y)
        res_mc2.append(np.sqrt(sigma_elec2) / evis / Y)
        res_mc3.append(np.sqrt(sigma_elec3) / evis / Y)


    res_mc1 = np.array(res_mc1)
    res_mc2 = np.array(res_mc2)
    res_mc3 = np.array(res_mc3)



    ymin3, ymax3 = [], []
    ymin4, ymax4 = [], []
    ymin5, ymax5 = [], []
    for i in range(len(elecE)):
        ymin3.append(1000000)
        ymax3.append(-100)
        ymin4.append(1000000)
        ymax4.append(-100)
        ymin5.append(1000000)
        ymax5.append(-100)


    for i in range(sampleSize):

        m_a1 = a1 + sigma1[0, i]
        m_b1 = b1 + sigma1[1, i]
        m_n1 = n1 + sigma1[2, i]

        m_a2 = a2 + sigma2[0, i]
        m_b2 = b2 + sigma2[1, i]
        m_n2 = n2 + sigma2[2, i]

        tmp_res1 = []
        tmp_res2 = []
        for j in elecNonl*elecE:
            Evistot = j*Y

            sigma_elec1 = sigma2FuncNew(m_a1, m_b1, m_n1, Y*j)
            sigma_elec2 = sigma2FuncNew(m_a2, m_b2, m_n2, Y*j)

            tmp_res1.append(np.sqrt(sigma_elec1) / Evistot )
            tmp_res2.append(np.sqrt(sigma_elec2) / Evistot )

        for j in range(len(elecE)):
            if tmp_res1[j] > ymax3[j]:
                ymax3[j] = tmp_res1[j]
            if tmp_res1[j] < ymin3[j]:
                ymin3[j] = tmp_res1[j]
            if tmp_res2[j] > ymax4[j]:
                ymax4[j] = tmp_res2[j]
            if tmp_res2[j] < ymin4[j]:
                ymin4[j] = tmp_res2[j]


    ymin3 = np.array(ymin3)
    ymax3 = np.array(ymax3)
    ymin4 = np.array(ymin4)
    ymax4 = np.array(ymax4)

    err3, err4 = [], []
    for i in range(len(elecE)):
        dmin3 = res_mc1[i] - ymin3[i]
        dmax3 = ymax3[i] - res_mc1[i]
        if dmax3 > dmin3 :
            err3.append(dmax3)
        else:
            err3.append(dmin3)

        dmin4 = res_mc2[i] - ymin4[i]
        dmax4 = ymax4[i] - res_mc2[i]
        if dmax4 > dmin4 :
            err4.append(dmax4)
        else:
            err4.append(dmin4)




    #ax2.plot(Evis_posi, diff2, "--", color="black")
    #ax2.fill_between(Evis_posi, (res2-err2-simData)/simData, (res2+err2-simData)/simData, color="slategray", alpha=0.5, label=r"$1\sigma$ errorbar")
    ax2.plot(elecNonl*elecE, (res_mc3-elecRes)/elecRes*100, "-.", lw=2, color="blue")
    ax2.plot(elecNonl*elecE, (res_mc1-elecRes)/elecRes*100, "-", lw=2, color="crimson")
    ax2.fill_between(elecNonl*elecE, (res_mc1-err3-elecRes)/elecRes*100, (res_mc1+err3-elecRes)/elecRes*100, alpha=0.3, color="crimson")
    ax2.plot(elecNonl*elecE, (res_mc2-elecRes)/elecRes*100, "--", lw=2, color="black")
    ax2.fill_between(elecNonl*elecE, (res_mc2-err4-elecRes)/elecRes*100, (res_mc2+err4-elecRes)/elecRes*100, alpha=0.5, color="slategray", label=r"$1\sigma$ errorbar")
    ax2.set_ylabel("Bias (%)", fontsize=15, color="black")
    ax2.tick_params(axis='both', which='major', labelsize=14, labelcolor="black")
    #ax2.set_ylim(-5, 5)
    ax2.legend(loc="upper left", prop={"size" : 15})
    ticknames = [ "-10", "-5", "0", "5", "10", "15", "20", "25"]
    ax2.set_yticks([-10, -5, 0, 5, 10, 15, 20, 25])
    ax2.set_yticklabels(ticknames)
    ax2.grid(True)

    ax3.errorbar(elecNonl*elecE, elecRes, yerr=elecReserr, fmt="o", color="blue", ms=6, mfc="w", label="Simulation", zorder=3)
    ax3.plot(Evis, res1, "-", lw=2, color="crimson", label=r"$\gamma$+$^{12}$B:$10\times 1\sigma$ errorbar", zorder=2)
    ax3.fill_between(Evis, res1-10*(res1-ymin1), res1+10*(ymax1-res1), color="crimson", alpha=0.3)
    ax3.plot(Evis, res2, "--", lw=2, color="black", label=r"$\gamma$+$^{12}$B+Michel $e^-$:" "\n" r"$10\times 1\sigma$ errorbar", zorder=2)
    ax3.fill_between(Evis, res2-10*(res2-ymin2), res2+10*(ymax2-res2), color="slategray", alpha=0.5)
    ax3.plot(Evis, res3, "-.", lw=1.5, color="blue", label=r"$\gamma$+$^{12}$B+Michel $e^-$: fix $n=2$", zorder=2)
    #ax3.set_xlim(-5, 70)
    ax3.legend(prop={"size" : 15})
    ax3.set_xlabel(r"Electron $E^{vis}$ [MeV]", fontsize=15)
    ax3.set_ylabel(r"$\sigma/E^{vis}$", fontsize=15, color="black")
    ax3.tick_params(axis='both', which='major', labelsize=14, labelcolor="black")
    ax3.grid(True)
    #ax.semilogy()
    ax3.set_title("(b)", y=-0.3, fontsize=15)
    ax3.semilogy()
    ax3.set_ylim(0.002, 0.2)
    



    plt.subplots_adjust(left=None, bottom=None, right=None, top=None , wspace=None, hspace=None)

    plt.tight_layout()
    plt.savefig("compare_electron_nonlres.pdf")
    plt.show()
















