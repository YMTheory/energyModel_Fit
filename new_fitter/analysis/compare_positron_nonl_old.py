import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import elecLoader as eloader
import random
import ROOT

# no correlation in fitter...
# no kB here :
Etrue = np.arange(0.1, 8, 0.5)
Etrue_posi = Etrue + 1.022
dnum = len(Etrue)

def cerpe_func(E, A2, A3, A4):
    E0 = 0.2
    A1 = 0
    x = np.log(1+E/E0)
    return (A1*x+A2*x**2+A3*x**3) * (1/E+A4) * E 


def cerpe_sim(E, kC):
    return kC * eloader.getCerNPE(E)


def sctpe_sim(E, kA, kB, scale):
    return kA * eloader.getQPE(E, kB, scale)


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

    gamEtrue_nonl = [0.511, 0.662, 0.835, 1.253, 1.461, 2.223, 4.43,  6.13]
    gamnonlData = [0.91675776, 0.93238476, 0.94729476, 0.97075413, 0.97924947, 1., 1.0236367, 1.03170508]
    gamnonlDataErr = [0.00037812,0.00048843,0.00044424,0.00027509,0.00037001,0.00032481,0.00024566,0.0002075]
    gamnonlCalc = [0.91667098,0.93275607,0.94690211,0.97106998,0.97902328,0.9997544,1.02360547, 1.03157114]

    gamEtrue_res = [0.662, 0.835, 1.022, 1.461, 2.223, 2.506, 4.43, 4.94, 6.13]
    gamresData = [0.03667723, 0.03270304, 0.02632371, 0.02269925,  0.01655166, 0.01399196]
    gamresDataErr = [0.00037148, 0.0003263, 0.00028073, 0.00023905,  0.00017784, 0.00014814]
    gamresCalc = [0.03637193, 0.03288746, 0.02644741, 0.02244664, 0.01647451, 0.01402209]



    filename1 = "../nonlcov_gam_kSimulationSct_kSimulationCer_kTotal.txt"
    filename2 = "../nonlcov_gam+B12_kSimulationSct_kSimulationCer_kTotal.txt"
    filename3 = "../nonlcov_gam+B12_kSimulationSct_kSimulationCer_kSeparate_fixeds0.txt"
    filename4 = "../tmp_results/gam+B12_kSimulation_kSimulationCer_kNPE_fixedq0_nonlcov1.txt"
    filename5 = "../tmp_results/gam+B12_kIntegralCalc_kAnalyticalCer_kNPE_fixedp0_nonlcov.txt"
    filename6 = "../gam_kIntegralCalc_kAnalyticalCer_kNPE_nonlcov.txt"


    sampleSize = 5000
    #sigma1 = sample_corelation(filename1, 4, sampleSize)
    sigma4 = sample_corelation(filename4, 3, sampleSize)
    sigma5 = sample_corelation(filename5, 5, sampleSize)
    #sigma6 = sample_corelation(filename6, 5, sampleSize)
    #sigma3 = sample_corelation(filename3, 4, sampleSize)

    kA1, kB1, scale1, kC1 = 0.99805, 6.31679e-3, 1409.35, 0.988685
    kA2, kB2, scale2, kC2 = 0.99879, 6.34823e-3, 1409.61, 0.984026
    kA3, kB3, scale3, kC3 = 0.99989, 6.35613e-3, 1409.68, 0.983110
    kA4, kB4, scale4, kC4 = 1, 6.29208e-03, 1.40845e+03, 9.90770e-01
    #kA4, kB4, scale4, kC4 = 1, 6.32175e-03, 1.40891e3, 9.86594e-01
    kA5, kB5, scale5, A25, A35, A45 = 1, 5.20018e-03, 1.41741e+03, -1.34659e+01 ,  1.46295e+01 , 1.08040e-02
    #kA6, kB6, scale6, A26, A36, A46 = 1, 5.06020e-03, 1.41203e+03, -1.05573e+01, 1.37141e+01, 2.72101e-02
    
    global_es = 3134.078/2.223

    #pars1 = sample_uncorelation(filename1, 6, sampleSize)

    ymin4, ymax4 = [], []
    ymin5, ymax5 = [], []
    ymin6, ymax6 = [], []
    for i in range(dnum):
        ymin4.append(100)
        ymax4.append(-100)
        ymin5.append(100)
        ymax5.append(-100)
        #ymin6.append(100)
        #ymax6.append(-100)
    

    for i in range(sampleSize):

        pe_arr1, pe_arr2, pe_arr3, pe_arr4, pe_arr5, pe_arr6 = [], [], [], [], [], []

        #m_kA1 = kA1 + sigma1[0, i]
        #m_kB1 = kB1 + sigma1[1, i]
        #m_scale1 = scale1 + sigma1[2, i]
        #m_kC1 = kC1 + sigma1[3, i]

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

        #m_kA6 = kA6
        #m_scale6 = scale6 + sigma6[0, i]
        #m_kB6 = kB6 + sigma6[1, i]
        #m_A26 = A26 + sigma6[2, i]
        #m_A36 = A36 + sigma6[3, i]
        #m_A46 = A46 + sigma6[4, i]
    
        #m_kA3 = kA3 + sigma3[0, i]
        #m_kB3 = kB3 + sigma3[1, i]
        #m_scale3 = scale3 + sigma3[2, i]
        #m_kC3 = kC3 + sigma3[3, i]

        for i, j in zip(Etrue,  Etrue_posi):
            #pe_arr1.append( (sctpe_sim(j, m_kA1, m_kB1, m_scale1) + cerpe_sim(j, m_kC1)) /j/global_es)
            #pe_arr2.append( (sctpe_sim(j, m_kA2, m_kB2, m_scale2) + cerpe_sim(j, m_kC2)) /j/global_es)
            #pe_arr3.append( (sctpe_sim(j, m_kA3, m_kB3, m_scale3) + cerpe_sim(j, m_kC3)) /j/global_es)
            pe_arr4.append( (sctpe_sim(i, m_kA4, m_kB4, m_scale4) + cerpe_sim(i, m_kC4) + 2*660.8 ) / j /global_es)
            pe_arr5.append( (sctpe_sim(i, m_kA5, m_kB5, m_scale5) + cerpe_func(i, m_A25, m_A35, m_A45) + 2*660.8) / j / global_es)
            #pe_arr6.append( (sctpe_sim(i, m_kA6, m_kB6, m_scale6) + cerpe_func(i, m_A26, m_A36, m_A46) + 2*660.8) / j / global_es)
            #pe_arr4.append( (sctpe_sim(j, m_kA4, m_kB4, m_scale4) + cerpe_sim(j, m_kC4)+2*660.8) /j/global_es)

        for n in range(dnum):
            if pe_arr4[n] < ymin4[n] :
                ymin4[n] = pe_arr4[n]
            if pe_arr4[n] > ymax4[n]:
                ymax4[n] = pe_arr4[n]
            if pe_arr5[n] < ymin5[n] :
                ymin5[n] = pe_arr5[n]
            if pe_arr5[n] > ymax5[n]:
                ymax5[n] = pe_arr5[n]
            #if pe_arr6[n] < ymin6[n] :
            #    ymin6[n] = pe_arr6[n]
            #if pe_arr6[n] > ymax6[n]:
            #    ymax6[n] = pe_arr6[n]


        #plt.plot(Etrue, pe_arr1, color="lightskyblue", alpha=0.01)
        #plt.plot(Etrue, pe_arr2, color="lightskyblue", alpha=0.01)
        #plt.plot(Etrue, pe_arr3, color="lightseagreen", alpha=0.01)


    #print(ymin2[9]*1*global_es, ymin2[19]*2*global_es, ymin2[29]*global_es*3, ymin2[49]*5*global_es, ymin2[79]*8*global_es)
    #print(ymax2[9]*1*global_es, ymax2[19]*2*global_es, ymax2[29]*3*global_es, ymax2[49]*5*global_es, ymax2[79]*8*global_es)

    #plt.fill_between(Etrue_posi, ymin6, ymax6, alpha=0.5, color="lightskyblue")

    kesim, musim, muerrsim, sigmasim, sigmaerrsim = loadSimTruth()
    gSim = ROOT.TGraph()
    for i in range(len(musim)):
        gSim.SetPoint(i, kesim[i], musim[i])


    best1, best2, best3 = [], [], []
    best4, best5, best6 = [], [], []
    diff, diff_err = [], []
    diffSim, diffSimerr = [], []
    diffAna, diffAnaerr = [], []
    nominal = []
    nn = 0
    for i, j in zip(Etrue, Etrue_posi):
        #best1.append( (sctpe_sim(i, kA1, kB1, scale1) + cerpe_sim(i, kC1) ) /i/global_es)
        #best2.append( (sctpe_sim(i, kA4, kB4, scale4) + cerpe_sim(i, kC4) ) /i/global_es)
        best4.append( (sctpe_sim(i, kA4, kB4, scale4) + cerpe_sim(i, kC4) + 2*660.8) /j/global_es)
        best5.append( (sctpe_sim(i, kA5, kB5, scale5) + cerpe_func(i, A25, A35, A45) + 2*660.8) / j / global_es)
        #best6.append( (sctpe_sim(i, kA6, kB6, scale6) + cerpe_func(i, A26, A36, A46) + 2*660.8) / j / global_es)
        diff.append( (best5[-1] - best4[-1]) / best4[-1] )
        err4 = ymax4[nn] - best4[nn]
        err5 = ymax5[nn] - best5[nn]
        diff_err.append( np.sqrt(err5**2/best4[nn]**2 + err4**2*best5[nn]**2/best4[nn]**2) )
        diffSim.append( (best4[-1] - gSim.Eval(i)/global_es/j) / (gSim.Eval(i)/global_es/j) )
        diffSimerr.append(err4/(gSim.Eval(i)/global_es/j))
        diffAna.append( (best5[-1] - gSim.Eval(i)/global_es/j) / (gSim.Eval(i)/global_es/j) )
        diffAnaerr.append(err5/(gSim.Eval(i)/global_es/j))
        nn += 1
        #best3.append( (sctpe_sim(i, kA3, kB3, scale3) + cerpe_sim(i, kC3) ) /i/global_es)
        #nominal.append( (sctpe_sim(i, 1, 0.0065, 1377.772/0.97940231) + cerpe_sim(i, 1)) / i/global_es)

    diff = np.array(diff)
    diff_err = np.array(diff_err)
    diffSim = np.array(diffSim)
    diffSimerr = np.array(diffSimerr)
    diffAna = np.array(diffAna)
    diffAnaerr = np.array(diffAnaerr)

    #plt.plot(Etrue, nominal, "-.", lw=1, color="darkviolet", label="nominal")
    ##plt.plot(Etrue, best1, "-", color="coral", label="only gamma")

    #plt.plot(Etrue, best2, "-", lw=1.5, color="#FF8D70", label="electron")
    #plt.plot(Etrue_posi, best4, "-.", color="#FF60F1", lw=1.5, label="positron")
    #plt.fill_between(Etrue_posi, ymin4, ymax4, alpha=1.0, color="lightskyblue", label=r"positron $1 \sigma$ zone")
    #plt.plot(Etrue, best3, "--", lw=1, color="peru",  label="best fit: gamma+B12")
    
    #plt.errorbar(gamEtrue_nonl, gamnonlData, yerr=gamnonlDataErr, fmt="o", color="#90F7C4", label="gamma: simulation")
    #plt.plot(gamEtrue_nonl, gamnonlCalc, "--", lw=1.5, color="#7054FF", label="gamma: fitting")
    

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
    #ax0.grid(True)
    ax0.set_ylabel("relaive difference", fontsize=15)
    ax0.tick_params(axis='both', which='major', labelsize=13)
    ax0.set_ylim(-0.01, 0.005)

    ax1.plot(Etrue_posi, best4, "-",  color="blue",   label="simulation-based")
    ax1.plot(Etrue_posi, best5, "--", color="#FF60F1",     label="analytical")
    ax1.plot(kesim+1.022, musim/global_es/(kesim+1.022), color="lightseagreen", label="simulation truth")
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
    #plt.savefig("Compare_positron_nonl.pdf")
    plt.show()


if __name__ == "__main__":
    main()
