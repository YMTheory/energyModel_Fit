import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import ROOT
from scipy.linalg import eigh, cholesky
from scipy.stats import norm

import elecLoader as el

#kA5, kB5, scale5, A25, A35, A45 = 1,  5.09593e-03, 1.40356e+03, -8.75277e+00, 1.41055e+01, 2.89894e-02
#kA5, kB5, scale5, A25, A35, A45 = 1,  5.086e-03, 1.40356e+03, -8.7533e+00, 1.41056e+01, 2.89873e-02
#kA4, scale4, kB4, p04, p14, p24, p34, p44 = 1, 1410.19, 5.45e-3, -0.407, 0.423, -0.0340, 3.40, 0.0192
#kA4, scale4, kB4, p04, p14, p24, p34, p44 = 1, 1410.71, 5.48e-3, -0.405, 0.422, -0.0393, 3.53, 0.0182
#kA4, scale4, kB4, p04, p14, p24, p34, p44 = 1, 1407.55, 6.18e-3, -0.409, 0.421, -0.0224, 3.02, 0.0189


def func(E, A2, A3, A4):
    E0 = 0.2
    A1 = 0
    x = np.log(1+E/E0)

    cerpe = (A1*x+A2*x**2+A3*x**3) * (1/E+A4) * E
    #print("%.2f, %.5f %.5f, %.5f, %.2f" %(E, A2, A3, A4, cerpe))
    return cerpe



def wenNcer(E, p0, p1, p2, p3, p4, E0):
    if E<=E0:
        return 0
    if E>E0:
        E = E - E0
        return p3*E**2/(p4+p0*E+p1*np.exp(-p2*E))


def load(filename):
    Etrue, cerpe = [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if float(data[0]) <= 15 :
                Etrue.append(float(data[0]))
                cerpe.append(float(data[1]))

    return Etrue, cerpe

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



def select(old1, old2):
    new1, new2 = [], []
    for i in range(900):
        if i%30 == 0 and i!=0:
            new1.append(old1[i])
            new2.append(old2[i])

        if i == 880 or i == 890 or i==900 or i==910 or i==872 or i==874 or i==876:
            new1.append(old1[i])
            new2.append(old2[i])


    new1 = np.array(new1)
    new2 = np.array(new2)
    return new1, new2


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

def main():
    Ac = 178.044 / 2.223
    
    Etrue, cerpe = load("../data/electron/cerPE4.txt")
    Etrue1 = np.array(Etrue)
    cerpe1 = np.array(cerpe)

    #Etrue1, cerpe1 = select(Etrue, cerpe) 


    filename1 = "/junofs/users/miaoyu/energy_model/fitter/energyModel_Fit/new_fitter/outputs/gamB12_freeE0_kSim_kAna_kNew/gamB12_freeE0_kSim_kAna_kNew_nonlcov.txt"

    npar = 8
    par, parerr = loadBestFit(filename1, npar)

    sampleSize = 5000
    sigma5 = sample_corelation(filename1, npar, sampleSize)

    ymin5, ymax5 = [], []
    dnum = len(Etrue1)
    for i in range(dnum):
        ymin5.append(1000000)
        ymax5.append(-1000)
        
    for i in range(sampleSize):

        pe_arr5 =  []

        m_p0 = par[2] + sigma5[2, i]
        m_p1 = par[3] + sigma5[3, i]
        m_p2 = par[4] + sigma5[4, i]
        m_p3 = par[5] + sigma5[5, i]
        m_p4 = par[6] + sigma5[6, i]
        m_E0 = par[7] + sigma5[7, i]

        for i in Etrue1:
            #pe_arr5.append( func(i, m_A25, m_A35, m_A45) )
            pe_arr5.append(wenNcer(i, m_p0, m_p1, m_p2, m_p3, m_p4, m_E0))


        for n in range(dnum):
            if pe_arr5[n] < ymin5[n] :
                ymin5[n] = pe_arr5[n]
            if pe_arr5[n] > ymax5[n]:
                ymax5[n] = pe_arr5[n]

    ymin5 = np.array(ymin5)
    ymax5 = np.array(ymax5)

    ymax6, ymin6 = [], []

    best = []
    for i in Etrue1:
        best.append(wenNcer(i, par[2], par[3], par[4], par[5], par[6], par[7]))
    best = np.array(best)

    for i in range(len(Etrue1)):
        dmin = best[i] - ymin5[i]
        ymin6.append(best[i] - 5*dmin)
        dmax = best[i] - ymax5[i]
        ymax6.append(best[i] - 5*dmax)


    ymin6 = np.array(ymin6)
    ymax6 = np.array(ymax6)

    fig, ax = plt.subplots()
    
    Ac_calc = Ac

    plt.plot(Etrue1, cerpe1/(Etrue1*Ac), "o", ms=6, color="blue", label="Simulation")
    #plt.plot(Etrue3, wenNcer(Etrue3, *popt)/Etrue3/(177.508/2.223), lw=2, color="black")
    plt.plot(Etrue1, best/(Ac_calc*Etrue1), "--", lw=2, color="crimson", label="Best fit")
    plt.fill_between(Etrue1, ymin6/Etrue1/Ac_calc, ymax6/Etrue1/Ac_calc, color="slategray", alpha=0.5)


    ax.set_xlabel("Electron $E$ [MeV]", fontsize=15)
    ax.set_ylabel(r"Normalized Cherenkov yield ($f_C$)", fontsize=14)

    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.legend(prop={"size":14})

    ax.grid(True)

    plt.tight_layout()
    plt.savefig("CerenkovPE_MeV.pdf")

    plt.show()



    """

    x1 = np.arange(0.1, 15, 0.1)
    #y1 = func(x1, *popt)


    y1, y2 = [], []
    for i in x1:
        y1.append( func(i, A25, A35, A45) )


    x1 = np.array(x1)
    y1 = np.array(y1)


    filename5 = "../output/articles/gam_kIntegral_kAnalytical_kNPE_fixedc0_nonlcov.txt"
    sampleSize = 5000
    sigma5 = sample_corelation(filename5, 5, sampleSize)
    ymin5, ymax5 = [], []
    dnum = len(Etrue1)
    for i in range(dnum):
        ymin5.append(1000000)
        ymax5.append(-1000)
        
    for i in range(sampleSize):

        pe_arr5 =  []

        m_A25 = A25 + sigma5[2, i]
        m_A35 = A35 + sigma5[3, i]
        m_A45 = A45 + sigma5[4, i]

        for i in Etrue1:
            pe_arr5.append( func(i, m_A25, m_A35, m_A45) )

        for n in range(dnum):
            if pe_arr5[n] < ymin5[n] :
                ymin5[n] = pe_arr5[n]
            if pe_arr5[n] > ymax5[n]:
                ymax5[n] = pe_arr5[n]

    ymin5 = np.array(ymin5)
    ymax5 = np.array(ymax5)


    fig, ax = plt.subplots()
    #ax.plot(x, y, "o", ms=1, color='blue')
    ax.plot(Etrue1, cerpe1, "o", ms=4, color='royalblue', label="Geant4")
    #ax.plot(Etrue2, cerpe2, "o", ms=4, color='red')
    plt.plot(x1, y1/x1/177.508*2.22, "-", color="red", label="Best Fit")
    plt.fill_between(x1, ymin5/x1/177.508*2.22, ymax5/x1/177.508*2.22, color="slategray", alpha=0.6)
    

    ax.set_xlabel("Electron $E_{dep}$ [MeV]", fontsize=15)
    ax.set_ylabel(r"$\frac{N_{Cer}/E_{dep}}{N_{Cer}^0/2.22}$", fontsize=21)

    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.legend(prop={"size":14})

    ax.grid(True)

    plt.tight_layout()
    plt.savefig("CerenkovPE_MeV.pdf")

    plt.show()

    """

if __name__ == "__main__":
    main()
            

