import numpy as np
import matplotlib.pyplot as plt
import elecLoader as eloader
import random

# no correlation in fitter...
# no kB here :
Etrue = np.arange(0.5, 7, 0.5)
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
    filename4 = "../output/draft/nonlCov.txt"


    sampleSize = 5000
    #sigma1 = sample_corelation(filename1, 4, sampleSize)
    sigma4 = sample_corelation(filename4, 3, sampleSize)
    #sigma3 = sample_corelation(filename3, 4, sampleSize)

    kA1, kB1, scale1, kC1 = 0.99805, 6.31679e-3, 1409.35, 0.988685
    kA2, kB2, scale2, kC2 = 0.99879, 6.34823e-3, 1409.61, 0.984026
    kA3, kB3, scale3, kC3 = 0.99989, 6.35613e-3, 1409.68, 0.983110
    kA4, kB4, scale4, kC4 = 1, 6.34928e-03, 1.40929e3, 9.83487e-01
    
    global_es = 3134.078/2.223

    #pars1 = sample_uncorelation(filename1, 6, sampleSize)

    ymin4, ymax4 = [], []
    #ymin3, ymax3 = [], []
    for i in range(dnum):
        ymin4.append(100)
        ymax4.append(-100)
        #ymin3.append(100)
        #ymax3.append(-100)
    
    for i in range(sampleSize):

        pe_arr1, pe_arr2, pe_arr3, pe_arr4 = [], [], [], []

        #m_kA1 = kA1 + sigma1[0, i]
        #m_kB1 = kB1 + sigma1[1, i]
        #m_scale1 = scale1 + sigma1[2, i]
        #m_kC1 = kC1 + sigma1[3, i]

        m_kA4 = kA4
        m_kB4 = kB4 + sigma4[0, i]
        m_kC4 = kC4 + sigma4[1, i]
        m_scale4 = scale4 + sigma4[2, i]

        #m_kA3 = kA3 + sigma3[0, i]
        #m_kB3 = kB3 + sigma3[1, i]
        #m_scale3 = scale3 + sigma3[2, i]
        #m_kC3 = kC3 + sigma3[3, i]

        for j in Etrue:
            #pe_arr1.append( (sctpe_sim(j, m_kA1, m_kB1, m_scale1) + cerpe_sim(j, m_kC1)) /j/global_es)
            #pe_arr2.append( (sctpe_sim(j, m_kA2, m_kB2, m_scale2) + cerpe_sim(j, m_kC2)) /j/global_es)
            #pe_arr3.append( (sctpe_sim(j, m_kA3, m_kB3, m_scale3) + cerpe_sim(j, m_kC3)) /j/global_es)
            pe_arr4.append( (sctpe_sim(j, m_kA4, m_kB4, m_scale4) + cerpe_sim(j, m_kC4)) /j/global_es)

        for n in range(dnum):
            if pe_arr4[n] < ymin4[n] :
                ymin4[n] = pe_arr4[n]
            if pe_arr4[n] > ymax4[n]:
                ymax4[n] = pe_arr4[n]
            #if pe_arr3[n] < ymin3[n] :
            #    ymin3[n] = pe_arr3[n]
            #if pe_arr3[n] > ymax3[n]:
            #    ymax3[n] = pe_arr3[n]


        #plt.plot(Etrue, pe_arr1, color="lightskyblue", alpha=0.01)
        #plt.plot(Etrue, pe_arr2, color="lightskyblue", alpha=0.01)
        #plt.plot(Etrue, pe_arr3, color="lightseagreen", alpha=0.01)


    #print(ymin2[9]*1*global_es, ymin2[19]*2*global_es, ymin2[29]*global_es*3, ymin2[49]*5*global_es, ymin2[79]*8*global_es)
    #print(ymax2[9]*1*global_es, ymax2[19]*2*global_es, ymax2[29]*3*global_es, ymax2[49]*5*global_es, ymax2[79]*8*global_es)

    #plt.fill_between(Etrue, ymin2, ymax2, alpha=1.0, color="lightskyblue", label="total resolution")
    #plt.fill_between(Etrue, ymin3, ymax3, alpha=0.5, color="pink", label="sct+cer resolution")
    plt.fill_between(Etrue, ymin4, ymax4, alpha=1.0, color="lightskyblue", label="total resolution")
    
    best1, best2, best3 = [], [], []
    nominal = []
    for i in Etrue:
        #best1.append( (sctpe_sim(i, kA1, kB1, scale1) + cerpe_sim(i, kC1) ) /i/global_es)
        best2.append( (sctpe_sim(i, kA2, kB2, scale2) + cerpe_sim(i, kC2) ) /i/global_es)
        #best3.append( (sctpe_sim(i, kA3, kB3, scale3) + cerpe_sim(i, kC3) ) /i/global_es)
        nominal.append( (sctpe_sim(i, 1, 0.0065, 1377.772/0.97940231) + cerpe_sim(i, 1)) / i/global_es)

    #plt.plot(Etrue, nominal, "-.", lw=1, color="darkviolet", label="nominal")
    ##plt.plot(Etrue, best1, "-", color="coral", label="only gamma")
    plt.plot(Etrue, best2, "-", lw=1, color="coral", label="electron")
    #plt.plot(Etrue, best3, "--", lw=1, color="peru",  label="best fit: gamma+B12")
    
    plt.errorbar(gamEtrue_nonl, gamnonlData, yerr=gamnonlDataErr, fmt="o", color="magenta", label="gamma: simulation")
    plt.plot(gamEtrue_nonl, gamnonlCalc, "--", color="blue", label="gamma: fitting")
    
    
    
    plt.legend(loc="lower right")

    plt.grid(True)
    plt.xlabel(r"$E_{dep}/MeV$")
    plt.ylabel(r"$E_{vis}/E_{dep}$")
    #plt.title("Electron NPE Nonlinearity")
    #plt.semilogy()
    plt.savefig("Compare_nonl.pdf")
    plt.show()


if __name__ == "__main__":
    main()
