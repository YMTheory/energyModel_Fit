import numpy as np
import matplotlib.pyplot as plt
import elecLoader as eloader
import random


# no correlation in fitter...

cerpe = np.arange(0, 1500, 10)
dnum = len(cerpe)


def cerFunc(x, a, b, c):
    s2 = a + b*x + c*x**2
    if s2<0:
        return 0
    else:
        return np.sqrt(s2)

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


def sctpe_sim(E, kA, kB, scale):
    return kA * eloader.getQPE(E, kB, scale)



def main():
    sampleSize = 1000
    sigma = sample_corelation("../rescov_gam+B12_kSimulationSct_kSimulationCer_kSeparate_fixeds0.txt", 3, sampleSize)

    c0, c1, c2 = -123.213, 3.39232, 2.53173e-3

    ymin, ymax = [], []
    for i in range(dnum):
        ymin.append(1000000)
        ymax.append(-100)

    for i in range(sampleSize):

        m_c0 = c0 + random.gauss(0, sigma[0, i])
        m_c1 = c1 + random.gauss(0, sigma[1, i])
        m_c2 = c2 + random.gauss(0, sigma[2, i])


        cersigma = []

        for i in cerpe:

            cersigma2 = cerFunc(i, m_c0, m_c1, m_c2)**2
            cersigma.append(np.sqrt(cersigma2))

        
        for n in range(dnum):
            if cersigma[n] < ymin[n]:
                ymin[n] = cersigma[n]
            if cersigma[n] > ymax[n]:
                ymax[n] = cersigma[n]
    
    ymin = np.array(ymin)
    ymax = np.array(ymax)

    best, nominal = [], []
    bc0, bc1, bc2 = -158.41, 4.333, 0.00181082
    for i in cerpe:
        cersigma2 = cerFunc(i, c0, c1, c2)**2
        best.append(np.sqrt(cersigma2)) 
        cersigma2 = cerFunc(i, bc0, bc1, bc2)**2
        nominal.append( np.sqrt(cersigma2))

    best = np.array(best)


    plt.fill_between(cerpe, ymin, ymax, color="lightskyblue", label=r"$1 \sigma$ zone")
    plt.plot(cerpe, nominal, "--", color="darkviolet", label="nominal")
    plt.plot(cerpe, best, "-", color="darkorange",  label="best fit")

    plt.legend(loc="upper left")

    plt.grid(True)
    plt.title("Electron NPE Sigma")
    plt.xlabel(r"cerenkov NPE")
    plt.ylabel(r"$\sigma_{NPE}$")

    plt.show()



if __name__ == "__main__":
    main()
