import numpy as np
import matplotlib.pyplot as plt
import elecLoader as eloader
import random


# no correlation in fitter...

Etrue = np.arange(0.1, 15, 0.1)
dnum = len(Etrue)

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
    sigma = sample_corelation("../rescov_gam+B12_kSimulationSct_kSimulationCer_kTotal.txt", 3, sampleSize)

    #kA, kB, scale, kC = 0.999805, 6.31677e-3, 1409.38, 0.987722
    #c0, c1, c2, s0 = -123.213, 3.39232, 2.53173e-3, 1
    ra, rb, rc = -1.10053e1, 1.48417e3, 1.21885e2

    ymin, ymax = [], []
    for i in range(dnum):
        ymin.append(1000000)
        ymax.append(-100)

    for i in range(sampleSize):

        #m_kA = kA + random.gauss(0, sigma[0, i])
        #m_kB = kB + random.gauss(0, sigma[1, i])
        #m_scale = scale + random.gauss(0, sigma[2, i])
        #m_kC = kC + random.gauss(0, sigma[3, i])
        #m_c0 = c0 + random.gauss(0, sigma[4, i])
        #m_c1 = c1 + random.gauss(0, sigma[5, i])
        #m_c2 = c2 + random.gauss(0, sigma[6, i])
        #m_s0 = 1

        m_ra = ra * random.gauss(0, sigma[0, i])
        m_rb = rb * random.gauss(0, sigma[1, i])
        m_rc = rc * random.gauss(0, sigma[2, i])

        totsigma = []

        for i in Etrue:

            #sctpe = sctpe_sim(i, m_kA, m_kB, m_scale)
            #cerpe = cerpe_sim(i, m_kC)
            #pp = sctpe / (sctpe + cerpe)
            #cersigma2 = cerFunc(cerpe, m_c0, m_c1, m_c2)**2
            #sctsigma2 = sctFunc(sctpe, m_s0)**2
            #totsigma2 = (sctsigma2+cersigma2) / (1-2*pp*(1-pp))
            totsigma.append( resFunc(i, m_ra, m_rb, m_rc) )
        
            #totsigma.append(totsigma)

        
        for n in range(dnum):
            if totsigma[n] < ymin[n]:
                ymin[n] = totsigma[n]
            if totsigma[n] > ymax[n]:
                ymax[n] = totsigma[n]
    
    ymin = np.array(ymin)
    ymax = np.array(ymax)

    print(ymin[9]**2, ymin[19]**2, ymin[29]**2, ymin[49]**2, ymin[79]**2)
    print(ymax[9]**2, ymax[19]**2, ymax[29]**2, ymax[49]**2, ymax[79]**2)

    best = []
    bc0, bc1, bc2 = -158.41, 4.333, 0.00181082
    for i in Etrue:
        #sctpe = sctpe_sim(i, kA, kB, scale)
        #cerpe = cerpe_sim(i, kC)
        #pp = sctpe / (sctpe + cerpe)
        #cersigma2 = cerFunc(cerpe, bc0, bc1, bc2)**2
        #sctsigma2 = sctFunc(sctpe, s0)**2
        #totsigma2 = (sctsigma2+cersigma2) / (1-2*pp*(1-pp))
      
        #best.append(np.sqrt(totsigma2))
        best.append(resFunc(i, ra, rb, rc)) 

    best = np.array(best)


    plt.fill_between(Etrue, ymin, ymax, color="lightskyblue", label=r"$1 \sigma$ zone")
    plt.plot(Etrue, nominal(), "--", color="darkviolet", label="nominal")
    plt.plot(Etrue, best, "-", color="darkorange",  label="best fit")

    plt.legend(loc="upper left")

    plt.grid(True)
    plt.title("Electron NPE Sigma")
    plt.xlabel(r"$E_{true}/MeV$")
    plt.ylabel(r"$\sigma_{NPE}$")

    plt.show()



if __name__ == "__main__":
    main()
