import numpy as np
import matplotlib.pyplot as plt
import elecLoader as eloader
import random

p01 = -8.31736
p0err1 = 8.88308e-1
p11 = 1.43134e3
p1err1 = 1.07852e1
p21 = 1.32963e2
p2err1 = 7.29333


p02 = -8.32822e0
p0err2 = 8.8812e-1
p12 = 1.43146e3
p1err2 = 1.07811e1
p22 = 1.32294e2
p2err2 = 7.70714


p03 = -1.09984e1
p13 = 1.48408e3
p23 = 1.21907e2


# no correlation in fitter...

Etrue = np.arange(0.1, 10, 0.1)

def resFunc(x, a, b, c):
    s2 = a + b*x + c*x**2
    if s2<0:
        return 0
    else:
        return np.sqrt(s2)

def samplePar(p0, p1, p2, p0err, p1err, p2err):
    m_p0 = random.gauss(p0, p0err)
    m_p1 = random.gauss(p1, p1err)
    m_p2 = random.gauss(p2, p2err)
    
    m_nonl = []
    for i in Etrue:
        m_nonl.append(resFunc(i, m_p0, m_p1, m_p2) )

    m_nonl = np.array(m_nonl)
    return m_nonl


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

def loadCov(filename):
    row, col = 0, 0
    cov_mat = np.ones((3, 3))
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


def sample_corelation():
    method = 'eigenvectors'
    
    num_sample = 50000

    cov = loadCov("../covar2.txt")
    print(cov)
    
    x = norm.rvs(size=(3, num_sample))

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
    """
    for i in range(50000):
        print(i)
        m_nonl1 = samplePar(p01, p11, p21, p0err1, p1err1, p2err1)
        m_nonl2 = samplePar(p02, p12, p22, p0err2, p1err2, p2err2)

        plt.plot(Etrue, m_nonl1, "-", color="lightskyblue", alpha=0.01)
        plt.plot(Etrue, m_nonl2, "-", color="pink", alpha=0.01)
    """
    
    sigma = sample_corelation()

    for i in range(50000):
        p0 = p03 + sigma[0, i]
        p1 = p13 + sigma[1, i]
        p2 = p23 + sigma[2, i]

        m_res = []
        for i in Etrue:
            m_res.append(resFunc(i, p0, p1, p2) )

        plt.plot(Etrue, m_res, "-", color="lightskyblue",  alpha=1)

    plt.plot(Etrue, nominal(), "-", color="coral", label="nominal")
    #plt.plot(Etrue, bestfit(p01, p11, p21), "--", color="blueviolet", label="best fit: only gamma")
    #plt.plot(Etrue, bestfit(p02, p12, p22), "--", color="deeppink", label="best fit: gamma+B12")
    plt.plot(Etrue, bestfit(p03, p13, p23), "--", color="deeppink", label="best fit: gamma+B12")

    plt.grid(True)
    plt.legend()
    plt.xlabel(r"$E_{true}$/MeV")
    plt.ylabel("p.e. resolution")
    plt.title("Electron P.E. Sigma")
    #plt.savefig("resCurve_fit.pdf")
    plt.show()


if __name__ == "__main__":
    main()
