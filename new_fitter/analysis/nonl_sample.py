import numpy as np
import matplotlib.pyplot as plt
import elecLoader as eloader
import random

kA1 = 0.999671
kAerr1 =  1.02112e-04
kC1 = 1.01579
kCerr1 = 2.24912e-03
es1 = 1.38648e+3
eserr1 = 1.41421

kA2 = 0.99644
kAerr2 = 1.02177e-4
kC2 = 1.01674
kCerr2 = 2.25917e-3
es2 = 1.38648e+3
eserr2 = 1.41421


kA3 = 9.9988e-1
kAerr3 = 6.61878e-4
kB3 = 6.34833e-3
kBerr3 = 1.15405e-4
kC3 = 9.84013e-1
kCerr3 = 1.42091e-2
es3 = 1.40961e3
eserr3 = 1.86591



# no correlation in fitter...
# no kB here :
Etrue = np.arange(0.1, 15, 0.1)


#def samplePar(kA, kC, es, kAerr, kCerr, eserr):
#    m_kA = random.gauss(kA, kAerr)
#    m_kC = random.gauss(kC, kCerr)
#    m_es = random.gauss(es, eserr)
#    
#    m_nonl = []
#    for i in Etrue:
#        m_nonl.append( (m_kA * eloader.getSctNPE(i) + m_kC * eloader.getCerNPE(i)) / i / m_es )
#
#    m_nonl = np.array(m_nonl)
#    return m_nonl
#
#
#def bestfit(kA, kC, es):
#    m_kA = kA
#    m_kC = kC
#    m_es = es
#
#    m_nonl = []
#    for i in Etrue:
#        m_nonl.append( (m_kA * eloader.getSctNPE(i) + m_kC * eloader.getCerNPE(i)) / i / m_es )
#
#    m_nonl = np.array(m_nonl)
#    return m_nonl
#
#def nominal(es):
#    m_kA = 1
#    m_kC = 1
#    #m_es = 3134.078/2.223
#    m_es = es
#
#    m_nonl = []
#    for i in Etrue:
#        m_nonl.append( (m_kA * eloader.getSctNPE(i) + m_kC * eloader.getCerNPE(i)) / i / m_es )
#
#    m_nonl = np.array(m_nonl)
#    return m_nonl
#
#
## Add kB fitting here :
#
#
#def samplePar_new(kA, kB, kC, es, kAerr, kBerr, kCerr, eserr):
#    m_kA = random.gauss(kA, kAerr)
#    m_kB = random.gauss(kB, kBerr)
#    m_kC = random.gauss(kC, kCerr)
#    m_es = -1
#    while m_es < 0:
#        m_es = random.gauss(es, eserr)
#    
#    m_nonl = []
#    for i in Etrue:
#        m_nonl.append( (m_kA * eloader.getQPE(i, m_kB, m_es) + m_kC * eloader.getCerNPE(i)) / i / m_es )
#
#    m_nonl = np.array(m_nonl)
#    return m_nonl
#
def bestfit_new(kA, kB, kC, es):
    m_kA = kA
    m_kB = kB
    m_kC = kC
    m_es = es

    m_nonl = []
    for i in Etrue:
        m_nonl.append( (m_kA * eloader.getQPE(i, m_kB, m_es) + m_kC * eloader.getCerNPE(i)) / i / m_es )

    m_nonl = np.array(m_nonl)
    return m_nonl

def nominal_new(es):
    m_kA = 1
    m_kC = 1
    m_kB = 0.0065
    m_es = es

    m_nonl = []
    for i in Etrue:
        m_nonl.append( (m_kA * eloader.getQPE(i, m_kB, m_es) + m_kC * eloader.getCerNPE(i)) / i / m_es )

    m_nonl = np.array(m_nonl)
    return m_nonl



def loadCov(filename):
    row, col = 0, 0
    cov_mat = np.ones((4, 4))
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
    
    num_sample = 5000

    cov = loadCov("../covar1.txt")
    print(cov)
    
    x = norm.rvs(size=(4, num_sample))

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

    #plt.subplot(2, 1, 1)
    #plt.plot(y[0], y[1], "b.")
    #plt.grid(True)

    #plt.subplot(2, 1, 2)
    #plt.plot(y[0], y[2], "b.")
    #plt.grid(True)

    #plt.show()


def main():

    #for i in range(5000):
    #    m_nonl3 = samplePar_new(kA3, kB3, kC3, es3, kAerr3, kBerr3, kCerr3, eserr3)
    #    #m_nonl2 = samplePar(kA2, kC2, es2, kAerr2, kCerr2, eserr2)

    #    plt.plot(Etrue, m_nonl3, "-", color="lightskyblue",  alpha=1)
    #    #plt.plot(Etrue, m_nonl2, "-", color="pink", alpha=0.01)

    #plt.plot(Etrue, nominal_new(es3), "-", color="coral", label="nominal")
    #plt.plot(Etrue, bestfit_new(kA3, kB3, kC3, es3), "--", color="blueviolet", label="best fit")
    ##plt.plot(Etrue, bestfit(kA2, kC2, es2), "--", color="deeppink", label="best fit: gamma+B12")
    #plt.legend()
    #plt.xlabel(r"$E_{true}$/MeV")
    #plt.ylabel(r"$E_{vis}/E_{true}$")
    ##plt.semilogy()
    ##plt.savefig("nonlCurve_fit_B12+gam.pdf")
    #plt.show()

    sigma = sample_corelation()

    for i in range(5000):
        kA = kA3 + sigma[0, i]
        kB = kB3 + sigma[1, i]
        kC = kC3 + sigma[2, i]
        es = es3 + sigma[3, i]

        m_nonl = []
        for i in Etrue:
            m_nonl.append( (kA * eloader.getQPE(i, kB, es) + kC * eloader.getCerNPE(i)) / i / es )
        plt.plot(Etrue, m_nonl, "-", color="lightskyblue",  alpha=1)

    plt.plot(Etrue, nominal_new(es3), "-", color="coral", label="nominal")
    plt.plot(Etrue, bestfit_new(kA3, kB3, kC3, es3), "--", color="blueviolet", label="best fit")
    plt.legend()
    plt.xlabel(r"$E_{true}$/MeV")
    plt.ylabel(r"$E_{vis}/E_{true}$")
    plt.grid(True)
    plt.title("Electron Nonlinearity")
    #plt.semilogx()
    plt.savefig("nonlCurve_fit_B12+gam_kBfit.pdf")
    plt.show()
    


if __name__ == "__main__":
    main()
