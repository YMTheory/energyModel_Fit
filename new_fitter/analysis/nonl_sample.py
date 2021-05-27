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


kA3 = 9.99805e-1
kAerr3 = 6.61878e-4
kB3 = 6.31679e-3
kBerr3 = 1.15405e-4
kC3 = 9.88685e-1
kCerr3 = 1.42091e-2
es3 = 1.40935e3
eserr3 = 1.86591


kA4 = 0.99988
kB4 = 6.34833e-3
kC4 = 0.984013
es4 = 1409.61


kA5 = 0.996993
kAerr5 = 5.1189e-5
kB5 = 5.445e-3
kBerr5 = 2.40804e-5
A25 = -7.34716
A2err5 = 5.41555e-2
A35 = 15.5519
A3err5 = 0.024233
A45 = 0.0267155
A4err5 = 6.043e-4
es5 = 1401.52
eserr5 = 0.143918

kA6 = 9.94149e-01
kAerr6 = 3.32354e-04
kB6 = 5.53304e-03
kBerr6 = 6.30502e-05
A26 = -7.29243e+00
A2err6 = 2.87015e-01
A36 = 1.54759e+01
A3err6 = 1.37904e-01
A46 = 2.69986e-02
A4err6 = 1.92680e-03
es6 = 

# no correlation in fitter...
# no kB here :
Etrue = np.arange(0.1, 15, 0.1)

def cerpe_func(E, A2, A3, A4):
    E0 = 0.2
    A1 = 0
    x = np.log(1+E/E0)
    print(E, A2, A3, A4, (A1*x+A2*x**2+A3*x**3)*(1/E+A4)*E)
    return (A1*x+A2*x**2+A3*x**3) * (1/E+A4) * E 

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
        #m_nonl.append( (m_kA * eloader.getQPE(i, m_kB, m_es) + m_kC * eloader.getCerNPE(i)) / i / m_es )
        m_nonl.append( (m_kA * eloader.getQPE(i, m_kB, m_es) + m_kC * eloader.getCerNPE(i)) )

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
    cov_mat = np.ones((6, 6))
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


def sample_corelation(filename):
    method = 'eigenvectors'
    
    num_sample = 5000

    cov = loadCov(filename)
    print(cov)
    
    x = norm.rvs(size=(6, num_sample))

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


    sigma1 = sample_corelation("../output/params/covar1_nonlres_gamB12.txt")
    #sigma2 = sample_corelation("../covar1.txt")

    for i in range(5000):
        kA = kA6 + sigma1[0, i]
        kB = kB6 + sigma1[1, i]
        A2 = A26 + sigma1[2, i]
        A3 = A36 + sigma1[3, i]
        A4 = A46 + sigma1[4, i]
        es = es6 + sigma1[5, i]

        m_nonl = []
        for i in Etrue:
            m_nonl.append( (kA * eloader.getQPE(i, kB, es) + cerpe_func(i, A2, A3, A4) ) )
        plt.plot(Etrue, m_nonl, "-", color="lightskyblue",  alpha=0.1)

    m_nonl_bestfit1 = []
    for i in Etrue:
        m_nonl_bestfit1.append( kA6 * eloader.getQPE(i, kB6, es6) + cerpe_func(i, A26, A36, A46)  )
        print(m_nonl_bestfit1[-1], m_nonl_bestfit2[-1])

    plt.plot(Etrue, m_nonl_bestfit1, "--", color="coral", label="best fit: simulation cerenkov")

    
    #plt.plot(Etrue, nominal_new(es3), "-", color="coral", label="nominal")
    plt.legend()
    plt.xlabel(r"$E_{true}$/MeV")
    plt.ylabel(r"$E_{vis}/E_{true}$")
    plt.grid(True)
    plt.title("Electron Nonlinearity")
    #plt.semilogx()
    plt.savefig("nonlCurve_fit_B12+gam_kBfit.pdf")
    plt.show()


    """
    #for i in range(100):
    #    kA = random.gauss(kA5, kAerr5)
    #    kB = random.gauss(kB5, kBerr5)
    #    A2 = random.gauss(A25, A2err5)
    #    A3 = random.gauss(A35, A3err5)
    #    A4 = random.gauss(A45, A4err5)
    #    es = random.gauss(es5, eserr5)
    #    m_nonl = []
    #    for i in Etrue:
    #        m_nonl.append( (kA * eloader.getQPE(i, kB, es) + cerpe_func(i, A2, A3, A4)) / i / es )
    #    plt.plot(Etrue, m_nonl, "-", color="lightseagreen",  alpha=0.1)


    m_nonl_bestfit1 = []
    m_nonl_bestfit2 = []
    for i in Etrue:
        m_nonl_bestfit2.append( kA5 * eloader.getQPE(i, kB5, es5) + cerpe_func(i, A25, A35, A45)  )
        m_nonl_bestfit1.append( kA4 * eloader.getQPE(i, kB4, es4) + kC4 * eloader.getCerNPE(i) )
        print(m_nonl_bestfit1[-1], m_nonl_bestfit2[-1])
        #m_nonl_bestfit1.append(kC4 * eloader.getCerNPE(i) / i)
        #m_nonl_bestfit2.append(cerpe_func(i, A25, A35, A45) / i )



    plt.plot(Etrue, m_nonl_bestfit2, "-", color="darkviolet", label="best fit: analytical cerenkov")
    plt.plot(Etrue, m_nonl_bestfit1, "--", color="coral", label="best fit: simulation cerenkov")
    #plt.plot(Etrue, bestfit_new(kA4, kB4, kC4, es4), "--", color="coral", label="best fit")
    #plt.plot(Etrue, nominal_new(es5), "-", color="coral", label="nominal")
    plt.legend()
    plt.xlabel(r"$E_{true}$/MeV")
    plt.ylabel("total NPE")
    plt.grid(True)
    #plt.title("Electron Nonlinearity")
    plt.semilogy()
    plt.show()
    """


if __name__ == "__main__":
    main()
