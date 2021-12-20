import numpy as np
import matplotlib.pyplot as plt

def func(p1, p2, N):
    #return np.sqrt(p1/N + p2) # resolution
    return p1*N + p2*N*N

def loadCovFile(filename):
    E, cerpe, sigma, sigma_err = [], [], [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            E.append(float(data[0]))
            cerpe.append(float(data[1]))
            sigma.append(float(data[2]))
            sigma_err.append(float(data[3]))

    E = np.array(E)
    cerpe = np.array(cerpe)
    sigma = np.array(sigma)
    sigma_err = np.array(sigma_err)

    return E, cerpe, sigma, sigma_err

import ROOT

if __name__ == "__main__" :

    """
    p1    = 0.972935
    p1err = 0.0142004
    p2    = 6.23307e-05   
    p2err = 1.07871e-06

    es = 3134.078/2.223
    
    N = np.arange(0, 20000, 10) 
    res1 = func(p1, p2, N)

    p1 = 1
    p2 = 0
    res2 = func(p1, p2, N)


    fig, ax = plt.subplots()
    #ax.plot(N/es, res1, "-", color="blue", label="Simulation data")
    #ax.plot(N/es, res2, "-", color="orange", label="Poisson statistics")
    ax.plot(N/es, (res1-res2)/res2, color="blue")

    plt.show()
    """

    """"
    E, cerpe, sigma, sigma_err = loadCovFile("../data/electron/elecCerPEResol1.txt")
    tmp_cerpe, tmp_cerpesigma, tmp_sigmaerr = [], [], []
    for i in np.arange(100, 800, 50):
        tmp_cerpe.append(cerpe[i])
        tmp_cerpesigma.append(sigma[i])
        tmp_sigmaerr.append(sigma_err[i])


    ge = ROOT.TGraphErrors()

    for i in range(len(tmp_cerpe)):
        ge.SetPoint(i, tmp_cerpe[i], tmp_cerpesigma[i]**2)
        ge.SetPointError(i, 0, 2*tmp_cerpesigma[i]*tmp_sigmaerr[i])

    f1 = ROOT.TF1("f1", "[0]*x + [1]*x*x", 0, 1000)
    ge.Fit(f1, "RE")

    """


    q01, q11, q21 = 0,  9.79299e-01, 6.12574e-05
    q03, q13, q23 = 0, 1.00710, 3.98696e-05

    elec, posi = [], []
    npe = np.arange(1000, 15000, 100)
    for i in npe:
        elec.append(np.sqrt(q11*i+q21*i**2)/i)
        posi.append(np.sqrt(q13*i+q23*i**2)/i)

    plt.plot(npe, elec)
    plt.plot(npe, posi)

    plt.show()


