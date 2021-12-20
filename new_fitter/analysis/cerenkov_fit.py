import numpy as np
import matplotlib.pyplot as plt
import elecLoader as el
from scipy.optimize import curve_fit

def readFile(filename):
    Etrue, totpe, sigma, sigmaErr = [], [], [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Etrue.append(float(data[0]))
            totpe.append(float(data[1]))
            sigma.append(float(data[2]))
            sigmaErr.append(float(data[3]))

    Etrue = np.array(Etrue)
    totpe = np.array(totpe)
    sigma = np.array(sigma)
    sigmaErr = np.array(sigmaErr)

    return Etrue, totpe, sigma, sigmaErr


def func(x, a, b, c):
    a = 0
    y = a*a+b*b*x+c*c*x**2
    if y.any() < 0:
        return 0
    else:
        return y



if __name__  == "__main__" :

    es = 3134.078/2.223

    covE, covPE, cov, covErr = readFile("../data/electron/elecPECov1.txt")
    sctE, sctPE, sctSigma, sctSigmaErr = readFile("../data/electron/elecSctPEResol1.txt")
    cerE, cerPE, cerSigma, cerSigmaErr = readFile("../data/electron/elecCerPEResol1.txt")

    E, correlation, correlation_err = [], [], []

    sctNPE, cerNPE, covArr, covErrArr = [], [], [], []
    sctSPE, cerSPE, sctSPEerr, cerSPEerr = [], [], [], []

    for i in range(55, len(covE), 10):
        if sctSigma[i] == 0 or cerSigma[i] == 0:
            continue
        E.append(covE[i])
        sctNPE.append(sctPE[i])
        cerNPE.append(cerPE[i])
        sctSPE.append(sctSigma[i])
        cerSPE.append(cerSigma[i])
        sctSPEerr.append(sctSigmaErr[i])
        cerSPEerr.append(cerSigmaErr[i])
        covArr.append(cov[i])
        covErrArr.append(covErr[i])
        correlation.append(cov[i]/sctSigma[i]/cerSigma[i])
        correlation_err.append( np.sqrt( covErr[i]**2/(sctSigma[i]*cerSigma[i])**2  \
            + (cov[i]/sctSigma[i])**2 * cerSigmaErr[i]**2/cerSigma[i]**4 +          \
            + (cov[i]/cerSigma[i])**2 * sctSigmaErr[i]**2/sctSigma[i]**4 ) )

    covArr = np.array(covArr)
    cerNPE = np.array(cerNPE)
    sctNPE = np.array(sctNPE)
    cerSPE = np.array(cerSPE)
    cerSPEerr = np.array(cerSPEerr)
    sctSPE = np.array(sctSPE)
    sctSPEerr = np.array(sctSPEerr)
    correlation = np.array(correlation)
    correlation_err = np.array(correlation_err)



    popt, pcov = curve_fit(func, (cerNPE), cerSPE**2)
    #popt, pcov = curve_fit(func, cerNPE+sctNPE, covArr, p0=[0, 0.203, 0.0025])
    print(popt)
    #popt = [0, 5.31798e-05, 9.524e-2]
    drawy = []
    drawx = np.arange(0, 1000,10)
    for i in drawx:
        drawy.append(func(i, *popt))


    plt.plot((cerNPE), cerSPE**2, "o")
    #plt.plot((sctNPE+cerNPE)/es, cerSPE**2/es**2, "o")
    plt.plot(drawx, drawy)

    plt.show()




