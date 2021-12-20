import uproot as up
import numpy as np
import matplotlib.pyplot as plt


def fileLoading(filename):
    tt = up.open(filename)["evt"]
    totpe = tt['totalPE'].array()

    return totpe


def resFunc(a, b, c, N):
    return np.sqrt(a+b*N+c*N**2)/N


if __name__ == "__main__" :

    global_es = 3134.078/2.223
    q01, q11, q21 = 0,  9.79299e-01, 6.12574e-05
    
    npe, spe = [], []
    micN, micS, micNerr, micSerr = 7.76588e+04, 8.19608e+02, 5.52924e+01, 2.13277e+02

    for i in [100, 200, 320, 400, 501, 583, 702, 789]:
        filename = "/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/electron/user-"+str(i)+".root"
        totpe = fileLoading(filename)
        npe.append(np.mean(totpe))
        spe.append(np.std(totpe))

    npe = np.array(npe)
    spe = np.array(spe)
    
    nn = np.arange(1000, micN+1000, 100)

    #plt.plot(npe, spe/npe, "o-")
    plt.plot(nn, resFunc(q01, q11, q21, nn), "--")
    plt.errorbar(micN, micS/micN, yerr=np.sqrt(micSerr**2/micN**2 + micNerr**2*micS**2/micN**4), fmt="o")
    plt.show()
