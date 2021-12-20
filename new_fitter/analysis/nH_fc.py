from singleGamma import singleGamma
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import prmBetaLoader as bloader
import elecLoader as el



def wenNcer(E, p0, p1, p2, p3, p4):
    if E<=0.2:
        return 0
    if E>0.2:
        E = E - 0.2
        return p3*E**2/(p4+p0*E+p1*np.exp(-p2*E))



def fitNsct(E, As, kB):
    if E > 15.9:
        return el.getFitNsct(15.9, kB, As, "Sim") / 15.9 * E
    else:
        return el.getFitNsct(E, kB, As, "Sim")

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


def sigma2FuncNew(a, b, n, x):
    s2 = a**2*x + b**2*np.power(x, n)
    return s2 if s2>0  else 0


if __name__ == "__main__" :

    nH = singleGamma("nH")
    Y = 3134.078 / 2.223
    filename = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamNewB12_kSimQ_kNewAnaCer_kNew/NewgamNewB12_kSimQ_kNewAnaCer_kNew_rescov.txt"
    par1, parerr = loadBestFit(filename, 3)
    a, b, n = par1

    filename = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamNewB12_kSimQ_kNewAnaCer_kNew/NewgamNewB12_kSimQ_kNewAnaCer_kNew_nonlcov.txt"
    par1, parerr1 = loadBestFit(filename, 7)
    scale, kB, p0, p1, p2, p3, p4 = par1


    filename = "../data/gamma/nH_J19.txt"
    prmBeta, prmAntiBeta = bloader.loadPrmBeta(filename)

    # calcSingleEvent :
    nSamples = 5000
    mu_arr, sigma_arr = [], []
    cer_arr = []
    for i in range(nSamples):
        tmpnpe, tmpspe = 0, 0
        cernpe = 0
        # Electrons :
        for j in prmBeta[i]:
            if j == 0:
                break
            onenpe  = (fitNsct(j, scale, kB) + wenNcer(j, p0, p1, p2, p3, p4))
            cernpe += wenNcer(j, p0, p1, p2, p3, p4)
            tmpnpe += onenpe
            tmpspe += sigma2FuncNew(a, b, n, onenpe)
        # Positrons:
        for j in prmAntiBeta[i]:
            if j == 0:
                break
            onenpe  = (fitNsct(j, scale, kB) + wenNcer(j, p0, p1, p2, p3, p4))
            tmpnpe += (onenpe + 2*660.8)
            cernpe += (wenNcer(j, p0, p1, p2, p3, p4) + 2*1.278)
            tmpspe += sigma2FuncNew(a, b, n, onenpe) 
            tmpspe += 2 * 27.07**2

        mu_arr.append(tmpnpe)
        cer_arr.append(cernpe)
        sigma_arr.append(np.sqrt(tmpspe))

    mu = np.array(mu_arr)
    cer_arr = np.array(cer_arr)
    sigma = np.array(sigma_arr)
    

    tmp_npe, tmp_spe = 0, 0
    tmp_cerpe = 0
    tmp_sigma_part1, tmp_sigma_part2 = 0, 0
    for i in mu:
        tmp_npe += i
    tmp_npe = tmp_npe/nSamples
    for i in cer_arr:
        tmp_cerpe += i
    tmp_cerpe /= nSamples
    for i, j in zip(mu, sigma):
        tmp_spe += (i-tmp_npe)**2 + j**2
        tmp_sigma_part1 += (i-tmp_npe)**2
        tmp_sigma_part2 += j**2

    tmp_spe = np.sqrt(tmp_spe/nSamples)
    tmp_sigma_part1 /= nSamples
    tmp_sigma_part2 /= nSamples

    
    print(tmp_npe, tmp_cerpe, tmp_spe)



