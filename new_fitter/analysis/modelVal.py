import numpy as np
import matplotlib.pyplot as plt
import ROOT
import prmBetaLoader as bloader
import elecLoader as eloader


def main():

    es = 3134.078/2.223

    for k in range(315, 380, 1):
        # Load PrmBeta
        filename = "../data/gamma/Livermore-"+str(k)+".txt"
        print(filename)
        prmBeta, prmAntiBeta = bloader.loadPrmBeta(filename)

        # calcSingleEvent :
        nSamples = 5000
        mu_arr, sigma_arr = [], []
        for i in range(nSamples):
            tmpnpe, tmpspe = 0, 0
            for j in prmBeta[i]:
                if j == 0:
                    break
                tmpnpe += (eloader.getNPE(j))
                tmpspe += (eloader.getSPE(j)**2)
            for j in prmAntiBeta[i]:
                if j == 0:
                    break
                tmpnpe += (eloader.getNPE(j) + 2*660.8)
                tmpspe += (eloader.getSPE(j)**2 + 27.07**2*2)

            mu_arr.append(tmpnpe)
            sigma_arr.append(np.sqrt(tmpspe))

        mu = np.array(mu_arr)
        sigma = np.array(sigma_arr)
    

        npe, spe = 0, 0
        sigma_part1, sigma_part2 = 0, 0
        for i in mu:
            npe += i
        npe = npe/nSamples
        for i, j in zip(mu, sigma):
            spe += (i-npe)**2 + j**2
            sigma_part1 += (i-npe)**2
            sigma_part2 += j**2

        spe = np.sqrt(spe/nSamples)
        sigma_part1 /= nSamples
        sigma_part2 /= nSamples

        print((k-300)/10., npe/es, spe/es)



if __name__ == "__main__" :
    main()



