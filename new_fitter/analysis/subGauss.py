import elecLoader as eloader
import prmBetaLoader as bloader
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import ROOT

def main():
    name = "AmC"
    secBetaArr, secAntiBetaArr = bloader.loadPrmBeta("../data/gamma/" + name + "_J19.txt")
    print(">>>>>>>>>> Loading primary beta distribution <<<<<<<<<<")
    
    mu_arr, sigma_arr = [], []

    for i in range(len(secBetaArr)):
        tmppe, tmpsigma = 0, 0
        for j in secBetaArr[i]:
            if j == 0:
                break
            tmppe += eloader.getNPE(j)
            tmpsigma += eloader.getSPE(j)*eloader.getSPE(j)
        for j in secAntiBetaArr[i]:
            if j == 0:
                break
            tmppe += eloader.getNPE(j) + 2*660.8
            tmpsigma += eloader.getSPE(j)*eloader.getSPE(j) + 2*27.02**2

        tmpsigma = np.sqrt(tmpsigma)

        mu_arr.append(tmppe)
        sigma_arr.append(tmpsigma)

    hist = ROOT.TH2D("hist", "", 100, 8300, 9300, 100, 70, 130)
    for i, j in zip(mu_arr, sigma_arr):
        hist.Fill(i, j)

    hist.SaveAs("AmC_model.root")

    #plt.plot(mu_arr, sigma_arr, "o", ms=0.3)
    #plt.xlabel("NPE")
    #plt.ylabel("PE sigma")
    #plt.grid(True)
    #plt.show()



if __name__ == "__main__":
    main()
