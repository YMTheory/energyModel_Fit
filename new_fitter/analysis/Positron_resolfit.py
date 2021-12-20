import ROOT
import numpy as np
import matplotlib.pyplot as plt

filename = '/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/positron/positronNPE.txt'
ke, mu, muerr, sigma, sigmaerr = [], [], [], [], []
with open(filename) as f:
    for lines in f.readlines():
        line = lines.strip("\n")
        data = line.split(" ")

        ke.append(float(data[0]))
        mu.append(float(data[1]))
        muerr.append(float(data[2]))
        sigma.append(float(data[3]))
        sigmaerr.append(float(data[4]))

ke = np.array(ke)
mu = np.array(mu)
muerr = np.array(muerr)
sigma = np.array(sigma)
sigmaerr = np.array(sigmaerr)

es = 3134.078/2.223

gr = ROOT.TGraphErrors()
for i in range(len(ke)):
    evis = mu[i]/es
    res = sigma[i] / mu[i]
    reserr = np.sqrt(sigmaerr[i]**2/mu[i]**2 + muerr[i]**2 * sigma[i]**2 / mu[i]**4)
    gr.SetPoint(i, mu[i], sigma[i])
    gr.SetPointError(i, 0, sigmaerr[i])



f1 = ROOT.TF1("f1", "[1]*x+[2]*x*x", 0, 15000)
gr.Fit(f1, "RE")






