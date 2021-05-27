import numpy as np
import matplotlib.pyplot as plt
import elecLoader as eLoader

import ROOT

ff = ROOT.TFile("../data/electron/Quench5.root", "read")
hist1 = ff.Get("kB55")
hist2 = ff.Get("kB65")
hist3 = ff.Get("kB75")

E1, n1, E2, n2, E3, n3 = [], [], [], [], [], []
E0, n0 = [], []

N = hist1.GetNbinsX()
for i in range(N):
    E1.append(hist1.GetBinCenter(i))
    n1.append(hist1.GetBinContent(i))
    E2.append(hist2.GetBinCenter(i))
    n2.append(hist2.GetBinContent(i))
    E3.append(hist3.GetBinCenter(i))
    n3.append(hist3.GetBinContent(i))

ff = ROOT.TFile("../data/electron/Quench_NumInt.root", "read")
hist4 = ff.Get("kB55")
hist5 = ff.Get("kB65")
hist6 = ff.Get("kB75")

E4, n4, E5, n5, E6, n6 = [], [], [], [], [], []

N = hist4.GetNbinsX()
for i in range(N):
    E4.append(hist4.GetBinCenter(i))
    n4.append(hist4.GetBinContent(i))
    E5.append(hist5.GetBinCenter(i))
    n5.append(hist5.GetBinContent(i))
    E6.append(hist6.GetBinCenter(i))
    n6.append(hist6.GetBinContent(i))

#E0 = np.arange(0.01, 10, 0.01)
#for i in E0:
#    n0.append(eLoader.Integral_BirkLaw(i))

#A1, A2, A3, A4, A5 = 1.019, 0.127, 6.6067e5, 0.117, 0.007
#A1, A2, A3, A4, A5, A6 = 1.31781e+00, 2.27665e-01, 1.72301e-01, 1.72301e-01, 3.93126e-01, 3.93126e-01
#Qeff = (A1+A2*np.log(E0)+A3*np.log(E0)**2 + A4*np.log(E0)**3 ) / (1+A5*np.log(E0)+A6*np.log(E0)**2+A4*np.log(E0)**3 )

#plt.plot(E0, Qeff, "-")

plt.plot(E1, n1, "-", color="royalblue", label=r"$k_B=0.0055$")
plt.plot(E2, n2, "-", color="peru", label=r"$k_B=0.0065$")
plt.plot(E3, n3, "-", color="orange", label=r"$k_B=0.0075$")
plt.plot(E4, n4, "--", color="royalblue", label=r"$k_B=0.0055$")
plt.plot(E5, n5, "--", color="peru", label=r"$k_B=0.0065$")
plt.plot(E6, n6, "--", color="orange", label=r"$k_B=0.0075$")
#plt.plot(E0, n0, "--", label=r"$k_A=0.0065$ Numerical Integral")

plt.legend()
plt.xlabel(r"$E_{true}/MeV$")
plt.ylabel(r"$E_{quenched}/E_{true}$")
#plt.xlim(0.1, 10)
#plt.ylim(0.8, 1.0)
plt.semilogx()

plt.grid(True)
plt.show()
