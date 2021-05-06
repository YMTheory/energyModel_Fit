import numpy as np
import matplotlib.pyplot as plt

import ROOT

ff = ROOT.TFile("../data/electron/Quench5.root", "read")
hist1 = ff.Get("kB55")
hist2 = ff.Get("kB65")
hist3 = ff.Get("kB75")

E1, n1, E2, n2, E3, n3 = [], [], [], [], [], []

N = hist1.GetNbinsX()
for i in range(N):
    E1.append(hist1.GetBinCenter(i))
    n1.append(hist1.GetBinContent(i))
    E2.append(hist2.GetBinCenter(i))
    n2.append(hist2.GetBinContent(i))
    E3.append(hist3.GetBinCenter(i))
    n3.append(hist3.GetBinContent(i))


plt.plot(E1, n1, "-", label=r"$k_B=0.0055$")
plt.plot(E2, n2, "-", label=r"$k_B=0.0065$")
plt.plot(E3, n3, "-", label=r"$k_B=0.0075$")

plt.legend()
plt.xlabel(r"$E_{true}/MeV$")
plt.ylabel(r"$E_{quenched}/E_{true}$")
plt.semilogx()

plt.grid(True)
plt.show()
