#!/usr/bin/env python
# coding=utf-8

import numpy as np, uproot as up
import matplotlib.pyplot as plt

file = up.open("/Users/yumiao/Documents/Works/Simulation/Nonlinearity/gamma/rootfiles/8MeV_onlyPrimarye-.root")
elecKE = file["evt"].array("ElectronKE")

cont, edge = np.histogram(elecKE.flatten(), bins=800, range=(0, 8.0))
Energy = np.arange(0.01, 8.01, 0.01)

print(len(Energy))

with open("../data/naked_gamma/primary_8MeV.txt", "w") as f:
    for i in range(len(cont)):
        f.write(str(round(Energy[i],2)) + " "+ str(round(cont[i], 2)) + "\n")

plt.hist(elecKE.flatten(), bins=800, range=(0, 8.0), histtype="step")
plt.yscale("log"); plt.xlabel("ElectronKE/MeV")
plt.show()
