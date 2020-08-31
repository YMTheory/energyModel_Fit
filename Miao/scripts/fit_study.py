#!/usr/bin/env python
# coding=utf-8
# study nonlinearity curve smeating with parameter


import matplotlib.pyplot as plt
import numpy as np
import uproot as up

scale = 3350/2.220

def plot_nonl(kA, kB, kC):
    quenchFile = up.open("../data/electron/Quench5.root")
    kbidx = "kB"+str(int(kB*1e4))
    quenchNonl = quenchFile[kbidx]
    value = quenchNonl.values
    etrue = np.arange(0.1, 8, 0.1)
    nonl = []
    for i in etrue:
        if i<0.1:
            idx = int(i/0.001)
        else:
            idx = int(i/0.01)+100
        idx1 = int(i/0.01)
        qNL = kA * value[idx];
        cNL = kC * cer[idx1]/i/scale
        nonl.append( qNL + cNL)
    return nonl

# read cerenkov :
cer = []
with open( "/Users/yumiao/Documents/Works/github/energyModel_Fit/Miao/data/electron/Cer.dat") as f:
    for lines in f.readlines():
        line = lines.strip("\n")
        data = line.split(" ")
        cer.append(float(data[2]))


etrue = np.arange(0.1, 8, 0.1)
nonl_nominal = plot_nonl(0.962, 6.5e-3, 1)
nonl1 = plot_nonl(0.962, 6.5e-3, 0.98)
nonl2 = plot_nonl(0.962, 6.5e-3, 0.99)
nonl3 = plot_nonl(0.962, 6.5e-3, 1.01)
nonl4 = plot_nonl(0.962, 6.5e-3, 1.02)

plt.plot(etrue, nonl_nominal, "-", label="nominal")
plt.plot(etrue, nonl1, "-", label="kC = 0.98")
plt.plot(etrue, nonl2, "-", label="kC = 0.99")
plt.plot(etrue, nonl3, "-", label="kC = 1.01")
plt.plot(etrue, nonl4, "-", label="kC = 1.02")

plt.legend()
plt.xlabel("etrue/MeV")
plt.ylabel("nonlinearity")

plt.show()
