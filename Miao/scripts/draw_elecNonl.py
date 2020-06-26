#!/usr/bin/env python
# coding=utf-8

import matplotlib.pyplot as plt
import numpy as np
import uproot as up

scale = 3350/2.220

# read quenching nonlinearity curve
quenchFile = up.open("../data/electron/Quench5.root")
quenchNonl = quenchFile["kB61"]
quenchNonl_up = quenchFile["kB53"]
quenchNonl_low = quenchFile["kB69"]
quenchNonl_nominal = quenchFile["kB65"]
quenchNonl_up_gamma = quenchFile["kB63"]
quenchNonl_low_gamma = quenchFile["kB58"]

edge = quenchNonl.edges
value = quenchNonl.values
value_low = quenchNonl_low.values
value_up = quenchNonl_up.values
value_nominal  = quenchNonl_nominal.values
value_low_gamma = quenchNonl_low_gamma.values
value_high_gamma = quenchNonl_up_gamma.values

binCenter = []
for i in range(edge.size-1):
    binCenter.append(edge[i])
binCenter = np.array(binCenter)


# read cerenkov :
cer = []
with open( "/Users/yumiao/Documents/Works/github/energyModel_Fit/Miao/data/electron/Cer.dat") as f:
    for lines in f.readlines():
        line = lines.strip("\n")
        data = line.split(" ")
        cer.append(float(data[2]))


kA_nominal = 0.962; kC_nominal = 1.0
kA = 0.961; kC = 1.0;
kA_up = 0.966; kC_up = 1.060;
kA_low = 0.955; kC_low = 0.939;
kA_gamma = 9.65736e-01; kC_gamma = 9.71846e-01;
kA_gamma_low = 9.65736e-01-3.45386e-03; kC_gamma_low = 9.71846e-01-6.90211e-02;
kA_gamma_high = 9.65736e-01+3.45386e-03; kC_gamma_high = 9.71846e-01+6.90211e-02;

etrue = np.arange(0.1, 8, 0.1)
nonl = []; nonl_low = []; nonl_up = []; nonl_nominal = [];
nonl_low_gamma = []; nonl_high_gamma = []; nonl_gamma = [];
for i in etrue:
    if i<0.1:
        idx = int(i/0.001)
    else:
        idx = int(i/0.01)+100
    idx1 = int(i/0.01)
    qNL = kA * value[idx];
    cNL = kC * cer[idx1]/i/scale
    nonl.append( qNL + cNL)
    #nonl_low.append(kA*value_low[idx]+kC*cer[idx1])
    nonl_low.append(kA_low*value_low[idx]+kC_low*cer[idx1]/i/scale)
    nonl_up.append(kA_up*value_up[idx]+kC_up*cer[idx1]/i/scale)
    nonl_nominal.append(kA_nominal*value_nominal[idx]+kC_nominal*cer[idx1]/i/scale)
    nonl_gamma.append(kA_gamma*value[idx]+kC_gamma*cer[idx1]/i/scale)
    nonl_low_gamma.append(kA_gamma_low*value_low_gamma[idx]+kC_gamma_low*cer[idx1]/i/scale)
    nonl_high_gamma.append(kA_gamma_high*value_high_gamma[idx]+kC_gamma_high*cer[idx1]/i/scale)

plt.plot(etrue, nonl_nominal, color="forestgreen", label="nominal")
plt.plot(etrue, nonl, "-.", color="chocolate", label="electron best fit") #, color='blue')
plt.plot(etrue, nonl_gamma, color="blue", label="gamma best fit")
#plt.plot(etrue, nonl_low, "-", color="blue", alpha=0.25) #, color='blue')
#plt.plot(etrue, nonl_up, "-",color="blue", alpha=0.25) #, color='blue')
#plt.text(3, 0.96, "best fit valuse:",fontsize=14)
#plt.text(3,0.952, 'kA=0.961+-0.0060')
#plt.text(3,0.944, "kB=0.00613+-0.000885")
#plt.text(3,0.936, "kC=1.000+-0.0605")
plt.text(3, 0.96, "best fit valuse:",fontsize=14)
plt.text(3,0.952, 'kA=9.65736e-01+-3.45386e-03')
plt.text(3,0.944, "kB=0.00612+-2.30433e-04")
plt.text(3,0.936, "kC=9.71846e-01+-6.90211e-02")
plt.text(0.2,1.01, "1$\sigma$ zone", fontsize=13, color="purple")

plt.legend(loc="center right", fontsize=9)
plt.xlabel("electron etrue/MeV")
plt.ylabel("nonlinearity")
a, b = 0, 8  # integral limits
xf = etrue[np.where((etrue>a)&(etrue<b))]
#plt.fill_between(xf, nonl_up, nonl_low, color='blue',alpha=0.25)
plt.fill_between(xf, nonl_high_gamma, nonl_low_gamma, color='purple',alpha=0.35)

plt.show()
