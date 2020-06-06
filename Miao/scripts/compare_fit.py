#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

kB_B12= [6.24499e-03]
kC_B12 = [1.01112e+00]
kBerr_B12 = [8.57618e-09]
kCerr_B12 = [2.18909e-06]

kB_gam = [7.94602e-03]
kC_gam = [1.19391e+00]
kBerr_gam = [1.38653e-03]
kCerr_gam  = [5.44589e-02]

kB_C11 = [6.10118e-03]
kC_C11 = [ 9.00125e-01]
kBerr_C11 = [1.30678e-03]
kCerr_C11 = [2.46137e-01]

kB_C10 = [6.65394e-03]
kC_C10 = [1.19106e+00]
kBerr_C10 = [1.01946e-05]
kCerr_C10 = [6.62637e-03]

kB = [0.0065]
kC = [1.0]

#kB_all = [6.56204e-03]
#kC_all = [1.02364e+00]
#kBerr_all = [5.80335e-09]
#kCerr_all = [1.06140e-06]

kB_all = [7.50885e-03]
kC_all = [1.24487e+00]
kBerr_all = [1.70448e-06]
kCerr_all = [1.73198e-04]


plt.errorbar(kB_B12, kC_B12 , xerr=kBerr_B12, yerr=kCerr_B12, fmt="o", label="B12")
plt.errorbar(kB_gam, kC_gam , xerr=kBerr_gam, yerr=kCerr_gam, fmt="s", label="gamma")
plt.errorbar(kB_C11, kC_C11 , xerr=kBerr_C11, yerr=kCerr_C11, fmt="^", label="C11")
plt.errorbar(kB_C10, kC_C10 , xerr=kBerr_C10, yerr=kCerr_C10, fmt="v", label="C10")
plt.errorbar(kB_all, kC_all , xerr=kBerr_all, yerr=kCerr_all, fmt="p", label="combined only gamma+B12",color="tab:pink")
plt.plot(kB, kC, "*", ms=10, label="sim(e-)", color="tab:purple")

plt.legend(loc="lower right"); plt.grid(True)
plt.xlim(2e-3, 14e-3)
plt.ylim(0.1,1.9)
plt.xlabel("kB")
plt.ylabel("kC")
plt.title("Fitting Parameters")
plt.show()
#plt.savefig("fitRes.png")