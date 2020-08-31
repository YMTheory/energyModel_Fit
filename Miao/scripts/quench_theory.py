#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

E_init = 0;
chi = 5.1e-29 # MeV*m2
mass = 0.511 # MeV
n = 2.7e29  #m-3, from DC LS
T = 56.4e-6 #MeV from DC LS
R = 10e-9 # m, should be a fitting parameters.
normalize = 0.00011;

energy = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2]
pred_NL = [];
binResol = 0.0001;
for E in energy:
    inteE = 0.0001
    sum = 0
    while inteE+binResol <= E:
        gamma_low = (inteE+mass)/mass; beta2_low = 1-1/gamma_low/gamma_low
        sigma_low = chi/beta2_low/T*(np.log(2*mass*mass*beta2_low*gamma_low**2/T)-beta2_low)
        dEdx_low = chi/beta2_low*n*(np.log(inteE**2*(1+gamma_low)/2/T**2)-beta2_low+1-(2*gamma_low-1)/gamma_low**2*np.log(2)+1/8*(gamma_low-1)**2/gamma_low**2)
        #print(n*R*sigma_low)
        exp_part_low = np.exp(-n*R*sigma_low)
        value_low = sigma_low/(dEdx_low)*exp_part_low
        #print(str(inteE) + " " + str(sigma_low) + " " +str(dEdx_low) + " " + str(exp_part_low))

        gamma_high = (inteE+mass+binResol)/mass; beta2_high = 1-1/gamma_high/gamma_high;
        sigma_high = chi/beta2_high/T*(np.log(2*mass*mass*beta2_high*gamma_high**2/T)-beta2_high)
        dEdx_high  = chi/beta2_high*n*(np.log((inteE+binResol)**2*(1+gamma_high)/2/T**2)-beta2_high+1-(2*gamma_high-1)/gamma_high**2*np.log(2)+1/8*(gamma_high-1)**2/gamma_high**2)
        exp_part_high = np.exp(-n*R*sigma_high)
        value_high = sigma_high/(dEdx_high)*exp_part_high

        area = (value_low + value_high) * binResol /2.;
        sum += area
        inteE = inteE + binResol
    pred_NL.append(sum/E*normalize*n)

# read nonl data :
etrue = []; data_NL = []; scale = 3350/2.22
with open("/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/Cerenkov/totPE_smearing.txt") as f:
    for lines in f.readlines():
        line = lines.strip("\n")
        data = line.split(" ")
        if(float(data[0])<=20):
            continue
        if (float(data[0]) > 200):
            break
        etrue.append(float(data[0])/1000.)
        data_NL.append(float(data[1])/scale/etrue[-1])



plt.plot(energy, pred_NL, "-",lw=2, label="theory")
plt.plot(etrue, data_NL, "o", ms=3,  label="data", color=".red")
plt.legend()
plt.show()
