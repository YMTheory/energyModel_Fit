#!/usr/bin/env python
# -*- coding=utf8 -*-
"""
# Author: MiaoYu ---> miaoyu@ihep.ac.cn
# Created Time : Thu Sep  3 16:04:42 2020
# File Name: double_gammaNL_fit.py
"""

import uproot as up
import numpy as np
import matplotlib.pyplot as plt

quench_file = "../data/electron/Quench5.root"
cer_file = "../data/electron/Cer.dat"

scale = 3350/2.22
cer_pe = []

def read_quench(kA, kB, E):
    infile = up.open(quench_file)
    if kB>74e-4:
        kB = 74e-4
    if kB<54e-4:
        kB=54e-4
    kbidx_low  = int(kB*1e4)
    kbidx_high =  kbidx_low+1
    kBResid = kbidx_high - kB*1e4
    kBResid = kbidx_high - kB*1e4
    conts_low, edges_low = infile["kB%d"%kbidx_low].numpy()
    conts_high, edges_high = infile["kB%d"%kbidx_high].numpy()

    idx = 0
    if E<0.1:
        idx = int(E/0.001)
    else:
        idx = int((E-0.1)/0.01)+100
    
    return kA*(kBResid*conts_low[idx] + (1-kBResid)*conts_high[idx])


def read_Cereknov():
    with open(cer_file) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            cer_pe.append(float(data[2]))

    

if __name__ == "__main__" : 
    
    read_Cereknov()    
    Elist1 = [i*0.1 for i in range(1, 79)]
    NLlist1, NLlist1_low, NLlist1_high  = [], [], []
    kA = 9.57520e-01
    kB = 7.19976e-03
    kC = 1.01999e+00
    kAerr = 1.86851e-03
    kBerr = 1.06004e-03
    kCerr = 3.67124e-02
    for Etru in Elist1:
        quenchpart = read_quench(kA, kB, Etru)
        cerpart = kC*cer_pe[int(Etru/0.01)]/Etru/scale

        quenchpart_low = read_quench(kA-kAerr, kB+kBerr, Etru)
        cerpart_low = (kC-kCerr)*cer_pe[int(Etru/0.01)]/Etru/scale

        quenchpart_high = read_quench(kA+kAerr, kB-kBerr, Etru)
        cerpart_high = (kC+kCerr)*cer_pe[int(Etru/0.01)]/Etru/scale

        totNL = quenchpart + cerpart
        NLlist1.append(totNL)
        NLlist1_low.append(quenchpart_low+cerpart_low)
        NLlist1_high.append(quenchpart_high+cerpart_high)



    kA = 9.62226e-01
    kB = 6.75940e-03
    kC = 9.80041e-01
    kAerr = 1.01195e-02
    kBerr = 1.56760e-04
    kCerr = 2.02534e-02
    Elist2 = [i*0.1 for i in range(1, 79)]
    NLlist2, NLlist2_low, NLlist2_high  = [], [], []
    for Etru in Elist2:
        quenchpart = read_quench(kA, kB, Etru)
        cerpart = kC*cer_pe[int(Etru/0.01)]/Etru/scale

        quenchpart_low = read_quench(kA-kAerr, kB+kBerr, Etru)
        cerpart_low = (kC-kCerr)*cer_pe[int(Etru/0.01)]/Etru/scale

        quenchpart_high = read_quench(kA+kAerr, kB-kBerr, Etru)
        cerpart_high = (kC+kCerr)*cer_pe[int(Etru/0.01)]/Etru/scale

        totNL = quenchpart + cerpart
        NLlist2.append(totNL)
        NLlist2_low.append(quenchpart_low+cerpart_low)
        NLlist2_high.append(quenchpart_high+cerpart_high)


    plt.plot(Elist1, NLlist1, "-", label="2-layer sampling")
    plt.fill_between(Elist1, NLlist1_high, NLlist1_low, color='blue',alpha=0.25)
    plt.plot(Elist2, NLlist2, "-", label="primElecDist")
    plt.fill_between(Elist2, NLlist2_high, NLlist2_low, color='orange',alpha=0.25)
    plt.legend()
    plt.xlabel("Etrue/MeV")
    plt.ylabel("Evis/Etrue")
    plt.show()
