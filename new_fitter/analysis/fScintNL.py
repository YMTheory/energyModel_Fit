# -*- coding: utf-8 -*-
"""
    scintillator nonlinearity
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2020 Xu Hangkun <xuhangkun@163.com>
    :license: MIT, see LICENSE for more details.
"""

from math import exp
import numpy as np
import pandas as pd
import ROOT
from scipy import interpolate
from scipy import integrate

def loadBinCenter():
    ff = ROOT.TFile("../data/electron/Quench5.root", "read")
    hist = ff.Get("kB65")
    binCenter = []
    for i in range(hist.GetNbinsX()):
        binCenter.append(hist.GetBinCenter(i))
    return binCenter



class ScintNL:

    def __init__(self,de_dx="/dybfs/users/xuhangkun/SimTAO/offline/NonLinearity/data/input/ESTAR_GdLS.txt"):
        self.de_d_dx = pd.read_csv(de_dx)
        self.f_de_d_dx = interpolate.interp1d(self.de_d_dx["Kinetic"].to_numpy(),self.de_d_dx["Total"].to_numpy())
        self.kb = 0.0158

    def __call__(self,es,kb=0.0065):
        """
        e should be a energy list in ascending order
        """
        self.kb = kb
        nes = np.concatenate((np.array([0],),es),axis=0)
        accum = 0.
        value = []
        for i in range(len(es)):
            accum += integrate.quad(self.dscintnl,nes[i],nes[i+1])[0]
            value.append(accum/nes[i+1])
        return np.array(value)

    def dscintnl(self,es):
        return 1./(1 + self.kb*self.f_de_d_dx(es))

class CherenkovNL:

    def __init__(self,cherenkov_file="/dybfs/users/xuhangkun/SimTAO/offline/NonLinearity/data/input/TAO_Cherenkov.csv"):
        self.cherenkov = pd.read_csv(cherenkov_file)
        self.f_cherenkov = interpolate.interp1d(self.cherenkov["energy"].to_numpy(),self.cherenkov["cov"].to_numpy())

    def __call__(self,e,kc):
        return self.f_cherenkov(e)*kc

class LSNonLin:

    def __init__(self,
            de_dx_file = "/dybfs/users/xuhangkun/SimTAO/offline/NonLinearity/data/input/ESTAR_GdLS.txt",
            cherenkov_file = "/dybfs/users/xuhangkun/SimTAO/offline/NonLinearity/data/input/TAO_Cherenkov.csv"
            ):
        self.cr_nl = CherenkovNL(cherenkov_file)
        self.scin_nl = ScintNL(de_dx_file)

    def __call__(self,es,p0,p1,p2):
        """
        p0 : A
        p1 : kb
        p2 : kc
        """
        return p0*(self.scin_nl(es, p1))
        #return p0*(self.scin_nl(es,p1)+self.cr_nl(es,p2))

def fnonlin(e,p0,p1,p2,p3):
    """
    nonlinearity model used in JUNO
    Args:
        e: energy of electron
        p0 ~ p3
    """
    #return np.asarray((p0 + p3*np.log(e) )/(1 + p1*np.exp(-1.0*p2*e)))
    return np.asarray((p0 + p3/e )/(1 + p1*np.exp(-1.0*p2*e)))

if __name__ == "__main__":
    # x = np.asarray([0.5,1,2,3,4,5,6])
    # print(fnonlin(x,10,10,-0.18,3.45))
    Etrue = loadBinCenter()
    scintnl = LSNonLin()
    print(scintnl([1, 2, 3],1,0.0065,0.019))
    for i in Etrue:
        print(scintnl([i], 1, 0.0065, 0.019))
















