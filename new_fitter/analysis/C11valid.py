import numpy as np
import matplotlib.pyplot as plt
import elecLoader as el
import ROOT
import uproot as up


def wenNcer(E, p0, p1, p2, p3, p4):
    if E<=0.2:
        return 0
    if E>0.2:
        E = E - 0.2
        return p3*E**2/(p4+p0*E+p1*np.exp(-p2*E))


def sigma2Func(a, b, c, A, x):
    s2 = a**2/A*x + b**2*x**2 + c**2/A**2
    return s2 if s2>0  else 0

def sigma2FuncNew(a, b, n, x):
    s2 = a**2*x + b**2*np.power(x, n)
    return s2 if s2>0  else 0


# C11 BETA+ DECAY
ff    = up.open("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/C11/C11.root")
totpe = ff["c11"]["totpe"].array()
ke  = ff["c11"]["ke"].array()

Y = 3134.078 / 2.223
a1, b1, c1 = 0.988, 7.89e-3, 0
a1err, b1err = 6.46e-3, 4.33e-4
p0, p1, p2, p3, p4 = 4.20, 3.64, 0.11, 405, -2.09


a2, b2, n2 = 0.939, 0.103, 1.439
a3, b3, n3 = 0.840, 0.286, 1.234
#a2, b2, n2 = 0.855, 0.262, 1.253


### calculation
Ecalc = []
for i in ke:
    Nsct = el.getFitNsct(i, 6.78e-3, 1417.77, "Sim")
    Ncer = wenNcer(i, p0, p1, p2, p3, p4)
    Ntot = Nsct + Ncer
    #k1 = ROOT.gRandom.Gaus(Ntot/Y, np.sqrt(sigma2Func(a1, b1, c1, Y, Ntot/Y)))
    k1 = ROOT.gRandom.Gaus(Ntot, np.sqrt(sigma2FuncNew(a2, b2, n2, Ntot))) / Y
    k2 = ROOT.gRandom.Gaus(Ntot, np.sqrt(sigma2FuncNew(a3, b3, n3, Ntot))) / Y

    E1 = ROOT.gRandom.Gaus(660.8, 27.07) / Y
    E2 = ROOT.gRandom.Gaus(660.8, 27.07) / Y

    #Ecalc.append(k)
    Ecalc.append(k2+E1+E2)

fig, ax = plt.subplots()
ax.hist(Ecalc,   bins=100, range=(0.5, 3), histtype="step", lw=2, color="black",   label="Model")
ax.hist(totpe/Y, bins=100, range=(0.5, 3), histtype="step", lw=2, color="crimson", label="Simulation truth")
ax.legend(loc="upper center", prop={'size':14}, ncol=2)
ax.set_xlabel(r"$E^{vis}$ [MeV]", fontsize=14)
ax.set_ylabel("A.U.", fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=14)

plt.tight_layout()
plt.show()






    
