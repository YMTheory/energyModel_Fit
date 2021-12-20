import numpy as np
import matplotlib.pyplot as plt
import ROOT
import elecLoader as el
import uproot as up

def wenNcer(E, p0, p1, p2, p3, p4):
    if E<=0.2:
        return 0
    if E>0.2:
        E = E - 0.2
        return p3*E**2/(p4+p0*E+p1*np.exp(-p2*E))




def fitNsct(E, As, kB):
    if E > 15:
        return el.getFitNsct(15.9, kB, As, "Sim") / 15.9 * i
    else:
        return el.getFitNsct(E, kB, As, "Sim")



### Parameters ###
As, kB = 1408.07, 0.0062
p0, p1, p2, p3, p4 = 4.20, 3.64, 0.11, 405, -2.09
A = 3134.078 / 2.223


Edep = np.arange(1, 55, 1)
Ntot = []
gNtot = ROOT.TGraph()
nn = 0
for i in Edep:
    Ntot.append(wenNcer(i, p0, p1, p2, p3, p4) + fitNsct(i, As, kB))
    gNtot.SetPoint(nn, i, Ntot[-1])
    nn+=1
Ntot = np.array(Ntot)



def sigma2Func(a, b, c, A, x):
    s2 = a**2/A*x + b**2*x**2 + c**2/A**2
    return s2 if s2>0  else 0

# Resolution
a1, b1, c1 = 0.988, 7.89e-3, 0
a1err, b1err = 6.46e-3, 4.33e-4



# Michel endpoint

Edep_max = 52.8 #MeV
Evis_max = gNtot.Eval(Edep_max)/A
print("Visible energy %.1f MeV" %Evis_max)


# Read kinetic energy spectum
edep  = up.open("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/michel/michel_totpe_LS_v2.root")["michel"]["edep"].array()
totpe = up.open("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/michel/michel_totpe_LS_v2.root")["michel"]["totpe"].array()
edep = np.array(edep)
edge = edep[edep >= 50]
#plt.hist(edge, bins=100)
#plt.show()

## Apply nonlinearity and resolution ...
applied = []
for i in edge:
    j = gNtot.Eval(i) / A
    k = ROOT.gRandom.Gaus(j, np.sqrt(sigma2Func(a1, b1, c1, A, j)) )
    applied.append(k)



plt.hist(applied, bins=30, range=(54, 57), histtype="step", label="Model")
plt.hist(totpe/A, bins=30, range=(54, 57), histtype="step", label="Simulation truth")

plt.legend()
plt.show()





















