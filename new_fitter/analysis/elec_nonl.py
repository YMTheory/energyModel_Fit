import numpy as np
import matplotlib.pyplot as plt
import elecLoader as el
import ROOT
from scipy.interpolate import make_interp_spline, BSpline
from matplotlib import gridspec


es = 3134.078/2.223
ceres = 177.508/2.223
Etrue = np.arange(0.1, 12.1, 0.5)

Nsct, Ncer, Ntot = [], [], []
for i in Etrue:
    Nsct.append(el.getSctNPE(i))
    Ncer.append(el.getCerNPE(i))
    Ntot.append(el.getNPE(i))

Nsct = np.array(Nsct)
Ncer = np.array(Ncer)
Ntot = np.array(Ntot)
#print(Ntot)


grin1 = ROOT.TGraph()
for i in range(len(Etrue)):
    if Etrue[i] > 1.5:
        grin1.SetPoint(i, Etrue[i], Ncer[i]/ceres/Etrue[i])

# smooth
gs = ROOT.TGraphSmooth("normal")
grout1 = gs.SmoothKern(grin1, "normal", 2.0)
x1, y1 = [], []
for i in range(grout1.GetN()):
    if grout1.GetPointX(i) > 1.5:
        x1.append(grout1.GetPointX(i))
        y1.append(grout1.GetPointY(i))

#print(x1)
#print(y1)


"""
grin2 = ROOT.TGraph()
for i in range(len(Etrue)):
    grin2.SetPoint(i, Etrue[i], Nsct[i]/Ntot[i])

# smooth
gs2 = ROOT.TGraphSmooth("normal")
grout2 = gs.SmoothKern(grin2, "normal", 2.0)
x2, y2 = [], []
for i in range(grout2.GetN()):
    x2.append(grout2.GetPointX(i))
    y2.append(grout2.GetPointY(i))

print(x2)
print(y2)
"""



fig = plt.figure(figsize=(12, 5))
spec = gridspec.GridSpec(ncols=2, nrows=1 )

ax0 = fig.add_subplot(spec[0])
ax1 = fig.add_subplot(spec[1])

bl = np.zeros(len(Ntot))
#ax0.plot(Etrue, Ntot/Etrue/es, "-", lw=3,  color="blue")
#ax0.plot(Etrue, Ncer/Ntot, "--", lw=3, color="seagreen")
ax0.plot(Etrue, Ntot, "-", lw=3, color="blue")
ax0.fill_between(Etrue, bl, Ncer, color="seagreen", alpha=0.3, label="Cherenkov")
ax0.fill_between(Etrue, Ncer, Ntot, color="orange", alpha=0.3, label="Scintillation")
#ax0.plot(Etrue, Nsct/Ntot, "--", lw=3, color="blue")
#ax.plot(x1, y1, "-", lw=3,  color="blue")
ax0.set_xlabel("Electron kinetic energy / MeV", fontsize=16)
#ax0.set_ylabel(r"$E_{vis}/E_{dep}$", fontsize=13, color="blue")
ax0.set_ylabel(r"$N_{tot}$", fontsize=16)
ax0.tick_params(axis='both', labelsize=14)
ax0.semilogy()
ax0.set_xlim(Etrue[0], Etrue[-1])
ax0.text(4, 1200, "Scintillation", fontsize=15)
ax0.text(6, 60, "Cherenkov", fontsize=15)

ax2 = ax0.twinx()
##ax2.plot(x2, y2, "--", lw=3, color="seagreen")
ax2.plot(Etrue, Nsct/Ntot, "--", lw=3, color="darkviolet")
ax2.set_ylabel(r"$N_{sct}/N_{tot}$", fontsize=16, color="darkviolet", labelpad=18, rotation=270)
ax2.tick_params(labelsize=12, labelcolor="darkviolet")


dEtrue, dNcer = [], []
for i, j in zip(Etrue, Ncer):
    if i >= 1.5:
        dEtrue.append(i)
        dNcer.append(j)

dEtrue = np.array(dEtrue)
dNcer = np.array(dNcer)


ax1.plot(Etrue, Ntot/Etrue/es, "-", lw=3,  color="blue", label="Total")
ax1.plot(x1, y1, "-.", lw=3,  color="seagreen", label="Cherenkov")
ax1.plot(Etrue, Nsct/Etrue/(es-ceres), "--", lw=3,  color="orange", label="Scintillation")
#ax1.plot(dEtrue, dNcer/dEtrue/ceres, "-", lw=3,  color="seagreen", label="Cherenkov")
#ax1.fill_between(Etrue, bl, Ncer/Etrue/es, color="orange", label="Cherenkov")
#ax1.fill_between(Etrue, Ncer/Etrue/es, Ntot/Etrue/es, color="seagreen", label="scintillation")
#ax1.semilogy()
#ax1.semilogx()
ax1.legend(prop={'size' : 13})
ax1.set_ylabel("NPE per MeV / NPE per MeV at nH point", fontsize=16)
ax1.set_xlabel("Electron kinetic energy / MeV", fontsize=16)
ax1.tick_params(axis='both', labelsize=12)
ax1.grid(True)

plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
            wspace=0.02, hspace=0.02)
plt.tight_layout()
plt.savefig("elecNonl+sctRatio.pdf")
plt.show()


