import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import elecLoader as el

es = 3134.078 / 2.223
es1 = (3134.078 - 111) / 2.223
es2 = 111 / 2.223

Edep = np.arange(0.1, 8, 0.1)
Evis = np.zeros((5, len(Edep)))
resFit = np.zeros((5, len(Edep)))

num = 0
for kC in [0.0, 0.2, 0.5, 0.8, 1.0]:  # cerPE be selected out
    nn = 0
    for i in Edep:
        NPE = el.getFitNsct(i) + el.getFitNcer(i) * kC + 660.8*2
        Evis[num, nn] = NPE / (es-kC*es2)

        cerSigma2 = el.getFitCerSigma(i)
        sctSigma2 = el.getFitSctSigma(i)
        cov = el.getFitCov(i)

        if cerSigma2 == 0 or cov == 0:
            resFit[num, nn] = np.sqrt(sctSigma2+2*27.07**2) / NPE
            print(num, nn, sctSigma2, NPE)
            nn += 1
            continue

        corr = cov / np.sqrt(cerSigma2 * sctSigma2)
        
        cerSigma2 *= (1-kC)**2

        tot = sctSigma2 + cerSigma2 + 2 * corr * np.sqrt(sctSigma2*cerSigma2) + 2*27.07**2

        resFit[num, nn] = np.sqrt(tot)/es/Evis[num, nn]
        nn += 1

    num += 1


fig = plt.figure(constrained_layout=True, figsize=(6, 6))
gs = fig.add_gridspec(2, 1, height_ratios=[1, 2])

ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[1, 0])

kc = [0, 20, 50, 80, 100]
cc = ['dimgray', 'royalblue', 'orange', 'green', 'crimson']
for i in range(4):
    ax0.plot(Evis[0, :], (resFit[i+1, :]-resFit[0,:])/resFit[0, :], "-", lw=2, color=cc[i+1])
for i in range(5):
    ax1.plot(Evis[0, :], resFit[i, :], lw=2, color=cc[i], label="%d%%"%kc[i])


ax0.grid(True)
ax0.set_xlabel(r"Positron $E_{vis}$[MeV]", fontsize=13)
ax0.set_ylabel(r"$\frac{R-R^0}{R^0}$", fontsize=15)
ax0.tick_params(axis='both', which='major', labelsize=11)

ax1.grid(True)
ax1.set_xlabel(r"Positron $E_{vis}$[MeV]", fontsize=13)
ax1.set_ylabel(r"$\sigma_E/E_{vis}$", fontsize=13)
ax1.legend(loc="upper right", prop={'size':13}, title=r"$N_{Cer}$ selection efficiency")
ax1.tick_params(axis='both', which='major', labelsize=11)


plt.tight_layout()

plt.savefig("CerOnResol.pdf")
plt.show()


