from singleGamma import singleGamma
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import prmBetaLoader as bloader
import elecLoader as el


def loadBestFit(filename, num):
    bestFit, err = [], []
    row = 0
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if row == num:
                for kk in range(num):
                    bestFit.append(float(data[kk]))
            if row == num+1:
                for kk in range(num):
                    err.append(float(data[kk]))
            row += 1

    return bestFit, err



def wenNcer(E, p0, p1, p2, p3, p4):
    if E<=0.2:
        return 0
    if E>0.2:
        E = E - 0.2
        return p3*E**2/(p4+p0*E+p1*np.exp(-p2*E))



def fitNsct(E, As, kB):
    if E > 15.9:
        return el.getFitNsct(15.9, kB, As, "Sim") / 15.9 * E
    else:
        return el.getFitNsct(E, kB, As, "Sim")


Y = 3134.078 / 2.223

filename = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamNewB12_kSimQ_kNewAnaCer_kNew/NewgamNewB12_kSimQ_kNewAnaCer_kNew_rescov.txt"
par1, parerr1 = loadBestFit(filename, 3)
a, b, n = par1
def elecResol(E, c, d):
    return np.sqrt(c**2/Y/E + d**2)
def sigma2FuncNew(a, b, n, x):
    s2 = a**2*x + b**2*np.power(x, n)
    return s2 if s2>0  else 0


# Electron
Evis_elec = np.arange(0.3, 6.5, 0.1)
sigma_elec = []
for i in Evis_elec:
    sigma_elec.append(np.sqrt(sigma2FuncNew(a, b, n, i*Y))/i/Y)
    #sigma_elec.append(elecResol(i, a, b))
    #sigma_elec.append( sGamArr[0].getFitSPE(i) )
sigma_elec = np.array(sigma_elec)


# Gamma Decomposition :

#filename = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/gamB12NewMic1whole/gamB12NewMic1whole_kSimQ_kCerNewAna_kNPE_nonlcov.txt"
filename = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12_kSimQ_kNewAnaCer_kNew/NewgamB12_kSimQ_kNewCerAna_kNew_nonlcov.txt"
par1, parerr1 = loadBestFit(filename, 7)
scale, kB, p0, p1, p2, p3, p4 = par1


es = 3134.078/2.223
npe, spe, sigma_part1, sigma_part2 = [], [], [], []


for i in range(305, 370, 5):
    mom = (i-300)/10
    filename = "../data/gamma/Livermore-"+str(i)+".txt"
    print(filename)
    prmBeta, prmAntiBeta = bloader.loadPrmBeta(filename)

    # calcSingleEvent :
    nSamples = 5000
    mu_arr, sigma_arr = [], []
    for i in range(nSamples):
        tmpnpe, tmpspe = 0, 0
        # Electrons :
        for j in prmBeta[i]:
            if j == 0:
                break
            onenpe  = (fitNsct(j, scale, kB) + wenNcer(j, p0, p1, p2, p3, p4))
            tmpnpe += onenpe
            #tmpspe += (elecResol(onenpe/Y, a, b) * onenpe )**2
            tmpspe += sigma2FuncNew(a, b, n, onenpe)
        # Positrons:
        for j in prmAntiBeta[i]:
            if j == 0:
                break
            onenpe  = (fitNsct(j, scale, kB) + wenNcer(j, p0, p1, p2, p3, p4))
            tmpnpe += (onenpe + 2*660.8)
            #tmpspe += ((elecResol(onenpe/Y, a, b) * onenpe )**2+ 27.07**2*2)
            tmpspe += sigma2FuncNew(a, b, n, onenpe) 

        mu_arr.append(tmpnpe)
        sigma_arr.append(np.sqrt(tmpspe))

    mu = np.array(mu_arr)
    sigma = np.array(sigma_arr)
    

    tmp_npe, tmp_spe = 0, 0
    tmp_sigma_part1, tmp_sigma_part2 = 0, 0
    for i in mu:
        tmp_npe += i
    tmp_npe = tmp_npe/nSamples
    for i, j in zip(mu, sigma):
        tmp_spe += (i-tmp_npe)**2 + j**2
        tmp_sigma_part1 += (i-tmp_npe)**2
        tmp_sigma_part2 += j**2

    tmp_spe = np.sqrt(tmp_spe/nSamples)
    tmp_sigma_part1 /= nSamples
    tmp_sigma_part2 /= nSamples

    npe.append(tmp_npe)
    spe.append(tmp_spe)
    sigma_part1.append(tmp_sigma_part1)
    sigma_part2.append(tmp_sigma_part2)


npe = np.array(npe)
spe = np.array(spe)
sigma_part1 = np.array(sigma_part1)
sigma_part2 = np.array(sigma_part2)
sigma = np.sqrt(sigma_part1 + sigma_part2)


#fig, ax = plt.subplots()
fig = plt.figure(figsize=(6, 6))
spec = gridspec.GridSpec(ncols=1, nrows=2, 
                     height_ratios=[1, 2])

ax0 = fig.add_subplot(spec[0])
ax = fig.add_subplot(spec[1])
ax.plot(npe/es, (spe/npe)**2, "o--", lw=2, ms=6, color="blue", label=r"$\gamma$ $\sigma^2$")
ax.plot(npe/es, np.sqrt(sigma_part2)**2/npe**2, "o--", lw=2, ms=6, zorder=2, color="dimgray", label=r"$\gamma$ $\sigma^2_{ave}$")
ax.fill_between(npe/es, (spe/npe)**2, sigma_part2/npe**2, color="slategray", alpha=0.5)

ax.plot(Evis_elec,sigma_elec**2, "-", lw=2, color="red", zorder=1, label=r"$e^- \sigma^2$")

ax.set_ylabel(r"$(\sigma/E^{vis})^{2}$", fontsize=17)
ax.tick_params(axis='x', which='major', labelsize=16)
ax.tick_params(axis='y', which='major', labelsize=16)
ax.set_xlabel(r"$E^{vis}$ [MeV]", fontsize=17)
ax.grid(True)
ax.ticklabel_format(style="sci", scilimits=(-3, 0), axis="y")
ax.set_xlim(0, 7)
ax.legend(prop={"size":16}, ncol=1)


#rect = [0.35,0.57,0.7,0.2]
#ax0 = add_subplot_axes(ax, rect)
ax0.plot(npe/es, (sigma_part1)/spe**2, "o-", ms=5, lw=2, color="black")
ax0.set_ylabel(r"$\sigma_{nonl}^2/\sigma^2$", fontsize=15)
ax0.tick_params(axis='x', which='major', labelsize=14)
ax0.tick_params(axis='y', which='major', labelsize=14)
#ax0.set_xlabel(r"$E^{vis}$ [MeV]", fontsize=15)
ax0.grid(True)
ax0.set_xlim(0, 7)
#ax0.set_ylim(0.1, 0.4)
#ax0.ticklabel_format(style="sci", scilimits=(-3, 0), axis="y")

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.1)


plt.tight_layout()
plt.savefig("gammaDecomp.pdf")
plt.show()




