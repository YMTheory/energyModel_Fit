import numpy as np
import matplotlib.pyplot as plt
import elecLoader as el

es = 3134.078 / 2.223

Edep = np.arange(0.1, 8, 0.1)
Evis = []
resFit, resConstFit = [], []
resCerFit, resSctFit, resCovFit = [], [], []
resCerConstFit, resCovConstFit = [], []
for i in Edep:
    NPE = el.getFitNsct(i) + el.getFitNcer(i)
    Evis.append(NPE/es)
    
    resFit.append(el.getFitResol(Evis[-1]))
    resConstFit.append(el.getFitConstResol(Evis[-1]))

    resCerFit.append(np.sqrt(el.getFitCerSigma(Evis[-1]))/NPE)
    resSctFit.append(np.sqrt(el.getFitSctSigma(Evis[-1]))/NPE)
    resCovFit.append(np.sqrt(el.getFitCov(Evis[-1]))/NPE)
    
    resCerConstFit.append(np.sqrt(el.getFitCerSigmaConst(Evis[-1]))/NPE)
    resCovConstFit.append(np.sqrt(el.getFitCovConst(Evis[-1]))/NPE)
    
    #print(el.getFitResol(Evis[-1]), el.getFitCerSigma(Evis[-1]), el.getFitSctSigma(Evis[-1]), el.getFitCov(Evis[-1]))


resFit         = np.array(resFit)
resConstFit    = np.array(resConstFit)
resCerFit      = np.array(resCerFit)
resSctFit      = np.array(resSctFit)
resCovFit      = np.array(resCovFit)
resCerConstFit = np.array(resCerConstFit)
resCovConstFit = np.array(resCovConstFit)



# Load Truth Points
NPE_truth, SPE_truth = [], []
for i in Edep:
    NPE_truth.append(el.getNPE(i))
    SPE_truth.append(el.getSPE(i))

NPE_truth = np.array(NPE_truth)
SPE_truth = np.array(SPE_truth)

fig, ax = plt.subplots()
#ax.plot(Evis, resFit**2, "-", lw=2, color="crimson", label=r"$\sigma_E^2$")
#ax.plot(Evis, resConstFit**2, "--", lw=2, color="crimson", label=r"$\sigma_{const}^2$")

ax.plot(Evis, resSctFit**2, "-",  lw=2, color="dimgray", label=r"$\sigma_{sct}^2$")
ax.plot(Evis, resSctFit**2+resCovFit**2+resCerFit**2, lw=2, color="blue", label=r"$\sigma_E^2$")
ax.plot(Evis, resCerFit**2, "--", lw=2, color="dimgray", label=r"$\sigma_C^2$")
ax.plot(Evis, resCovFit**2, "-.", lw=2, color="dimgray", label=r"$cov$")
#ax.plot(Evis, resCerFit**2+resCovFit**2, ".", lw=2, color="dimgray", label=r"$\sigma_C^2$")
#ax.plot(Evis, resCerConstFit**2+resCovConstFit**2, "-.", lw=2, color="dimgray", label=r"$\sigma_C^2$")

#ax.plot(NPE_truth/es, SPE_truth**2/NPE_truth**2, "o", ms=4, color="orange")

ax.semilogy()

plt.grid(True)

plt.show()







