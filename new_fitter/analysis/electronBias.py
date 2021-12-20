import numpy as np
import matplotlib.pyplot as plt
import elecLoader as el
from matplotlib import gridspec

A, kB, kC, p1, p2 = 1408.45, 6.29e-3, 0.991, 0.979, 6.13e-5

if __name__ == "__main__" :

    globalES = 3134.078/2.223
    Edep = np.arange(0.50, 8.10, 0.10)
    Edep1 = np.arange(0.50, 8.10, 0.50)

    # load simulation
    npeSim = []
    speSim = []
    for i in Edep:
        npeSim.append(el.getNPE(i))
        speSim.append(el.getSPE(i))
    npeSim = np.array(npeSim)
    speSim = np.array(speSim)
    nonlSim = npeSim / Edep / globalES 
    resSim = speSim/npeSim

    # fitting results
    npeCalc = []
    npeCalc1 = []
    speCalc = []
    for i in Edep:
        sctpe = el.getQPE(i, kB, A)
        cerpe = kC * el.getCerNPE(i)
        npeCalc.append((sctpe+cerpe))
        speCalc.append(np.sqrt(p1*(sctpe+cerpe) + p2*(sctpe+cerpe)**2))
    for i in Edep1:
        sctpe = el.getQPE(i, kB, A)
        cerpe = kC * el.getCerNPE(i)
        npeCalc1.append((sctpe+cerpe))
    npeCalc = np.array(npeCalc)
    npeCalc1 = np.array(npeCalc1)
    speCalc = np.array(speCalc)
    nonlCalc = npeCalc / Edep / globalES
    nonlCalc1 = npeCalc1 / Edep1 / globalES
    resCalc = speCalc/npeCalc


    nonlBias = (nonlCalc - nonlSim) / nonlSim
    resBias = (resCalc - resSim) / resSim

    fig = plt.figure(figsize=(9, 4))
    spec = gridspec.GridSpec(ncols=2, nrows=1)

    ax0 = fig.add_subplot(spec[0])
    ax1 = fig.add_subplot(spec[1])


    ax0.plot(Edep, nonlSim, "o", ms=3, color='blue')
    ax0.plot(Edep1, nonlCalc1, color='orange')
    ax0.set_xlabel(r"$E_{dep}$/MeV")
    ax0.set_ylabel(r"$E_{vis}/E_{dep}$")


    ax1.plot(Edep, resSim, "o", ms=3, color='blue')
    ax1.plot(Edep, resCalc, color='orange')
    ax1.set_xlabel(r"$E_{vis}$/MeV")
    ax1.set_ylabel(r"$\sigma/N_{tot}$)

    plt.show()






