import numpy as np
import elecLoader as el
import matplotlib.pyplot as plt

kB1, es1, kC1 = 6.29e-3, 1408.45, 0.991
kB2, es2, kC2 = 6.03479e-3, 1404.65, 1.02747

Etrue = np.arange(0.1, 15, 0.1)
npe_sim = []
cerpe_sim, sctpe_sim = [], []
npe1, npe2 = [], []
for i in Etrue:
    npe1.append(el.getQPE(i, kB1, es1) + kC1*el.getCerNPE(i))
    npe2.append(el.getQPE(i, kB2, es2) + kC2*el.getCerNPE(i))
    npe_sim.append(el.getNPE(i))
    cerpe_sim.append(el.getCerNPE(i))
    sctpe_sim.append(el.getSctNPE(i))

sctpe_sim = np.array(sctpe_sim)
cerpe_sim = np.array(cerpe_sim)
npe1 = np.array(npe1)
npe2 = np.array(npe2)
npe_sim = np.array(npe_sim)
scale = 1400

plt.plot(Etrue, npe1/Etrue/scale, label="previous")
#plt.plot(Etrue, npe2/Etrue/scale, "--", label="w/15MeV gamma")
plt.plot(Etrue, npe_sim/Etrue/scale, label="sim")
plt.semilogx()
#plt.semilogy()
plt.legend()

plt.show()
