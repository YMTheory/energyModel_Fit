import numpy as np
import matplotlib.pyplot as plt

def func(p0, p1, p2, x):
    s2 = p0 + p1*x + p2 * x**2
    if s2 < 0:
        return 0
    else:
        #return np.sqrt(s2)
        return s2

c0 = -158.415
c1 = 4.3332
c2 = 0.0018108

s0 = 0
s1 = 1
s2 = 0

global_es = 3134.078/2.223
ratio = np.arange(0.8, 0.99, 0.001)

resol1, resol2 = [], []
for p in ratio:
    sctsigma2 = func(s0, s1, s2, p*global_es)
    cersigma2 = func(c0, c1, c2, global_es*(1-p))
    sigma1 = np.sqrt( (2-p)/p*sctsigma2 + cersigma2 )
    resol1.append(sigma1/global_es)
    sigma2 = np.sqrt( (sctsigma2+cersigma2)/1-2*p*(1-p) )
    resol2.append(sigma2)

plt.plot(ratio, resol1, "-", color="royalblue")
plt.plot(ratio, resol2, "-", color="orange")
plt.xlabel(r"$N_{sct}/N_{tot}$")
plt.ylabel(r"$\sigma_{N_{tot}}/N_{tot}$")
plt.grid(True)

plt.show()