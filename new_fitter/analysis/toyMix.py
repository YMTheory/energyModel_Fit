import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import norm
#from scipy.optimizer import curve_fit

def gauss(x, A, m, s):
    return A * np.exp(-(x-m)**2/2/s**2)

X = np.arange(50, 150, 1)
p1 = norm.pdf(X, 95, 8)
p2 = norm.pdf(X, 100, 20)
p3 = norm.pdf(X, 105, 9)

p4 = []
d1 = np.random.normal(95, 8,  size=10000)
d2 = np.random.normal(100, 20, size=10000)
d3 = np.random.normal(105, 9, size=10000)

for i, j, k in zip(d1, d2, d3):
    p4.append(i)
    p4.append(j)
    p4.append(k)

p4 = np.array(p4)
#popt, pcov = curve_fit(gauss, p4,)

p5 = norm.pdf(X, 100, np.sqrt((50+8**2+9**2+20**2)/3.) )


plt.plot(X, p1, "--", label="gauss 1")
plt.plot(X, p2, "--", label="gauss 2")
plt.plot(X, p3, "--", label="gauss 3")
plt.hist(p4, range=(50, 150), bins=100, density=True, alpha=0.3, label="sampling")

plt.text(55, 0.04, r"sampling dist.$\mu=%.2f, \sigma=%.3f$ " %(np.mean(p4), np.std(p4)))
plt.text(55, 0.03, r"calculated dist.$\mu=%.2f, \sigma=%.3f$ " %(100, np.sqrt((50+64+81+400)/3)))

print(np.mean(p4), np.std(p4))

plt.legend()

plt.show()
