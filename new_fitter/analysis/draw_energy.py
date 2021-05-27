
import numpy as np
import matplotlib.pyplot as plt

name = ["Cs137", "Mn54", "K40", "nH", "AmBe", "nC12", "AmC"]
Etrue = [0.662, 0.835,  1.461, 2.223, 4.43, 4.94, 6.13]
name1 = ["Ge68", "Co60"]
Etrue1 = [1.022, 2.506]
Etrue_nonl = [0.662, 0.835, 0.511, 1.461, 2.223, 1.253, 4.43, 4.94, 6.13]

plt.plot([0, 1, 3, 4, 6, 7, 8], Etrue, "o-")
plt.plot([2, 5], Etrue1, "s")
plt.xticks([0, 1,3,4,6,7,8], name)
plt.xticks([2, 5], name1)
plt.ylabel(r"$E_{true}/MeV$")
plt.show()
