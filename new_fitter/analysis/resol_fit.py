import numpy as np
import matplotlib.pyplot as plt
import elecLoader as eloader

def func(p0, p1, p2, x):
    s2 = p0 + p1*x + p2 * x**2
    if s2 < 0:
        return 0
    else:
        return np.sqrt(s2)

p0 = -2.17203
p1 = 1.31498e3
p2 = 1.60508e2

def main():
    resolE, resolData = eloader.loadResFile("../data/electron/elecResol1.txt")

    Etrue = np.arange(0.1, 8, 0.1)

    resolFit =  []
    for i in Etrue:
        resolFit.append( func(p0, p1, p2, i) )

    plt.plot(resolE, resolData, "o", ms=1.5, label="Simulation", color="darkviolet", alpha=0.3)
    plt.plot(Etrue, resolFit, "-", label="Parameterization", color="seagreen")

    plt.xlabel(r"$E_{true}/MeV$")
    plt.ylabel(r"$\sigma_{PE}$")

    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()




