import numpy as np
import matplotlib.pyplot as plt
import elecLoader as eloader

def func(p0, p1, p2, x):
    s2 = p0 + p1*x + p2 * x**2
    if s2 < 0:
        return 0
    else:
        #return np.sqrt(s2)
        return s2


p0 = -2.17203
p1 = 1.31498e3
p2 = 1.60508e2

a0 = -1.10053e1
a1 = 1484.17
a2 = 121.885
#p0 = -158.41
#p1 = 4.333
#p2 = 0.00181082

p0 = 22.333       
p1 = 0.982622     
p2 = 6.12389e-05  


def loadCerResFile(filename):
    E, cerpe, sigma = [], [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            E.append(float(data[0]))
            cerpe.append(float(data[1]))
            sigma.append(float(data[2]))

    E = np.array(E)
    cerpe = np.array(cerpe)
    sigma = np.array(sigma)

    return E, cerpe, sigma
            


def main():
    #E, cerpe, sigma = loadCerResFile("../data/electron/elecCerPEResol.txt")

    #petrue = np.arange(0, 1600, 10)
    etrue = np.arange(0, 26000, 10)

    resolE, totpe, resolData, resolerr = eloader.getResolArray()
    resolE1, resolData1, resolerr1 = [], [], []
    for i in [0, 90, 120, 150, 180, 210, 240, 270, 300, 350, 400, 450, 500, 550,  600, 620, \
              650, 700, 720, 750, 770, 800, 820, 840, 850, 860, 870, 880, 890, 900, 905, 910]:
        resolE1.append(totpe[i])
        resolData1.append(resolData[i])
        resolerr1.append(resolerr[i])
    resolFit  = []
    resolModel = []
    #resolFit1 =  []
    for i in etrue:
        resolFit.append( np.sqrt(func(p0, p1, p2, i)) )
        resolModel.append( np.sqrt(func(a0, a1, a2, i)) )
    #    resolFit1.append( np.sqrt(func(-123.213, 3.39232, 2.53173e-3, i)) )

    #plt.plot(cerpe, sigma, "o", ms=1.5, label="Simulation", color="darkviolet", alpha=0.3)
    plt.errorbar(resolE1, resolData1, yerr=resolerr1, fmt="o", ms=5, color="violet",  label="Simulation")
    plt.plot(etrue, resolFit, "-", label="Parameterization", color="royalblue")
    #plt.plot(etrue, resolModel, "--", label="Model Fitting", color="orange")
    #plt.plot(petrue, resolFit1, "-", label="Parameterization", color="hotpink")

    #plt.xlabel(r"$E_{dep}$/MeV")
    plt.xlabel(r"$N_{tot}$")
    plt.ylabel(r"$\sigma_{N_{tot}}$")

    plt.grid(True)

    plt.legend()
    plt.savefig("sigmavsNPE_electron.pdf")

    plt.show()


if __name__ == "__main__":
    main()




