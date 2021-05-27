import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import elecLoader as el

def func(E, A2, A3, A4):
    E0 = 0.2
    A1 = 0
    x = np.log(1+E/E0)

    cerpe = (A1*x+A2*x**2+A3*x**3) * (1/E+A4) * E
    print("%.2f, %.5f %.5f, %.5f, %.2f" %(E, A2, A3, A4, cerpe))
    return cerpe


def load(filename):
    Etrue, cerpe = [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Etrue.append(float(data[0]))
            cerpe.append(float(data[1]))

    return Etrue, cerpe


def select(old1, old2):
    new1, new2 = [], []
    for i in range(900):
        if i%30 == 0 and i!=0:
            new1.append(old1[i])
            new2.append(old2[i])

        if i == 880 or i == 890 or i==900 or i==910 or i==872 or i==874 or i==876:
            new1.append(old1[i])
            new2.append(old2[i])


    new1 = np.array(new1)
    new2 = np.array(new2)
    return new1, new2


def main():
    
    Etrue, cerpe = load("../data/electron/cerPE1.txt")
    Etrue = np.array(Etrue)
    cerpe = np.array(cerpe)
    Etrue1, cerpe1 = [], []
    for i in range(0, 915, 30):
        Etrue1.append(Etrue[i])
        cerpe1.append(cerpe[i]/55.226/Etrue[i])
    #Etrue1, cerpe1 = select(Etrue, cerpe) 

    #popt, pcov = curve_fit(func, Etrue1, cerpe1/Etrue1)
    #print(popt)

    #x1 = np.arange(0.1, 15, 0.1)
    #y1 = func(x1, *popt)
    #y2 = []
    #for i in x1:
    #    y2.append( func(i, -7.34716, 15.5519, 0.0267155) )

    plt.plot(Etrue1, cerpe1, "o-", color='blue')
    #plt.plot(x1, y1, "-")
    #plt.plot(x1, y2, "-")

    plt.xlabel(r"$E_{dep}$/MeV")
    plt.ylabel(r"$N_{cer}$/MeV")
    plt.grid(True)
    plt.savefig("CerenkovPE_MeV.pdf")

    plt.show()


if __name__ == "__main__":
    main()
            

