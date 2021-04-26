import numpy as np
import matplotlib.pyplot as plt
import elecLoader as eloader
import random

p0 = -4.27498
p0err = 0.698514
p1 = 1360
p1err = 50.927
p2 = 154.691
p2err = 7.53067

# no correlation in fitter...

Etrue = np.arange(0.1, 8, 0.1)

def resFunc(x, a, b, c):
    s2 = a + b*x + c*x**2
    if s2<0:
        return 0
    else:
        return np.sqrt(s2)

def samplePar():
    m_p0 = random.gauss(p0, p0err)
    m_p1 = random.gauss(p1, p1err)
    m_p2 = random.gauss(p2, p2err)
    
    m_nonl = []
    for i in Etrue:
        m_nonl.append(resFunc(i, m_p0, m_p1, m_p2) )

    m_nonl = np.array(m_nonl)
    return m_nonl


def bestfit():
    m_p0 = p0
    m_p1 = p1
    m_p2 = p2

    m_nonl = []
    for i in Etrue:
        m_nonl.append(resFunc(i, m_p0, m_p1, m_p2) )

    m_nonl = np.array(m_nonl)
    return m_nonl

def nominal():
    m_p0 = -2.17203
    m_p1 = 1.31498e3
    #m_p2 = 3134.078/2.223
    m_p2 = 1.60508e2

    m_nonl = []
    for i in Etrue:
        m_nonl.append(resFunc(i, m_p0, m_p1, m_p2) )

    m_nonl = np.array(m_nonl)
    return m_nonl

def main():
    for i in range(5000):
        print(i)
        m_nonl1 = samplePar()

        plt.plot(Etrue, m_nonl1, "-", color="lightskyblue", alpha=0.05)

    plt.plot(Etrue, nominal(), "-", color="coral", label="nominal")
    plt.plot(Etrue, bestfit(), "--", color="blueviolet", label="best fit")
    plt.grid(True)
    plt.legend()
    plt.xlabel(r"$E_{true}$/MeV")
    plt.ylabel("p.e. resolution")
    plt.savefig("resCurve_fit.pdf")
    plt.show()


if __name__ == "__main__":
    main()
