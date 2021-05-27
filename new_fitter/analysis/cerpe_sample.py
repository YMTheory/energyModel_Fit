import matplotlib.pyplot as plt
import numpy as np


def func(E, A2, A3, A4):
    E0 = 0.2
    A1 = 0
    x = np.log(1+E/E0)
    return (A1*x+A2*x**2+A3*x**3) * (1/E+A4) 

A2 = 9.99975e+00
A2err = 1.78279e+01
A3 = 1.98147e+01
A3err = 4.72424e-01
A4 = 1.00000e-02
A4err = 3.33673e-02

def sample():
    m_A2 = np.random.normal(A2, A2err, size=1000)
    m_A3 = np.random.normal(A3, A3err, size=1000)
    m_A4 = np.random.normal(A4, A4err, size=1000)

    E0 = np.arange(0.1, 16, 0.1)
    for i, j, k in zip(m_A2, m_A3, m_A4):
        cerpe = func(E0, i, j, k)

        plt.plot(E0, cerpe, "-", color="royalblue", alpha=0.1)
    
    plt.show()



if __name__ == "__main__":
    sample()
