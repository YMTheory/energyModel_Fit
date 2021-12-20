import numpy as np
import matplotlib.pyplot as plt
import elecLoader as el

def func(E, p0, p1, p2, p3, p4):
    if E <= 0.2 :
        return 0
    else:
        E = E - 0.2
        N = p3*E**2/(p4+p0*E+p1*np.exp(-p2*E))
        return N



if __name__ == "__main__":
    
    x = np.arange(0.1, 40, 0.1)
    
    par = [4, 0.1, 3, 410, 1]

    Ac = 177.508 / 2.223

    Edep = np.arange(0.1, 15, 0.1)
    Edep1 = np.arange(0.1, 55, 0.1)
    truth = []
    shape1, shape2, shape3, shape4 = [], [], [], []
    for i in Edep:
        truth.append(el.getCerNPE(i))
    for i in Edep1:
        shape1.append(func(i, *par))
    truth  = np.array(truth)
    shape1 = np.array(shape1)
    Ac = 177.508/2.223

    Edep2 = [20, 30, 40, 50, 60]
    truth2 = [1985, 2995, 3968, 4955, 5952]
    Edep2 = np.array(Edep2)
    truth2 = np.array(truth2)

    #plt.plot(Edep, truth , "o", ms=3)
    #plt.plot(Edep1, shape1, label="All")
    #plt.plot(Edep2, truth2, "s")
    plt.plot(Edep,  truth/Ac/Edep , "o", ms=3)
    plt.plot(Edep1, shape1/Ac/Edep1, label="All")
    plt.plot(Edep2, truth2/Ac/Edep2, "s")
    #plt.plot(Edep1, shape2/Ac/Edep1, label="No p1")
    #plt.plot(Edep1, shape3/Ac/Edep1, label="No p0")
    #plt.plot(Edep1, shape4/Ac/Edep1, label="No p4")
    plt.legend()

    plt.show()

