import numpy as np
import sys

def sigma2Func(a, b, c, A, x):
    s2 = a**2/A*x + b**2*x**2 + c**2/A**2
    return s2 if s2>0  else 0

# Resolution
a1, b1, c1 = 0.988, 7.89e-3, 0
a1err, b1err = 6.46e-3, 4.33e-4
A = 3134.078/2.223



if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Wrong input number !!!")

    else:
        E = float(sys.argv[1])
        sigma2 = sigma2Func(a1, b1, c1, A, E)
        resol = np.sqrt(sigma2) / E

        print("Energy resolution for %.1f MeV electron is %.3f" %(E, resol))



