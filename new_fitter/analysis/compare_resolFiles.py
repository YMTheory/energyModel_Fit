import numpy as np
import matplotlib.pyplot as plt

def read(ff):
    npem, npes = [], []
    with open(ff) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            npem.append(float(data[1]))
            npes.append(float(data[3]))

    npem = np.array(npem)
    npes = np.array(npes)



    return npem, npes




def main():
    npem1, npes1 = read("../data/electron/elecResol1.txt")
    npem2, npes2 = read("../data/electron/elecResol2.txt")
    npem3, npes3 = read("../data/electron/elecResol3.txt")
    
    plt.plot(npem1, npes1, "o", ms=2)
    plt.plot(npem2, npes2, "o", ms=2)
    plt.plot(npem3, npes3, "o", ms=2)
    plt.xlim(0, 1000)
    plt.ylim(0, 50)


    plt.show()



if __name__ == "__main__":
    main()
