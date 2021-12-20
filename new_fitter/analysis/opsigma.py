import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def readCerSigma(name):
    data = pd.read_csv(name, sep=" ")
    ke = data["ke"].to_numpy()
    k = data["k"].to_numpy()
    ave_cerpe = data["ave_cerpe"].to_numpy()
    std_cerpe = data["std_cerpe"].to_numpy()
    ave_cerop = data["ave_cerop"].to_numpy()
    std_cerop = data["std_cerop"].to_numpy()

    return ke, k, ave_cerpe, std_cerpe, ave_cerop, std_cerop


if __name__ == "__main__" :

    filename = "/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/electron/op_recorder/cerSigma.txt"
    ke, k, ave_cerpe, std_cerpe, ave_cerop, std_cerop = readCerSigma(filename)

    plt.plot(ke, std_cerpe**2, "o-")
    plt.plot(ke, k**2*std_cerop**2+ave_cerpe, "v-")
    plt.plot(ke, std_cerop**2+(1-k)*ave_cerop-2*np.sqrt((1-k)*ave_cerop)*std_cerop)

    plt.show()


    






