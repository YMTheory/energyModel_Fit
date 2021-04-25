import numpy as np
import matplotlib.pyplot as plt

name = ["Cs137", "Mn54", "Ge68", "K40", "nH", "Co60", "AmBe", "nC12", "AmC"]
Etrue = [0.662, 0.835, 1.022, 1.461, 2.223, 2.506, 4.43, 4.94, 6.13]
Etrue_nonl = [0.662, 0.835, 0.511, 1.461, 2.223, 1.253, 4.43, 4.94, 6.13]

def ReadModelPrediction(filename):
    simnonl, simnonl_err, simres, simres_err = [], [], [], []
    calcnonl, calcnonl_err, calcres, calcres_err = [], [], [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            simnonl.append(float(data[0]))
            simnonl_err.append(float(data[1]))
            #simnonl_err.append(0.001)
            calcnonl.append(float(data[2]))
            calcnonl_err.append(float(data[3]))
            #calcnonl_err.append(0.001)
            simres.append(float(data[4]))
            simres_err.append(float(data[5]))
            calcres.append(float(data[6]))
            calcres_err.append(float(data[7]))


    # change nonlinearity data order :
    new_id = [2, 0, 1, 5, 3, 4, 6, 7, 8]
    simnonl_new, simnonl_new_err, calcnonl_new, calcnonl_new_err, Etrue_nonl_new = [], [], [], [], []
    for i in new_id:
        Etrue_nonl_new.append(Etrue_nonl[i])
        simnonl_new.append(simnonl[i]*Etrue[i]/Etrue_nonl[i])
        simnonl_new_err.append(simnonl_err[i]*Etrue[i]/Etrue_nonl[i])
        calcnonl_new.append(calcnonl[i]*Etrue[i]/Etrue_nonl[i])
        calcnonl_new_err.append(calcnonl_err[i]*Etrue[i]/Etrue_nonl[i])


    nonl_diff = (np.array(calcnonl_new) - np.array(simnonl_new)) / np.array(simnonl_new)
    nonl_diff_err = np.sqrt(np.array(calcnonl_new_err)**2/np.array(simnonl_new) + np.array(simnonl_new_err)**2 * np.array(calcnonl_new)**2 / np.array(simnonl_new)**4)
    res_diff = (np.array(calcres) - np.array(simres)) / np.array(simres)
    res_diff_err = np.sqrt(np.array(calcres_err)**2/np.array(simres) + np.array(simres_err)**2 * np.array(calcres)**2 / np.array(simres)**4)

    return Etrue_nonl_new, simnonl_new, simnonl_new_err, calcnonl_new, calcnonl_new_err, Etrue, simres, simres_err, calcres, calcres_err, res_diff, res_diff_err
