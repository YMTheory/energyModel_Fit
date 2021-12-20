import numpy as np
import matplotlib.pyplot as plt


Ndim = 5

coeff_matrix = np.ones((Ndim, Ndim))

row, column = 0, 0
with open("/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/gam+B12_kSimulation_kSimulationCer_kNPE_fixedq0_correlation.txt") as f:
    for lines in f.readlines():
        line = lines.strip("\n")
        data = line.split(" ")
        for column in range(Ndim):
            coeff_matrix[row, column] = float(data[column])
        row += 1

print(coeff_matrix)

#coeff_matrix_bk = np.ones((Ndim, Ndim))
#for i in range(Ndim):
#    for j in range(Ndim):
#        coeff_matrix_bk[i ,j] = coeff_matrix[i, j]
#
#coeff_matrix_bk[0, 1] = coeff_matrix[1, 0]
#coeff_matrix_bk[1, 0] = coeff_matrix[0, 1]
#for i in range(Ndim):
#    coeff_matrix_bk[i, 0] = coeff_matrix[i, 1]
#    coeff_matrix_bk[i ,1] = coeff_matrix[i, 0]
#
#for i in range(Ndim):
#    coeff_matrix_bk[0, i] = coeff_matrix[1, i]
#    coeff_matrix_bk[1, i] = coeff_matrix[0, i]
#coeff_matrix_bk[0, 0] = coeff_matrix[0, 0]
#coeff_matrix_bk[1, 1] = coeff_matrix[1, 1]


#fig, ax = plt.subplots(figsize=(8, 8))
#im = ax.imshow(coeff_matrix, extent=[0, Ndim, 0, Ndim], aspect="auto", cmap="coolwarm")
#cbar = ax.figure.colorbar(im, ax=ax)

import seaborn as sns
ax = sns.heatmap(coeff_matrix, annot=True, fmt=".3f", cmap="coolwarm")


par_name = [r"$A$", r"$k_B$", r"$k_C$", r"$p_1$", r"$p_2$"]
#par_name = [r"$k_A$", r"$k_B$", r"$k_C$, "r"$A$", r"$p_0$", r"$p_1$", r"$p_2$"]
revert_par = list(reversed(par_name))
ax.set_xticks(np.arange(0, Ndim, 1))
ax.set_yticks(np.arange(0, Ndim, 1))
ax.set_xticklabels(par_name, fontsize=15)
ax.set_yticklabels(par_name, fontsize=15)
plt.xticks(rotation=45)

#plt.title("Correlation Coefficients")

plt.savefig("CovMat.pdf")

plt.show()

