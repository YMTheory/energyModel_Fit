import numpy as np
import matplotlib.pyplot as plt

coeff_matrix = np.ones((15, 15))

row, column = 0, 0
with open("../coeff.txt") as f:
    for lines in f.readlines():
        line = lines.strip("\n")
        data = line.split(" ")
        print(data)
        for column in range(15):
            coeff_matrix[row, column] = float(data[column])
        row += 1


fig, ax = plt.subplots()
im = ax.imshow(coeff_matrix, extent=[0, 16, 0, 16], aspect="auto", cmap="GnBu")
cbar = ax.figure.colorbar(im, ax=ax)


par_name = ["kA", "kB1", "kC", "Scale", "Cs137", "Mn54", "Ge68", "K40", "nH", "Co60", \
            "AmBe", "nC12", "AmC", "ra", "rb", "rc"]
revert_par = list(reversed(par_name))
ax.set_xticks(np.arange(0, 16, 1))
ax.set_yticks(np.arange(0, 16, 1))
ax.set_xticklabels(par_name)
ax.set_yticklabels(revert_par)
plt.xticks(rotation=45)

plt.title("Correlation Coefficients")

plt.show()

