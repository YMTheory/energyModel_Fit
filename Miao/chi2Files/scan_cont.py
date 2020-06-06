#!/usr/bin/env python
# coding=utf-8

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set()

tmpkB = []; tmpkC = []; tmpChi2 = [];
with open("./Gamma.txt") as f:
    for lines in f.readlines():
        line = lines.strip("\n")
        data = line.split(" ")
        tmpkB.append(float(data[0]))
        tmpkC.append(float(data[1]))
        tmpChi2.append(float(data[2]))

scanChi2 = np.zeros((20, 20))
index = 0
for i in range(20):
    for j in range(20):
        scanChi2[i][j] = tmpChi2[index]
        index += 1

f, ax = plt.subplots(figsize=(16, 9))
sns.heatmap(scanChi2, annot=False, fmt=".2f",linewidths=.0, ax=ax,cmap="viridis")
plt.show()
    
