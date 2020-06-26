#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

etrue = []; nonl_data = []; nonl_model = []
with open("./nonl_compare.txt") as f:
    for lines in f.readlines():
        line = lines.strip("\n")
        data = line.split(" ")
        etrue.append(float(data[0]))
        nonl_data.append(float(data[1]))
        nonl_model.append(float(data[2]))

plt.plot(etrue, nonl_data, "o", ms=0.1, label="data")
plt.plot(etrue, nonl_model, "o",ms=0.1, label="model")
plt.legend()
plt.show()
