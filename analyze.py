import os
import glob
import re
import shutil
import pathlib
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')

file = r"data/path_data/path_data_9.csv"
data = np.loadtxt(file, delimiter=",", dtype=np.complex_)
y1 = data[:, 1]
y2 = data[:, 2]
dy = y2 - y1
j = 0
for i in range(len(dy)):
    if abs(dy[i]) < 0.01:
        if j + 1 == i:
            j = i
        else:
            print(i, data[i])
            j = i
