import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

i = 0
with open(f'data/intersection_data/unsuccessful_{i}.csv', 'r') as f:
    lines = f.readlines()

with open('data/intersection_data/test.csv', 'r') as f:
    inter_lines = f.readlines()

inter_pts = np.array([complex(n) for n in inter_lines[0].split(",")])

fig = plt.figure()
# plt.scatter(inter_pts[i].real, inter_pts[i].imag)
for l in lines:
    intersection_data = np.array([complex(n) for n in l.split(",")])
    print(intersection_data)
    plt.plot(intersection_data.real, intersection_data.imag)
    
# intersection_data_0 = np.loadtxt(intersection_file, userows=0, delimiter=",", dtype=np.complex_)
# intersection_data_1 = np.loadtxt(intersection_file, userows=1, delimiter=",", dtype=np.complex_)
# print(intersection_data_0)
# print(intersection_data_1)

# 
#
# plt.plot(intersection_data_1.real, intersection_data_1.imag)
plt.savefig("graphics/debug.png", dpi = fig.dpi)
