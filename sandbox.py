import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

s1 = int(sys.argv[1])
fig = plt.figure()
# for i in range(6):
filename= f'data/path_data/path_data_{s1}.csv'
data = np.loadtxt(filename,delimiter=",", dtype=np.complex_)
x_data = data[:, 0]
plt.plot(x_data.real, x_data.imag)
plt.savefig('graphics/test_graphic.png', dpi=fig.dpi)