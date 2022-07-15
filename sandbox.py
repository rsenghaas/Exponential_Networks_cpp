import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

intersection_file = 'data/intersection_data/test.csv'
intersection_data = np.loadtxt(intersection_file, delimiter=",", dtype=np.complex_)
print(intersection_data)

if sys.argv[1] == "all":
    s1 = 0
    fig = plt.figure()
    plt.plot(intersection_data.real, intersection_data.imag, "ro")
    filename= f'data/path_data/path_data_{s1}.csv'
    path = Path(filename)
    while path.is_file():
        s1 += 1
        data = np.loadtxt(filename,delimiter=",", dtype=np.complex_)
        x_data = data[:, 0]
        plt.plot(x_data.real, x_data.imag)
        filename= f'data/path_data/path_data_{s1}.csv'
        path = Path(filename)
    plt.savefig('graphics/test_graphics.png', dpi=fig.dpi)

else:
    s1 = int(sys.argv[1])
    fig = plt.figure()
    # for i in range(6):
    filename= f'data/path_data/path_data_{s1}.csv'
    data = np.loadtxt(filename,delimiter=",", dtype=np.complex_)
    x_data = data[:, 0]
    plt.plot(x_data.real, x_data.imag)
    plt.savefig('graphics/test_graphic.png', dpi=fig.dpi)

# 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
