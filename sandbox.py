import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob, os
import pathlib

intersection_file = 'data/intersection_data/test.csv'
intersection_data = np.array([]) #np.loadtxt(intersection_file, delimiter=",", dtype=np.complex_)

def transform(z):
    # return z
    return z / ( 1/4 - z)

if sys.argv[1] == "all":
    fig = plt.figure(dpi=2000)
    intersection_data = transform(intersection_data);
    plt.plot(intersection_data.real, intersection_data.imag, "ro", markersize=1, markeredgewidth=0, markerfacecolor="red")
    current_path = pathlib.Path(__file__).parent.resolve()
    os.chdir('data/path_data')
    for file in glob.glob("*.csv"):
        data = np.loadtxt(file ,delimiter=",", dtype=np.complex_)
        data = transform(data)
        x_data = data[:, 0]
        plt.plot(x_data.real, x_data.imag, linewidth=0.1)
    os.chdir(current_path)
    # plt.axis([-3.0, 3.0, -3.0, 3.0])
    plt.axis([-1.8, 0.1, -0.8, 0.8])
    plt.savefig('graphics/test_graphic.png', dpi=fig.dpi)
    # plt.savefig('graphics/n=8/(1,2,2,1,1,1).png', dpi=fig.dpi)

else:
    s1 = int(sys.argv[1])
    fig = plt.figure(dpi=1000)
    # for i in range(6):
    filename= f'data/path_data/path_data_{s1}.csv'
    data = np.loadtxt(filename,delimiter=",", dtype=np.complex_)
    x_data = data[:, 0]
    plt.plot(x_data.real, x_data.imag)
    plt.savefig('graphics/test_graphic.png', dpi=fig.dpi)

