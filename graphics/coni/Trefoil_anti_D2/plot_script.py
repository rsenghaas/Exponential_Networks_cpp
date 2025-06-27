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


# WC at -0.0290095


def transform(z):
    return z / (0.3 - z)


green = '#11ff1f'
green_dark = '#0abb10'
black = '#222222'
grey = '#aaaaaa'
blue = '#1515dd'
pink = '#ff66bb'

red = '#ff1011'
light_red = '#ff7070'
dark_red = '#880505'
orange = '#ffa500'
yellow = '#ffd020'
light_yellow = '#ffee66'


path_colors = {3: black, 21: black, 
               23: black, 5: black, 18: black, 19: black,
               25: black, 10: black, 9: black,
            }
               # path_colors = {0: green, 1: grey, 2: blue, 3:red, 4: green}

partition = 'Trefoil_anti_D2'
# TODO: Need to change directories, so we can actually run that from the subdirectories.
if os.path.exists('./plot_script.py'):
    output_dir = '.'
else:
    output_dir = 'graphics/coni/' + partition
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    if os.path.exists(output_dir + '/data'):
        shutil.rmtree(output_dir + '/data')
    shutil.copytree('data', output_dir + '/data')
    shutil.copy('networks/coni.cpp', output_dir + '/coni.cpp')
    shutil.copy('./coni.py', output_dir + '/plot_script.py')


# sing_tf = np.array([0, -1])
sing_tf = np.array([0, -1])
branch = np.loadtxt("data/map_data/branch_points.csv", dtype=np.complex_)
branch_tf = transform(branch)
print(branch)
print(branch_tf)

if sys.argv[1] == "all":
    fig = plt.figure(dpi=800)
    current_path = pathlib.Path(__file__).parent.resolve()
    os.chdir('data/path_data')
    for file in glob.glob("*.csv"):
        i = int(re.search('path_data_(.+?).csv', file).group(1))
        if ( 
            i in [] + [] # [18,19,20,21,22,23]
        ): # [3, 6, 21, 25]                                             
            continue
        print(i)
        data = np.loadtxt(file, delimiter=",", dtype=np.complex_)
        data = transform(data)
        x_data = data[:, 0]
        if i in path_colors:
            order = 3
            color = path_colors[i]
        else:
            order = 2
            color = grey
        if i in []:
            order = 1
        if i in [5]:
            order = 4
        plt.plot(x_data.real, x_data.imag, linewidth=0.2, color=color,
                 zorder=order)
    os.chdir(current_path)
    # plt.axis([-5.0, 5.0, -5.0, 5.0])

    plt.plot(sing_tf.real, sing_tf.imag, color='white', marker='o',
             markersize=0.6,
           fillstyle='full', linestyle='none', mew=0.4, zorder=5)

    plt.plot(sing_tf.real, sing_tf.imag, color=blue, marker='o', markersize=0.6,
             fillstyle='none', linestyle='none', mew=0.4, zorder=5)

    plt.plot(branch_tf.real, branch_tf.imag, color=orange, marker='x',
             markersize=3,
             fillstyle='none', linestyle='none', mew=0.5, zorder=5)

    plt.axis([-2.5, 2.5, -2.5, 2.5])
    ax = plt.gca()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.axis('off')
    fig.tight_layout()
    plt.savefig('graphics/test_graphic.png', dpi=fig.dpi)
    plt.savefig(output_dir + '/network.png', dpi=fig.dpi)
    plt.savefig(output_dir + '/network.pdf', dpi=fig.dpi)


else:
    print("Rendering single path")
    s1 = int(sys.argv[1])
    fig = plt.figure(dpi=1000)
    # for i in range(6):
    filename = f'data/path_data/path_data_{s1}.csv'
    data = np.loadtxt(filename, delimiter=",", dtype=np.complex_)
    mask = [i for i in range(min(len(data),3000))]
    x_data = transform(data[:, 0])[mask]
    plt.plot(x_data.real, x_data.imag, linewidth=0.8, color=black,
             zorder=2)
    plt.plot(sing_tf.real, sing_tf.imag, color='white', marker='o', markersize=2,
             fillstyle='full', linestyle='none', mew=0.4)

    plt.plot(sing_tf.real, sing_tf.imag, color=blue, marker='o', markersize=2,
             fillstyle='none', linestyle='none', mew=0.4)

    plt.plot(branch_tf.real, branch_tf.imag, color=orange, marker='x',
             markersize=5,
             fillstyle='none', linestyle='none', mew=2)
    plt.axis([-0.5, 0.5, -1, 1])
    ax = plt.gca()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.axis('off')
    fig.tight_layout()
    plt.savefig('graphics/test_graphic.png', dpi=fig.dpi)
    plt.savefig(output_dir + '/network.png', dpi=fig.dpi)
    plt.savefig(output_dir + '/network.pdf', dpi=fig.dpi)
