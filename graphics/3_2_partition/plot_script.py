import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob, os
import pathlib
import shutil
import re

def transform(z):
    # return z
    return z / ( 1/4 - z)



green = '#11ff1f'
green_dark = '#0abb10'
red = '#ff1011'
black = '#222222'
grey = '#aaaaaa'
blue = '#1515dd'
orange = '#ffa500'

path_colors = {3: green, 5 :red, 7: green, 9: green, 11: red}

partition = '3_2_partition'
#TODO: Need to change directories, so we can actually run that from the subdirectories.
if os.path.exists('./plot_script.py'):
    output_dir = '.'
else:
    output_dir = 'graphics/' + partition

sing_tf = np.array([0, -1])
branch = np.array([-0.25])
branch_tf = transform(branch) 

if True:
    fig = plt.figure(dpi=800)
    current_path = pathlib.Path(__file__).parent.resolve()
    os.chdir('data/path_data')
    for file in glob.glob("*.csv"):
        i = int(re.search('path_data_(.+?).csv', file).group(1))
        data = np.loadtxt(file ,delimiter=",", dtype=np.complex_)
        data = transform(data)
        x_data = data[:, 0]
        if i in path_colors:
            order = 1
            color = path_colors[i]
        else:
            order = 2
            color = black
        plt.plot(x_data.real, x_data.imag, linewidth=0.8, color=color,
                 zorder=order)
    os.chdir(current_path)
    # plt.axis([-5.0, 5.0, -5.0, 5.0])
    
    
    plt.plot(sing_tf.real, sing_tf.imag, color='white', marker='o', markersize=4,
                fillstyle='full', linestyle='none', mew=0.4);

    plt.plot(sing_tf.real, sing_tf.imag, color=blue, marker='o', markersize=4,
                fillstyle='none', linestyle='none', mew=0.4);

    plt.plot(branch_tf.real, branch_tf.imag, color=orange, marker='x',
             markersize=5,
                fillstyle='none', linestyle='none', mew=2);

    
    plt.axis([-1.8, 0.1, -0.8, 0.8])
    ax =plt.gca()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.axis('off')
    fig.tight_layout()
    plt.savefig(output_dir + '/network.png', dpi=fig.dpi)
    plt.savefig(output_dir + '/network.pdf', dpi=fig.dpi)


