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

path_colors = {3: red, 5 :red}

partition = '3_2_1_partition'
#TODO: Need to change directories, so we can actually run that from the subdirectories.
if os.path.exists('./plot_script.py'):
    output_dir = './'
else:
    output_dir = 'graphics/' + partition

sing_tf = np.array([0, -1])
branch = np.array([-0.25])
branch_tf = transform(branch) 

if sys.argv[1] == "all":
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    if os.path.exists(output_dir + '/data'):
        shutil.rmtree(output_dir + '/data')
    shutil.copytree('data', output_dir + '/data')
    shutil.copy('./sandbox.py', output_dir + '/plot_script.py');
    intersection_file = 'data/intersection_data/test.csv'
    intersection_data = np.array([]) #np.loadtxt(intersection_file, delimiter=",", dtype=np.complex_)


    fig = plt.figure(dpi=800)
    intersection_data = transform(intersection_data);
    plt.plot(intersection_data.real, intersection_data.imag, "ro", markersize=1, markeredgewidth=0, markerfacecolor="red")
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
        plt.plot(x_data.real, x_data.imag, linewidth=0.2, color=color,
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
    plt.savefig('graphics/test_graphic.png', dpi=fig.dpi)
    plt.savefig(output_dir + '/network.png', dpi=fig.dpi)
    plt.savefig(output_dir + '/network.pdf', dpi=fig.dpi)


else:
    print("Rendering single path")
    s1 = int(sys.argv[1])
    fig = plt.figure(dpi=1000)
    # for i in range(6):
    filename= f'data/path_data/path_data_{s1}.csv'
    data = np.loadtxt(filename,delimiter=",", dtype=np.complex_)
    data = transform(data)
    x_data = data[:, 0]
    plt.plot(x_data.real, x_data.imag)
    plt.savefig('graphics/test_graphic.png', dpi=fig.dpi)

