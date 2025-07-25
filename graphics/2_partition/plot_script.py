import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob, os
import pathlib
import shutil
import re

plt.rcParams["font.family"] = 'serif'
plt.rcParams["mathtext.fontset"] = 'cm'
plt.rc('font', size=5);

def transform(z):
    # return z
    return z / ( 1/4 - z)

green = '#11ff1f'
green_dark = '#0abb10'
red = '#ff1011'
black = '#222222'
grey = '#aaaaaa'
gray_dark = '#444444'
blue = '#1515dd'
orange = '#ffa500'

path_colors = {3: green, 5 :green, 7: green, 10: grey, 11: red}

partition = '2_partition'
#TODO: Need to change directories, so we can actually run that from the subdirectories.
if os.path.exists('./plot_script.py'):
    output_dir = '.'
    data_dir = './data/path_data'
else:
    output_dir = 'graphics/' + partition
    data_dir = 'data/path_data'
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    if os.path.exists(output_dir + '/data'):
        shutil.rmtree(output_dir + '/data')
    shutil.copytree('data', output_dir + '/data')
    shutil.copy('./sandbox.py', output_dir + '/plot_script.py') 
    intersection_file = 'data/intersection_data/test.csv'
    intersection_data = np.array([]) #np.loadtxt(intersection_file, delimiter=",", dtype=np.complex_)


sing_tf = np.array([0, -1])
branch = np.array([-0.25])
branch_tf = transform(branch) 

fig = plt.figure(dpi=800, figsize=(5,2.5))
current_path = pathlib.Path(__file__).parent.resolve()
os.chdir(data_dir)
if True:
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
        plt.plot(x_data.real, x_data.imag, linewidth=0.6, color=color,
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

    cut_x = np.array([- 1.0/2 - i / 5000.0 for i in range(2500)])
    cut_y = 0.01 * np.sin(np.pi * cut_x * 50)
    
    # plt.plot(cut_x, cut_y, color=orange, linewidth=0.1,zorder=0)
    log_cut_minus = np.array([-1, -10]);
    log_cut_plus = np.array([0, 10]);
    '''plt.plot(log_cut_minus, 0*log_cut_minus, 
             color=gray_dark, 
             linestyle='dashed',
             dashes=(40, 25),
             zorder=0, linewidth=0.1)
    plt.plot(log_cut_plus, 0*log_cut_plus, 
             color=gray_dark, 
             linestyle='dashed',
             dashes=(40, 25),
             zorder=0, linewidth=0.1)
    ''' 
    # plt.text(-0.65, 0.25, r"$(+-)_0$")
    # plt.text(-0.7, -0.2, r"$(+-)_0$")
    
    plt.text(-0.65, 0.25, "2")
    plt.text(-0.65, -0.27, "1") 
    plt.plot([-0.865, -0.8], [0, 0.15], color=black, linewidth=0.2)
    plt.plot([-0.855, -0.77], [0, 0.15], color=black, linewidth=0.2)
    plt.text(-0.795, 0.14, "1")
    plt.text(-0.765, 0.14, "1")
    plt.text(-1.38, -0.04, "1")
    plt.text(-1.45, -0.04, "1")

    plt.axis([-1.8, 0.1, -0.55, 0.45])
    ax =plt.gca()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.axis('off')
    fig.tight_layout()
    if output_dir != '.':
        plt.savefig('graphics/test_graphic.png', dpi=fig.dpi)
    plt.savefig(output_dir + '/network_with_mult.png', dpi=fig.dpi)
    plt.savefig(output_dir + '/network_with_mult.pdf', dpi=fig.dpi)

