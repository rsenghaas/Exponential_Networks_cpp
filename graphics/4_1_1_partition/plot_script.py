import os
import glob
import re
import shutil
import pathlib
import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('Agg')

plt.rcParams["font.family"] = 'serif'
plt.rcParams["mathtext.fontset"] = 'cm'
plt.rc('font', size=5)

def transform(z):
    return z / (1/4 - z)


green = '#11ff1f'
green_dark = '#0abb10'
red = '#ff1011'
black = '#222222'
grey = '#aaaaaa'
light_grey = '#eeeeee'
gray_dark = '#444444'
medium_gray = '#888888'
blue = '#1515dd'
orange = '#ffa500'
pink = '#ee0aaa'
pastel_blue = '#c1c1f4'
pastel_red = '#f4c1c1'
pastel_green = '#c8ffcc'
lime = '#caff00'
purple = '#9a80f4'


path_colors = {0: grey, 3: green, 6: green, 
               8: pastel_green, 9: red, 11: grey, 12:grey, 13: green, 14: grey}
# path_colors = {0: light_grey, 4: green, 7: red, 8: green}

partition = '4_1_1_partition'
# TODO: Need to change directories, so we can actually run that from the subdirectories.
if os.path.exists('./plot_script.py'):
    output_dir = '.'
else:
    output_dir = 'graphics/' + partition
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    if os.path.exists(output_dir + '/data'):
        shutil.rmtree(output_dir + '/data')
    shutil.copytree('data', output_dir + '/data')
    shutil.copy('networks/adhm.cpp', output_dir + '/adhm.cpp')
    shutil.copy('./sandbox.py', output_dir + '/plot_script.py')


sing_tf = np.array([0, -1])
branch = np.array([-0.25])
branch_tf = transform(branch)

if True:
    fig = plt.figure(dpi=800, figsize=(5,3))
    current_path = pathlib.Path(__file__).parent.resolve()
    os.chdir('data/path_data')
    for file in glob.glob("*.csv"):
        i = int(re.search('path_data_(.+?).csv', file).group(1))
        print(i)
        if i in [0]:
            continue
        data = np.loadtxt(file, delimiter=",", dtype=np.complex_)
        data = transform(data)
        x_data = data[:, 0]
        if i in path_colors:
            order = 1
            color = path_colors[i]
            if path_colors[i] == green:
                order = 0
            if path_colors[i] == grey:
                order = 2
        else:
            order = 3
            color = black
        plt.plot(x_data.real, x_data.imag, linewidth=0.6, color=color,
                 zorder=order)
    os.chdir(current_path)
    # plt.axis([-5.0, 5.0, -5.0, 5.0])

    plt.plot(sing_tf.real, sing_tf.imag, color='white', marker='o', markersize=4,
             fillstyle='full', linestyle='none', mew=0.4, zorder=3)

    plt.plot(sing_tf.real, sing_tf.imag, color=blue, marker='o', markersize=4,
             fillstyle='none', linestyle='none', mew=0.4, zorder=3)

    plt.plot(branch_tf.real, branch_tf.imag, color=orange, marker='x',
             markersize=5,
             fillstyle='none', linestyle='none', mew=2, zorder=3)


    plt.axis([-1.8, 0.1, -0.65, 0.55])
    ax = plt.gca()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.axis('off')


    plt.plot([-1.37, -1.7], [0.05, 0.5], color=black, linewidth=0.2)
    plt.plot([-1.36, -1.5], [0.01, 0.5], color=black, linewidth=0.2)
    plt.plot([-1.343, -1.3], [-0.03, 0.5], color=black, linewidth=0.2)
    plt.text(-1.75, 0.52, "$a_1 b_2 | b_1 b_2$")
    plt.text(-1.55, 0.52, "$b_1 | c_1, c_3 | c_3$")
    plt.text(-1.35, 0.52, "$c_1 | d_1$")

    plt.plot([-1.38, -1.7], [-0.07, -0.5], color=black, linewidth=0.2)
    plt.plot([-1.356, -1.5], [-0.05, -0.5], color=black, linewidth=0.2)
    
    plt.text(-1.75, -0.55, "$a_1 | b_1 b_2$")
    plt.text(-1.55, -0.55, "$b_1 | c_1, b_2 | c_3$")

    plt.plot([-1.367, -1.6], [-0.0, 0.1], color=black, linewidth=0.2)
    plt.plot([-1.349, -1.6], [-0.04, -0.1], color=black, linewidth=0.2)
    
    plt.text(-1.62, 0.1, "$b_2$")
    plt.text(-1.62, -0.1, "$c_3$")


    fig.tight_layout()
    plt.savefig(output_dir + '/network.png', dpi=fig.dpi)
    plt.savefig(output_dir + '/network.pdf', dpi=fig.dpi)
