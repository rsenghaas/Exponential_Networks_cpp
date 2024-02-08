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


path_colors = {3: blue, 6: pastel_green, 7: red, 9: pastel_green , 10: red, 12:
               pastel_green, 13: red}
# path_colors = {0: light_grey, 4: green, 7: red, 8: green}

partition = '1_partition_wallcrossing'
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
    fig = plt.figure(dpi=800, figsize=(5,2.5))
    current_path = pathlib.Path(__file__).parent.resolve()
    os.chdir('data/path_data')
    for file in glob.glob("*.csv"):
        i = int(re.search('path_data_(.+?).csv', file).group(1))
        print(i)
        data = np.loadtxt(file, delimiter=",", dtype=np.complex_)
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
             fillstyle='full', linestyle='none', mew=0.4)

    plt.plot(sing_tf.real, sing_tf.imag, color=blue, marker='o', markersize=4,
             fillstyle='none', linestyle='none', mew=0.4)

    plt.plot(branch_tf.real, branch_tf.imag, color=orange, marker='x',
             markersize=5,
             fillstyle='none', linestyle='none', mew=2)

    plt.axis([-1.8, 0.1, -0.5, 0.5])
    ax = plt.gca()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.axis('off')
    fig.tight_layout()
    plt.savefig(output_dir + '/network.png', dpi=fig.dpi)
    plt.savefig(output_dir + '/network.pdf', dpi=fig.dpi)
