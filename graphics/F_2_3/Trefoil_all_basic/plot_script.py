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
import random
import matplotlib.patches as mpatchespoint
# matplotlib.use('Agg')


def transform(z):
    return z
    # return z / (-4 - z)



def random_hex_colors(n: int) -> dict[int, str]:
    """Return dict {i: '#RRGGBB'} for i = 1..n with random RGB colors."""
    rng = random.Random(5) 
    def rand_hex():
        return "#{:02X}{:02X}{:02X}".format(
            rng.randint(50, 255),
            rng.randint(50, 255),
            rng.randint(50, 255),
        )
    
    return {i: rand_hex() for i in range(0, n)}

green = '#11ff1f'
green_dark = '#0abb10'
red = '#ff1011'
black = '#222222'
grey = '#aaaaaa'
grey_dark = '#555555'
blue = '#1515dd'
orange = '#ffa500'

Mblue = '#5E81B5'
Morange ='#E19C24'
Mgreen = '#8FB032'


# path_colors = random_hex_colors(25);
path_colors = {12: Morange,16: Morange,18: Morange,23: Morange,29: Morange,25:
               Morange, 
               0: Mblue, 1: Mblue, 2: Mblue, 3: Mblue, 4: Mblue, 5: Mblue, 6:
               Mblue, 7: Mblue, 8: Mblue, 9: Mblue,}

# 12 15

partition = 'Trefoil_D0_D2_full'
# TODO: Need to change directories, so we can actually run that from the subdirectories.
if os.path.exists('./plot_script.py'):
    output_dir = '.'
else:
    output_dir = 'graphics/F_2_3/' + partition
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    if os.path.exists(output_dir + '/data'):
        shutil.rmtree(output_dir + '/data')
    shutil.copytree('data', output_dir + '/data')
    shutil.copy('./F_2_3.py', output_dir + '/plot_script.py')


sing_tf = np.array([0])
branch = np.loadtxt("data/map_data/branch_points.csv", dtype=np.complex128)
print(branch)
branch_tf = transform(branch)

if sys.argv[1] == "all":
    current_path = pathlib.Path(__file__).parent.resolve()
    for k in range(30, 31):
        print(k)
        os.chdir('data/path_data')
        fig = plt.figure(dpi=1000)
        for file in glob.glob("*.csv"):
            highlight = [4,5,6,7, 8, 9,18]
            # highlight = [12,15,16,17,19,20,22, 24,29,79, 83]
            i = int(re.search('path_data_(.+?).csv', file).group(1))
            if (i not in highlight):
                continue
            if i > 300:
                continue
            data = np.loadtxt(file, delimiter=",", dtype=np.complex128)
            x_data = transform(data[:, 0])
            if i == 12:
                x_data = x_data[:-30]
            if i == 16:
                x_data = x_data[:-54]
            if True:
                if i == 6:
                    x_data = x_data[:-236]
                if i == 7:
                    x_data = x_data[:-82] 
                if i == 4:
                    x_data = x_data[:596]
                if i == 5:
                    x_data = x_data[:442]

            if i in path_colors:
                order = 1
                color = path_colors[i]
                linesize = 1.3
            elif i in highlight:
                order = 2
                color = black
                linesize = 1.3
            else:
                order = 1
                color = grey
                linesize = 0.3
            plt.plot(x_data.real, x_data.imag, linewidth=linesize, color=color,
                     zorder=order)
        os.chdir(current_path)
        # plt.axis([-5.0, 5.0, -5.0, 5.0])
        
        plt.plot(sing_tf.real, sing_tf.imag, color='white', marker='o', markersize=4,
                 fillstyle='full', linestyle='none', mew=0.4)

        plt.plot(sing_tf.real, sing_tf.imag, color=blue, marker='o', markersize=4,
                 fillstyle='none', linestyle='none', mew=0.4)

        plt.plot(branch_tf.real, branch_tf.imag, color=orange, marker='x',
                 markersize=5,
                 fillstyle='none', linestyle='none', mew=0.4)

        plt.axis([-4, 5, -5,4])
        ax = plt.gca()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        plt.axis('off')
        fig.tight_layout()
        """ ax.legend(
            handles=patches,
            loc='upper left',
            fontsize=8,              # smaller text
            frameon=True,
            title="Colors",
            title_fontsize=9,
            handlelength=1,          # smaller color box width
            handleheight=0.8,        # smaller color box height
            labelspacing=0.3,        # vertical spacing
            borderpad=0.4
        )   """
        plt.savefig(output_dir + f'/network.png', dpi=fig.dpi)
        plt.savefig(output_dir + f'/network.pdf', dpi=fig.dpi)
        plt.show()
        matplotlib.pyplot.close()


else:
    print("Rendering single path")
    fig = plt.figure(dpi=1000)
    try:
        intersection_file = 'data/map_data/intersections.csv'
        inter_data = np.loadtxt(
            intersection_file, delimiter=",", dtype=np.complex128)
        for inter in inter_data:
            if str(int(inter[3])) in sys.argv[1:] and str(int(inter[4])) in sys.argv[1:]:
                pt = transform(inter[0])
                plt.plot(pt.real, pt.imag, color=red, marker='o', markersize=2,
                         fillstyle='full', linestyle='none', mew=0.4)
    except:
        print("No intersections found")
    for s in sys.argv[1:]:
        s1 = int(s)
        print(s1)
        # for i in range(6):
        filename = f'data/path_data/path_data_{s1}.csv'
        data = np.loadtxt(filename, delimiter=",", dtype=np.complex128)
        data = data
        x_data = transform(data[:, 0])
        plt.plot(x_data.real, x_data.imag, linewidth=0.5, color=black,
                 zorder=2)
    plt.plot(sing_tf.real, sing_tf.imag, color='white', marker='o', markersize=4,
             fillstyle='full', linestyle='none', mew=0.2)

    plt.plot(sing_tf.real, sing_tf.imag, color=blue, marker='o', markersize=4,
             fillstyle='none', linestyle='none', mew=0.2)

    plt.plot(branch_tf.real, branch_tf.imag, color=orange, marker='x',
             markersize=3,
             fillstyle='none', linestyle='none', mew=0.2)

    plt.axis([-5, 5, -5, 5])
    ax = plt.gca()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.axis('off')
    fig.tight_layout()
    plt.savefig(output_dir + '/network.png', dpi=fig.dpi)
    plt.savefig(output_dir + '/network.pdf', dpi=fig.dpi)
