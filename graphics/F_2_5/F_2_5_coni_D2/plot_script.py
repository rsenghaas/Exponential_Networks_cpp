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
from collections import defaultdict, deque
import csv
# matplotlib.use('Agg')

partition ="F_2_5_coni_D2"

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

def build_truncate(intersection_file, highlight):
    lookup = defaultdict(list)

    with open(intersection_file, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            id_new = int(row[0])
            id_A = int(row[1])
            cutoff_A = int(row[2])
            id_B = int(row[3])
            cutoff_B = int(row[4])

            lookup[id_new].append((id_A, cutoff_A + 1,  id_B, cutoff_B + 1))

# ----------------------
# Step 2: Recursive propagation (BFS)
# ----------------------
    truncate = {}
    queue = deque(highlight)

    while queue:
        current = queue.popleft()

        if current not in lookup:
            continue

        for id_A, cutoff_A, id_B, cutoff_B in lookup[current]:

            # --- process id_A ---
            old_A = truncate.get(id_A, -1)
            if cutoff_A > old_A:
                truncate[id_A] = cutoff_A
                if id_A >= 30:
                    queue.append(id_A)

            # --- process id_B ---
            old_B = truncate.get(id_B, -1)
            if cutoff_B > old_B:
                truncate[id_B] = cutoff_B
                if id_B >= 30:
                    queue.append(id_B)

    return truncate

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
path_colors = {}

# 12 15

# TODO: Need to change directories, so we can actually run that from the subdirectories.
if os.path.exists('./plot_script.py'):
    output_dir = '.'
else:
    output_dir = 'graphics/F_2_5/' + partition
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    if os.path.exists(output_dir + '/data'):
        shutil.rmtree(output_dir + '/data')
    shutil.copytree('data', output_dir + '/data')
    shutil.copy('./F_2_3.py', output_dir + '/plot_script.py')

try:
    intersection_file = r"/home/spamdoodler/Uni/Exponential_Networks_cpp/data/intersection_data/database.csv"
except:
    pass

sing_tf = np.array([0])
branch = np.loadtxt("data/map_data/branch_points.csv", dtype=np.complex128)
print(branch)
branch_tf = transform(branch)

if sys.argv[1] == "all":
    current_path = pathlib.Path(__file__).parent.resolve()
    show = False
    for k in range(30, 31):
        os.chdir('data/path_data')
        fig = plt.figure(dpi=1000)
        for file in glob.glob("*.csv"):
            highlight = []
            if len(sys.argv[2:]) == 0:
                pass #highlight = []
            else:
                for s in sys.argv[2:]:
                    highlight.append(int(s))
            truncate = {}
            if True:
                truncate = build_truncate(intersection_file, highlight)
            # path_colors[3] = red
            # path_colors[201] = green
            # truncate[3] += 60
            # truncate[201] += 60
            # path_colors[19] = blue
            # truncate[18] = -85
            # truncate[1] = -100
            # path_colors[3] = red
            # path_colors[17] = green
            # truncate[17] = -100

            i = int(re.search('path_data_(.+?).csv', file).group(1))
            if (i not in highlight and i not in truncate):
                if len(highlight) > 0:
                    continue
            print(i)
            data = np.loadtxt(file, delimiter=",", dtype=np.complex128)
            if len(data.shape) < 2:
                continue
            x_data = transform(data[:, 0])  
            for k in truncate:
                if i == k:
                    x_data = x_data[:truncate[k]]
            if i in path_colors:
                order = 3
                color = path_colors[i]
                linesize = 0.3
                if show:
                    linesize = 0.3
            elif i in highlight or i in truncate:
                order = 2
                color = black
                linesize = 1.3
                if show:
                    linesize = 0.3
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
        plt.savefig('graphics/test_graphic.png', dpi=fig.dpi)
        plt.savefig(output_dir + f'/network.png', dpi=fig.dpi)
        plt.savefig(output_dir + f'/network.pdf', dpi=fig.dpi)
        plt.savefig(output_dir + f'/network_{k}.pdf', dpi=fig.dpi)
        if show:
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
    plt.savefig('graphics/test_graphic.png', dpi=fig.dpi)
    plt.savefig(output_dir + '/network.png', dpi=fig.dpi)
    plt.savefig(output_dir + '/network.pdf', dpi=fig.dpi)
