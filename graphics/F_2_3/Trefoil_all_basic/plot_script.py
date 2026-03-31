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
from scipy.interpolate import CubicSpline
from collections import defaultdict, deque
import matplotlib as mpl

mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 4,          # base font size
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
})


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

def smooth_curve_through(points, n_per_segment=200):
    """
    points: list/array of complex numbers
    returns smooth complex-valued curve passing through all points
    """
    points = np.asarray(points, dtype=np.complex128)

    # Parameter along curve (arc-length style)
    distances = np.abs(np.diff(points))
    t = np.concatenate(([0], np.cumsum(distances)))
    t /= t[-1]

    # Separate real/imag interpolation
    cs_real = CubicSpline(t, points.real, bc_type='natural')
    cs_imag = CubicSpline(t, points.imag, bc_type='natural')

    t_fine = np.linspace(0, 1, n_per_segment * (len(points)-1))

    curve = cs_real(t_fine) + 1j * cs_imag(t_fine)
    return curve

def add_wave(curve, amplitude=0.05, wavelength=0.5):
    """
    Adds sinusoidal displacement along the normal direction.
    """
    # Tangent
    dcurve = np.gradient(curve)
    normal = 1j * dcurve
    normal /= np.abs(normal)

    # Arc-length parameter
    s = np.linspace(0, 1, len(curve))

    wave = amplitude * np.sin(2*np.pi * s / wavelength)

    return curve + wave * normal



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
            A1 = [12, 23,18]
            A2 = [18]
            A3 = [16, 23, 29, 25]
            B1 = [0,1]
            B2 = [4,5,6,7,8,9]
            B3 = [2,3]
            for i in B2:
                path_colors[i] = Mgreen
            highlight = []
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
            if i == 18:
                pass
                # x_data = x_data[:-340]
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
                linesize = 0.6
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
                 fillstyle='full', linestyle='none', mew=0.4, zorder=4)

        plt.plot(sing_tf.real, sing_tf.imag, color=blue, marker='o', markersize=4,
                 fillstyle='none', linestyle='none', mew=0.4, zorder =5)

        plt.plot(branch_tf.real, branch_tf.imag, color=orange, marker='x',
                 markersize=5,
                 fillstyle='none', linestyle='none', mew=0.4)
        
        draw_branch_cuts = True
        if draw_branch_cuts:
            branch_cuts = [
                [branch_tf[6],0.5 + 0.2j, 1 - 0.5j, branch_tf[4]],
                [branch_tf[7],-0.0-0.3j,0.1,
                 -0.2+0.4j, -0.5+1j, branch_tf[8]],
                [branch_tf[5],2.7+2j,branch_tf[9]],
                [0, -0.5+0.4j, -3 - 0.1j, -5, -15],
            ]
            for pts in branch_cuts:
                arc = 0
                for i in range(1, len(pts)):
                    arc += abs(pts[i] - pts[i - 1])
                wavelength = 0.15 / arc
                n = int(200 * arc / (len(pts) - 1))
                curve = smooth_curve_through(pts, n_per_segment=n)
                wave_curve = add_wave(curve, amplitude=0.03, wavelength=wavelength)
                plt.plot(wave_curve.real, wave_curve.imag, linewidth=0.3, color=orange)
            
            plt.plot([0,-5], [0,-0.5], linewidth=0.3, linestyle='--')

            branch_label = {
                4: "$[13]$",
                5: "$[02]$",
                6: "$[13]$",
                7: "$[23]$",
                8: "$[23]$",
                9: "$[02]$",
            }
            label_offset = {
                4: (-0.08, 0.1),
                5: (-0.08,0.1),
                6: (-0.08,0.1),
                7: (0.05,0),
                8: (0.1,0),
                9: (0.1,0),
            }
            for i in range(4,10):
                offset_x, offset_y = label_offset[i]
                plt.text(branch_tf[i].real + offset_x, 
                         branch_tf[i].imag + offset_y, 
                         branch_label[i])

            plt.text(-0.25,
                     0, 
                     "$[01]$")

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
        # plt.show()
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
