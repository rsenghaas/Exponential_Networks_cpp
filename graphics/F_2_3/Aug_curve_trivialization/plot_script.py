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
from scipy.interpolate import CubicSpline
from collections import defaultdict, deque
import csv
# matplotlib.use('Agg')

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
path_colors = {}

# 12 15

partition = 'Aug_curve_D0'
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


intersection_file = r"/home/spamdoodler/Uni/Exponential_Networks_cpp/data/intersection_data/database.csv";
sing_tf = np.array([0])
branch = np.loadtxt("data/map_data/branch_points.csv", dtype=np.complex128)
print(branch)
branch_tf = transform(branch)

if sys.argv[1] == "all":
    current_path = pathlib.Path(__file__).parent.resolve()
    fig = plt.figure(dpi=1000)
    os.chdir(current_path)
    # plt.axis([-5.0, 5.0, -5.0, 5.0])
    
    plt.plot(sing_tf.real, sing_tf.imag, color='white', marker='o', markersize=4,
             fillstyle='full', linestyle='none', mew=0.4)

    plt.plot(sing_tf.real, sing_tf.imag, color=blue, marker='o', markersize=4,
             fillstyle='none', linestyle='none', mew=0.4)

    plt.plot(branch_tf.real, branch_tf.imag, color=orange, marker='x',
             markersize=5,
             fillstyle='none', linestyle='none', mew=0.4)
    branch_cuts = [
            [branch_tf[6],0.5 + 0.2j, 1 - 0.5j, branch_tf[4]],
            [branch_tf[7],-0.3-0.5j,-0.7 + 0.5j, branch_tf[8]],
            [branch_tf[5],2.7+2j,branch_tf[9]],
            [0, 0.5+0.4j, 3 - 0.1j, 5, 15]
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
    
    branch_label = {
        4: "$[13]$",
        5: "$[02]$",
        6: "$[13]$",
        7: "$[23]$",
        8: "$[23]$",
        9: "$[02]$",
    }
    label_offset = {
        4: (-0.1, 0.1),
        5: (-0.1,0.1),
        6: (-0.25,0),
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
