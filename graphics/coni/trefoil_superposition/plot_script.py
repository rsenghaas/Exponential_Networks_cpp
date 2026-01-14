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
plt.rcParams["font.family"] = 'serif'
plt.rcParams["mathtext.fontset"] = 'cm'
plt.rc('font', size=3)


# WC at -0.0290095


def transform(z):
    return z / (0.3 - z)

def deform_path(x, param):
    r = param[0]
    v = param[1]
    a = param[2]
    x = np.array(
            [x[i] + v * (r * (i / (len(x) - 1))**(a)
                         * (1 - i / (len(x) - 1)))**4 
             for i in range(len(x))])
    return x


def draw_cut(points):
# Parameter values in the interval [0, 1] (you can choose others)
    t = np.linspace(0, 1, len(points))
    t_vals = np.linspace(0,1, 1000)
    points = np.array(points)
# Separate real and imaginary parts
    x = t
    y_real = points.real
    y_imag = points.imag

# Polynomial regression of chosen degree
    degree = len(points) - 1
    coef_real = np.polyfit(x, y_real, degree)
    coef_imag = np.polyfit(x, y_imag, degree)

# Construct polynomial functions
    poly_real = np.poly1d(coef_real)
    poly_imag = np.poly1d(coef_imag)

    dpoly_real = np.polyder(poly_real)
    dpoly_imag = np.polyder(poly_imag)

    vals = poly_real(t_vals) + 1j * poly_imag(t_vals)
    normals = 1j * dpoly_real(t_vals) - dpoly_imag(t_vals)
    


    vals += 0.005 * normals * np.sin(50 * 2 * np.pi * t_vals) / np.abs(normals)

    plt.plot(vals.real, vals.imag, linewidth=0.1, color=orange)

green = '#11ff1f'
green_dark = '#0abb10'
red = '#ff1011'
black = '#222222'
grey = '#aaaaaa'
blue = '#1515dd'
orange = '#ffa500'
pink = '#ff66bb'


path_colors = {
         0: blue, 
         # 1: orange, 
         # 2: green_dark, 
         3: red, 
         # 4 :pink, 
         # 5: green
        }
deform = {1: (1.6, 1 + 0.5j, 1), 3: (3, -2 + 2j, 3)}

aux_points = transform(np.array([]))

partition = 'trefoil_superposition'
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

branch_cuts = [
        [branch_tf[3], 0.11 + 0.32j, -0.2 + 0.4j,-0.7 + 0.35j, branch_tf[0]],
        [branch_tf[1], -0.7 + 0.1j,-1],
        [branch_tf[2], -1.05 - 0.4j ,-1],
        [branch_tf[8], -0.7 - 0.2j,  -0.3 ,0],
        ]

if sys.argv[1] == "all":
    fig = plt.figure(dpi=1000)
    current_path = pathlib.Path(__file__).parent.resolve()
    for j in range(100):
        try:
            os.chdir(f'data/path_data_{j}')
        except:
            break
        for file in glob.glob("*.csv"):
            i = int(re.search('path_data_(.+?).csv', file).group(1))
            if ( 
                j == 2 and i in [] or
                j == 5 and i in [] 
            ): # [3, 6, 21, 25]                              
                continue
            data = np.loadtxt(file, delimiter=",", dtype=np.complex_) 
            if False:
                data[:, 0] = (data[:,0] - data[0,0])*np.exp((0.066 - 0.732)*1j) + data[0,0]
                data[:, 1] = (data[:,1] - data[0,1])*np.exp((0.066 - 0.732)*1j / 2) + data[0,1]
                data[:, 2] = (data[:,2] - data[0,2])*np.exp((0.066 - 0.732)*1j / 2) + data[0,2]
                print(data[230,:])
            if j == 3: 
                print(data[-30,:])
                print(branch[8])
            data = transform(data)
            x_data = data[:, 0]
            if j in deform:
                x_data = deform_path(x_data, deform[j])
            if j in path_colors:
                order = 3
                color = path_colors[j]
            else:
                order = 2
                color = grey
            if i in []:
                order =4
            plt.plot(x_data.real, x_data.imag, linewidth=0.3, color=color,
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
    for cut in branch_cuts:
        draw_cut(cut)
    branch_labels = {
             0: "$(jl)$",
             1: "$(kl)$",
             2: "$(ik)$",
             3: "$(jl)$",
             4: "$(kl)$",
             6: "$(ij)$",
             8: "$(ij)$", 
        }
    for i in branch_labels.keys():
        b = branch_tf[i]
        plt.text(b.real + 0.01, b.imag + 0.01, branch_labels[i], zorder = 5)
        plt.axis('off')
    plt.plot(aux_points.real, aux_points.imag, markersize=0.7, 
             linestyle=" ", marker="o")
    fig.tight_layout()
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
    plt.savefig(output_dir + '/network.png', dpi=fig.dpi)
    plt.savefig(output_dir + '/network.pdf', dpi=fig.dpi)
