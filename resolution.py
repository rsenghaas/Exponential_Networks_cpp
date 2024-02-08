import matplotlib.pyplot as plt
import numpy as np

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
plt.rc('font', size=5)
plt.rcParams["pdf.use14corefonts"] = True
output_dir = "graphics"
fig = plt.figure(dpi=800, figsize = (5,2))

phase = 0.5
l = 0.06
linewidth = 0.6
rot = np.matrix([[np.cos(phase), np.sin(phase)], [-np.sin(phase), np.cos(phase)]])

def plot_edge(pts, inter, color=black):
    v0 = np.array( [pts[inter[0]][0], pts[inter[0]][1]] )
    v1 = np.array( [pts[inter[1]][0], pts[inter[1]][1]] )
    v = (v1 + v0) / 2.0
    print(v)
    u_unit = (v1 - v0) / np.linalg.norm(v1 - v0)
    print(u_unit)
    epsilon = -l*(np.matmul(u_unit, rot)[0])
    print(np.linalg.norm(epsilon))
    print(epsilon)
    print(epsilon[0,0])
    plt.plot([v[0], v[0] + epsilon[0,0]], [v[1], (v[1] + epsilon[0,1])], color =
             color, linewidth=linewidth)
    epsilon = -l*(np.matmul(u_unit, np.transpose(rot)))
    print(np.linalg.norm(epsilon))
    plt.plot([v[0], v[0] + epsilon[0,0]], [v[1], v[1] + epsilon[0,1]], color =
                                         color, linewidth=linewidth)
    plt.plot(*zip(*[pts[i] for i in inter]), color=color,
             linewidth=linewidth)




A = [(-2, -0.75), (-0.5, -0.75), (-1.25, 0), (-1.25, 0.75)]
B = [(0.5, -0.75), (2, -0.75), (1.25, -0.25), (2, -0.5), 
     (1.25, 0),  (0.5, -0.25),
     (1.25, 0.25), (2, 0), (1.25, 0.5), (1.25, 0.75)]

plot_edge(A, [0,2])
plot_edge(A, [1,2])
plot_edge(A, [2,3])

plot_edge(B, [0,2])
plot_edge(B, [1,2])
plot_edge(B, [2,4])
plot_edge(B, [3,4])
plot_edge(B, [4,6])
plot_edge(B, [5,6])
plot_edge(B, [6,8])
plot_edge(B, [7,8])
plot_edge(B, [8,9])

plt.text(-2.1, -0.85, r"$2(ij)_{n_1}$")
plt.text(-0.6, -0.85, r"$3(ji)_{n_2}$")
plt.text(-1.35, 0.8, r"$(ji)_{2(n_1 + n_2) + n_2}$")
plt.rc('font', size=10)
plt.text(-0.2, 0, r"$= \lim$")
plt.rc('font', size=5)
plt.text(0.4, -0.85, r"$(ij)_{n_1}$")
plt.text(0.4, -0.35, r"$(ij)_{n_1}$")
plt.text(0.93, -0.14, r"$(ii/jj)_{n_1 + n_2}$")
plt.text(0.87, 0.36, r"$(ii/jj)_{2(n_1 + n_2)}$")
plt.text(1.3, 0.125, r"$(ji)_{(n_1 + n_2) + n_2}$")
plt.text(1.15, 0.8, r"$(ji)_{2(n_1 + n_2) + n_2}$")
plt.text(1.9, -0.6, r"$(ji)_{n_2}$")
plt.text(1.9, -0.1, r"$(ji)_{n_2}$")
plt.text(1.9, -0.85, r"$(ji)_{n_2}$")


plt.axis([-2.2, 2.2, -0.875, 0.875])
ax =plt.gca()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.axis('off')
fig.tight_layout()
plt.savefig(output_dir + '/resolution.png', dpi=fig.dpi)
plt.savefig(output_dir + '/resolution.pdf', dpi=fig.dpi)
