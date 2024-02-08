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
label_offset = 0.03
stretch = 0.55
linewidth = 1
l = 0.1
phase = 0.5
rot = np.matrix([[np.cos(phase), np.sin(phase)], [-np.sin(phase), np.cos(phase)]])
fig_size = (5,2)
ax_size = [3, -2., -1.6, 1.6]
x_stretch = fig_size[0] / (ax_size[0] - ax_size[1])
y_stretch = fig_size[1] / (ax_size[3] - ax_size[2])

if True:
  x = np.array([-1.8, 2.4])
  gens = [1, 2, 3, 8]

  fig = plt.figure(dpi=800, figsize=fig_size)
  for slope in [1, -1]:
    plt.plot([x[0], x[1]],[ x[0] * slope * stretch, x[1] * slope * stretch], color=black, linewidth =
             linewidth, zorder=2)
    plt.arrow(x[0], x[0] * slope * stretch, 0.25*(x[1] - x[0]), 
              0.25* (x[1] - x[0]) * slope * stretch, 
              color=black, linewidth = linewidth, zorder=2, head_width = 0.05)
    #plot_edge(A, [0,1], color=black, zorder=2, factor=0.25)
  plt.text(2.9, x[1] * stretch- label_offset, "$(+ -)_{0}$")
  plt.text(2.9, -x[1] * stretch - label_offset, "$(-+)_{-1}$")



  x = np.array([0, 2.4])
  plt.plot(x, 0*x, linewidth=1, color=green, zorder=1)
  plt.text(2.9, -label_offset , "$(\pm \pm)_{-w}$")

  for g in gens:
      slope = 1.7 / (g + 1) * stretch
      plt.plot(x, slope*x, linewidth=1, color=red, zorder=1) 
      label = f"$(+-)_{{ -{g} }}$"
      if g == 1:
        label = f"$(+-)_{{ -1 }}$"
      if g == gens[-1]:
          label =r"$(+-)_{-w}$"
      plt.text(2.9, 2.4* slope-label_offset, label)

  for g in gens:
      slope = -1.7 / (g + 1) * stretch
      plt.plot(x, slope*x, linewidth=1, color=blue, zorder=1) 
      label = f"$(-+)_{{ -{g + 1} }}$"
      if g == 1:
        label = f"$(-+)_{{ -2 }}$"
      if g == gens[-1]:
          label =r"$(-+)_{-(w + 1)}$"
      plt.text(2.9, 2.4* slope - label_offset, label)


  plt.text(-1.9, -1.1-label_offset, r"$(+-)_{0}$")
  plt.text(-1.9, 1.1-label_offset, r"$(-+)_{-1}$")
  plt.axis(ax_size)
  plt.text(2.3,1.7/gens[-1] *0.55, "$\dots$")
  plt.text(2.3, - 1.7/gens[-1]*0.55, "$\dots$")
  plt.text(2.3,1.7/(3 + 1)*0.55 + 0.1 , "$\dots$")
  plt.text(2.3, - 1.7/(3 + 1)*0.55 - 0.1, "$\dots$")

  ax = plt.gca()
  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)
  plt.axis('off')
  fig.tight_layout()
  plt.savefig('graphics/colored_cone.png', dpi=fig.dpi)
  plt.savefig('graphics/colored_cone.pdf', dpi=fig.dpi)

