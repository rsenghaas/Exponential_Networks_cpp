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

plt.rc('font', size=6)
label_offset = 0.03


if True:
  x = np.array([-2, 1.5])
  gens = 7

  fig = plt.figure(dpi=800, figsize=(5,3))
  for slope in [1, -1]:
    plt.plot(x, slope*x, linewidth=1, color=black, 
             zorder=2)
  plt.text(1.55, 1.5 - label_offset, "$(ij)_{n_1}$")
  plt.text(1.55, -1.52, "$(ji)_{n_2}$")



  x = np.array([0, 1.5])
  plt.plot(x, 0*x, linewidth=1, color=green, zorder=1)
  plt.text(1.55, -label_offset , "$(ii/jj)_{k(n_1 + n_2}$")

  for g in range(1, gens + 1):
      slope = 1.7 / (g + 1)
      plt.plot(x, slope*x, linewidth=1, color=red, zorder=1) 
      label = f"$(ij)_{{ {g + 1}n_1 + {g}n_2 }}$"
      plt.text(1.55, 1.5* slope-label_offset, label)

  for g in range(1, gens + 1):
      slope = -1.7 / (g + 1)
      plt.plot(x, slope*x, linewidth=1, color=blue, zorder=1) 
      label = f"$(ji)_{{ {g}n_1 + {g + 1}n_2 }}$"
      plt.text(1.55, 1.5* slope - label_offset, label)





  plt.axis([-0.8, 1.6, -1.6, 1.6])
  plt.text(1.2,1.7/(gens + 1)*0.75, "$\dots$")
  plt.text(1.2, - 1.7/(gens + 1)*0.75, "$\dots$")
  ax = plt.gca()
  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)
  plt.axis('off')
  fig.tight_layout()
  plt.savefig('graphics/test_cone.png', dpi=fig.dpi)
  plt.savefig('graphics/test_cone.pdf', dpi=fig.dpi)

