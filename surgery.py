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

black = '#222222'

fig = plt.figure(dpi=800, figsize=(5,5))

plt.axis([-1.6, 1.6, -1.6, 1.6])
# plt.text(1.2,1.7/(gens + 1)*0.75, "$\dots$")
# plt.text(1.2, - 1.7/(gens + 1)*0.75, "$\dots$")
ax = plt.gca()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.axis('off')
fig.tight_layout()
x = np.array([-1.5, 1.5])
y = np.array([1.5, 1.5])
plt.plot(x, -y, color=black, mew=0.5)
plt.plot(x, y, color=black, mew=0.5)
plt.plot(y, x, color=black, mew=0.5)
plt.plot(-y, x, color=black, mew=0.5)
plt.plot(x, 0*y, color=black, mew=1)
plt.plot(0*y, x, color=black, mew=1)

fig.savefig('graphics/surgery_before.pdf', dpi=fig.dpi)

fig2 = plt.figure(dpi=800, figsize=(5,4))
plt.axis([-1.6, 1.6, -1.6, 1.6])
# plt.text(1.2,1.7/(gens + 1)*0.75, "$\dots$")
# plt.text(1.2, - 1.7/(gens + 1)*0.75, "$\dots$")
ax = plt.gca()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.axis('off')
fig.tight_layout()

plt.plot(x, -y, color=black, mew=0.5)
plt.plot(x, y, color=black, mew=0.5)
plt.plot(y, x, color=black, mew=0.5)
plt.plot(-y, x, color=black, mew=0.5)

x = np.linspace(0.001/1.5, 1.5, 10007)
y = 0.001 * 1 / x

plt.plot(x, y, color=black, mew=1, zorder = -1)
plt.plot(-x, -y, color=black, mew=1, zorder = -1)

fig2.tight_layout()
fig2.savefig('graphics/surgery_after.pdf', dpi=fig.dpi)


