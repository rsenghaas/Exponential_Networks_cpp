import numpy as np
import matplotlib.pyplot as plt

# ------------------------------
# Torus knot parameters
# ------------------------------
n, m = 2, 3     # (n, m) torus knot
N = 2000          # number of points for smooth curve

# S^3 embedding of the torus knot
phi = np.linspace(0, 2*np.pi, N+1)
z1 = np.exp(1j * n * phi) / 2
z2 = np.exp(1j * m * phi) / 2

# Coordinates in R^4
X = np.real(z1)
Y = np.imag(z1)
Z = np.real(z2)
W = np.imag(z2)

# ------------------------------
# Stereographic projection from north pole (0,0,0,1)
# ------------------------------
x = X / (1 - W)
y = Y / (1 - W)
z = Z / (1 - W)   # use z to determine over/under crossings

# ------------------------------
# Plot with gaps at undercrossings
# ------------------------------
plt.figure(figsize=(6,6))
eps = 0.02 # no additional gap; we just skip undercross segments
for i in range(N):
    # if z[i] < 0, this strand is "under", so skip plotting a tiny segment
    alpha = (x[i] + 1j*y[i])**m
    if z[i] < 0 and abs(np.imag(alpha)) < eps and np.real(alpha)> 0:
        continue
    plt.plot(x[i:i+2], y[i:i+2], 'k', lw=3)

plt.axis('equal')
plt.axis('off')
plt.savefig("knot_curve.pdf", dpi=600)
