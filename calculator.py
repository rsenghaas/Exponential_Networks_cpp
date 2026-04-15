import cmath
import numpy as np


k = 16
Q_coni = 11 / 8 * np.exp(2 * cmath.pi * 1j * 0)
print(np.abs(Q_coni))
Q_4_loc = 11.089582 * cmath.exp(-0.4365 * 1j)
Q_4 = 201.9 * cmath.exp(-0.01 * 1j)
Q_04 = 237.5416 * cmath.exp(-0.008364474849797432 * 1j)
Q_2 = -((cmath.log(Q_coni)) * 2*cmath.pi * 1j)
Q_0 = 4*cmath.pi**2
n = 1
D0_factor = 1
n = 0
m = 1
print(cmath.polar(Q_2))
print(cmath.polar(m * (4 * cmath.pi**2 - Q_2) + n * Q_2))

# print(cmath.phase(n*Q_0 + Q_2))
# print(abs((n*Q_0 + Q_2) * D0_factor))
# print(cmath.phase(Q_04-Q_0))
# print(abs(Q_04 + 2*Q_0))
# print(cmath.phase(6.3 * cmath.exp(0.0084*1j) + Q_0))

# print(cmath.phase(Q_0 + Q_4_loc))

# print(x)

print(cmath.polar(10.28 * cmath.exp(-0.306 * 1j) + 3*18.001 * cmath.exp(-0.075 * 1j)))
