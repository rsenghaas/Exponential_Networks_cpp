import cmath
import numpy as np

print(np.pi - 0.4365)
Q_4_loc = 11.089582 * cmath.exp(-0.4365 * 1j)
Q_4 = 201.9 * cmath.exp(-0.01 * 1j)
Q_04 = 237.5416 * cmath.exp(-0.008364474849797432 * 1j)
Q_2 = cmath.log(9.0/8.0) * 2*cmath.pi*1j
Q_0 = 4*cmath.pi**2
n = 1
D0_factor = 1
print(Q_2)

print(cmath.phase(n*Q_0 + Q_2))
print(abs((n*Q_0 + Q_2) * D0_factor))
print(cmath.phase(Q_04-Q_0))
print(abs(Q_04 - Q_0))
# print(cmath.phase(6.3 * cmath.exp(0.0084*1j) + Q_0))

print(cmath.phase(Q_0 + Q_4_loc))
