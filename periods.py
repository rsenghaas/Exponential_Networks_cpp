#!/usr/bin/env python3
"""
Seiberg-Witten periods for curve:
    y^2 = x + 1/x - 2u

Computes electric period a(u), magnetic period a_D(u), and dyon central charge
Z = n_e a + n_m a_D.

Outputs both Cartesian and polar form so phases can be read directly.

Requires:
    pip install mpmath

Usage:
    python seiberg_witten_periods.py 2j
    python seiberg_witten_periods.py 0.5+1.2j
"""

import os
import sys
import cmath
import mpmath as mp
import numpy as np

mp.mp.dps = 10  # precision

# -----------------------------
# Period formulas
# -----------------------------

def electric_period(u):
    """
    Electric period a(u)
    Solution of Picard-Fuchs equation regular at infinity
    a(u) = sqrt(u+1) * 2F1(1/2,1/2;1; 2/(u+1))
    """
    return mp.sqrt(u) * mp.hyper([-0.25, 0.25], [1.0], 1 / u**2) * 17.772


def magnetic_period(u):
    """
    Magnetic period a_D(u)
    Second independent solution
    a_D(u) = sqrt(u-1) * 2F1(1/2,1/2;1; 2/(1-u))
    """
    # return 1j/ 4 * (u**2 - 1) * mp.hyper([0.75, 0.75], [2], 1 -  1 / u**2)
    return 1j* (u-1) * mp.hyper([0.5,0.5],[2], (1-u)/2) * 2 * mp.pi


def dyon_charge(u, n_e=1, n_m=1):
    a = electric_period(u)
    aD = magnetic_period(u)
    return n_e * a + n_m * aD


# -----------------------------
# Polar representation
# -----------------------------

def polar(z):
    r = abs(z)
    theta = cmath.phase(complex(z))
    return r, theta


# -----------------------------
# Pretty printing
# -----------------------------

def print_quantity(name, z):
    z = z
    r, theta = polar(z)

    print(f"{name}:")
    print(f"  cartesian = {z}")
    print(f"  magnitude = {r}")
    print(f"  phase(rad) = {theta}")
    print()


# -----------------------------
# Main
# -----------------------------

def main():
     
    try:
        alpha = float(sys.argv[1])
    except:
        alpha = 0
    u = 2*np.exp(2*np.pi*1j * alpha)

    print("Seiberg-Witten periods for y^2 = x + 1/x - 2u")
    print("u =", u)
    print()

    a = electric_period(u)
    aD = magnetic_period(u)
    Zdyon = dyon_charge(u, 1, 1)
    Zdyon1 = dyon_charge(u, 2, -1)
    Zdyon2 = dyon_charge(u, 3, 1)
    Zdyon3 = dyon_charge(u, 4, 1)

    print_quantity("a (electric)", a)
    print_quantity("a_D (monopole)", aD)
    print_quantity("Z_dyon (1,1)", Zdyon)
    print_quantity("Z_dyon (2,1)", Zdyon1)
    print_quantity("Z_dyon (3,1)", Zdyon2)
    print_quantity("Z_dyon (4,1)", Zdyon3)
    
    r, theta = polar(dyon_charge(u,-1,1))
    sign_str = "+"
    if u.imag < 0:
        sign_str = ""
    os.system(f"bazelnr run networks:pureYM_networks -- {theta} 80 {u.real:.5f}{sign_str}{u.imag:.5f}i")
    os.system("python spectral.py all")

    print_quantity("a (electric)", a)
    print_quantity("a_D (monopole)", aD)
    print_quantity("Z_dyon (1,1)", Zdyon2)



if __name__ == "__main__":
    main()
