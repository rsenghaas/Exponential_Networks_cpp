import pandas as pd 
import matplotlib.pyplot as plt 
import matplotlib
import numpy as np
import sympy as sym

I = sym.Matrix([np.sqrt(2),0,0,0])

a,b,c,d,e,f,g,h,i,j,k,l = sym.symbols("a,b,c,d,e,f,g,h,i,j,k,l")
B1 = sym.Matrix([[0, 0, 0, 0], 
                 [1, 0, 0, 0],
                 [a, 1, 0, 0],
                 [b, c, d, 0]])

B2 = j*sym.Matrix([[0, 0, 0, 0], 
                 [e, 0, 0, 0],
                 [f, g, 0, 0],
                 [1, h, i, 0]])

rels = [(g, e), (i, d*e)]
B1 = B1.subs(rels)
B2 = B2.subs(rels)

print(B1 * B2 - B2 * B1)

print(B1 * B1.T + B2 * B2.T + I* I.T)
print(B1, B2)


