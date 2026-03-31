import sympy as sp
import matplotlib.pyplot as plt
import numpy as np

# define symbols
x, y, Q = sp.symbols('x y Q')

# your polynomial
F_2_3 = y - y**2 - x + x*y - 2*x*y**2 + 2*Q*x*y**2 + Q*x*y**3 - x**2*y**2 - Q**2*x*y**4 + Q*x**2*y**3

F_4_1 = (
    x
    - x**2
    - 2*Q*x*y
    + 2*x**2*y
    - Q**2*y**2
    + x**3*y**2
    + Q**2*y**3
    - Q*x**3*y**3
    + 2*Q**2*x*y**4
    - 2*Q**2*x**2*y**4
    - Q**3*x*y**5
    + Q**2*x**2*y**5
)

F_5_1 = (
    - Q**2
    + x
    - Q**3*y
    + Q*x*y
    - 2*Q*x*y**2
    + 2*Q**2*x*y**2

    - 2*Q**2*x*y**3
    + 2*Q**3*x*y**3
    + Q**2*x*y**4
    - 4*Q**3*x*y**4
    + 3*Q**4*x*y**4

    + Q**3*x*y**5
    + Q**4*x*y**5
    - 2*Q**3*x**2*y**5
    + 2*Q**4*x*y**6

    - Q**3*x**2*y**6
    - Q**4*x**2*y**6
    - Q**3*x**2*y**7
    + 4*Q**4*x**2*y**7

    - 3*Q**5*x**2*y**7
    + 2*Q**4*x**2*y**8
    - 2*Q**5*x**2*y**8

    + 2*Q**4*x**2*y**9
    - 2*Q**5*x**2*y**9
    - Q**5*x**2*y**10

    + Q**6*x**3*y**10
    - Q**5*x**2*y**11
    + Q**6*x**3*y**11
)

F_framed = F_5_1.subs(x, x*y**(-1)) * y
F_framed = sp.simplify(F_framed)
F_framed = sp.together(F_framed)
F = F_framed

def tropical_linear_forms(F, variables):
    F = sp.expand(F)   # expand but DO NOT use together()

    monomials = F.as_ordered_terms()
    tropical_forms = []

    for m in monomials:
        coeff = m.as_coeff_Mul()[0]   # coefficient
        powers = m.as_powers_dict()   # works for negative powers too

        # build linear tropical form
        lin = 0
        for v in variables:
            lin += powers.get(v, 0) * v

        tropical_forms.append(sp.simplify(lin))

    return tropical_forms

tropical_forms = tropical_linear_forms(F, [x,y,Q])

# show tropical linear forms
for f in tropical_forms:
    print(f)

x, y, t = sp.symbols('x y t')

Q_num = -1
# Example linear forms
linear_forms = [
        f.subs({Q: Q_num}) for f in tropical_forms
]

linear_forms = [sp.sympify(f) for f in linear_forms]


plt.tight_layout()

# Step 1: Compute parametric lines for f_i = f_j
lines = []  # store tuples (i, j, p, v)
for i in range(len(linear_forms)):
    for j in range(i+1, len(linear_forms)):
        f_i = linear_forms[i]
        f_j = linear_forms[j]
        # compute coefficients
        a_i, b_i, c_i = sp.linear_eq_to_matrix([f_i], [x, y])[0][0], sp.linear_eq_to_matrix([f_i], [x, y])[0][1], f_i.subs({x:0, y:0})
        a_j, b_j, c_j = sp.linear_eq_to_matrix([f_j], [x, y])[0][0], sp.linear_eq_to_matrix([f_j], [x, y])[0][1], f_j.subs({x:0, y:0})
        # difference
        a = a_i - a_j
        b = b_i - b_j
        c = c_i - c_j
        if a == 0 and b == 0:
            # parallel or identical, skip
            continue
        # parametrize line: choose t = x if b != 0 else t = y
        if b != 0:
            y_expr = (-a*x - c)/b
            slope = sp.diff(y_expr, x)
            p = sp.Matrix([0, y_expr.subs(x,0)])
            v = sp.Matrix([1, slope])
        else:
            # vertical line, parametrize by y
            x_expr = (-b*y - c)/a
            slope = sp.diff(x_expr, y)
            p = sp.Matrix([x_expr.subs(y,0), 0])
            v = sp.Matrix([slope, 1])
        lines.append((i, j, p, v))

print(f"Total parametric lines: {len(lines)}")

# Step 2: Compute all intersections algebraically
intersections = {k: [] for k in range(len(lines))}
for i in range(len(lines)):
    def is_numeric(expr):
        return expr.free_symbols == set()
    p1, v1 = lines[i][2], lines[i][3]
    for j in range(i+1, len(lines)):
        p2, v2 = lines[j][2], lines[j][3]
        # solve p1 + t1*v1 = p2 + t2*v2
        t1, t2 = sp.symbols('t1 t2')
        sol_list = sp.solve([p1[0] + t1*v1[0] - (p2[0] + t2*v2[0]),
                             p1[1] + t1*v1[1] - (p2[1] + t2*v2[1])],
                            (t1, t2), dict=True)
        for sol in sol_list:
            if t1 in sol and is_numeric(sol[t1]):
                intersections[i].append(sol[t1].evalf())
            if t2 in sol and is_numeric(sol[t2]):
                intersections[j].append(sol[t2].evalf())

# Step 3: Check segments for maximality
segments_to_plot = []
for idx, (i, j, p, v) in enumerate(lines):
    t_values = intersections[idx]
    t_values = sorted(t_values)
    # add -inf and +inf if you want full line segments
    t_values = [-10] + t_values + [10]
    
    for t0, t1 in zip(t_values[:-1], t_values[1:]):
        t_mid = (t0 + t1)/2
        pt = p + t_mid*v
        x_val, y_val = pt
        # check if f_i and f_j are maximal
        vals = [f.subs({x: x_val, y: y_val}) for f in linear_forms]
        max_val = max(vals)
        if vals[i] == max_val or vals[j] == max_val:
            segments_to_plot.append((p, v, t0, t1))

# Step 4: Plot segments
plt.figure(figsize=(8,8))
for p, v, t0, t1 in segments_to_plot:
    t_vals = np.linspace(float(t0), float(t1), 50)
    x_vals = [float((p + t*v)[0]) for t in t_vals]
    y_vals = [float((p + t*v)[1]) for t in t_vals]
    plt.plot(x_vals, y_vals, 'b')

plt.xlabel('x')
plt.ylabel('y')
# plt.title('Algebraic tropical diagram')
plt.axis('off')
plt.tight_layout()
plt.savefig("graphics/toric/toric.pdf", dpi=1000);
plt.show()
