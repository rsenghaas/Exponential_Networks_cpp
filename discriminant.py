from sympy import Rational, I, symbols, discriminant, expand

x, y = symbols('x y')

H = (
    x**3*y +
    (Rational(9,8) - Rational(1,4)*I)*x**2*y**3 +
    (Rational(1,4) - Rational(1,2)*I)*x**2*y**2 +
    x**2*y -
    (Rational(1545,1280) - Rational(9,16)*I)*y**4 -
    x**3*y**2 -
    x**4 -
    x*y**2 +
    (Rational(9,8) - Rational(1,4)*I)*x*y**3
)

Delta_x = expand(discriminant(H, y))
print(Delta_x)
