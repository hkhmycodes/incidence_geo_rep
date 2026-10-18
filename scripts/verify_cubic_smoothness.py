#!/usr/bin/env python3
"""Exact finite-field smoothness certificate for the interpolated cubic in the paper.

Requires SymPy. It verifies that the affine singularity condition has no solution
by taking the gcd of two resultants modulo p, and checks that the binary cubic at
infinity has nonzero discriminant. This certifies smoothness over the algebraic
closure of F_p for this particular cubic only; it does not verify dynamical
invariance under either generator.
"""
import sympy as sp

P = 1_000_003
x, y, X, Y = sp.symbols("x y X Y")
f = (
    x**3 + 967124*x**2*y + 791970*x*y**2 + 660105*y**3
    + 517859*x**2 + 530013*x*y + 631848*y**2
    + 468845*x + 643754*y + 66533
)
fx, fy = sp.diff(f, x), sp.diff(f, y)
r_fx = sp.Poly(sp.resultant(f, fx, y), x, modulus=P)
r_fy = sp.Poly(sp.resultant(f, fy, y), x, modulus=P)
g = r_fx.gcd(r_fy)
if g.degree() != 0:
    raise SystemExit(f"Affine singularity not ruled out: resultant gcd = {g.as_expr()}")

leading = X**3 + 967124*X**2*Y + 791970*X*Y**2 + 660105*Y**3
discriminant = int(sp.discriminant(leading.subs(Y, 1), X)) % P
if discriminant == 0:
    raise SystemExit("The cubic at infinity has a repeated root; smoothness not certified.")

print(f"field characteristic: {P}")
print(f"deg Res_y(f,f_x): {r_fx.degree()}")
print(f"deg Res_y(f,f_y): {r_fy.degree()}")
print(f"gcd of affine resultant polynomials: 1 (degree {g.degree()})")
print(f"discriminant of binary cubic at infinity: {discriminant} (nonzero mod {P})")
print("certificate: this projective plane cubic is smooth over the algebraic closure")
