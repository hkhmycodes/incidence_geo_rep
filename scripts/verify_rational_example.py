#!/usr/bin/env python3
"""Exact characteristic-zero verification of the rational cubic and a 3-step orbit."""
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp

OUT = Path(__file__).resolve().parents[1] / "results" / "rational_example.json"
x, y = sp.symbols("x y")
A2, A4 = sp.Integer(-7), sp.Integer(15)
f = x**3 - 7*x**2 + 15*x
poly_disc = sp.discriminant(f, x)
assert poly_disc == -2475


def P(xv, yv):
    return sp.Matrix([sp.Rational(xv), sp.Rational(yv), sp.Integer(1)])


def cross(a, b):
    return sp.Matrix(a).cross(sp.Matrix(b))


def normalize(v):
    v = sp.Matrix([sp.factor(c) for c in v])
    if v[2] == 0:
        raise AssertionError("orbit test expected a finite affine intersection")
    return sp.Matrix([sp.factor(v[0]/v[2]), sp.factor(v[1]/v[2]), sp.Integer(1)])


def meet(p1, p2, p3, p4):
    return normalize(cross(cross(p1, p2), cross(p3, p4)))


def sigma(seed, state):
    p, pbar = seed
    q, qbar = state
    r = meet(p, q, pbar, qbar)
    rbar = meet(p, qbar, pbar, q)
    return r, rbar


def on_curve(pt):
    xv, yv = sp.Rational(pt[0]), sp.Rational(pt[1])
    return sp.simplify(yv**2 - (xv**3 - 7*xv**2 + 15*xv)) == 0


def add_T_xy(pt):
    xv, yv = map(sp.Rational, pt)
    if xv == 0:
        raise AssertionError("the affine translation-by-T formula is undefined at x=0")
    return (sp.factor(15/xv), sp.factor(-15*yv/xv**2))


def ec_neg(pt):
    if pt is None:
        return None
    return (pt[0], -pt[1])


def ec_add(P1, P2):
    """Group addition on y^2=x^3-7x^2+15x, using O=None."""
    if P1 is None:
        return P2
    if P2 is None:
        return P1
    x1, y1 = map(sp.Rational, P1)
    x2, y2 = map(sp.Rational, P2)
    if x1 == x2 and y1 == -y2:
        return None
    if P1 != P2:
        lam = (y2-y1)/(x2-x1)
        nu = (y1*x2-y2*x1)/(x2-x1)
    else:
        if y1 == 0:
            return None
        lam = (3*x1**2 + 2*A2*x1 + A4)/(2*y1)
        nu = y1 - lam*x1
    x3 = sp.factor(lam**2 - A2 - x1 - x2)
    y3 = sp.factor(-lam*x3 - nu)
    return (x3, y3)


# Base incidence configuration.
Pxy = (sp.Rational(1), sp.Rational(3))
Qxy = (sp.Rational(3), sp.Rational(-3))
P0, Pbar = P(*Pxy), P(*add_T_xy(Pxy))
Q0, Qbar = P(*Qxy), P(*add_T_xy(Qxy))
D = meet(P0, Q0, Pbar, Qbar)
Dbar = meet(P0, Qbar, Pbar, Q0)

assert add_T_xy(Pxy) == (sp.Integer(15), sp.Integer(-45))
assert add_T_xy(Qxy) == (sp.Integer(5), sp.Integer(5))
assert D == P(12, -30)
assert Dbar == P(sp.Rational(5, 4), sp.Rational(25, 8))
for point in [P0, Pbar, Q0, Qbar, D, Dbar]:
    assert on_curve((point[0], point[1])), f"point not on cubic: {point}"
assert add_T_xy((D[0], D[1])) == (Dbar[0], Dbar[1])

# Exact 3-step orbit for g = sigma_B o sigma_A.
# On the group, 2Q=P, so P-Q=Q and q_{n+1}=q_n+P-Q=q_n+Q.
assert ec_add(Qxy, Qxy) == Pxy
delta = ec_add(Pxy, ec_neg(Qxy))
assert delta == Qxy
q = ec_add(Pxy, Pxy)  # q_0 = 2P
orbit = []
for n in range(4):
    qbar = add_T_xy(q)
    assert on_curve(q) and on_curve(qbar)
    orbit.append((q, qbar))
    if n < 3:
        q = ec_add(q, delta)

# Independently apply the homogeneous projective incidence construction.
seed_A = (P0, Pbar)
seed_B = (Q0, Qbar)
state = (P(*orbit[0][0]), P(*orbit[0][1]))
for n in range(3):
    state = sigma(seed_A, state)
    state = sigma(seed_B, state)
    expected = (P(*orbit[n+1][0]), P(*orbit[n+1][1]))
    assert state == expected, f"projective incidence mismatch at step {n}"
    # Check the translation relation q_{n+1}+b = q_n+a in the group law.
    assert ec_add(orbit[n+1][0], Qxy) == ec_add(orbit[n][0], Pxy)

def pair_json(pair):
    return [[str(v) for v in pt] for pt in pair]

result = {
    "field": "Q",
    "projective_cubic": "Y^2 Z = X^3 - 7 X^2 Z + 15 X Z^2",
    "cubic_polynomial_discriminant": int(poly_disc),
    "smooth": True,
    "two_torsion_point": [0, 0],
    "P": [1, 3],
    "P_bar": [15, -45],
    "Q": [3, -3],
    "Q_bar": [5, 5],
    "D": [12, -30],
    "D_bar": ["5/4", "25/8"],
    "all_six_construction_points_on_cubic": True,
    "pair_differences_are_T": True,
    "orbit_initial_state": "q_0=2P=(49/9,-161/27)",
    "orbit_states_and_partners": [
        {"n": n, "q": [str(v) for v in pair[0]], "q_plus_T": [str(v) for v in pair[1]]}
        for n, pair in enumerate(orbit)
    ],
    "steps_checked_by_exact_projective_cross_products": 3,
    "translation_relation_q_next_plus_Q_equals_q_plus_P": True,
}
OUT.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
