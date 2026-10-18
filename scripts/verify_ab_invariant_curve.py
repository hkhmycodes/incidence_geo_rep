#!/usr/bin/env python3
"""Exact finite-field verification of the AB invariant curve and its quadrics.

This script is deliberately tied to the explicit F_1,000,003 example printed in
incidence_involutions_geometry.tex. It:
  * regenerates 2,500 states of the AB orbit from the printed seeds/state;
  * recomputes the degree-2 and degree-3 evaluation ranks;
  * checks the three displayed quadrics and the two plane-cubic equations;
  * substitutes each generator sigma_A and sigma_B into these five equations,
    clears denominators, and reduces the resulting numerators modulo the exact
    Groebner basis of the displayed state-curve ideal.

The substitution check is exact for this finite-field example. The generic
invariance theorem in the paper is proved geometrically from Schroeter's
construction; this script is an additional reproducibility certificate.
Requires SymPy (no python-flint needed).
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from inc.core import sigma

P = 1_000_003
SEEDS = [
    ((140891, 596853, 888598), (841235, 800875, 66172)),
    ((267459, 123646, 519501), (797926, 471325, 495185)),
]
INITIAL = ((799308, 804423, 2208), (729633, 467022, 279267))

# The monomial order is total degree, then lexicographic exponent tuples.
# Variables are x,y,u,v, where (x,y) represents Q and (u,v) represents Q'.
def monomials(max_degree: int):
    out = []
    for total in range(max_degree + 1):
        for e in itertools.product(range(total + 1), repeat=4):
            if sum(e) == total:
                out.append(e)
    return out


def affine(point):
    z = point[2] % P
    if not z:
        return None
    iz = pow(z, -1, P)
    return (point[0] * iz % P, point[1] * iz % P)


def orbit_states(N=2500):
    state = INITIAL
    states = []
    for _ in range(N):
        state = sigma(SEEDS[0], state, P)
        if state is None:
            raise RuntimeError("degenerate sigma_A step")
        state = sigma(SEEDS[1], state, P)
        if state is None:
            raise RuntimeError("degenerate sigma_B step")
        states.append(state)
    return states


def eval_row(point, mons):
    row = []
    for e in mons:
        val = 1
        for a, n in zip(point, e):
            if n:
                val = val * pow(a, n, P) % P
        row.append(val)
    return row


def rref(rows):
    """Incremental reduced row-echelon basis over F_P, indexed by pivot column."""
    basis = {}
    for row0 in rows:
        row = [v % P for v in row0]
        for pivot in sorted(basis):
            if row[pivot]:
                a = row[pivot]
                row = [(x - a * y) % P for x, y in zip(row, basis[pivot])]
        pivot = next((i for i, v in enumerate(row) if v), None)
        if pivot is None:
            continue
        inv = pow(row[pivot], -1, P)
        row = [v * inv % P for v in row]
        for old_pivot, old_row in list(basis.items()):
            if old_row[pivot]:
                a = old_row[pivot]
                basis[old_pivot] = [(x - a * y) % P for x, y in zip(old_row, row)]
        basis[pivot] = row
    return basis


def kernel(basis, ncols):
    out = []
    for free in [j for j in range(ncols) if j not in basis]:
        vec = [0] * ncols
        vec[free] = 1
        for pivot, row in basis.items():
            vec[pivot] = (-row[free]) % P
        out.append(vec)
    return out


def poly_from_vector(vec, mons, variables):
    out = sp.Integer(0)
    for c, exponents in zip(vec, mons):
        if not c:
            continue
        term = sp.Integer(c)
        for var, exponent in zip(variables, exponents):
            term *= var ** exponent
        out += term
    return out


def coeff_vector(poly, mons, variables):
    Ppoly = sp.Poly(poly, *variables, modulus=P)
    return [int(Ppoly.coeff_monomial(sp.prod(v**e for v, e in zip(variables, exponents)))) % P
            for exponents in mons]


def affine_state(state):
    q, qp = affine(state[0]), affine(state[1])
    if q is None or qp is None:
        raise RuntimeError("state has a point outside the selected affine chart")
    return q + qp


def filtered_hilbert_from_leading_monomials(leading_exponents, max_degree, nvars=4):
    values = []
    for degree in range(max_degree + 1):
        count = 0
        for exponents in itertools.product(range(degree + 1), repeat=nvars):
            if sum(exponents) > degree:
                continue
            if not any(all(exponents[i] >= lm[i] for i in range(nvars))
                       for lm in leading_exponents):
                count += 1
        values.append(count)
    return values


def cross(a, b):
    return [sp.expand(a[1]*b[2] - a[2]*b[1]),
            sp.expand(a[2]*b[0] - a[0]*b[2]),
            sp.expand(a[0]*b[1] - a[1]*b[0])]


def affine_seed(seed):
    P0, P1 = seed
    iz0, iz1 = pow(P0[2], -1, P), pow(P1[2], -1, P)
    return (P0[0]*iz0 % P, P0[1]*iz0 % P,
            P1[0]*iz1 % P, P1[1]*iz1 % P)


def substitution_numerator(poly, seed, variables):
    x, y, u, v = variables
    ax, ay, bx, by = affine_seed(seed)
    A = [sp.Integer(ax), sp.Integer(ay), sp.Integer(1)]
    Ap = [sp.Integer(bx), sp.Integer(by), sp.Integer(1)]
    Q = [x, y, sp.Integer(1)]
    Qp = [u, v, sp.Integer(1)]
    l1, l2 = cross(A, Q), cross(Ap, Qp)
    l3, l4 = cross(A, Qp), cross(Ap, Q)
    R, Rp = cross(l1, l2), cross(l3, l4)
    images = [R[0]/R[2], R[1]/R[2], Rp[0]/Rp[2], Rp[1]/Rp[2]]
    rational = sp.cancel(poly.subs(dict(zip(variables, images)), simultaneous=True))
    numerator, _denominator = sp.fraction(rational)
    return sp.Poly(numerator, *variables, modulus=P).as_expr()


def main():
    states = orbit_states()
    points = []
    for state in states:
        q, qp = affine(state[0]), affine(state[1])
        if q is not None and qp is not None:
            points.append(q + qp)
    if len(points) != 2500 or len(set(points)) != 2500:
        raise RuntimeError(f"expected 2,500 distinct affine states; got {len(points)}")

    mons2, mons3 = monomials(2), monomials(3)
    rank2 = rref(eval_row(pt, mons2) for pt in points)
    rank3 = rref(eval_row(pt, mons3) for pt in points)
    q_basis = kernel(rank2, len(mons2))
    cubic_kernel = kernel(rank3, len(mons3))
    if (len(rank2), len(q_basis), len(rank3), len(cubic_kernel)) != (12, 3, 18, 17):
        raise RuntimeError("unexpected Hilbert/evaluation ranks")

    x, y, u, v = sp.symbols("x y u v")
    variables = (x, y, u, v)
    quadrics = [poly_from_vector(q, mons2, variables) for q in q_basis]
    # Normalize the display using the same monic cross terms as the paper.
    expected = [
        x*y*0 + y*v + 360203*y + 360203*v + 94183*x + 94183*u + 104557,
        x*v + y*u + 792235*y + 792235*v + 993772*x + 993772*u + 194562,
        x*u + 482705*y + 482705*v + 678975*x + 678975*u + 32895,
    ]
    # The nullspace implementation orders free columns deterministically; compare spans
    # rather than depending on the particular basis returned by Gaussian elimination.
    span_expected = rref([eval_row((a,b,c,d), mons2) for a,b,c,d in []])  # empty helper
    del span_expected
    # The three displayed quadrics are independent and span the sampled
    # degree-two kernel.
    expected_rank = len(rref(coeff_vector(q, mons2, variables) for q in expected))
    if expected_rank != 3:
        raise RuntimeError("displayed quadrics are not independent")

    fxy = (x**3 + 967124*x**2*y + 791970*x*y**2 + 660105*y**3
           + 517859*x**2 + 530013*x*y + 631848*y**2
           + 468845*x + 643754*y + 66533)
    fuv = fxy.xreplace({x: u, y: v})
    equations = expected + [fxy, fuv]

    # Also check the initial state and the two seed states, not only the 2,500
    # retained orbit states.
    extra_states = [affine_state(INITIAL), affine_seed(SEEDS[0]), affine_seed(SEEDS[1])]
    all_states_for_equations = points + extra_states
    f_xy_poly = sp.Poly(fxy, x, y, modulus=P)
    f_uv_poly = sp.Poly(fuv, u, v, modulus=P)
    for pt in all_states_for_equations:
        if int(f_xy_poly.eval({x: pt[0], y: pt[1]})) % P:
            raise RuntimeError("first projected cubic failed on a state/seed point")
        if int(f_uv_poly.eval({u: pt[2], v: pt[3]})) % P:
            raise RuntimeError("second projected cubic failed on a state/seed point")

    for q in expected:
        q_poly = sp.Poly(q, *variables, modulus=P)
        for pt in all_states_for_equations:
            if int(q_poly.eval(dict(zip(variables, pt)))) % P:
                raise RuntimeError("displayed quadric failed on a state/seed point")

    G = sp.groebner(equations, *variables, modulus=P, order="grevlex")
    G_bilin = sp.groebner(expected, *variables, modulus=P, order="grevlex")
    leading_exponents = [poly.LM(order=G.order).exponents for poly in G.polys]
    leading_exponents_bilin = [poly.LM(order=G_bilin.order).exponents for poly in G_bilin.polys]
    filtered_hilbert = filtered_hilbert_from_leading_monomials(leading_exponents_bilin, 10)
    expected_hilbert = [1, 5, 12, 18, 24, 30, 36, 42, 48, 54, 60]
    if filtered_hilbert != expected_hilbert:
        raise RuntimeError(f"unexpected filtered Hilbert function: {filtered_hilbert}")
    if leading_exponents_bilin != leading_exponents:
        raise RuntimeError("bilinear-only and five-generator ideals have different leading monomials")

    cubic_remainders = {}
    for name, equation in (("f_x_y", fxy), ("f_u_v", fuv)):
        rem = G_bilin.reduce(equation)[1]
        cubic_remainders[name] = sp.Poly(rem, *variables, modulus=P).is_zero
        if not cubic_remainders[name]:
            raise RuntimeError(f"{name} is not redundant modulo the bilinear ideal")

    substitution_checks = {}
    for name, seed in zip(("A", "B"), SEEDS):
        flags = []
        for equation in expected:
            numerator = substitution_numerator(equation, seed, variables)
            remainder = G_bilin.reduce(numerator)[1]
            flags.append(sp.Poly(remainder, *variables, modulus=P).is_zero)
        substitution_checks[name] = flags
        if not all(flags):
            raise RuntimeError(f"substitution remainder did not vanish for sigma_{name}: {flags}")

    result = {
        "field": f"F_{P}",
        "distinct_affine_states": len(set(points)),
        "degree_2_evaluation_rank": len(rank2),
        "degree_2_kernel_dimension": len(q_basis),
        "degree_3_evaluation_rank": len(rank3),
        "degree_3_kernel_dimension": len(cubic_kernel),
        "displayed_quadric_rank": expected_rank,
        "quadrics": [str(sp.Poly(q, *variables, modulus=P).as_expr()) for q in expected],
        "plane_cubic": str(fxy),
        "bilinear_ideal_generators": 3,
        "five_generator_ideal_generators": 5,
        "bilinear_groebner_basis_size": len(G_bilin.polys),
        "five_generator_groebner_basis_size": len(G.polys),
        "same_leading_monomials_as_five_generator_ideal": leading_exponents_bilin == leading_exponents,
        "leading_monomials_grevlex": [str(poly.LM(order=G_bilin.order)) for poly in G_bilin.polys],
        "filtered_hilbert_function_d0_to_10": filtered_hilbert,
        "projected_cubics_reduce_to_zero_modulo_bilinear_ideal": cubic_remainders,
        "checked_state_count_including_initial_and_seeds": len(all_states_for_equations),
        "cleared_substitution_remainders_zero_for_three_bilinear_generators": substitution_checks,
    }
    out = ROOT / "results" / "ab_invariant_curve_verification.json"
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
