#!/usr/bin/env python3
"""Recompute the four sampled cubic coefficient vectors printed in the paper's appendix.

Uses the pure-Python incidence map and Gaussian elimination modulo p, so this exact
fixed-input check does not require python-flint. It is intended to reproduce the
finite-sample interpolation, not to prove dynamical invariance.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from inc.core import sigma

P = 1_000_003
SEEDS = [
    ((140891, 596853, 888598), (841235, 800875, 66172)),
    ((267459, 123646, 519501), (797926, 471325, 495185)),
]
STATES = [
    (((799308, 804423, 2208), (729633, 467022, 279267)),
     (66533, 643754, 631848, 660105, 468845, 530013, 791970, 517859, 967124, 1)),
    (((905035, 993869, 890298), (59298, 96033, 88994)),
     (500557, 460022, 54257, 213337, 826484, 978972, 89414, 109945, 72725, 1)),
    (((249523, 621429, 570665), (136758, 387926, 960437)),
     (610463, 52030, 550045, 45730, 63475, 569353, 15186, 369450, 109185, 1)),
    (((247514, 318031, 108177), (756250, 415297, 502140)),
     (391879, 636114, 972060, 631647, 147514, 811934, 935913, 818900, 486500, 1)),
]
MONOMIALS = [(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(1,2),(2,0),(2,1),(3,0)]


def affine(point):
    z = point[2] % P
    if z == 0:
        return None
    iz = pow(z, -1, P)
    return (point[0] * iz % P, point[1] * iz % P)


def orbit_ab(state, N=2500):
    current = state
    result = []
    for _ in range(N):
        current = sigma(SEEDS[0], current, P)
        if current is None:
            raise RuntimeError("degenerate sigma_A step")
        current = sigma(SEEDS[1], current, P)
        if current is None:
            raise RuntimeError("degenerate sigma_B step")
        result.append(current)
    return result


def rref_kernel(rows, ncols):
    """Return pivot-column dictionary in RREF over F_P using incremental elimination."""
    basis = {}
    for row0 in rows:
        row = [v % P for v in row0]
        for col in sorted(basis):
            a = row[col]
            if a:
                pivot = basis[col]
                row = [(v - a*w) % P for v, w in zip(row, pivot)]
        pivot_col = next((j for j, v in enumerate(row) if v), None)
        if pivot_col is None:
            continue
        inv = pow(row[pivot_col], -1, P)
        row = [(v * inv) % P for v in row]
        for old_col, old in list(basis.items()):
            a = old[pivot_col]
            if a:
                basis[old_col] = [(v - a*w) % P for v, w in zip(old, row)]
        basis[pivot_col] = row
    return basis


def interpolation(points, degree):
    mons = [(a,b) for a in range(degree+1) for b in range(degree+1) if a+b <= degree]
    rows = []
    for x, y in points:
        xp = [1]
        yp = [1]
        for _ in range(degree):
            xp.append(xp[-1]*x % P)
            yp.append(yp[-1]*y % P)
        rows.append([(xp[a] * yp[b]) % P for a,b in mons])
    basis = rref_kernel(rows, len(mons))
    free_cols = [j for j in range(len(mons)) if j not in basis]
    if not free_cols:
        return mons, 0, None
    # A one-dimensional kernel has one free column. Construct its generator.
    if len(free_cols) != 1:
        return mons, len(free_cols), None
    free = free_cols[0]
    coeffs = [0]*len(mons)
    coeffs[free] = 1
    for col, pivot_row in basis.items():
        coeffs[col] = (-pivot_row[free]) % P
    return mons, 1, tuple(coeffs)


def fit(initial_state):
    states = orbit_ab(initial_state)
    projected = []
    for state in states:
        projected.extend([affine(state[0]), affine(state[1])])
    projected.extend(affine(q) for pair in SEEDS for q in pair)
    projected = [pt for pt in projected if pt is not None]
    # No equation of degree <=2, and a unique one-dimensional kernel in degree 3.
    for degree in (1, 2):
        _, nul, _ = interpolation(projected, degree)
        if nul:
            raise RuntimeError(f"unexpected degree-{degree} equation: nullity={nul}")
    mons, nullity, coeffs = interpolation(projected, 3)
    if mons != MONOMIALS or nullity != 1 or coeffs is None:
        raise RuntimeError(f"expected unique cubic kernel; got {nullity}")
    x3_index = MONOMIALS.index((3,0))
    inv = pow(coeffs[x3_index], -1, P)
    return tuple(c*inv % P for c in coeffs)


if __name__ == "__main__":
    checks = []
    for i, (state, expected) in enumerate(STATES, 1):
        actual = fit(state)
        match = actual == expected
        checks.append({"run": i, "coefficients": list(actual), "matches_paper": match})
        print(f"run {i}: {actual} {'OK' if match else 'MISMATCH'}")
        if not match:
            raise SystemExit(1)
    out = ROOT / "results" / "app_cubic_reproduction.json"
    out.write_text(json.dumps(checks, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out.relative_to(ROOT)}")
