"""S2: finite-sample orbit Hilbert ranks and plane-cubic interpolation diagnostics.

The affine Hilbert-function and plane-curve ranks are diagnostics of finite samples.
A rank pattern may suggest a low-dimensional closure but does not prove that a generic
orbit lies on that variety. A one-dimensional degree-three kernel gives a unique
interpolating cubic for the tested points; it does not prove invariance.

usage: python -m experiments.s2_orbit_dim --word AB --p 1000003 --N 2500 --D 10
"""
import argparse, itertools, json, math, random
from flint import nmod_mat
from inc.core import rand_seeds, rand_state, sigma, normalize, Language

L2I = {c: i for i, c in enumerate("ABCDEFGHIJKL")}


def orbit(word, seeds, state, p, N):
    pts = []
    while len(pts) < N:
        for ch in word:
            state = sigma(seeds[L2I[ch]], state, p)
            if state is None:
                return pts, True
        pts.append(state)
    return pts, False


def aff(pt, p):
    """(x,y) chart of a projective point; None if on the line at infinity."""
    if pt[2] % p == 0:
        return None
    inv = pow(pt[2], -1, p)
    return (pt[0] * inv % p, pt[1] * inv % p)


def monomials(nv, d):
    return [e for e in itertools.product(range(d + 1), repeat=nv) if sum(e) <= d]


def eval_matrix(points, nv, d, p):
    mons = monomials(nv, d)
    maxd = d
    pw = []
    for pt in points:
        row = []
        for x in pt:
            r, c = [1], 1
            for _ in range(maxd):
                c = c * x % p; r.append(c)
            row.append(r)
        pw.append(row)
    entries = []
    for i in range(len(points)):
        for e in mons:
            v = 1
            for j, ej in enumerate(e):
                v = v * pw[i][j][ej] % p
            entries.append(v)
    return nmod_mat(len(points), len(mons), entries, p), mons


def hilbert_function(points, nv, D, p):
    hf = []
    for d in range(D + 1):
        M, mons = eval_matrix(points, nv, d, p)
        hf.append(M.rank())
    return hf


def plane_curve_degree(points, p, dmax=8):
    for d in range(1, dmax + 1):
        M, mons = eval_matrix(points, 2, d, p)
        r = M.rank()
        if r < len(mons):
            X, nullity = M.nullspace()
            return d, nullity, mons, X
    return None, 0, None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--word", default="AB")
    ap.add_argument("--p", type=int, default=1000003)
    ap.add_argument("--N", type=int, default=2500)
    ap.add_argument("--D", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="results/s2")
    a = ap.parse_args()
    rng = random.Random(a.seed); p = a.p
    seeds = rand_seeds(max(4, 1 + max(L2I[c] for c in a.word)), p, rng)
    st = rand_state(p, rng)
    pts, dead = orbit(a.word, seeds, st, p, a.N)
    print(f"word={a.word} p={p}: {len(pts)} orbit points (degenerate stop: {dead})")
    res = {"word": a.word, "p": p, "N": len(pts)}

    # (a) Hilbert function in the 4-variable chart
    chart = []
    for (Q, Qp) in pts:
        u, v = aff(Q, p), aff(Qp, p)
        if u and v: chart.append(u + v)
    hf = hilbert_function(chart, 4, a.D, p)
    full = [math.comb(d + 4, 4) for d in range(a.D + 1)]
    diff = [hf[0]] + [hf[i] - hf[i - 1] for i in range(1, len(hf))]
    d2 = [diff[i] - diff[i - 1] for i in range(1, len(diff))]
    print("d      :", list(range(a.D + 1)))
    print("HF(d)  :", hf)
    print("full   :", [min(f, len(chart)) for f in full])
    print("diff   :", diff, "\ndiff^2 :", d2)
    verdict = "generic (no low-degree variety)"
    if hf[-1] < min(full[-1], len(chart)):
        tail = diff[-3:]
        if len(set(tail)) == 1:
            m = tail[0]; g = 1 + m * a.D - hf[-1]
            verdict = f"CURVE of degree m={m}, arithmetic genus ~ {g}"
        elif len(set(d2[-3:])) == 1:
            verdict = "SURFACE-like growth (quadratic HF): lies on a 2-dim variety"
        else:
            verdict = "proper subvariety, growth not yet stabilised (raise D and N)"
    print("verdict:", verdict)
    res.update(HF=hf, full=full, verdict=verdict)

    # (b) plane-curve test on projections (+ seeds)
    sp = [pt for s in seeds[:2] for pt in s]  # A,A',B,B' as extra points
    for label, pl in (("Q-orbit", [aff(q, p) for q, _ in pts]),
                      ("Q'-orbit", [aff(q, p) for _, q in pts]),
                      ("Q+Q' orbits", [aff(q, p) for pr in pts for q in pr]),
                      ("Q+Q' orbits + seeds A,A',B,B'", [aff(q, p) for pr in pts for q in pr] + [aff(q, p) for q in sp])):
        pl = [x for x in pl if x]
        d, nul, mons, X = plane_curve_degree(pl, p)
        msg = f"  {label:34s}: lowest plane-curve degree = {d}, kernel dim = {nul}"
        print(msg)
        res.setdefault("plane", {})[label] = dict(degree=d, kernel_dim=nul)
        if label.startswith("Q+Q' orbits +") and d == 3 and nul == 1:
            coef = {str(e): int(X[i, 0]) for i, e in enumerate(mons)}
            json.dump(dict(p=p, cubic=coef, orbit_pt=pl[0]), open(a.out + "_cubic.json", "w"))
            print("  -> unique cubic through sampled orbit points AND seeds; saved", a.out + "_cubic.json")
    json.dump(res, open(a.out + f"_{a.word}.json", "w"), indent=1)


if __name__ == "__main__":
    main()

