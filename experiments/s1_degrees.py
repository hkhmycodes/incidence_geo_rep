"""S1: specialized degree growth for words acting on a random line in state space.

For each word, the six homogeneous state coordinates are polynomials in a parameter
t, with the initial pair state specialized to a random projective line over F_p.
Common factors are removed after each step. The resulting degrees are diagnostics:
they need not equal global degrees, and growth patterns alone do not prove or rule out
an invariant fibration.

usage: python -m experiments.s1_degrees --p 1000003 --L 10 --lang reduced3 --seed 1
"""
import argparse, csv, json, math, random, statistics, time
from flint import nmod_poly
from inc.core import Language, rand_seeds, cross, is_zero


def reduce_vec(v):
    g = v[0].gcd(v[1]).gcd(v[2])
    if g.degree() > 0:
        v = tuple(x // g for x in v)
    return v


def sigma_poly(seed, st):
    (P, Pp), (Q, Qp) = seed, st
    l1, l2 = reduce_vec(cross(P, Q)), reduce_vec(cross(Pp, Qp))
    l3, l4 = reduce_vec(cross(P, Qp)), reduce_vec(cross(Pp, Q))
    if any(is_zero(l) for l in (l1, l2, l3, l4)):
        return None
    R, Rp = reduce_vec(cross(l1, l2)), reduce_vec(cross(l3, l4))
    if is_zero(R) or is_zero(Rp):
        return None
    return R, Rp


def deg(st):
    return max(x.degree() for v in st for x in v)


def longest_alt_run(w):
    """longest run where w[i]==w[i+2] (ABAB...), counted in letters."""
    best = cur = 2 if len(w) >= 2 else len(w)
    for i in range(2, len(w)):
        cur = cur + 1 if w[i] == w[i - 2] else 2
        best = max(best, cur)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=int, default=1000003)
    ap.add_argument("--L", type=int, default=10)
    ap.add_argument("--lang", default="reduced3", choices=["reduced3", "cyclic"])
    ap.add_argument("--m", type=int, default=3, help="pairs for cyclic language")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="results/s1")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    p, lang = a.p, Language(a.lang, a.m)
    seeds = rand_seeds(lang.n_seeds, p, rng)
    # constant polynomials for the seeds
    seeds_poly = [tuple(tuple(nmod_poly([c], p) for c in pt) for pt in s) for s in seeds]
    lin = lambda: nmod_poly([rng.randrange(p), rng.randrange(p)], p)
    Q = tuple(lin() for _ in range(3)); Qp = tuple(lin() for _ in range(3))

    rows = []  # (word, degree)
    t0 = time.time()

    def rec(word, st):
        if word:
            rows.append(("".join("ABCDEFGHIJKL"[i] for i in word), deg(st)))
        if len(word) == a.L:
            return
        for o in lang.options(len(word), word[-1] if word else None):
            nxt = sigma_poly(seeds_poly[o], st)
            if nxt is not None:
                rec(word + [o], nxt)

    rec([], (Q, Qp))
    print(f"{len(rows)} words in {time.time()-t0:.1f}s")

    with open(a.out + "_words.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["word", "length", "degree", "alt_run"])
        for word, d in rows:
            w.writerow([word, len(word), d, longest_alt_run(word)])

    summary = {"p": p, "L": a.L, "lang": a.lang, "by_length": {}, "weak_words": []}
    print(f"{'len':>3} {'#':>6} {'min':>8} {'median':>8} {'max':>8}  min/2^k  max/2^k")
    for k in range(1, a.L + 1):
        ds = [d for word, d in rows if len(word) == k]
        if not ds: continue
        summary["by_length"][k] = dict(n=len(ds), min=min(ds), median=statistics.median(ds), max=max(ds))
        print(f"{k:>3} {len(ds):>6} {min(ds):>8} {statistics.median(ds):>8.0f} {max(ds):>8}"
              f"  {min(ds)/2**k:7.3f}  {max(ds)/2**k:7.3f}")
    # which words are 'weak' (degree far below the generic ~2^k)?
    Lk = [(w, d) for w, d in rows if len(w) == a.L]
    Lk.sort(key=lambda x: x[1])
    summary["weak_words"] = Lk[:15]
    print("\nlowest-degree words at max length (candidates for structured / weak subwords):")
    for w, d in Lk[:10]:
        print(f"  {w}  deg={d}  ({d/2**a.L:.3f} of 2^L)")
    # alternation-run feature
    by_run = {}
    for w, d in rows:
        if len(w) == a.L:
            by_run.setdefault(longest_alt_run(w), []).append(d)
    print("\nlongest alternating run -> median degree ratio to 2^L")
    summary["by_alt_run"] = {}
    for r in sorted(by_run):
        m = statistics.median(by_run[r]) / 2 ** a.L
        summary["by_alt_run"][r] = dict(n=len(by_run[r]), median_ratio=round(m, 4))
        print(f"  run {r:>2}: n={len(by_run[r]):>5}  {m:.3f}")
    # quadratic-family regression test: (AB)^k should be 3k^2+1 for reduced3
    if a.lang == "reduced3":
        fam = {len(w): d for w, d in rows if set(w) <= {"A", "B"} and all(w[i] != w[i+1] for i in range(len(w)-1))}
        summary["AB_family"] = fam
        print("\n(AB)* family degrees:", fam)
    json.dump(summary, open(a.out + "_summary.json", "w"), indent=1)


if __name__ == "__main__":
    main()

