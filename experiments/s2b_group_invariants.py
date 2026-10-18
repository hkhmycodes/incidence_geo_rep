"""S2b: bounded-degree vanishing-polynomial diagnostics for sampled word orbits.

Apply random reduced words with n seed involutions to one starting state and evaluate
a bounded-degree monomial basis at the sampled endpoints. The matrix nullity is the
dimension of the space of polynomials vanishing on that finite sample. It does not
measure group invariants and does not prove a generic orbit-closure statement.

usage: python -m experiments.s2b_group_invariants --nseeds 3 4 5 6 --D 2 3
"""
import argparse, json, random
from inc.core import rand_seeds, rand_state, sigma
from experiments.s2_orbit_dim import aff, eval_matrix


def run(nseeds, d, p, N, steps, rng):
    seeds = rand_seeds(nseeds, p, rng); st = rand_state(p, rng); pts = []
    for _ in range(N):
        cur, prev = st, None
        for i in range(steps):
            o = rng.choice([x for x in range(nseeds) if x != prev]); prev = o
            cur = sigma(seeds[o], cur, p)
        pts.append(cur)
    ch = [aff(q, p) + aff(qq, p) for q, qq in pts if aff(q, p) and aff(qq, p)]
    M, mons = eval_matrix(ch, 4, d, p)
    return len(mons) - M.rank(), len(mons)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=int, default=1000003)
    ap.add_argument("--nseeds", type=int, nargs="+", default=[3, 4, 5, 6])
    ap.add_argument("--D", type=int, nargs="+", default=[2, 3])
    ap.add_argument("--steps", type=int, default=14)
    ap.add_argument("--N", type=int, default=600)
    ap.add_argument("--trials", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="results/s2b")
    a = ap.parse_args()
    rng = random.Random(a.seed); res = {}
    print("n_seeds  deg  sample kernel dimension (min..max over trials)  / #monomials")
    for n in a.nseeds:
        for d in a.D:
            ks = [run(n, d, a.p, a.N, a.steps, rng) for _ in range(a.trials)]
            kd = [k for k, _ in ks]; nm = ks[0][1]
            res[f"{n}_{d}"] = dict(min=min(kd), max=max(kd), monomials=nm)
            print(f"{n:>7} {d:>4}   {min(kd)}..{max(kd)}   / {nm}")
    json.dump(res, open(a.out + ".json", "w"), indent=1)


if __name__ == "__main__":
    main()

