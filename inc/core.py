"""Core incidence operator Phi and the involutions sigma_s it induces.

Phi((P,P'),(Q,Q')) = (PQ ^ P'Q',  PQ' ^ P'Q)      (^ = intersection of lines)

With homogeneous coordinates, the line through U,V is U x V and the meet of
lines l,m is l x m, so Phi is six cross products.  Everything works over ints
mod p (pass p) and over any ring with +,-,* (pass p=None), e.g. flint nmod_poly
for symbolic degree tracking.

For a fixed first argument s=(P,P'), sigma_s(Q,Q') = Phi(s,(Q,Q')) is an
involution (generically): the structural fact the experiments rely on.
"""
from __future__ import annotations
import random


# ----------------------------------------------------------------- vectors
def cross(u, v, p=None):
    a = u[1] * v[2] - u[2] * v[1]
    b = u[2] * v[0] - u[0] * v[2]
    c = u[0] * v[1] - u[1] * v[0]
    if p is None:
        return (a, b, c)
    return (a % p, b % p, c % p)


def is_zero(v):
    return not (v[0] or v[1] or v[2])


def normalize(v, p):
    """Scale so the first nonzero coordinate is 1 (canonical projective rep)."""
    for x in v:
        if x % p:
            inv = pow(x, -1, p)
            return tuple((y * inv) % p for y in v)
    raise ValueError("zero vector")


def proportional(u, v, p):
    return is_zero(cross(u, v, p))


# ------------------------------------------------------------------- Phi
def phi_raw(P, Pp, Q, Qp, p=None):
    """Raw (unnormalised) Phi. Returns None on a degenerate configuration."""
    l1, l2 = cross(P, Q, p), cross(Pp, Qp, p)
    l3, l4 = cross(P, Qp, p), cross(Pp, Q, p)
    if is_zero(l1) or is_zero(l2) or is_zero(l3) or is_zero(l4):
        return None
    R, Rp = cross(l1, l2, p), cross(l3, l4, p)
    if is_zero(R) or is_zero(Rp):
        return None
    return R, Rp


def sigma(seed, state, p):
    """sigma_seed(state) = Phi(seed, state);  None if degenerate."""
    return phi_raw(seed[0], seed[1], state[0], state[1], p)


def canon_state(state, p):
    return (normalize(state[0], p), normalize(state[1], p))


# ------------------------------------------------------------- sampling
def rand_point(p, rng=random):
    while True:
        v = (rng.randrange(p), rng.randrange(p), rng.randrange(p))
        if any(v):
            return v


def rand_seeds(n, p, rng=random):
    return [(rand_point(p, rng), rand_point(p, rng)) for _ in range(n)]


def rand_state(p, rng=random):
    return (rand_point(p, rng), rand_point(p, rng))


# ------------------------------------------------------------- languages
class Language:
    """Which seed may be applied at step i, given the previous letter.

    'reduced3' : 3 seeds, no immediate repeat (a repeat cancels: sigma^2 = id).
                 3 choices at step 0, then 2.
    'cyclic'   : 2m seeds grouped in m pairs; step i chooses one of the two seeds
                 of pair (i mod m).  Consecutive steps use different pairs, so the
                 word is always reduced, there are exactly 2 choices per step, and
                 the choice set is public (needed by the zero-knowledge circuit).
                 m=2 is the 'bipartite' 4-seed language; the experiments in
                 s2b_group_invariants say use m >= 3 (>= 5 seeds) so the group has
                 no low-degree invariant.
    """

    def __init__(self, kind="reduced3", m=3):
        assert kind in ("reduced3", "cyclic")
        self.kind, self.m = kind, m
        self.n_seeds = 3 if kind == "reduced3" else 2 * m

    def options(self, step, prev, binary=False):
        if self.kind == "cyclic":
            j = step % self.m
            return [2 * j, 2 * j + 1]
        if prev is None:
            return [0, 1] if binary else [0, 1, 2]
        return [i for i in range(3) if i != prev]

    def word_from_bits(self, bits):
        out, prev = [], None
        for i, b in enumerate(bits):
            prev = self.options(i, prev, binary=True)[b]
            out.append(prev)
        return out

    def iter_words(self, L):
        """All words of length exactly L (depth-first)."""
        def rec(w):
            if len(w) == L:
                yield tuple(w)
                return
            for o in self.options(len(w), w[-1] if w else None):
                w.append(o)
                yield from rec(w)
                w.pop()
        yield from rec([])


def apply_word(word, seeds, state, p):
    for s in word:
        state = sigma(seeds[s], state, p)
        if state is None:
            return None
    return state

