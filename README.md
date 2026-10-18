# Reproducibility archive: Paired-Point Incidence Dynamics

This archive accompanies **Bilinear Geometry of Paired-Point Incidence Dynamics: Elliptic Fibres, Translations, and Doubling**.

## Scope

The archive contains the exact rational example, its three-step incidence orbit, and the finite-field certificate used in the manuscript. The main algebraic-geometric results are proved in the paper; the finite-field material is a reproducible example and is not used in place of a proof.

## Environment

Python 3.10 or newer is recommended. Install dependencies from this directory:

```sh
python -m pip install -r requirements.txt
```

## Reproduce the principal calculations

Run from this directory:

```sh
# Characteristic-zero cubic example and exact three-step paired-point orbit
python scripts/verify_rational_example.py

# Exact finite-field invariant-curve, bilinear-ideal, and Hilbert-function certificate
python scripts/verify_ab_invariant_curve.py

# Exact smoothness certificate for the finite-field plane cubic
python scripts/verify_cubic_smoothness.py

# Recreate the four exploratory cubic coefficient vectors in the appendix
python scripts/reproduce_app_cubics.py
```

The rational script checks the six base construction points, the nonzero cubic discriminant, the partner formula, and three successive applications of the projective cross-product incidence maps. It verifies the elliptic translation identity at each step.

The finite-field script uses \(\mathbb F_{1,000,003}\). It regenerates 2,500 distinct orbit states; recomputes degree-two and degree-three evaluation ranks; calculates Gröbner bases for the three bilinear equations and the five-generator ideal; verifies that both projected cubic equations reduce to zero modulo the bilinear ideal; compares leading monomials and filtered Hilbert functions; and verifies invariance of the three-generator bilinear ideal under both rational maps by exact reduction.

## Contents

- `incidence_involutions_geometry.tex`: the LaTeX manuscript source; the compiled PDF is distributed separately.
- `incidence_rational_example.pdf`: the figure used by the manuscript.
- `inc/core.py`: projective cross-product incidence map and fixed-seed involutions.
- `scripts/verify_rational_example.py`: exact rational example and orbit verification.
- `scripts/verify_ab_invariant_curve.py`: finite-field invariant-curve and bilinear-ideal certificate.
- `scripts/verify_cubic_smoothness.py`: exact smoothness check for the finite-field plane cubic.
- `scripts/reproduce_app_cubics.py`: reproduction of the exploratory cubic coefficient vectors.
- `experiments/`: exploratory orbit-rank and specialized-degree experiments retained for provenance.
- `results/`: machine-readable outputs, including the rational orbit and finite-field checks.

## Persistent archive record

This ZIP is a local supplementary archive and has not yet been deposited in a public repository; it therefore has no DOI. Before journal submission, deposit this exact archive and the manuscript source in a repository such as Zenodo, then replace this note with the permanent DOI and citation metadata.
