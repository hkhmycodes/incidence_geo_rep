# Interpretation of legacy output files

The source snapshot's `s2_AB.json` and related files contain program-generated `verdict` fields such as `CURVE` or `SURFACE-like growth`. These are heuristic classifications from a finite-sample Hilbert-function pattern. In the manuscript, they are described only as consistency checks and are not asserted as proofs about the full orbit closure.

Likewise, the `s2b.json` nullities are dimensions of spaces of bounded-degree polynomials vanishing on finite samples. They are not dimensions of a space of rational first integrals and do not prove absence of invariants of any degree.

The one-parameter degrees in the `s1_*.json` and `s1_*_words.csv` files are not global multidegrees. All mathematical claims in the manuscript are limited accordingly.
