# R7N-040 — Exact quartic phase census for the fixed negative-z6 fixture

## Theorem
For the exact decimal amplitude fixture

- r3 = 0.1131879146,
- r4 = 0.1698528641,
- r5 = 0.2269339093,
- z6 = -0.3380663037,

on the full phase torus `(R/2πZ)^3`, the quartic phase function has exactly **60** critical points.

Their negative-Hessian-index histogram is:

- index 0: 12,
- index 1: 24,
- index 2: 18,
- index 3: 6.

The fixed-sign symmetry group used for orbit accounting has order 12: translations `a` must be even because `z6 -> (-1)^a z6`, while both reflection signs are allowed. Under this subgroup the 60 roots form nine orbits: eight of size 6 and one of size 12.

## Proof objects
1. 60 disjoint radius-1e-7 local interval root boxes with certified existence, uniqueness and Hessian index.
2. 60 pairwise-disjoint radius-0.05 injectivity collars; each contains at most one quartic root and contains its corresponding local root box.
3. Complete normalized-torus complement cover: 27,272 terminal leaves excluded by a nonzero gradient component and 640 terminal leaves contained in a certified root collar; zero unresolved leaves.
4. Exact derivation of the fixed-sign subgroup plus verified action on the 60-root catalog.

## Nonconclusions
This is a fixed-amplitude quartic theorem only. It is not a full log-mgf census, not an amplitude-robustness theorem, not a full-D12 quotient inside one fixed-sign fixture, and not a statement about full X7 stationary points.
