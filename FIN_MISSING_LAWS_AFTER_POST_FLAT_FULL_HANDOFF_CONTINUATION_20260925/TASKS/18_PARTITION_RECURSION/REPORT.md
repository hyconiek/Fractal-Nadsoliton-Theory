# GEOMETRY-FREE-ENERGY-04 — exact colored hierarchy partition recursion

Status: **PROOF_GRADE_RECURSION**

## Main result
For a fixed phase-count vector, the full Gibbs partition sum over unordered balanced colored hierarchies can be computed recursively without explicit tree enumeration.  Child energies are multiplied by s under one depth shift, so the recursive subpartition is evaluated at `s beta`.  Equal child count-vectors require an unordered-pair correction.

## Core formulas
For `a!=b`: \[
Z_n(eta)\supset e^{-eta\Delta(a,b)}Z_a(seta)Z_b(seta).
\]
For `a=b`: \[
Z_n(eta)\supsetrac{e^{-eta\Delta}}2\left[Z_a(seta)^2+Z_a(2seta)ight].
\]

## Evidence / reproduction
Implemented in `REPLAYS/replay_four_phase_partition.py` for scalable symmetric subalphabets.

## Caveats
Assumes perfect balanced binary hierarchies and the declared scale-weighted Ward energy.

## Next question
Use symmetry/composition compression to approach 8/12 phases at larger M.
