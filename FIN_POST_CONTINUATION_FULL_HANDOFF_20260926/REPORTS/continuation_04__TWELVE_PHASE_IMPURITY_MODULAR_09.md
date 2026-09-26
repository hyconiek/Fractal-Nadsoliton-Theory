# TWELVE-PHASE-IMPURITY-MODULAR-09

Status: **EXHAUSTIVE FINITE M=32 CERTIFICATE FOR THE DECLARED IMPURITY CLASS**.

Start from the 495 M=16 near-equal compositions in which four of the twelve
phases are doubled.  Replicate once to M=32, so those four phases have count 4
and the remaining eight count 2.  The D12 action reduces the 495 four-subsets
to 29 orbits.  The impurity dynamic program was evaluated for all 29 orbit
representatives.

Only two impurity gaps occur:

- 417/495 compositions: Delta_imp = 12.401035516734...;
- 78/495 compositions:  Delta_imp = 12.721916246139....

For every orbit the result satisfies, to <2.5e-14,

Delta_imp = 16 E_cherry(d),
E_cherry(d)=1/2 ||mu_i-mu_(i+d)||^2,

with d=3 in the first class and d=4 in the second.

The exceptional 78 compositions have an exact modular characterization.  Let
S be the four doubled phase labels.  Then d=4 occurs iff:

1. S contains exactly one representative from each residue class modulo 4;
2. S is not one of the three cosets {r,r+3,r+6,r+9}.

There are 3^4 choices of a modulo-4 transversal and exactly 3 excluded cosets,
so the count is 3^4-3=78.  Exhaustive expansion of the 29 D12 orbits back to all
495 subsets gives zero classifier mismatches.

Equivalent Fourier wording: the exceptional class is a modulo-4 transversal,
so its indicator has vanishing k=3 and k=6 character sums; the three excluded
cosets are the special cases with maximal k=4 character amplitude.

This is a theorem only for the declared localized impurity-defect class at
finite M=32.  It is not the full second spectral gap of the Gibbs hierarchy and
not a thermodynamic statement.
