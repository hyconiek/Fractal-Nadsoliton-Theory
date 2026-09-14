# R7P-081--090 phase/cumulant foundation

## R7P-081 amplitude convention
The imported decimal coexistence fixture is
`(|z3|,|z4|,|z5|,z6)=(0.1131879146,0.1698528641,0.2269339093,-0.3380663037)`.
It is not a raw Cartesian `theta7` vector.  With the accepted feature normalization,

`||theta||^2 = 2|z3|^2/lambda3 + 2|z4|^2/lambda4 + 2|z5|^2/lambda5 + z6^2/lambda6`.

The fixture reconstructs `||theta||=0.364555570112...`, matching the handoff's
angular coexistence radius `0.364555570128...` to the supplied decimal precision.
The corresponding field is

`h_j = sum_{k=3,4,5} Re[z_k exp(2 pi i k j/12)]/sqrt(3) + z6 (-1)^j/sqrt(12)`.

This is a fixed decimal fixture; no upstream interval amplitude certificate is inferred.

## R7P-082 exact cumulants
For the twelve-label uniform average, `E h=0` and
`K4 = kappa2/2 + kappa3/6 + kappa4/24`, with
`kappa2=E h^2`, `kappa3=E h^3`, and `kappa4=E h^4-3(E h^2)^2`.
`src/phase_cumulants.py` independently evaluates direct finite sums and the
mod-12 Fourier-resonance expansion. Random regression fixtures agree to roundoff.
The disconnected `3(Eh^2)^2` term is explicitly subtracted.

## R7P-083 cubic locks
The only phase-sensitive cubic resonances are proportional to
`s cos(2 phi3)`, `cos(phi3+phi4+phi5)`, and `cos(3 phi4)` with positive
amplitude prefactors.  Hence every global cubic phase maximum simultaneously
sets all three factors to +1.  This gives six locks per sign of `z6` and twelve
total.  The explicit D12 action reconstructs one orbit of size 12 with a
representative reflection stabilizer of order 2.

## R7P-084 quartic compatibility
Finite resonance enumeration gives exactly the five reported phase-sensitive
quartic factors.  At every cubic lock all five signed cosine factors equal +1.
Because their amplitude prefactors are nonnegative for positive amplitudes,
the cubic locks simultaneously maximize every quartic phase-sensitive term.
This is a fixed-positive-amplitude statement, not a radial or global theorem.

## R7P-085 derivatives and Sobol diagnostic
`src/phase_cumulants.py` contains analytic value/gradient/Hessian formulas for
both full log-mgf and K4 in the same three phase coordinates.  Finite-difference
regressions test the chain-rule `E[h'']` contribution. A fixed-seed scrambled
Sobol run of 65,536 points reproduces the imported C2 scales; every sample and
its three error diagnostics is stored in `results/R7P-085_sobol_samples.npz`.
Sample maxima are diagnostics, not uniform proof bounds.

## R7P-086 bounded failed uniform-C2 attempt
A deliberately assumption-light amplitude-only triangle bound gives
`||H_full-H_K4|| <= 0.163701140638...`, roughly 2738 times the smallest
quartic stationary Hessian margin. It is valid only as a demonstration that
bounds discarding resonance/phase cancellations are useless for structural
stability. A sharp uniform C2 remainder remains unresolved; the needed missing
atom is a dependency-aware interval/analytic bound preserving those cancellations.

## R7P-087 phase representation
The proof representation is the full periodic angle cube with lifted real boxes
and identification modulo `2 pi`. Boxes may cross a conventional seam because
the equations are periodic. This avoids tangent-half-angle infinity faces and
is complete for the real three-torus.

## R7P-088 route selection
The chosen route is local interval isolation plus an independent complement
cover. Numerical smallness of the remainder is not used to declare topology.
For the full model, direct local isolation near quartic roots is preferred until
a rigorous uniform remainder is available.

## R7P-089 numerical quartic catalog
Four thousand fixed-seed random solves return 60 distinct roots with zero failed
residual checks, without forcing a target count: 6 maxima, 42 saddles, 12 minima.
The minimum absolute Hessian eigenvalue is `5.978843480...e-5`.

## R7P-090 local interval isolation
All 60 quartic candidates are isolated by Krawczyk interval boxes. Radius 0.05
still certifies existence, local uniqueness, pairwise distinctness modulo the
torus, and Hessian inertia for every root. The portable certificate retains a
smaller canonical radius where appropriate; the larger-radius replay is useful
for complement-cover planning.

No statement here proves quartic complement exhaustion or exclusion of extra
full-function roots.
