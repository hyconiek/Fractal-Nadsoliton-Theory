# NL-07 — Fisher kinetic candidate audit

Status: **CONDITIONAL_CANONICAL_STATE_METRIC__NOT_A_KINETIC_SOURCE**

At the accepted positive-amplitude branch, the pullback Fisher metric is
well-defined.  After removing the translation direction by the metric Schur
quotient, the `(beta,gamma)` block is

```
 0.159237363481  -0.076525210522
-0.076525210522   0.041894296289
```

with eigenvalues `0.004137324245` and `0.196994335525`.

This is not proportional to the FIN phase-potential Hessian, so if one *adds*
the rule `M=Fisher`, the two transverse normal modes receive a definite
dimensionless relative spectrum.  The generalized squared frequencies are
`4.208149781086` and `16.216261009204`, frequency ratio
`1.963042798820`.  This remains a conditional
kinetic construction, not a sourced clock.

The metric is stable under the already-certified large-q localized branch:
at `q=190` the maximum entrywise defect from the continuum quotient Fisher
matrix is only `3.070e-12`.

It fails exactly where a phase ceases to be observable.  Setting any one of the
three active amplitudes to zero produces one zero eigenvalue in the two-
dimensional quotient Fisher metric.  This matches the earlier warning that
amplitude-zero strata cannot support the same reduced phase dynamics.

There is an important *conditional* uniqueness route: the Čencov/Campbell
characterization makes Fisher unique up to scale when one requires the relevant
family of probability-simplex metrics to respect the appropriate Markov
morphisms/congruent embeddings.  The current FIN q-refinement has **not** been
shown to be such a morphism.  As a diagnostic, the most natural adjacent-pair
aggregation from q=384 to q=192 has L1 defect `0.795967` on the full
parity branch (`0.065561` even with parity removed).  This does not
prove that no suitable Markov refinement exists; it shows it has not been
supplied by the obvious sampling refinement.

**Disposition:** Fisher is now a serious candidate for an *internal state
metric*, but it still does not supply an independent momentum, symplectic
structure, reversible evolution law or physical time normalization.
