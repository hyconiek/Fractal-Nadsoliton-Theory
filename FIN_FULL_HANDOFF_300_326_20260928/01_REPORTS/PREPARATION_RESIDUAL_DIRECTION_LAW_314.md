# PREPARATION-RESIDUAL-DIRECTION-LAW-314
## The cross-N preparation residual is almost two-dimensional, but its amplitude is not yet a reliable standalone law

Date: 2026-09-27

Status: grouped-N finite-state analysis on N=3,...,8; frozen before any N=9 microscopic opening.

## Question

After reports 311-313 localized the dominant large-N error in the preparation/prior transfer, can the residual

\[
r_{N,\kappa}=\beta^{\rm micro}_{N,\kappa}-\beta^{\rm baseline}_{N,\kappa}
\in\mathbb R^6
\]

be described by a small number of preparation-dependent directions?

Here \(\beta_k\) are the six cosine Fourier coefficients of the 12-state localized prior at burn time t=24.

## Data provenance

- N=3,...,7: honest grouped leave-one-N-out residuals already frozen before N=8 in task 312.
- N=8: direct microscopic residual opened only after the 312 freeze.
- preparation coordinate: \(m=E[n_0/N]\).

No N=9 microscopic data were used.

## 1. Low-dimensional geometry

Centered PCA of all N=3,...,8 residuals gives

- PC1: 81.47998% of variance;
- PC1+PC2: 98.73992%;
- PC1+PC2+PC3: 99.80522%.

Thus the residual is genuinely close to a two-dimensional surface.

## 2. Model selection

Candidate laws were fitted only in grouped leave-one-N-out fashion. For each training fold the PCA basis was recomputed from the training N-groups, so the held-out group did not leak into the directions.

The best maximum LOO error is obtained by rank 3 with cubic dependence on m:

\[
\max TV=0.0169880.
\]

Rank 2 with the same cubic dependence gives

\[
\boxed{\max TV=0.0170175}
\]

and mean TV

\[
0.00824785.
\]

The absolute gain from the third direction is only

\[
2.95\times10^{-5}\;TV,
\]

so task 314 selects the rank-2 model by parsimony.

The selected form is

\[
r(N,m)=\mu+U_2 z(m),
\]

with two fixed directions and two cubic score laws

\[
z_a(m)=c_{a0}+c_{a1}m+c_{a2}m^2+c_{a3}m^3.
\]

Notably, adding explicit 1/N or N terms worsened the worst grouped-LOO error. The best simple law depends on the operational preparation coordinate m, not explicitly on N.

## 3. Important negative result

The baseline uncorrected residual already has the following honest groupwise max TV magnitudes:

- N=3: 0.4367%;
- N=4: 0.7472%;
- N=5: 0.2970%;
- N=6: 0.7417%;
- N=7: 1.9231%;
- N=8: 1.5099%.

The rank-2 law reduces the global worst case only from about 1.92% to 1.70%, while worsening several small-N folds.

Therefore the low-dimensional geometry is real, but a point correction is not yet a uniformly superior cross-N law.

## 4. Frozen use for N=9

Before opening N=9 the following were frozen:

- rank: 2;
- cubic-in-m score law;
- grouped-LOO empirical residual tube:

\[
\boxed{B_{314}=0.0170175382\;TV};
\]

- preparation-mean model selected independently by grouped LOO;
- predicted N=9 residual curve for kappa=0,...,12.

The frozen artifact hash is

`32f099cdf028d5e9c2ac610d831a10c1bac4e87f99a522f273fde0621b9c92b9`.

## Verdict

**PASS** for the statement that the preparation residual is strongly low-dimensional.

**PARTIAL / NOT PROMOTED** for a point cross-N correction law. The improvement is too modest and not Pareto-uniform across N.

Canonical use: retain the report-313 baseline as primary; use the rank-2 law only as a frozen diagnostic/optional central correction with the 1.70% TV residual tube carried explicitly.
