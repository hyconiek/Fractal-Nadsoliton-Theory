# RECURSIVE-FIBER-SOURCE-NOGO-119
## Scale covariance, symmetry, positivity, and an unused refinement level still leave a continuum of fiber laws

Date: 2026-09-26

Status:
- exact constructive nonuniqueness theorem;
- directly addresses the proposed “unused refinement level” test.

## 1. Stronger source requirements

Suppose we demand that the fiber rate be generated from the current coarse
operator itself.

Let

    delta(L)

be the spectral gap of a positive connected coarse Laplacian L.

For any dimensionless constant

    c>0

define the target-blind rule

    boxed:
    mu(L)=c delta(L).

This rule is:
- permutation/equivariant;
- positive;
- constructed only from L;
- homogeneous under clock scaling
      L -> rho L
  because
      mu -> rho mu;
- independent of any selected vertex or target label.

Use the ST231 local graph refinement

    L^+ =
      L tensor I_2
      +mu(L) I tensor L_2.

## 2. Exact next-level recursion

The spectrum of the Kronecker sum gives

    delta(L^+)
      =
      min[
        delta(L),
        2 mu(L)
      ]

      =
      delta(L) min(1,2c).

Therefore:

### if c >= 1/2

    delta_n = delta_0

at every refinement level, and

    mu_n = c delta_0

at every level.

### if 0<c<1/2

    delta_n
      =
      (2c)^n delta_0,

and

    mu_n
      =
      c(2c)^n delta_0.

In both cases the same source rule can be applied indefinitely.

## 3. Unused-level test fails to select c

Choose, for example,

    c=1/2,
    c=1,
    c=2.

All three rules:
- are sourced only from the current coarse Laplacian;
- preserve the complete embedded coarse semigroup exactly;
- satisfy positivity;
- satisfy fiber-swap symmetry;
- respect global clock covariance;
- survive a second unused refinement level;
- survive arbitrarily many further levels.

Yet they give different fiber spectra.

Hence:

    boxed:
    recursive consistency does not select a unique fiber rate.

This is stronger than ST231's original observation that mu can simply be
chosen arbitrarily at one level.

## 4. Generalization

The spectral gap is only one example.

Any positive homogeneous symmetry-invariant functional F(L), such as suitable
combinations of:
- spectral gap;
- mean positive eigenvalue;
- trace per state;
- operator norm;
- effective-conductance averages

can generate

    mu=F(L).

Without an additional law choosing F, “source it from the coarse dynamics” is
not a unique prescription.

## 5. Verdict

The ST231 nonuniqueness survives:
- dynamic sourcing;
- scale covariance;
- symmetry;
- positivity;
- and repeated refinement.

So the missing ingredient is not another refinement replay.

It is a NEW SELECTION PRINCIPLE on the fiber-resolving sector.
