# FINITE-N-CONVENTION-ROBUSTNESS-20 — robustness of the quartic closure tensor

Status: **FULL-TENSOR NUMERICAL INVARIANCE + STRUCTURAL THEOREM CANDIDATE**.

## Tested one-parameter family

Interpolate finite-N target conventions by

    q^(theta)_{j|i}(p) proportional to
      exp[g(Ap)_j - theta (g/N) A_{ji}].

`theta=0` is the empirical-refresh rule and `theta=1` is the leave-one-out
mean-field Gibbs refresh.

The complete 210-monomial degree-four tensor of

    lim N (L_full^3-L_ME7^3) F_N(u)

was recomputed for `theta=-1,0,0.5,1,2`.  Relative to theta=0 the largest
absolute changes are at the `1e-14` level; the `g^2` tensor remains zero.

Because the h^2 coefficient can depend on theta only polynomially through
quadratic order, the repeated zero is strong evidence of exact theta
independence.  A line-by-line algebraic proof for the whole family is still to
be exported.

## Broader correction tests

Replacing the leave-one-out matrix `A` by the retained projector `P_V` or the
hidden projector `P_H` also leaves the full tensor unchanged to about `1e-14`.

A deterministic zero-mean **circulant** correction with unrelated spectral
weights leaves

- pure k5 unchanged to roundoff;
- an independently chosen mixed 7D probe unchanged to about `4e-15`.

As a negative control, a deterministic zero-mean **noncirculant** correction
changes the full tensor substantially: the maximum change in the coefficient
linear in g is about `0.08827`.

## Current structural interpretation

These checks indicate that the leading quartic closure defect is insensitive to
finite-N target corrections that preserve the translation/Fourier block
structure, but not to arbitrary label-dependent corrections that mix those
blocks.

This is currently a **theorem candidate**, not yet a proved classification of
all admissible correction matrices.

## Consequence

The k=5/general-quartic `1/N` fingerprint is not an artifact of choosing between
the empirical-refresh and leave-one-out Gibbs conventions, even though those
two chains have different exact finite-N stationary measures.
