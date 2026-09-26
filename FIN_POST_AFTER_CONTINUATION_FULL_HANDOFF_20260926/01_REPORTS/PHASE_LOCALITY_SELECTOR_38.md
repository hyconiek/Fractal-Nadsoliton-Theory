# PHASE-LOCALITY-SELECTOR-38
## Refinement alone does not select point evaluation; zero-order internal locality does

Date: 2026-09-26

Status:
- exact no-go for selection by Fourier/carrier refinement alone;
- exact uniqueness theorem after adding zero-order point locality;
- conditional physical interpretation only.

## 1. Starting family

From `PHASE-LOCKED-HIDDEN-SCALAR-37`, the most general reflection-even
D12-invariant scalar linear in the hidden sector and relative to the intrinsic
cell phase is

    F_{a,b}[h,theta]
      = a h_1(theta) + b h_2(theta),

where h_1 and h_2 are the k=1 and k=2 components of the hidden trigonometric
field evaluated at the cell's own phase.

Equivalently this is

    a Re[z1 conj(chi)] + b Re[z2 conj(chi)^2].

## 2. Exact refinement no-go

Embed the same trigonometric polynomial into any finer cyclic sampling q' that
resolves k=1,2.

Fourier restriction/prolongation preserves the two low harmonics separately:

    h_1 -> h_1,
    h_2 -> h_2.

Therefore for every constants a,b,

    F_{a,b}^{fine} = F_{a,b}^{coarse}

on the common low-mode subspace.

So exact phase/carrier refinement commutation leaves the full two-parameter
family untouched.

Hence:

    refinement covariance alone cannot force a=b.

This is an exact no-go.

## 3. Why "no new spectral scale" is still insufficient

The coefficients a,b are dimensionless relative weights between already
present k=1 and k=2 sectors. No new dimensional length or frequency scale is
needed to choose a != b.

Therefore banning new dimensional scales does not collapse the family either.

A stronger notion of locality is required.

## 4. Zero-order internal point locality

Assume now that the cell-phase coupling is zero-order local:

> If two hidden fields h and g have the same value at the cell's intrinsic
> phase theta, then the scalar coupling must be the same:
>
>     h(theta)=g(theta)  =>  F[h,theta]=F[g,theta].

For the retained hidden space,

    h(theta)=h_1(theta)+h_2(theta).

Since F is linear, this assumption means F factors through the one-dimensional
evaluation map

    ev_theta : h -> h(theta).

Therefore there exists one scalar beta such that

    F[h,theta]=beta h(theta).

Comparing separately a pure k=1 field and a pure k=2 field gives

    a=b=beta

(up to the fixed normalization convention used for the two orthonormal Fourier
planes).

Thus zero-order point locality uniquely selects point evaluation up to one
overall coupling.

## 5. Minimality of the premise

If first derivatives are allowed, another local functional exists:

    F = a h(theta) + b d_theta h(theta).

Reflection parity constrains scalar/pseudoscalar combinations but does not
generically remove all derivative freedom.

If second derivatives are allowed,

    h''(theta)

distinguishes k=1 and k=2 by eigenvalues -1 and -4 and therefore reproduces
nontrivial spectral filtering locally.

Hence the uniqueness is specifically a **zero-order locality theorem**.

It must not be restated as "all local couplings are point evaluation."

## 6. Combined bridge

With the zero-order locality premise:

    u_i = beta_h h_i(theta_i),

and with arbitrary-split metric refinement:

    delta c_ij/c_ij
      = beta_g (u_i+u_j)/2.

The two constants can be combined into one overall dimensionless coupling

    beta = beta_h beta_g

until a normalization or separate physical calibration is supplied.

Thus the complete functional form becomes

    delta c_ij/c_ij
      = beta/2 [
          h_i(theta_i)+h_j(theta_j)
        ].

No shell-specific or hidden-mode-specific coefficient remains.

## 7. Source status

Derived from:
- the existing H4 Fourier representation;
- the existing intrinsic localized phase;
- zero-order internal point locality;
- the existing metric-edge arbitrary-split refinement law.

Still supplied:
- zero-order locality itself;
- the global coupling beta;
- any kinetic term for the hidden scalar;
- physical units.

## 8. Consequence

The bridge-source problem has been compressed from:

    six shell coefficients
        ->
    two state-relative hidden-mode coefficients
        ->
    one global coefficient

once two explicit naturality principles are imposed:

1. point-local sampling in internal phase;
2. arbitrary-split metric refinement.

The remaining question is whether the resulting one-parameter bridge can arise
from one common variational/balance law rather than being a one-way response.

## 9. Next atom

`RECIPROCAL-HIDDEN-GEOMETRY-ACTION-39`

Construct the first-order interaction energy corresponding to the
refinement-rigid conductance modulation and test:

- whether varying the visible/intercell field reproduces the bridge J(u);
- whether varying u produces a reciprocal local source;
- positivity/boundedness around the declared small-coupling regime;
- what additional kinetic/storage object is still needed.

Acceptance:
one common conditional interaction functional with exact reciprocal variations,
or a no-go.
