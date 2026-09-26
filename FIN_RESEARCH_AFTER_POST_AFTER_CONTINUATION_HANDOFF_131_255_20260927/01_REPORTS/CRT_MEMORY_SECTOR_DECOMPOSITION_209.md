# CRT-MEMORY-SECTOR-DECOMPOSITION-209
## The three-sector Mori-Zwanzig memory has only k=0,4,8 symmetry content

Date: 2026-09-26

Status:
- exact representation-selection theorem from Z12 equivariance;
- exact finite-state N=6 numerical decomposition.

Projection:
three localized Z3 sectors plus one residual transition class.

## 1. Resolved representation

Under one-label rotation j->j+1:
- the three localized sectors cycle among themselves;
- the residual class is invariant.

Therefore the resolved four-dimensional representation decomposes as

    boxed:
    k=0
      plus
    k=0
      plus
    k=4
      plus
    k=8.

The k=4 and k=8 components are the two real degrees of freedom of the
nontrivial Z3 quotient.

## 2. Equivariance theorem

Let:
- S be the self-adjoint reversible microscopic generator;
- P the equilibrium orthogonal projection onto the resolved subspace;
- Q=I-P.

Both S and P commute with Z12 rotations.

Therefore the memory source

    C=QSP

is an equivariant map.

An equivariant map cannot create irreducible representation sectors absent
from its domain.

Hence:

    boxed:
    image(C)
      contains only
    k=0,4,8.

The hidden memory propagator

    QSQ

also preserves these sectors.

Thus the exact Mori-Zwanzig kernel

    K(t)=C^* exp(tQSQ) C

cannot receive contributions from k=1,2 hidden Fourier sectors.

This is an all-N symmetry statement for the equivariant projection.

## 3. Exact N=6 decomposition

The Frobenius power of C splits as:

    k=0:
      0.806268596850

    k=4:
      0.096865701575

    k=8:
      0.096865701575.

All other k sectors are zero to numerical precision.

The M0 sector reconstruction residual is approximately

    2.3e-14.

## 4. Slow quotient versus residual mode

The large k=0 contribution belongs to the invariant combination:
- total localized-sector occupancy;
- residual/transition occupancy.

It affects the fast leakage/residual mode.

The actual slow three-sector relaxation lives entirely in the conjugate pair

    k=4,8.

So raw memory norm is not a good measure of the memory relevant to the slow
clock.

One must resolve the memory by symmetry sector.
