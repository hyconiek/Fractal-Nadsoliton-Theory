# FIN: exact conditional Z12 relational-memory object

Date: 2026-10-19. Scope: **mathematical construction from supplied strict spectral weights, conditional on the already declared diffusive evolution**. This is not an intrinsic source law, a newly adopted FIN kinetic rule, or evidence for physical spacetime.

## Frozen input

Let `k1,...,k6` be the six positive strict radial weights listed in `RESULTS.json`, `W_ij=k_{min(|i-j|,12-|i-j|)}` for `i!=j`, `W_ii=0`, `s=2 sum_{r=1}^5 k_r+k6`, and `A=sI-W` on real Z12 coordinates. `A` is symmetric PSD with `A 1=0`. The finite matrix and numerical values are explicit, but the printed decimal weights are rounded source values; arithmetic tests are consequently numerical, while the symmetry proofs below hold symbolically for arbitrary radial weights.

Choose the **already available** model semigroup `dot z=-Az`. This choice is marked as a kinetic premise: `A` alone also supports `exp(-itA)` and a wave propagator, with different physical semantics.

Partition `z=(x,y)` into even coordinates `(0,2,...,10)` and odd coordinates `(1,3,...,11)`. By cyclic symmetry:

```
A = [[D,C], [C^T,D]],
D = circulant(s,-k2,-k4,-k6,-k4,-k2),
C row 0 = -(k1,k3,k5,k5,k3,k1).
```

Other rows are cyclic shifts. Each eigenvalue of `D` is strictly positive because it is a proper principal block of the connected graph Laplacian; here `D>0` numerically.

## Exact memory derivation

The two equations are `x'=-Dx-Cy`, `y'=-C^Tx-Dy`. Direct variation of constants gives

```
y(t) = exp(-tD)y(0) - integral_0^t exp(-(t-s)D) C^T x(s) ds,
x'(t)= -D x(t) - C exp(-tD)y(0)
      + integral_0^t M(t-s)x(s) ds,
M(t) = C exp(-tD) C^T.
```

This is an **exact identity**, not a Markovian closure on `x`. Retained initial data alone do not fix the `-C exp(-tD)y(0)` term. A prepared conditional distribution of `y(0)` is additional input. In Laplace space the operator is

`A_eff(z) = D-C(zI+D)^(-1)C^T` (Re z>0),

with static Schur `D-CD^{-1}C^T` at z=0. The frequency dependence is an explicit memory diagnostic.

## Positivity and the exactly dark channel

For t>=0, `M(t)` is symmetric PSD because `e^{-tD}` is PSD. Every Z6 Fourier vector `v_l(j)=exp(2 pi i l j/6)/sqrt(6)` diagonalizes `D` and `M(t)` (the latter because the factors are circulant), with

`D v_l=d_l v_l`, `M(t) v_l=m_l exp(-d_l t) v_l`, `m_l=|c_l|^2>=0`,

where `c_l` is the Z6 Fourier eigenvalue of `C`. For `l=3`, `v_3(j)=(-1)^j/sqrt(6)` and the paired odd couplings cancel **exactly**:

`C^T v_3=0` because `k1-k3+k5-k5+k3-k1=0`.

Therefore `M(t)v_3=0` for every t and the projected component `v_3^T x` follows `d/dt(v_3^Tx)=-d_3(v_3^Tx)` independent of hidden initial data. This is an exact decoupling of **one representation channel** in the chosen even/odd partition. It does NOT mean all coordinates in either half evolve autonomously, nor does it describe a physical isolated particle or a lower informational substrate.

Numeric mode data in `RESULTS.json`, in Z6 Fourier order:

| l | d_l | m_l |
|---|---:|---:|
| 0 | 1.171091020573 | 1.371454178467 |
| 1,5 | 1.526363713143 | 0.596358569832 |
| 2,4 | 1.888309181880 | 0.096882580583 |
| 3 | 1.961406861976 | **0 exactly by symmetry** |

The nonzero kernel has three distinct exponential rates and symmetry-degenerate channels. This supplies an exact conditional multiscale **memory function**, not physical fractal spacetime or a derived law of evolution.

## Status and source boundary

**CONDITIONAL_CLOSED** algebraically, given the frozen symmetric finite Z12 weights and explicit choice `z'=-Az`.

**OPEN_SOURCE** for why nature or the FIN nadsoliton must choose this semigroup, these even/odd readouts, this initial hidden preparation, this time variable, or these as physically complete state coordinates. Historical FIN Green–Schur and Mori–Zwanzig-style reduction results are predecessors; this calculation is a concrete partition-specific instantiation, not a claim that memory reduction was unknown in FIN.
