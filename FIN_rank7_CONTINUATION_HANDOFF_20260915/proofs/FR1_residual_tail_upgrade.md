# Frontier 1 continuation — projected-tail reduction of the unresolved 4D core

This note starts **after** the completed R7P-001--128 campaign.  It does not
rewrite the historical task ledger.  Its purpose is to shrink the first ranked
frontier: the compact residual core of the positive-orthant four-amplitude
ceiling.

Write

- `a=exp(-J3)`,
- `s=exp(-3J4/2)`,
- `t=exp(-J5/2)`,
- `y=exp(-2J6)`.

The old R7P-069 tail used a full anchored trace and therefore had to count the
slow `j=5,7` weights proportional to `t^(2-sqrt(3))`.  Those two labels have the
**same C4 feature vector**.  This continuation removes the corresponding
one-dimensional line by Courant--Fischer before taking a trace bound.

## 1. General projection lemma used three times

For any unit vector `v`, Courant--Fischer gives

`lambda2(M) <= lambda_max(P M P | v^perp) <= tr(P M P)`.

For any fixed anchor feature `F0`,

`tr(P M P) <= E ||P(F-F0)||^2`.

In the compact charts the anchor weight is `w0=1`, hence every normalized
probability satisfies `p_j<=w_j`.  If all states surviving on a compact face
lie on the removed affine line, only weights vanishing away from that face
remain in the compressed trace.

## 2. Large-J3 tail

At `a=0` only the two distinct feature values `F0` and `F4=F8` survive.  Remove
the line `F4-F0`.  The six odd labels have weight at most `a`, while labels
`2,6,10` have weight at most `a^2`.  Exact summation of squared projected
feature distances gives

`C1 = lambda3 + 2 lambda6 + lambda4 lambda5/(lambda4+lambda5)`,

`C2 = 2 lambda3 + lambda4 lambda5/(lambda4+lambda5)`.

Therefore

`lambda2(M4) <= C1 a + C2 a^2`.

Using the accepted strict spectral intervals, at `a=1/30` the upper bound is
strictly below `sigma_*` with a positive gap of about `2.84e-3`.  Thus the full
domain `exp(-J3)<=1/30` is certified, independently of `J4,J5,J6`.

## 3. Large-J4 tail

At `s=0`, labels `0,3,6,9` lie in the plane spanned by the `c3/c5` line `u` and
the alternating coordinate `e6`.  Compress to `u^perp`.  The alternating
coordinate has variance at most `lambda6/12` for every distribution.  The
other two compressed coordinates are constant on the `s=0` support; every
state changing them carries a factor `s`.

The exact sum of their squared projected distances is

`Cs=(3 lambda3 lambda4 + 2 lambda3 lambda5 + 3 lambda4 lambda5)/(lambda3+lambda5)`.

Hence

`lambda2(M4) <= lambda6/12 + Cs s`.

At `s=1/128` the strict spectral intervals leave a positive gap about
`4.17e-3`.  Therefore `exp(-3J4/2)<=1/128` is globally safe.

## 4. Strong large-J5 tail

At `t=0` the slow states `j=5,7` approach the anchor more slowly than all other
non-anchor weights, but `F5=F7`.  Remove the affine line through `F0,F5`.
Then those `t^(2-sqrt(3))` states contribute exactly zero.  The remaining
weights are bounded by groups `t`, `t^2`, `t^3`, `t^4`; for labels `1,11` we
use only `2+sqrt(3)>=3`.

Exact projected-distance sums produce positive spectral coefficients
`C1,...,C4` recorded in `results/FR1_residual_tail_upgrade.json`, with

`lambda2(M4) <= C1 t + C2 t^2 + C3 t^3 + C4 t^4`.

At `t=1/9` the interval upper bound is about `0.262573`, while
`sigma_*>0.267443`, leaving a strict gap about `4.87e-3`.  Thus

**`exp(-J5/2)<=1/9` is globally certified for arbitrary nonnegative J3,J4,J6.**

This replaces the much smaller R7P-069 tail `t<=2^-11`.

## 5. Strengthened physical q constraint

The parity partition sums also imply the exact inequality

`Z_even >= cosh(J3) Z_odd`.

After factoring the `J4` dependence, the only nontrivial term is

`A0(x)=cosh(x)+2cosh(x/2)-1-2cosh(sqrt(3)x/2)`.

Its `x^2` and `x^4` coefficients vanish; for order `x^(2n)`, `n>=3`, the sign
is that of `h_n=4^n+2-2*3^n`.  Since `h_3=12` and
`h_(n+1)=4h_n+2*3^n-6>0`, all remaining coefficients are positive.  The odd
part is

`sinh(x)-2sinh(x/2)=2sinh(x/2)(cosh(x/2)-1)>=0`.

Thus for all nonnegative fields

`q/(1-q)>=cosh(J3)`.

With `a=exp(-J3)` this is the polynomial constraint

**`q(1+a)^2 >= 1+a^2`.**

This is intended for the next dependency-aware cover.

## 6. New residual hull

Any still-unresolved physical point must now satisfy simultaneously

`a>1/30`, `s>1/128`, `t>1/9`,

in addition to lying outside the previously certified boundary, extreme-face
and local-cone regions.  This is a strict shrinkage of the post-handoff
residual core; it is not yet a complete global 4D ceiling theorem.
