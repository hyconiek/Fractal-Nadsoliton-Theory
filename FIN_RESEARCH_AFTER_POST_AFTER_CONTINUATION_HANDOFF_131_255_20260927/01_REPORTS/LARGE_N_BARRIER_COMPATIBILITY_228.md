# LARGE-N-BARRIER-COMPATIBILITY-228
## B4 remains compatible with the finite-N rate once polynomial prefactors are allowed, but is not yet certified as the asymptotic exponent

Date: 2026-09-26

Status:
finite-N compatibility analysis, not an Eyring-Kramers proof.

The mapped communication hierarchy gives

    B3=0.644851587327878
    B4=0.662219137127460

with

    B4-B3=0.017367549799582.

The observed effective Z3 exit rate over N=3,...,8 has the simple log-linear
fit

    r_eff
      ~
      0.5006 exp(-0.58435 N).

The exponent 0.58435 is below B4=0.66222.

That alone does NOT refute B4 because Eyring-Kramers / capacity asymptotics
generically contain N-dependent prefactors.

## 1. Fix the exponent to the mapped communication barrier

Fit only

    r_N
      =
      A N^alpha exp(-B4 N).

For the effective exit rate:

    A≈0.401251347
    alpha≈0.393198.

The log-RMSE is only

    0.023403,

and the maximum absolute log residual is

    0.033181.

Observed/predicted ratios stay within a few percent over all six finite-N
points.

For the deep-core capacity:

    alpha≈0.587513

with log-RMSE

    0.031338.

Therefore B4 is quantitatively COMPATIBLE with the current finite-N data.

## 2. But the data do not select B4 uniquely

If both beta and alpha are fitted freely,

    r_N
      =
      A N^alpha exp(-beta N),

the effective-exit fit prefers approximately

    beta≈0.569261,
    alpha≈-0.078044,

with a smaller residual.

For the deep-core rate the free fit prefers

    beta≈0.539109.

Thus N<=8 cannot distinguish:
- a pre-asymptotic exponent below B4;
from
- an eventual B4 exponent with polynomial corrections.

## 3. Current proof boundary

The barrier filtration is exact CONDITIONAL on the mapped d3/d4 saddle graph.

What is still missing for proof-grade large-N metastability:

1. a global lower bound excluding an unseen communication path below B4;
2. a potential-theory capacity asymptotic for the declared leave-one-out chain;
3. control of prefactors and finite-N correction terms.

So:

    boxed:
    B4 is a strong asymptotic candidate,
    not a certified switching exponent.
