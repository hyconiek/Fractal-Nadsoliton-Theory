# FINITE-N-HEAT-BATH-CORRECTION-62
## The exact Gibbs sampler gives the first 1/N correction to the self-consistent heat-bath drift

Date: 2026-09-26

Status:
- exact first-order Taylor formula for the finite-N count drift;
- identifies the finite-size deterministic bias separately from stochastic
  fluctuations;
- no uniform theorem at the fold singularity is claimed.

## 1. Exact source-conditioned target

For current empirical p and old label i,

    q^{(i)}
      =
      softmax[
        h-(g/N)A7 e_i
      ],

where

    h=gA7 p,
    q=softmax(h).

Let

    S_q=diag(q)-q q^T.

Taylor expansion gives

    q^{(i)}
      =
      q
      -(g/N)S_q A7 e_i
      +O(N^-2).

## 2. Averaged finite-N drift

Average over the old-label distribution p:

    sum_i p_i q^{(i)}
      =
      q
      -(g/N)S_q A7 p
      +O(N^-2).

Therefore

    boxed:
    dot p_N
      =
      q(p)-p
      -(g/N)S_q A7 p
      +O(N^-2).

The first finite-N deterministic correction is thus fully fixed by the same
rank-seven operator and categorical covariance.

## 3. Meaning

There are two distinct finite-N effects:

1. deterministic bias:
       O(1/N);

2. stochastic empirical noise:
       O(N^-1/2).

Away from singular points the stochastic fluctuation is therefore parametrically
larger than the deterministic finite-N drift correction.

## 4. Near the simple fold

The soft deterministic restoring rate scales as

    lambda_soft ~ delta^(1/2).

The O(1/N) bias is amplified by the inverse soft rate, giving a stationary
branch displacement of rough order

    N^-1 delta^-1/2.

The branch separation is

    delta^(1/2).

Their ratio is therefore

    ~1/(N delta).

So deterministic finite-N branch displacement becomes comparable only near

    delta~N^-1.

By contrast, the stochastic saddle-node crossover from report 60 is

    delta~N^-2/3.

For large N,

    N^-1 << N^-2/3.

Thus stochastic broadening/barrier loss is expected before the deterministic
1/N shift becomes comparable with the branch separation.

This is a scaling statement, not a uniform remainder theorem.

## 5. Consistency with fluctuation-dissipation

The equilibrium soft variance scales as

    Var_soft ~ 1/[N lambda_static],

and the fold static curvature behaves as

    lambda_static~delta^(1/2).

Hence

    std_soft
      ~ N^-1/2 delta^-1/4.

Comparing with branch separation delta^(1/2),

    std_soft / separation
      ~ N^-1/2 delta^-3/4.

This becomes O(1) precisely when

    boxed:
    delta~N^-2/3.

So the barrier argument and the fluctuation-width argument give the same
finite-size exponent.

## 6. Next research atom

`FOLD-FINITE-N-METASTABILITY-63`:

Use the exact count Gibbs sampler of report 61 to determine whether its
quasipotential is exactly the R7P-111 V_g and whether a controlled local
mean-first-passage asymptotic can be obtained.

Acceptance:
- a reversible-chain quasipotential theorem plus the exponential escape
  exponent;
- keep the prefactor open unless a genuine Eyring-Kramers theorem is verified.
