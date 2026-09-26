# FINITE-N-FOLD-BARRIER-60
## Certified saddle-node barrier and the N^{-2/3} crossover scale

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Inputs:
- interval-certified fold `R7P-031`;
- exact dual/primal bridge `MP7-011`;
- conditional finite-copy Gibbs realization `R7P-111`.

Status:
- local asymptotic theorem with certified leading coefficient;
- finite-N crossover is an exponential-scale statement for the conditional Gibbs model;
- no Eyring-Kramers prefactor or physical escape time is claimed.

## 1. Fold normal form

Let

    delta = g-g_f,

and use the normalized fold coordinate y along the certified kernel vector v.

The dual stationary equation has

    F_y
      = a delta + (b/2)y^2
        + O(delta |y| + |y|^3 + delta^2),

with certified signs

    a<0,
    b>0.

Integrating in y gives the local dual potential

    Phi
      = Phi_f
        + a delta y
        + (b/6)y^3
        + higher order.

For delta>0, the two stationary branches are

    y_+/- =
      +/- sqrt(-2a/b) delta^(1/2)
      + O(delta).

The + branch is locally stable in the declared maximum-entropy heat-bath flow;
the - branch is the neighboring saddle branch.

## 2. Barrier exponent

At a stationary branch, use

    a delta = -(b/2)y^2.

Then the leading potential value is

    Phi(y)
      = Phi_f -(b/3)y^3 + ...

Therefore the saddle-minus-minimum barrier is

    Delta Phi
      = (2b/3)
        (-2a delta/b)^(3/2)
        + O(delta^2)

      = C_fold delta^(3/2)
        + O(delta^2),

where

    C_fold
      = (2b/3)(-2a/b)^(3/2).

Propagating the accepted R7P-031 intervals gives

    boxed:
    0.5357594331509
      < C_fold <
    0.5357598007359.

Thus the 3/2 barrier exponent is fixed locally by the already certified fold.

## 3. Why the same barrier belongs to V_g

At every stationary dual/primal pair,

    theta = g X^T p,
    p = softmax(X theta),

the exact MP7-011 bridge gives

    V_g(p)=Phi_g(theta).

Therefore the difference between the stable and saddle stationary values is
the same in the primal finite-copy rate function and the dual fold potential.

Hence

    boxed:
    Delta V_g
      = C_fold (g-g_f)^(3/2)
        + O((g-g_f)^2).

No identification of the full local p-space and theta-space potentials is
needed for this stationary-value statement.

## 4. Finite-copy Gibbs exponent

R7P-111 supplies the exact labelled-copy measure

    Pi_N(x)
      proportional to
      12^(-N)
      exp[
        (g/(2N))
        sum_{a,b} A7[x_a,x_b]
      ].

At empirical p its N-speed large-deviation rate is V_g.

Therefore the exponential suppression associated with the local saddle barrier
is

    exp[-N Delta V_g]

up to subexponential factors.

The exponent is

    N C_fold (g-g_f)^(3/2)
    + lower-order terms.

## 5. Crossover scale

A finite-N crossover occurs when the exponential barrier is O(1):

    N Delta V_g ~ 1.

Hence

    boxed:
    g-g_f = O(N^(-2/3)).

If the convention `N Delta V_g=1` is used for the leading normal form, then

    g-g_f
      ~ C_fold^(-2/3) N^(-2/3),

with certified coefficient

    boxed:
    1.515955956
      < C_fold^(-2/3) <
    1.515956650.

This number is not a universal physical threshold; it belongs to the declared
conditional Gibbs normalization.

## 6. Associated branch scale

Since

    |y|~delta^(1/2),

the crossover gives

    |y|=O(N^(-1/3)).

So the natural saddle-node finite-size scaling pair is

    delta ~ N^(-2/3),
    y     ~ N^(-1/3).

This is also the canonical scaling obtained by balancing

    N[a delta y + (b/6)y^3]=O(1).

## 7. Numerical replay

Direct stationary solves of the exact reflection-even C4 equations give

    Delta Phi / delta^(3/2)
        -> 0.5357596...

as delta -> 0+.

Representative values:

    delta=1e-6   : 0.53575899
    delta=1e-5   : 0.53575826
    delta=1e-4   : 0.53574624
    delta=1e-3   : 0.53562593.

The deviations at larger delta are the expected higher-order fold corrections.

## 8. What is NOT yet proved

This report does not prove an Eyring-Kramers mean first-passage formula.

To promote

    escape time ~ exp(N Delta V_g)

with a controlled prefactor, one needs the exact finite-N dynamics and a
metastability theorem for that chain.

The next report constructs the exact finite-N Gibbs sampler corresponding to
R7P-111.
