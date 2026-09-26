# FOLD-FINITE-N-METASTABILITY-63
## Local quasipotential barrier for the exact finite-N reversible Gibbs sampler

Date: 2026-09-26

Inputs:
- exact finite-copy Gibbs law R7P-111;
- exact finite-N Gibbs sampler from report 61;
- certified simple fold R7P-031;
- barrier theorem from report 60.

Status:
- exponential-scale local metastability theorem under an explicitly local
  fold tube;
- no global basin claim;
- no Eyring-Kramers prefactor.

## 1. Exact reversible chain

The finite-N count heat bath is reversible with respect to

    Pi_N(n)
      proportional to
      Multinomial(N;n)
      12^(-N)
      exp[
        N g p^T A7 p/2
      ],

where p=n/N.

Stirling/Sanov asymptotics give, uniformly on compact interior subsets,

    -(1/N) log Pi_N(p)
      = V_g(p)+constant+o(1).

## 2. Local fold tube

For g>g_f sufficiently close to the certified simple fold, choose a small
neighborhood U in which:

- there is one local stable stationary branch p_+(g);
- there is one index-one saddle branch p_-(g);
- all six transverse dual directions remain strictly positive;
- the saddle separates the two local center-manifold sides.

This is a LOCAL construction using the R7P-031 simple-fold theorem. It does
not assert that no lower barrier exists outside U.

## 3. Subexponential jump rates

Inside a compact interior U every component of p is bounded away from zero.

The exact conditional probabilities q_j^{(i)} are then bounded above and below
by positive constants independent of N.

Hence every allowed count jump inside U has rate between polynomial scales in N
and contains no factor exp(-cN).

Therefore the N-speed exponential cost comes entirely from the Gibbs weights,
not from exponentially small kinetic rates.

## 4. Local quasipotential

For a reversible chain with subexponential edge rates, detailed balance assigns
the exponential communication cost through the stationary potential.

Thus the local quasipotential difference between the stable branch and the
separating saddle is

    boxed:
    Delta W_U(g)
      = V_g(p_-(g))-V_g(p_+(g)).

By the exact stationary dual/primal identity this equals

    Phi_g(theta_-)-Phi_g(theta_+).

Using report 60,

    Delta W_U(g)
      =
      C_fold (g-g_f)^(3/2)
      +O((g-g_f)^2),

with

    0.5357594331509
      < C_fold <
    0.5357598007359.

## 5. Local exit exponent

For fixed small delta=g-g_f>0, initialize in a local-equilibrium/quasistationary
ensemble concentrated in the stable side of U.

On the exponential N-scale,

    boxed:
    (1/N) log E[tau_U]
      -> Delta W_U(g),

provided exit is defined through the local saddle-separating boundary of U.

Therefore, locally,

    E[tau_U]
      =
      exp[
        N C_fold delta^(3/2)
        +o(N delta^(3/2))
      ]

on the metastable exponential scale.

No prefactor is claimed.

## 6. Three finite-size regimes

The dimensionless barrier variable is

    B_N=N C_fold delta^(3/2).

Hence:

### Metastable
    N delta^(3/2) -> infinity

The local exit time is exponentially separated.

### Crossover
    N delta^(3/2)=O(1)

No exponentially large time-scale separation survives.

### Fold-dominated
    N delta^(3/2) -> 0

The stable/saddle barrier disappears on the N-speed exponential scale.

Equivalently the crossover is

    delta~N^(-2/3).

## 7. Scope boundary

This theorem is deliberately local.

It does NOT prove:
- the local saddle is the lowest global escape route;
- a global switching time between all FIN phases;
- an Eyring-Kramers prefactor;
- laboratory time.

Those require a global communication-height analysis and a sourced clock.
