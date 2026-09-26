# INFORMATION-SURPRISAL-CLOCK-195
## Transition information supplies a second internal dimensionless clock that agrees with refresh depth asymptotically

Date: 2026-09-26

Status:
exact ergodic theorem for the finite leave-one-out refresh chain.

Let K_N be the discrete one-refresh transition kernel and pi_N its stationary
law.

For one observed transition

    X_t -> X_(t+1),

define the transition surprisal

    s_t
      =
      -log K_N(X_t,X_(t+1)).

The stationary entropy rate is

    h_N
      =
      E_pi[s_t].

Define the information clock

    boxed:
    T_info(T)
      =
      (1/h_N)
      sum_(t=0)^(T-1) s_t.

## 1. Asymptotic agreement

The finite leave-one-out chain is irreducible and aperiodic because every
softmax transition probability is strictly positive and self-refreshes occur.

Therefore the Markov ergodic theorem gives

    (1/T) sum s_t
      ->
      h_N

almost surely.

Hence

    boxed:
    T_info(T)/T
      ->
      1.

So:
- raw refresh count T;
- accumulated path information divided by h_N;

are independently constructed clocks that agree asymptotically.

## 2. Exact current entropy rates

At g=5.145228719489142 the stationary count-process entropy rates are:


    N=3:
      h_N=1.654059170165 bits/refresh

    N=4:
      h_N=1.389660322536 bits/refresh

    N=5:
      h_N=1.171214165819 bits/refresh

    N=6:
      h_N=0.991254126943 bits/refresh

    N=7:
      h_N=0.848943363111 bits/refresh

    N=8:
      h_N=0.740451994778 bits/refresh


The value of h_N changes with internal copy number N, but the normalized clock
remains asymptotically one refresh-depth unit per refresh.

## 3. Why this is not a physical second

K_N describes the embedded one-refresh chain.

A global continuous-time rate rescaling changes how fast refreshes occur in an
external time coordinate but does not change K_N or h_N.

Therefore the information clock measures:

    transformation / refresh depth,

not absolute duration.

## 4. Importance

FIN now has three mutually distinct internal constructions that can agree on
the same dimensionless history parameter:

1. reversible shift depth;
2. normalized Poisson attempt count;
3. normalized accumulated transition surprisal.

Their agreement is nontrivial but all share one residual global rate gauge.
