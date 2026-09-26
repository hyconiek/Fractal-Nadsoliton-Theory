# LOCAL-CAPACITY-GATE-134
## The exact Gibbs chain turns the k6 isolation question into two local capacity exponents

Date: 2026-09-26

Status:
- exact reversible-network capacity identities;
- local exponential barriers use the same compact-tube logic as report 63;
- global capacity comparison is deliberately left open.

## 1. Exact conductance network

For the leave-one-out count chain, let

    pi_N(n)

be the exact stationary count Gibbs law and

    r_N(n,n')

the exact jump rate.

Detailed balance gives symmetric edge conductance

    c_N(n,n')
      =
      pi_N(n) r_N(n,n').

For disjoint sets A,B, the capacity is the Dirichlet minimum

    cap_N(A,B)
      =
      inf_h
      (1/2) sum_(n,n')
        c_N(n,n')
        [h(n)-h(n')]^2,

with h=1 on A and h=0 on B.

The dual Thomson principle gives the equivalent minimum-resistance flow
description.

This is the correct object for transition rates in the chosen microscopic
process.

## 2. Explicit-path lower bound

Keep only the edges of any count-lattice path gamma_N from A to B.

Rayleigh monotonicity gives

    cap_N(A,B)
      >=
      [
        sum_(e in gamma_N)
          1/c_N(e)
      ]^(-1).

For q=12 the number of count states and a smooth lattice path length are only
polynomial in N.

Inside a compact interior tube the leave-one-out jump probabilities are
bounded away from exponential zero.

Therefore the N-speed exponent of the path resistance is controlled entirely
by the maximum Gibbs rate function V_g along that path.

This reproduces, now as a capacity bound rather than a bare barrier analogy,

    local path exponent
      =
      max_gamma V_g
      -V_g(start).

## 3. Local cut upper bound

Conversely, choose a local basin tube U whose boundary is separated from the
stable point by one nondegenerate index-one saddle and for which all boundary
crossing states have

    V_g >= V_saddle-o(1).

The cut-set bound and the polynomial number of count states give

    cap_N(U,U^c)
      <=
      exp[
        -N(V_saddle-constant)
        +o(N)
      ].

Together with a lattice path through the saddle this yields the local
capacity exponent

    -(1/N) log cap_N
      =
      V_saddle-constant+o(1),

equivalent to the local exit theorem of report 63.

## 4. Application to the two k6 channels

### Direct pair channel

In a local pitchfork tube containing +k6, uniform, -k6:

    B_pair
      =
      Phi_uniform-Phi_k6.

### Localized escape channel

In a local tube around the main index-one boundary saddle connecting k6 and
the localized basin:

    B_out
      =
      Phi_main_saddle-Phi_k6.

So the two local capacities have exponential orders

    cap_pair
      ~ exp[-N B_pair]

and

    cap_out
      ~ exp[-N B_out]

up to common basin-weight conventions and subexponential factors.

This converts report 132's energy comparison into a local
potential-theoretic comparison.

## 5. What is still missing for a global two-state theorem

A true autonomous binary reduction requires the GLOBAL capacity relation

    cap(pair-internal)
      >>
    cap(pair -> outside)

in the appropriate quasistationary normalization.

Local tubes are insufficient because another route outside the selected tubes
may have lower communication height.

The remaining proof obligation is therefore:

1. define the +k6 and -k6 metastable basin sets in the finite-N count space;
2. construct a global trial equilibrium potential for an upper bound;
3. construct matching Thomson flows for lower bounds;
4. show that all other exits are exponentially subdominant in a nonempty
   parameter interval.

The result may be positive only on a subwindow of
`(g6,g_x)`.

This is a bounded capacity problem; no new stationary atlas is required.
