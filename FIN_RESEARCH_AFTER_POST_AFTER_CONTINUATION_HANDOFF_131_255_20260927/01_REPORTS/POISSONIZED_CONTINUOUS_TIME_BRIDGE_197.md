# POISSONIZED-CONTINUOUS-TIME-BRIDGE-197
## The exact leave-one-out continuous-time process is the Poissonization of refresh-depth dynamics

Date: 2026-09-26

Status:
exact identity for the declared homogeneous leave-one-out heat-bath process.

Let K_N be the one-refresh count-state kernel:
- choose one of N copies uniformly;
- redraw its label from the leave-one-out Gibbs conditional.

Let each copy attempt refresh at common rate rho.

Then total refresh attempts occur at rate

    N rho.

The continuous-time count generator is exactly

    boxed:
    L_N
      =
      N rho (K_N-I).

## 1. Exact semigroup

Uniformization gives

    boxed:
    exp(t L_N)
      =
      exp(-N rho t)
      sum_(m=0)^infinity
      [(N rho t)^m/m!]
      K_N^m.

So the continuous-time process is:
1. draw a Poisson number of discrete refresh-depth steps;
2. evolve the embedded refresh chain that many steps.

No limiting approximation is required.

## 2. History depth as a continuous-time clock

Let M_t be the number of refresh attempts by time t.

Then

    M_t ~ Poisson(N rho t),

so

    E[M_t]=N rho t,
    Var(M_t)=N rho t.

Define

    tau_depth
      =
      M_t/N.

Then

    E[tau_depth]
      =
      rho t,

and

    relative standard deviation
      =
      1/sqrt(N rho t).

Thus discrete transformation depth converges statistically to a continuous
dimensionless duration as the number of observed refreshes grows.

## 3. Additive independent increments

For disjoint intervals:

    M_(t+s)-M_t

is independent of the previous count and has mean N rho s.

Therefore the Poissonized history clock is:
- additive;
- composition-consistent;
- compatible with the reversible shift depth.

## 4. Residual gauge

The transformation

    rho -> c rho,
    t -> t/c

leaves:
- the embedded kernel K_N;
- all refresh-depth records;
- the state law exp(tL_N)

invariant after reparametrization.

So Poissonization constructs a consistent continuous clock coordinate, but not
an absolute time unit.

This closes the discrete-depth -> continuous-time bridge within the declared
kinetic model.
