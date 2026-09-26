# INTERNAL-POISSON-CLOCK-AGREEMENT-194
## Homogeneous microscopic refresh gives asymptotically agreeing independent internal clocks up to one global rate gauge

Date: 2026-09-26

Status:
exact conditional theorem inside the declared continuous-time leave-one-out
heat-bath contract.

Assume every microscopic copy carries the SAME Poisson attempt rate rho.

Take two disjoint groups A,B with sizes

    m_A,
    m_B.

Let

    N_A(t),
    N_B(t)

be their attempt counts.

Because Poisson clocks superpose,

    N_A(t) ~ Poisson(m_A rho t),
    N_B(t) ~ Poisson(m_B rho t).

Define normalized clocks

    tau_A=N_A/m_A,
    tau_B=N_B/m_B.

Then

    E[tau_A]
      =
      E[tau_B]
      =
      rho t.

Their difference has

    boxed:
    Var(tau_A-tau_B)
      =
      rho t(1/m_A+1/m_B).

Hence

    RMS(tau_A-tau_B)/(rho t)
      =
      sqrt[
        (1/m_A+1/m_B)/(rho t)
      ]
      ->
      0.

So independently constructed internal attempt clocks synchronize by a law of
large numbers.

## 1. What is non-arbitrary

The relative calibration between equal microscopic copies is fixed:

    one attempt per copy

is the same clock construction everywhere.

No separate rate is fitted for A and B.

## 2. What remains arbitrary

A global replacement

    rho -> c rho,
    t -> t/c

leaves all dimensionless event records unchanged.

Therefore internal clock AGREEMENT is obtained, but absolute physical duration
is not.

## 3. Failure mode

If different copies are allowed independent attempt rates rho_x, the normalized
clocks converge to different rates and synchronization fails.

Thus clock agreement is a genuine consequence of the HOMOGENEOUS heat-bath
kinetic contract, not of the static FIN potential.

## 4. Importance

This is stronger than report 193.

FIN now has, conditionally on one common microscopic refresh law:
- additive transformation depth;
- multiple internal clocks;
- asymptotic agreement between those clocks;
- one residual global scale gauge.

That is close to the strongest clock statement possible before external
calibration or a new internally sourced dimensional anchor.
