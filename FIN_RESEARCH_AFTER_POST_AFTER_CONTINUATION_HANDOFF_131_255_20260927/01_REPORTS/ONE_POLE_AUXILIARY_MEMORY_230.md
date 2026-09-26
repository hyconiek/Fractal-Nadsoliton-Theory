# ONE-POLE-AUXILIARY-MEMORY-230
## A two-variable model determined only by A,M0,M1 reproduces the full N=6 transient from t=0 at sub-percent accuracy

Date: 2026-09-26

Status:
moment-matched dynamical approximation tested against the exact microscopic
semigroup.

Approximate one symmetry-sector memory kernel by

    K(t)
      =
      a exp(-gamma t).

Match its first two integrated moments:

    M0 = a/gamma,
    M1 = a/gamma^2.

Therefore

    boxed:
    gamma=M0/M1,
    a=M0^2/M1.

No time-domain fitting is used.

Introduce one auxiliary memory variable z:

    u_dot = A u + z
    z_dot = a u - gamma z

with

    u(0)=1,
    z(0)=0.

This reproduces:
- the exact instantaneous slope A at t=0;
- the measured M0 and M1;
- the slow renormalized pole.

## 1. Memory scale

For N=6 the six Fourier sectors give

    gamma≈2.94 ... 3.28

or

    tau_mem≈0.305 ... 0.340

microscopic clock units.

## 2. Exact-semigroup comparison

Maximum error over all six modes:


    t=0.1:
      absolute=0.0006068
      relative=0.0615 %

    t=0.25:
      absolute=0.0023447
      relative=0.2411 %

    t=0.5:
      absolute=0.0043997
      relative=0.4601 %

    t=1:
      absolute=0.0044737
      relative=0.4795 %

    t=2:
      absolute=0.0016645
      relative=0.1856 %

    t=4:
      absolute=0.0002256
      relative=0.0270 %

    t=8:
      absolute=0.0004460
      relative=0.0615 %

    t=16:
      absolute=0.0002894
      relative=0.0528 %

    t=32:
      absolute=0.0001073
      relative=0.0330 %

    t=64:
      absolute=0.0000220
      relative=0.0089 %

The worst sampled relative error is below 0.5%.

From t=2 onward it is below about 0.19%, and from t=4 onward below about
0.062%.

Unlike the slip-only exponential, this model is accurate already from t=0.

## 3. Interpretation

The short memory can be represented, to high accuracy, by ONE auxiliary
relaxing degree of freedom per symmetry sector.

That auxiliary variable is an effective bookkeeping state, not a newly
identified physical particle or field.

The important point is structural:

    non-Markov projected dynamics
      ->
    one extra short-lived state
      ->
    nearly Markov enlarged dynamics.

This gives a concrete realization of the general principle that hidden memory
can be restored as explicit state.
