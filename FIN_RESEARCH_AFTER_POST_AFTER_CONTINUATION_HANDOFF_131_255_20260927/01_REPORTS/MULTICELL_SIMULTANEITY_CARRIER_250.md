# MULTICELL-SIMULTANEITY-CARRIER-250
## A finite cyclic register gives an explicit simultaneous reversible carrier that reproduces the Z3 reset for a controlled time window

Date: 2026-09-26

Status:
exact finite-state construction.

Let there be L simultaneous Z3 slots

    z_0,z_1,...,z_(L-1).

Slot 0 is the currently observed subsystem.

The other L-1 slots form a finite environment.

Define the global deterministic update

    boxed:
    T(
      z_0,z_1,...,z_(L-1)
    )
      =
    (
      z_1,z_2,...,z_(L-1),z_0
    ).

This is one cyclic permutation of the SLOT CONTENTS.

It is bijective and preserves global information exactly.

## 1. Exact discrete reset window

Prepare:
- z_0 arbitrarily;
- z_1,...,z_(L-1) independently uniform on Z3.

After m update events,

    z_0(m)
      =
      z_m(0)

for

    0<=m<L.

Therefore for every

    m=1,...,L-1,

the observed output is:
- uniform;
- independent of the original subsystem state;
- independent of all earlier observed outputs.

Thus the reduced subsystem experiences EXACT full MaxEnt resets for

    boxed:
    L-1 consecutive events.

At event L:

    z_0(L)=z_0(0).

The first exact recurrence appears.

This is the finite-environment version of the all-time history dilation.

## 2. Poissonized clock

Let the cyclic update occur at rate rho and set

    tau=rho t.

The event count M_t is Poisson(tau).

The subsystem remembers its initial label exactly when

    M_t=0 mod L.

Define

    p_0^(L)(tau)
      =
      exp(-tau)
      sum_(q=0)^infinity
      tau^(qL)/(qL)!.

Then the exact tagged-slot kernel is

    boxed:
    K_L(t)
      =
      p_0^(L)(tau) I
      +
      [1-p_0^(L)(tau)] U.

The ideal infinite fresh-ancilla heat bath is

    K_inf(t)
      =
      exp(-tau) I
      +
      [1-exp(-tau)] U.

Therefore the recurrence error is completely explicit.

For a point initial state the total-variation error is

    boxed:
    TV
      =
      (2/3)[
        p_0^(L)(tau)-exp(-tau)
      ]

and is bounded by

    TV
      <=
      (2/3)
      P[
        Poisson(tau)>=L
      ].

## 3. High-order agreement

The first recurrence contribution occurs at exactly L events.

Hence

    K_L(t)-K_inf(t)
      =
      O(tau^L/L!).

The finite register matches the ideal heat-bath expansion through order L-1.

## 4. Examples

For L=12:

    tau=1:
      TV≈5.12e-10

    tau=5:
      TV≈2.29e-3.

For L=24:

    tau=10:
      TV≈4.88e-5

    tau=12:
      TV≈5.25e-4.

So a modest finite simultaneous environment can reproduce the open Z3
heat-bath extremely accurately over many refresh times before recurrence.

## Why this is different from report 188

Report 188 used an infinite two-sided HISTORY carrier.

Here all L trit slots coexist at one instant.

The transformation acts on SLOT IDENTITY inside one global state.

This is the first current construction that explicitly separates:

    simultaneous slot structure

from

    history/time ordering.
