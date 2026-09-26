# CRT-BARRIER-SEPARATION-206
## The base/fiber quotient is dynamically favored because Z4 moves have a lower barrier than Z3 moves

Date: 2026-09-26

Status:
numerical saddle-level result at the current working gain, combined with exact
CRT typing.

At

    g=5.145228719489142,

the relevant localized saddle barriers are

    d=3:
      Delta V3≈0.644851906971

    d=4:
      Delta V4≈0.662219456770

    d=5:
      Delta V5≈0.782655...

The CRT interpretation is:

    d3 = fiber-only Z4 move,
    d4 = base-only Z3 move,
    d5 = mixed move.

Thus

    boxed:
    Delta V_fiber
      <
    Delta V_base
      <
    Delta V_mixed

for these identified transition families.

The base-fiber barrier split is

    Delta
      =
      Delta V4-Delta V3
      ≈0.017367549800.

## 1. Large-N implication

At exponential level and ignoring prefactor differences,

    tau_base / tau_fiber
      ~
      exp[N Delta].

Illustrative values:


    N=8:
      exp(N Delta)≈1.14906

    N=20:
      exp(N Delta)≈1.41531

    N=50:
      exp(N Delta)≈2.38304

    N=100:
      exp(N Delta)≈5.67889

    N=200:
      exp(N Delta)≈32.2497

    N=300:
      exp(N Delta)≈183.143

    N=500:
      exp(N Delta)≈5906.3

    N=1000:
      exp(N Delta)≈3.48844e+07


At the currently computable N<=8 this separation is mild.

But if the same barriers control the large-N asymptotics, the separation grows
exponentially.

This gives a concrete mechanism for:

    fast equilibration in the Z4 fiber
      ->
    slower motion of the Z3 base coordinate.

## 2. Connection to reports 139-144

The exact finite-N reduction already found that:
- the three Z3 sectors carry the slow observable;
- eliminated dynamics creates short memory;
- M0+M1 reproduces the slow Z3 rate with improving accuracy.

The barrier ordering now supplies a structural explanation for why that
particular quotient is dynamically useful.

## 3. Caveat

This is not yet a certified large-N capacity theorem.

The exp[N Delta] formula is a saddle-level asymptotic candidate until:
- metastable cores;
- prefactors;
- competing paths

are controlled.
