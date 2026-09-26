# PHASE-ROBUSTNESS-THEOREM-31 — analytic full-7D phase reduction for the drift-robust quartic probe

Status: **PROOF-GRADE CANDIDATE, CONDITIONAL ON THE DECLARED FINITE-N HEAT-BATH + ME7 LANE AND THE EXISTING CO-PHASED INTERVAL CERTIFICATE**.

Baseline: `FIN_POST_CONTINUATION_FULL_HANDOFF_20260926`, commit `fe14a6f4e436815635df54429102f23c22862296`.

This result targets the open item `P0 — PHASE-ROBUSTNESS-THEOREM` from `NEXT_RESEARCH_PROGRAM.md`.
It does not change any foundational physical guardrail.

## 1. Input already accepted in the baseline

For a retained real probe `phi` on the coefficient sphere, let

    H = ||P_H(phi^2)||_u^2,
    C = -12 H + 2 g <P_H(phi^2), P_H(phi A7 phi)>_u,
    rho(phi) = 4 C / (g sqrt(H)).

The co-phased four-sector optimization is already interval-certified at

    rho_* in [0.71753326806113993, 0.71753326806121198].

The only missing step was to exclude a larger value at non-co-phased points in the full seven-real-coordinate retained space.

Declared constants:

    g  = 3.7183448981203875
    l3 = 1.9614068619764455
    l4 = 2.1995688493332102
    l5 = 2.298606272079097
    l6 = 2.3421820411462999

## 2. Complex-channel reduction

Write the retained Fourier amplitudes as `b3,b4,b5` and real `b6` and define

    t1 = b4 conj(b3)
    t2 = b5 conj(b4)
    t3 = b6 conj(b5)

    t4 = b5 conj(b3)
    t5 = b6 conj(b4)
    t6 = conj(b5)^2 / 2.

Then

    H1 = 2(t1+t2+t3),
    H2 = 2(t4+t5+t6).

Define

    S1=t1+t2+t3,  S2=t4+t5+t6.

For the weighted channel define six positive constants

    d1 = g(l3+l4)-12
    d2 = g(l4+l5)-12
    d3 = g(l5+l6)-12
    d4 = g(l3+l5)-12
    d5 = g(l4+l6)-12
    d6 = 2 g l5-12

and

    U1=d1 t1+d2 t2+d3 t3,
    U2=d4 t4+d5 t5+d6 t6.

The general quartic theorem becomes exactly

    H = 8 (|S1|^2+|S2|^2),
    C = 8 Re[conj(S1)U1 + conj(S2)U2].

There is also an exact one-parameter phase gauge

    theta3 -> theta3+3a,
    theta4 -> theta4+2a,
    theta5 -> theta5+a,

under which channel 1 only acquires a common phase and channel 2 twice that phase. Therefore only two relative phase variables are physical for this functional.

## 3. Co-phased comparison at fixed amplitudes

For fixed sector magnitudes, let subscript `0` denote the co-phased point, i.e. every `t_i` in a given channel is positive real. Then

    H0 >= H,
    C0 >= C,
    ||U0|| >= ||U||.

More precisely, pair by pair,

    H0-H = sum 16 |ti tj| (1-cos Delta_ij),

while

    C0-C = sum 8 (di+dj)|ti tj| (1-cos Delta_ij).

Hence, with

    d = H0-H,
    e = C0-C,

we have

    e >= r_min d,

where the smallest pair-average is

    r_min = min_pairs (di+dj)/2
          = 4.098854660453307861... > 4.098.

The six constituent weights satisfy

    d_min = 3.471942807351108633... > 3.471,
    d_max = 5.256051547738373390... < 5.257.

Since `C0/H0` is a positive weighted average of the six `di` and their pair averages,

    R := C0/H0 <= d_max < 5.257.

## 4. A candidate beating rho_* would have to lose very little H

Cauchy-Schwarz gives

    C/sqrt(H) <= sqrt(8) ||U||.

Therefore any hypothetical point with `rho > rho_U`, where

    rho_U = 0.71753326806121198,

must satisfy

    ||U|| > T := rho_U g/(4 sqrt(8))

with

    T = 0.235823308225280253... > 0.23582.

The loss in the weighted vector controls the loss in H:

    ||U0||^2-||U||^2 >= (d_min^2/8)(H0-H).

Also `||U0|| <= d_max sqrt(H0/8)`.

It remains to bound `H0`. On the coefficient sphere `c3^2+c4^2+c5^2+c6^2=1`, write

    S10 = c3 c4/24 + c4 c5/24 + c5 c6/(12 sqrt(2)),
    S20 = c3 c5/24 + c4 c6/(12 sqrt(2)) + c5^2/48.

These are quadratic forms `c^T M1 c` and `c^T M2 c`. Their exact spectral norms are

    ||M1||^2 = (2+sqrt(2))/2304,
    ||M2||   = (1+sqrt(5))/96.

Thus

    H0 = 8(S10^2+S20^2)
       <= (7+sqrt(5)+2sqrt(2))/576
       = 0.020945303996954826... < 0.020946.

Combining the previous inequalities, every hypothetical point above `rho_U` must obey

    x := (H0-H)/H0 < 0.530906.

Using the non-rounded declared constants gives the stronger value

    x < 0.529681780828920.

## 5. Phase loss below this threshold cannot improve rho

From `e >= r_min d` and `R<=d_max`,

    e/C0 >= (r_min/d_max) x.

Let

    a = r_min/d_max > 4.098/5.257.

For `0<=x<=1`, the condition

    a x >= 1-sqrt(1-x)

is equivalent (for `a>1/2`) to

    x <= (2a-1)/a^2.

The conservative bounds give

    (2a-1)/a^2 > 0.920012,

and the declared constants give

    (2a-1)/a^2 = 0.92029428255273958...

But a hypothetical counterexample was already forced to have

    x < 0.530906 < 0.920012.

Therefore

    e/C0 >= 1-sqrt(1-x),

so

    C = C0-e <= C0 sqrt(1-x),
    sqrt(H) = sqrt(H0) sqrt(1-x),

and hence

    C/sqrt(H) <= C0/sqrt(H0).

Thus every point capable of competing with the certified optimum is dominated by the co-phased point with the same sector magnitudes.

## 6. Full-7D conclusion

Assume for contradiction that a non-co-phased retained probe satisfies

    rho(phi) > rho_U.

Sections 4–5 imply that its co-phased probe with the same sector magnitudes has

    rho(phi0) >= rho(phi) > rho_U.

This contradicts the already accepted interval branch-and-bound + Krawczyk certificate that `rho_U` is a global upper bound in the co-phased four-sector class.

Therefore the co-phased global certificate extends to the full seven-real-coordinate retained probe space:

    rho_global in [0.71753326806113993, 0.71753326806121198].

No non-co-phased probe can exceed the certified co-phased optimum.

## 7. Audit status and boundaries

This closes the mathematical gap targeted by `P0 — PHASE-ROBUSTNESS-THEOREM`, conditional on the already accepted GENERAL-QUARTIC-THEOREM-13 and the existing co-phased interval certificate.

Before promotion in the repository claim ledger, independently replay:
1. the six `di` values and conservative inequalities;
2. the exact spectral-norm formulas for `M1,M2`;
3. the algebraic identities `H=8||S||^2` and `C=8 Re<S,U>`;
4. the existing co-phased interval upper bound `rho_U`.

No physical clock, apparatus, activity law, scale/r source, QW-2191, legacy-to-strict bridge/role transfer, role-bearing `L_total`, SM/GR or ToE conclusion follows.
