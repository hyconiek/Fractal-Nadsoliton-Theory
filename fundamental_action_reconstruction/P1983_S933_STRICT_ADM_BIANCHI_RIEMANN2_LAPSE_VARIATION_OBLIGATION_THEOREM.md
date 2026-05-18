# P1983 / S933 — Strict ADM/Bianchi-I Riemann² Lapse-Variation Obligation Theorem Packet

Status: `OPEN_OBSTRUCTION_WITH_TRACE`  
Local theorem result: `PASS_RIEMANN2_LAPSE_VARIATION_OBLIGATION_EXPORTED`  
Route: `strict_only`  
Date anchor: `2026-05-18`

## Professor-level decision

`P1982` exported the Ricci-squared lapse obligation and showed that the
curvature-squared sector contains higher-derivative anisotropic structures.  The
next honest step is the remaining independent quadratic invariant before
Gauss-Bonnet combination: `Riemann_{mu nu rho sigma}R^{mu nu rho sigma}`.

`P1983` performs that calculation in the same strict ADM/Bianchi-I lane with
lapse `N(t)` retained.  It does not claim Gauss-Bonnet, coefficient-level,
spatial-EOM, background-independence, or ToE closure.

## Strict inputs

This packet continues the strict-only ADM/Bianchi-I route after `P1980/P1981/P1982`.
It uses no legacy-kernel physical role transfer.

## Riemann block setup

For diagonal Bianchi-I, define

```text
h_i = dot(a_i)/a_i,
sigma1 + sigma2 + sigma3 = 0.
```

The lapse-retaining orthonormal curvature blocks are represented as

```text
E_i = [dot(h_i) + h_i^2 - h_i*Ndot/N]/N^2,
F_ij = h_i*h_j/N^2.
```

The Kretschmann scalar is then

```text
Riemann^2 = 4*sum_i E_i^2 + 4*sum_{i<j} F_ij^2.
```

The FRW check reduces to

```text
12*((Hdot + H^2 - H*Ndot/N)^2/N^4 + H^4/N^4),
```

which is the standard flat-FRW Kretschmann form with lapse retained.

## Riemann² lapse Euler operator

The density difference is

```text
Delta L_Riemann2 = N*V*(Riemann2_BianchiI - Riemann2_FRW).
```

Because `E_i` contains `Ndot`, the correct lapse operator is

```text
E_N(Delta L_Riemann2)
= partial_N Delta L_Riemann2 - d/dt(partial_Ndot Delta L_Riemann2).
```

The full expression is exported machine-readably in the JSON witness.  The
gatekeeper checks prove:

```text
E_N(Delta L_Riemann2) -> 0 in the isotropic limit,
E_N(Delta L_Riemann2) contains d2sigma terms,
E_N(Delta L_Riemann2) contains Nddot terms,
E_N(Delta L_Riemann2) contains shear-velocity-square terms,
Riemann2_BianchiI - Riemann2_FRW contains cubic shear/velocity mixing.
```

## Machine-check export

Execution script:

```text
fundamental_action_reconstruction/p1983_s933_strict_adm_bianchi_riemann2_lapse_variation_obligation.py
```

Generated witness:

```text
fundamental_action_reconstruction/generated/p1983_s933_strict_adm_bianchi_riemann2_lapse_variation_obligation.json
```

The witness exports:

1. the lapse-retaining orthonormal Riemann blocks,
2. `Riemann2_BianchiI`, `Riemann2_FRW`, and their difference,
3. the full lapse Euler operator,
4. structural checks for isotropic vanishing and higher-derivative terms,
5. rational numeric replay samples,
6. explicit no-false-pass boundaries.

## Theorem statement actually licensed

For the diagonal trace-free Bianchi-I minisuperspace reduction, the
Riemann-squared lapse Euler operator difference is an exactly exported symbolic
expression that vanishes in the isotropic limit but contains higher-derivative
shear/lapse structures.

Therefore the strict curvature-squared background-independence problem still
requires a Gauss-Bonnet combination and strict coefficient-level comparison; it
is not closed by the individual `R^2`, `Ricci^2`, or `Riemann^2` packets.

## False-pass boundary

This packet does **not** prove:

1. Gauss-Bonnet ADM variation,
2. strict coefficient-level curvature-squared cancellation,
3. spatial Bianchi-I EOM cancellation,
4. global background-independence closure,
5. PO2/PO3 closure,
6. BRST or Cutkosky closure,
7. `QW-2191` selector discharge,
8. ToE closure.

It only exports the Riemann-squared lapse-equation obligation.

## Progress toward ToE

This is progress because all three elementary curvature-squared lapse operators
needed for the Gauss-Bonnet combination are now available as machine-checkable
objects: `R^2`, `Ricci^2`, and `Riemann^2`.  The next ToE-relevant question is
whether their higher-derivative shear terms cancel in the Gauss-Bonnet
combination and then in the strict coefficient-level gravity sector.

## Next honest step

Build the Gauss-Bonnet ADM/Bianchi-I lapse combination

```text
G_GB = Riemann^2 - 4*Ricci^2 + R^2
```

from `P1981/P1982/P1983`, test cancellation of the higher-derivative shear
structures, and only then combine with the strict coefficients from
`P1972/P1853`.

## Lay explanation

Policzyliśmy trzeci składnik kwadratowy: `Riemann^2`.  Tak jak `Ricci^2`, ma on
trudne wyrazy z wyższymi pochodnymi.  Dobra wiadomość jest taka, że mamy już
trzy potrzebne klocki do sprawdzenia kombinacji Gaussa-Bonneta.  Zła wiadomość:
bez tej kombinacji i bez współczynników strict nadal nie wolno mówić o domknięciu
ToE.
