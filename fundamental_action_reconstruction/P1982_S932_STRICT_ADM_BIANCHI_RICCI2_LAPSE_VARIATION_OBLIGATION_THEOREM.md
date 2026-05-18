# P1982 / S932 — Strict ADM/Bianchi-I Ricci² Lapse-Variation Obligation Theorem Packet

Status: `OPEN_OBSTRUCTION_WITH_TRACE`  
Local theorem result: `PASS_RICCI2_LAPSE_VARIATION_OBLIGATION_EXPORTED`  
Route: `strict_only`  
Date anchor: `2026-05-18`

## Professor-level decision

`P1981` showed that the `R^2` curvature-squared lapse equation already contains
new shear structures beyond the Einstein-Hilbert `-Q_shear` term.  The next
honest step is not to guess cancellation.  It is to compute the next independent
curvature-squared operator: `Ricci_{mu nu} Ricci^{mu nu}`.

`P1982` performs that calculation in the same diagonal Bianchi-I minisuperspace
with lapse `N(t)` retained.  It does not claim Riemann-squared, Gauss-Bonnet,
spatial-EOM, coefficient-level, or background-independence closure.

## Strict inputs

This packet continues the strict-only ADM/Bianchi-I route after `P1980/P1981`.
It uses no legacy-kernel physical role transfer.

## Ricci-component setup

For

```text
ds^2 = -N(t)^2 dt^2 + a1(t)^2 dx^2 + a2(t)^2 dy^2 + a3(t)^2 dz^2,
h_i = dot(a_i)/a_i,
sigma1 + sigma2 + sigma3 = 0,
```

the witness uses the lapse-retaining diagonal Ricci components

```text
R00 = -sum_i [dot(h_i) + h_i^2 - h_i*Ndot/N],
Rii/a_i^2 = [dot(h_i) + h_i*sum_j h_j - h_i*Ndot/N]/N^2.
```

Then

```text
Ricci^2 = R00^2/N^4 + sum_i (Rii/a_i^2)^2.
```

The FRW expression is subtracted before the lapse Euler operator is applied.

## Ricci² lapse Euler operator

The density difference is

```text
Delta L_Ricci2 = N*V*(Ricci2_BianchiI - Ricci2_FRW).
```

Because the Ricci components contain `Ndot`, the correct lapse operator is

```text
E_N(Delta L_Ricci2)
= partial_N Delta L_Ricci2 - d/dt(partial_Ndot Delta L_Ricci2).
```

The full expression is long and is therefore exported machine-readably in the
JSON witness.  The gatekeeper checks prove the following structural facts:

```text
E_N(Delta L_Ricci2) -> 0 in the isotropic limit,
E_N(Delta L_Ricci2) contains d2sigma terms,
E_N(Delta L_Ricci2) contains Nddot terms,
E_N(Delta L_Ricci2) contains shear-velocity-square terms.
```

These structures are absent from the simple Einstein-Hilbert lapse witness and
cannot be silently discarded.

## Machine-check export

Execution script:

```text
fundamental_action_reconstruction/p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.py
```

Generated witness:

```text
fundamental_action_reconstruction/generated/p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.json
```

The witness exports:

1. the lapse-retaining Ricci component formulas,
2. `Ricci2_BianchiI`, `Ricci2_FRW`, and their difference,
3. the full lapse Euler operator,
4. structural checks for isotropic vanishing, shear acceleration, lapse second
   derivative, and shear-velocity-square dependence,
5. rational numeric replay samples,
6. explicit no-false-pass boundaries.

## Theorem statement actually licensed

For the diagonal trace-free Bianchi-I minisuperspace reduction, the Ricci-squared
lapse Euler operator difference is an exactly exported symbolic expression that
vanishes in the isotropic limit but contains higher derivative shear/lapse
structures:

```text
d2sigma_i, Nddot, (dsigma_i)^2.
```

Therefore the strict curvature-squared background-independence problem cannot be
closed by the `R^2` result alone.  A coefficient-level combination with
Riemann-squared and Gauss-Bonnet remains necessary.

## False-pass boundary

This packet does **not** prove:

1. Riemann-squared ADM variation,
2. Gauss-Bonnet ADM variation,
3. strict coefficient-level curvature-squared cancellation,
4. spatial Bianchi-I EOM cancellation,
5. global background-independence closure,
6. PO2/PO3 closure,
7. BRST or Cutkosky closure,
8. `QW-2191` selector discharge,
9. ToE closure.

It only exports the Ricci-squared lapse-equation obligation.

## Progress toward ToE

This is progress because a second independent curvature-squared operator has
now been executed symbolically.  The calculation reveals that the strict
curvature-squared sector is genuinely higher-derivative in anisotropic
backgrounds.  A ToE-level closure must therefore show real cancellations or
consistent transport among all curvature-squared operators; it cannot rely on a
single simplified shear scalar.

## Next honest step

Compute the ADM/Bianchi-I lapse Euler operator for

```text
Riemann^2
```

with `N(t)` retained, then compute the Gauss-Bonnet combination and compare all
curvature-squared shear structures using the strict coefficients from
`P1972/P1853`.

## Lay explanation

Policzyliśmy kolejny trudny składnik: `Ricci^2`.  Wynik jest jeszcze bardziej
złożony niż dla `R^2`: pojawiają się przyspieszenia ścinania oraz druga pochodna
funkcji lapse.  To znaczy, że pełnej zgodności tła nie da się uzyskać prostym
„dopasowaniem znaku”.  Trzeba policzyć cały zestaw składników krzywiznowych i
sprawdzić, czy razem tworzą spójne równanie.
