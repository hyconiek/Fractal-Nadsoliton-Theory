# P1981 / S931 — Strict ADM/Bianchi-I R² Lapse-Variation Obligation Theorem Packet

Status: `OPEN_OBSTRUCTION_WITH_TRACE`  
Local theorem result: `PASS_R2_LAPSE_VARIATION_OBLIGATION_EXPORTED`  
Route: `strict_only`  
Date anchor: `2026-05-18`

## Professor-level decision

`P1980` gave a positive Einstein-Hilbert lapse result, but stopping there would
be a false pass.  The strict gravity sector also contains curvature-squared
terms.  The next honest step is therefore to keep the lapse `N(t)` and execute
the first curvature-squared ADM/Bianchi-I lapse variation.

`P1981` performs that calculation for the `R^2` density only.  It does not claim
Ricci-squared, Riemann-squared, Gauss-Bonnet, spatial-EOM, or full
background-independence closure.

## Strict inputs

This packet continues the strict-only ADM/Bianchi-I route after `P1980`.  It uses
no legacy-kernel physical role transfer.

## Ricci scalar with lapse retained

For diagonal Bianchi-I,

```text
ds^2 = -N(t)^2 dt^2 + a1(t)^2 dx^2 + a2(t)^2 dy^2 + a3(t)^2 dz^2,
h_i = dot(a_i)/a_i.
```

The lapse-retaining Ricci scalar used in the minisuperspace witness is

```text
R = 2/N^2 * [sum_i dot(h_i) + sum_i h_i^2 + sum_{i<j} h_i*h_j
             - (Ndot/N)*sum_i h_i].
```

Under trace-free shear,

```text
h1 = H + sigma1,
h2 = H + sigma2,
h3 = H - sigma1 - sigma2,
Q_shear = sigma1^2 + sigma1*sigma2 + sigma2^2.
```

The script proves

```text
R_BianchiI - R_FRW = 2*Q_shear/N^2.
```

## R² lapse Euler operator

The density difference is

```text
Delta L_R2 = N*V*(R_BianchiI^2 - R_FRW^2),
```

where `V = a1*a2*a3`.  Because `R` contains `Ndot`, the correct lapse equation is
not just `partial Delta L_R2 / partial N`; it is the Euler-Lagrange operator

```text
E_N(Delta L_R2) = partial_N Delta L_R2 - d/dt(partial_Ndot Delta L_R2).
```

The exported symbolic result is

```text
E_N(Delta L_R2)
= -12*V*(6*H^2*Q_shear - 2*H*Qdot_shear + 4*Hdot*Q_shear + Q_shear^2)/N^4.
```

It vanishes in the isotropic limit `Q_shear = Qdot_shear = 0`, but it is not a
mere constant multiple of the Einstein-Hilbert `-Q_shear` term.

## Machine-check export

Execution script:

```text
fundamental_action_reconstruction/p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.py
```

Generated witness:

```text
fundamental_action_reconstruction/generated/p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.json
```

The witness exports:

1. the lapse-retaining Ricci scalar specialization,
2. `R_BianchiI - R_FRW = 2*Q_shear/N^2`,
3. the `R^2` density difference,
4. the full lapse Euler operator including the `Ndot` term,
5. exact checks for `Qdot_shear` and `Q_shear^2` dependence,
6. rational numeric replay samples,
7. explicit no-false-pass boundaries.

## Theorem statement actually licensed

For the diagonal Bianchi-I trace-free shear minisuperspace reduction, the `R^2`
lapse equation contributes the exact anisotropic difference

```text
-12*V*(6*H^2*Q_shear - 2*H*Qdot_shear + 4*Hdot*Q_shear + Q_shear^2)/N^4.
```

Therefore curvature-squared closure cannot be inferred from the Einstein-Hilbert
`-Q_shear` success alone.  The strict coefficient-level combination with the
remaining curvature-squared operators must still be computed.

## False-pass boundary

This packet does **not** prove:

1. Ricci-squared ADM variation,
2. Riemann-squared ADM variation,
3. Gauss-Bonnet ADM variation,
4. spatial Bianchi-I EOM cancellation,
5. global background-independence closure,
6. PO2/PO3 closure,
7. BRST or Cutkosky closure,
8. `QW-2191` selector discharge,
9. ToE closure.

It only exports the `R^2` lapse-equation obligation.

## Progress toward ToE

This is progress because it prevents a premature ToE claim after the positive EH
result.  The first curvature-squared calculation shows new shear structures:
`Qdot_shear` and `Q_shear^2`.  A genuine ToE closure must show how all strict
curvature-squared terms combine, not just that Einstein-Hilbert has the right
sign.

## Next honest step

Compute the ADM/Bianchi-I lapse Euler operators for

```text
Ricci^2, Riemann^2, G_GB
```

with `N(t)` retained, then combine all curvature-squared shear terms using the
strict coefficients from `P1972/P1853`.

## Lay explanation

Po dobrym wyniku dla zwykłej grawitacji Einsteina sprawdziliśmy pierwszy
trudniejszy składnik: `R^2`.  On również widzi anizotropię, ale daje bardziej
skomplikowane wyrażenie — pojawia się pochodna ścinania i kwadrat ścinania.  To
znaczy, że pełnej teorii nie wolno zamknąć na samym składniku Einsteina; trzeba
policzyć wszystkie składniki kwadratowe i zobaczyć, czy razem tworzą spójny
wynik.
