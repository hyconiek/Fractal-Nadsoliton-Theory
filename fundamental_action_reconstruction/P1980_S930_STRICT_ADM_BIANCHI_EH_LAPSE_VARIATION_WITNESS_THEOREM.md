# P1980 / S930 — Strict ADM/Bianchi-I Einstein-Hilbert Lapse Variation Witness Theorem Packet

Status: `OPEN_OBSTRUCTION_WITH_TRACE`  
Local theorem result: `PASS_EH_LAPSE_SHEAR_TERM_DERIVED`  
Route: `strict_only`  
Date anchor: `2026-05-18`

## Professor-level decision

`P1979` was too weak to be the final word: it was only a current-export registry
audit.  It correctly said that no ADM/Bianchi-I lapse/shear certificate was
already exported, but a token audit cannot decide whether the covariant gravity
term actually contains the needed shear-energy contribution.

The next honest step is therefore to stop auditing strings and perform the first
real ADM/Bianchi-I lapse variation.  `P1980` does this for the Einstein-Hilbert
bulk term only.

This is a positive local result, not a global closure result.

## Strict inputs

This packet uses the strict-side gravity sector already present in the
`P1846/P1907` chain and the post-`P1979` obligation to keep the lapse `N(t)`
until variation.  It uses no legacy-kernel physical role transfer.

## ADM/Bianchi-I setup

The diagonal Bianchi-I metric is

```text
ds^2 = -N(t)^2 dt^2 + a1(t)^2 dx^2 + a2(t)^2 dy^2 + a3(t)^2 dz^2.
```

For Bianchi-I, the spatial curvature is zero.  Let

```text
v_i = dot(a_i)/a_i,
K^i_i = v_i/N,
V = a1*a2*a3.
```

The Einstein-Hilbert ADM bulk density, with the overall constant removed, is

```text
L_EH_ADM = N*V*(K_ij*K^ij - K^2).
```

For the diagonal ansatz, SymPy reduces it to

```text
L_EH_ADM = -2*V*(v1*v2 + v1*v3 + v2*v3)/N.
```

Varying with respect to the lapse gives

```text
dL_EH_ADM/dN = 2*V*(v1*v2 + v1*v3 + v2*v3)/N^2.
```

Removing the positive conventional factor `2*V/N^2`, the normalized lapse
constraint is

```text
v1*v2 + v1*v3 + v2*v3.
```

## Trace-free shear specialization

Set

```text
H1 = H + sigma1,
H2 = H + sigma2,
H3 = H - sigma1 - sigma2.
```

Then

```text
H1*H2 + H1*H3 + H2*H3
= 3*H^2 - sigma1^2 - sigma1*sigma2 - sigma2^2.
```

Define

```text
Q_shear = sigma1^2 + sigma1*sigma2 + sigma2^2.
```

Therefore

```text
G_nn^EH(Bianchi-I) - G_nn^EH(FRW) = -Q_shear.
```

This exactly matches the `P1974` `G00`/energy residual sign that `P1978/P1979`
identified as necessary.

## Machine-check export

Execution script:

```text
fundamental_action_reconstruction/p1980_s930_strict_adm_bianchi_eh_lapse_variation_witness.py
```

Generated witness:

```text
fundamental_action_reconstruction/generated/p1980_s930_strict_adm_bianchi_eh_lapse_variation_witness.json
```

The witness exports:

1. the ADM lapse-retaining Einstein-Hilbert bulk density,
2. the symbolic lapse derivative,
3. the normalized constraint `v1*v2 + v1*v3 + v2*v3`,
4. the trace-free shear specialization,
5. exact and SciPy eigenvalue checks for positive-definiteness of `Q_shear`,
6. rational numeric replay samples,
7. explicit no-false-pass boundaries.

## Theorem statement actually licensed

For diagonal Bianchi-I in the ADM Einstein-Hilbert bulk term, lapse variation
derives the normalized energy constraint

```text
H1*H2 + H1*H3 + H2*H3.
```

Under trace-free shear this equals

```text
3*H^2 - Q_shear.
```

Hence the Einstein-Hilbert lapse equation derives the Bianchi-I minus FRW shear
energy residual

```text
-Q_shear.
```

## Correction to the P1979 reading

`P1979` remains valid only as a current-export certificate audit.  `P1980` shows
that the absence of a textual provider certificate in the registry does **not**
mean the Einstein-Hilbert covariant term lacks the shear-energy contribution.
After actual ADM variation, the EH term does carry the correct `-Q_shear` sign.

## False-pass boundary

This packet does **not** prove:

1. full curvature-squared ADM variation,
2. spatial Bianchi-I EOM transport,
3. global background-independence closure,
4. PO2/PO3 closure,
5. BRST or Cutkosky closure,
6. `QW-2191` selector discharge,
7. ToE closure.

It is one positive local lapse-equation witness for the Einstein-Hilbert term.

## Progress toward ToE

This is stronger progress than `P1979`: we moved from a superficial export audit
to an actual symbolic variational calculation.  The result is encouraging for
ToE because the basic strict gravity term already produces the required shear
energy sign.  But ToE remains open because the strict program still needs the
curvature-squared terms and the spatial equations, not only the EH lapse
constraint.

## Next honest step

Extend the ADM/Bianchi-I lapse variation from the Einstein-Hilbert bulk term to
strict curvature-squared terms:

```text
R^2, Ricci^2, Riemann^2, G_GB.
```

Keep `N(t)` until variation, export the resulting shear contributions, and test
whether they preserve or obstruct the `P1974/P1975` component obligations.

## Lay explanation

W poprzednim kroku tylko sprawdzaliśmy, czy w repozytorium jest gotowy opis
brakującego składnika.  Teraz zrobiliśmy prawdziwy rachunek dla podstawowego
składnika grawitacji Einsteina.  Wynik jest dobry: sama geometria anizotropowa
daje dokładnie brakujący znak energii `-Q_shear`.  To jednak nie kończy teorii,
bo trudniejsze składniki kwadratowe w krzywiźnie i równania przestrzenne nadal
czekają na obliczenie.
