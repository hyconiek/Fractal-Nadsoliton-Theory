# P1984 / S934 — Strict ADM/Bianchi-I Gauss-Bonnet Lapse-Cancellation Witness Theorem Packet

Status: `OPEN_OBSTRUCTION_WITH_TRACE`  
Local theorem result: `PASS_GB_LAPSE_CANCELLATION_WITNESS`  
Route: `strict_only`  
Date anchor: `2026-05-18`

## Professor-level decision

`P1981`, `P1982`, and `P1983` exported the three elementary curvature-squared
lapse obligations: `R^2`, `Ricci^2`, and `Riemann^2`.  The next honest step is
not to jump directly to ToE closure.  The next honest step is to test the known
geometric combination

```text
G_GB = Riemann^2 - 4*Ricci^2 + R^2
```

in the same ADM/Bianchi-I lapse framework.

`P1984` performs that combination and proves exact cancellation of the lapse
Euler operator in this minisuperspace witness.  It does not claim spatial-EOM,
strict coefficient-level, background-independence, or ToE closure.

## Strict inputs

This packet continues the strict-only ADM/Bianchi-I route after
`P1981/P1982/P1983`.  It uses no legacy-kernel physical role transfer.

## Combination being tested

The witness recomputes the three lapse-retaining invariant differences:

```text
Delta R^2,
Delta Ricci^2,
Delta Riemann^2,
```

where `Delta` means Bianchi-I minus FRW in the trace-free diagonal shear class.
It then forms

```text
Delta G_GB = Delta Riemann^2 - 4*Delta Ricci^2 + Delta R^2.
```

The density difference is

```text
Delta L_GB = N*V*Delta G_GB.
```

Because the individual invariants contain `Ndot`, the correct lapse operator is
again

```text
E_N(Delta L_GB) = partial_N Delta L_GB - d/dt(partial_Ndot Delta L_GB).
```

## Result

The Gauss-Bonnet density difference is not the zero polynomial before variation.
However, the lapse Euler operator cancels exactly:

```text
E_N(Delta L_GB) = 0.
```

Thus the higher-derivative shear/lapse structures exposed separately by
`P1981/P1982/P1983` cancel in the Gauss-Bonnet lapse equation, as expected for a
four-dimensional topological combination at this minisuperspace lapse level.

## Machine-check export

Execution script:

```text
fundamental_action_reconstruction/p1984_s934_strict_adm_bianchi_gauss_bonnet_lapse_cancellation_witness.py
```

Generated witness:

```text
fundamental_action_reconstruction/generated/p1984_s934_strict_adm_bianchi_gauss_bonnet_lapse_cancellation_witness.json
```

The witness exports:

1. the recomputed `R^2`, `Ricci^2`, and `Riemann^2` differences,
2. the Gauss-Bonnet density difference,
3. `dL/dN`, `dL/dNdot`, and `d/dt(dL/dNdot)`,
4. the exact identity `E_N(Delta L_GB)=0`,
5. numeric replay samples confirming zero lapse Euler output,
6. explicit no-false-pass boundaries.

## Theorem statement actually licensed

For the diagonal trace-free Bianchi-I ADM minisuperspace witness,

```text
E_N[N*V*(Delta Riemann^2 - 4*Delta Ricci^2 + Delta R^2)] = 0.
```

This proves Gauss-Bonnet lapse-Euler cancellation in the current minisuperspace
calculation.

## False-pass boundary

This packet does **not** prove:

1. strict coefficient-level cancellation for non-GB `R^2/Ricci^2/Riemann^2`
   terms,
2. spatial Bianchi-I EOM cancellation,
3. global background-independence closure,
4. PO2/PO3 closure,
5. BRST or Cutkosky closure,
6. `QW-2191` selector discharge,
7. ToE closure.

It only proves the Gauss-Bonnet lapse-Euler cancellation in this ADM/Bianchi-I
minisuperspace witness.

## Progress toward ToE

This is genuine progress: a hard consistency expectation was checked rather than
assumed.  The difficult higher-derivative shear terms from `R^2`, `Ricci^2`, and
`Riemann^2` cancel in the Gauss-Bonnet lapse combination.  But ToE remains open
because the non-GB strict curvature-squared coefficients still have to be
combined and spatial equations remain unclosed.

## Next honest step

Use `P1981/P1982/P1983/P1984` to form the strict non-GB curvature-squared lapse
operator with coefficients

```text
a_R2, a_Ric2, a_Riem2
```

from `P1972/P1853`, then test whether the remaining shear structures cancel or
persist.

## Lay explanation

Z trzech trudnych składników zbudowaliśmy specjalną kombinację Gaussa-Bonneta.
W równaniu lapse wszystkie trudne wyrazy z wyższymi pochodnymi kasują się do
zera.  To bardzo dobry test spójności geometrii.  Ale to nie kończy teorii:
pozostałe kombinacje nie-Gauss-Bonnetowe ze ścisłymi współczynnikami nadal mogą
zostawić resztki i trzeba je policzyć osobno.
