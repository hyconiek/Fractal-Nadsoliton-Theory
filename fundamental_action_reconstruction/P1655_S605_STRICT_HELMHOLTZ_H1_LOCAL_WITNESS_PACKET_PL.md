# P1655 / S605 — Strict Helmholtz H1 local witness for reverse-chain start

## Cel
Wykonać pierwszy twardy krok z P1654: zbudować lokalny witness H1 (Helmholtz)
dla fragmentu pełnego układu sprzężonego, jako start toru odwrotnego `EOM -> L_total`.

## Zakres
- strict-only (`legacy_bridge_used=false`),
- bez deklaracji domknięcia H2/H3/H4,
- bez false-pass: eksport tylko lokalnego kroku H1.

## Konstrukcja
- Rozpatrujemy sektor `phi + H` z mieszaniem `lambda_{phiH} phi^2(H^†H)`.
- Liczymy lokalnie symetrię pochodnych krzyżowych operatora Eulera-Lagrange'a
  jako warunek integrabilności typu Helmholtza H1.
- Wynik H1 lokalny zapisujemy jako `PARTIAL`, bo pełny układ (gauge+fermion+metryka)
  nadal wymaga osobnych witnessów theorem-level.

## Wyjście
- `generated/p1655_s605_strict_helmholtz_h1_local_witness_summary.json`
