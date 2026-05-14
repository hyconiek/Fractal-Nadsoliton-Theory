# P1651 / S601 — Strict det(J) lower-bound certificate

## Cel
Wykonać pierwszy discharge dla O1 z P1650: dostarczyć lokalny certyfikat
niezerowości Jacobianu mapy projektorowej na badanym regionie.

## Zakres
- obliczenie `eps_hat = min |detJ|` z regionalnego witnessu,
- jawny status `PARTIAL` (certyfikat obliczeniowy, nie pełny dowód analityczny),
- utrzymanie `OPEN` dla global closure.

## Wyjście
- `generated/p1651_s601_strict_detj_lower_bound_certificate_summary.json`
