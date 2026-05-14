# P1667 / S617 Gauge-Fixed Quadratic Kernel Diagonalization Proxy Packet

Status: `P1667_EXECUTED_QUADRATIC_KERNEL_DIAGONALIZATION_PROXY_NO_FALSE_PASS`
As of: `2026-05-14`

## Cel

Przejść z `S616` do jawnego proxy-diagonalization kroku dla sektora kwadratowego
(gauge-fixed), aby dalej domykać tor strict:

`K_strict -> coefficients -> full L -> quadratic kernel -> poles/residues -> EOM consistency`.

## Zakres

- budowa proxy-macierzy jądra dla sektorów spin-2/spin-0,
- diagonalizacja numeryczna (`2x2` operator proxy),
- eksport wartości własnych i znaku residuów,
- kontrola zgodności z warunkami nontachyonic.

## Rygor

To nadal poziom operacyjny, nie pełny dowód theorem-level.
Status globalny pozostaje `OPEN_OBLIGATION` do czasu pełnej diagonalizacji operatorowej
w formalizmie projekcyjnym na tełach zakrzywionych.

## Rekomendowany następny uczciwy krok

`S618`: połączyć diagonalizację proxy z jawnie wyprowadzonymi równaniami ruchu metryki i skalaru oraz wykonać test dwukierunkowy: `kernel->EOM` i `EOM->kernel constraints` na tej samej klasie teł.

## Omówienie dla laika

To jak rozkład drgań złożonego układu na proste tryby.
Jeśli każdy tryb ma sensowne parametry, system jest bardziej wiarygodny,
ale pełna pewność wymaga jeszcze formalnego dowodu dla całego modelu.
