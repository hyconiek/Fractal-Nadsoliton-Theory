# P1704 S654 Strict Nonproxy Metric+Spinor EOM Bundle Export Scaffold Packet (PL)

Status: `P1704_EXECUTED_STRICT_NONPROXY_METRIC_SPINOR_SCAFFOLD_NO_FALSE_PASS`  
As of: `2026-05-15`

## Cel

Wykonać pierwszy realny krok wykonawczy po `P1703`: przygotować jawny scaffold
nieproxy eksportu równania metrycznego i spinorowego z pełnego lagranżianu.

## Co wyeksportowano

1. Schemat pól i celów EOM dla sektora metrycznego (`g_{μν}`) i spinorowego (`psi, psibar, tetrada, spin-connection`).
2. Listy wymaganych obiektów wariacyjnych do policzenia w kolejnym kroku.
3. Jawne sprzężenie do filarów `Helmholtz_nonproxy` i `BRST_nonproxy`.

## Rygor

- strict-only,
- bez bridge do legacy,
- bez fałszywego pass,
- status globalny nadal `KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED`.

## Dla laika

To etap przygotowujący „duże równania” w pełnej wersji bez uproszczeń. Dzięki temu
kolejny krok będzie już bezpośrednio liczył właściwe wzory dla grawitacji i fermionów,
a nie tylko planował, co zrobić.

## Następny uczciwy krok (rekomendacja)

Wyeksportować jawne równania wariacyjne metryczne i spinorowe (z połączeniem spinowym)
w postaci obliczalnej, aby zasilić dowody Helmholtza i BRST.
