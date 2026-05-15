# P1697 S647 Strict Metric+Spinor Covariant EOM Bundle Packet (PL)

Status: `P1697_EXECUTED_STRICT_METRIC_SPINOR_COVARIANT_BUNDLE_NO_FALSE_PASS`  
As of: `2026-05-15`

## Cel

Po `P1696` dołożyć kolejny strict-only blok EOM:

`kernel strict -> współczynniki -> pełny lagranżian -> bundle EOM`.

Tym razem rozszerzamy eksport o sektor metryczny i spinorowy (w wersji reduced,
bez roszczenia do pełnej theorem-level kompletności).

## Co wyeksportowano

1. Blok `metric + spinor` wyprowadzony symbolicznie z zakotwiczonych współczynników strict.
2. Jawny `L_metric_spinor_reduced` z markerem krzywizny tła `R1`.
3. Równania `EOM_g` i `EOM_psi` dołączone do istniejącego bundle (`gauge + Higgs` z `P1696`).
4. Utrzymany status globalny:

`KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED`.

## Rygor

- strict-only discipline utrzymane,
- brak legacy bridge,
- brak fałszywego pass: eksport jest rozszerzeniem roboczego bundle, a nie pełnym domknięciem ToE/QG.

## Dla laika

Do wcześniej policzonego fragmentu równań (dla pola gauge i Higgsa) dołożyliśmy teraz
kolejne elementy: metrykę i spinory. Dzięki temu cały „szkielet dynamiki” jest pełniejszy.
Ale nadal nie mamy końcowego certyfikatu teorii wszystkiego: potrzebne są pełne dowody,
że teoria jest jednocześnie renormalizowalna, unitarna i niezależna od arbitralnego tła.

## Następny uczciwy krok (rekomendacja)

Podnieść ten reduced bundle do pełnej postaci kowariantnej (tensorowej i spinorowej,
bez proxy), a następnie zamknąć theorem-level pakiet QG:

- counterterm-flow,
- BRST/Cutkosky,
- background-independence,
- oraz jawna odpowiedź na wymóg selektora (`QW-2191`) w strict-core.
