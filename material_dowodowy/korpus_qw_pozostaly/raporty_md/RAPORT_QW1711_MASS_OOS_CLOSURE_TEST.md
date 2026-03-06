# RAPORT QW-1711: MASS OOS CLOSURE TEST

- Data UTC: 2026-03-02T15:33:20.945273+00:00
- Werdykt: **MASS_SECTOR_NOT_OOS_CLOSED**

## 1) Model korekcyjny
- Delta = l1*(Q/24) + l2*sector + l3*(gen-2)
- l1=-0.096670, l2=-0.031035, l3=+0.023793

## 2) Wynik out-of-sample
- Train mean error: 0.00%
- Test mean error: 9.00%
- Generalization gap: 9.00 p.p.

## 3) Test holdout
- Top: base=0.00%, corrected=2.41%
- Tau: base=15.02%, corrected=13.42%
- Charm: base=18.90%, corrected=11.17%

## Interpretacja
- Jeśli test nie przechodzi, obecny sektor mas nie jest domknięty predykcyjnie i wymaga bogatszych niezmienników topologicznych.

## Artefakty
- JSON szczegółowy: `report_qw1711_mass_oos_closure_test.json`
