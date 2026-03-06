# RAPORT QW-1708: REPAIR CONSTANTS ESTIMATOR

- Data UTC: 2026-03-02T15:26:05.202349+00:00

## 1) Poprawka Weinberga
- sin²θ_W(base)=4ln2/12 = 0.231049060
- sin²θ_W(exp) = 0.231220000
- wymagane δ_W = +7.398420630e-04

## 2) Poprawka α_EM
- α⁻¹(base)=4ln2/(2β)*(1-β) = 137.243141751
- α⁻¹(exp) = 137.035999000
- wymagane δ_vac = -1.509312219e-03

## 3) Poprawki masowe Δ_a
- Top: Q=0, pred=173000.0000 MeV, exp=173000.0000 MeV, Δ_a=+0.000000, błąd_base=0.00%
- Bottom: Q=7, pred=4330.7840 MeV, exp=4180.0000 MeV, Δ_a=-0.035437, błąd_base=3.61%
- Tau: Q=9, pred=1510.0834 MeV, exp=1776.9000 MeV, Δ_a=+0.162705, błąd_base=15.02%
- Charm: Q=9, pred=1510.0834 MeV, exp=1270.0000 MeV, Δ_a=-0.173148, błąd_base=18.90%
- Muon: Q=14, pred=108.4144 MeV, exp=105.7000 MeV, Δ_a=-0.025356, błąd_base=2.57%
- Electron: Q=24, pred=0.5588 MeV, exp=0.5110 MeV, Δ_a=-0.089428, błąd_base=9.35%

## 4) Agregat
- średnia |Δ_a| (bez top): 0.097215
- max |Δ_a| (bez top): 0.173148

## Artefakty
- JSON szczegółowy: `report_qw1708_repair_constants_estimator.json`
