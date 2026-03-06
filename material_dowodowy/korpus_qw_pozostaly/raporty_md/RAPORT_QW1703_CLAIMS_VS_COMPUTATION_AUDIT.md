# RAPORT QW-1703: CLAIMS VS COMPUTATION AUDIT

- Data UTC: 2026-03-02T15:11:43.743660+00:00
- Źródła: `README.md`, `FINAL_SUMMARY_TEORIA_WSZYSTKIEGO.md`
- Daty źródeł (mtime UTC): 2025-12-30T01:45:04.924618+00:00 ; 2025-11-19T03:25:14.085528+00:00

## 1) Recomputacja kluczowych wzorów
- sin²θ_W = 4ln2/12 = 0.231049, exp=0.231220, błąd=0.0739%
- α⁻¹ (4ln2/(2β), β=0.01) = 138.629436, exp=137.035999, błąd=1.1628%
- α⁻¹ (wariant verify_values) = 137.243142, exp=137.035999, błąd=0.1512%
- 4ln2 vs φ√3: 2.772589 vs 2.802517, luka=1.0679%

## 2) Test wzoru masowego (Q, γ=1.52)
- Top: Q=0, pred=173000.0000 MeV, exp=173000.0000 MeV, błąd=0.00%
- Bottom: Q=7, pred=4330.7840 MeV, exp=4180.0000 MeV, błąd=3.61%
- Tau: Q=9, pred=1510.0834 MeV, exp=1776.9000 MeV, błąd=15.02%
- Charm: Q=9, pred=1510.0834 MeV, exp=1270.0000 MeV, błąd=18.90%
- Muon: Q=14, pred=108.4144 MeV, exp=105.7000 MeV, błąd=2.57%
- Electron: Q=24, pred=0.5588 MeV, exp=0.5110 MeV, błąd=9.35%

## 3) Kontrast deklaracji vs obliczeń
- Wzmianki `EXACT/0.00%` (README + FINAL_SUMMARY): 11
- Wzmianki `fitting/calibration` (README + FINAL_SUMMARY): 5

## 4) Artefakty
- JSON szczegółowy: `report_qw1703_claims_vs_computation_audit.json`
