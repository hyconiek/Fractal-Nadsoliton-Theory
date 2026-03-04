# PLAN: KERNEL-DERIVED T2 NON-ANCHOR INPUTS (QW-2088..QW-2089)

Date: 2026-03-04

## 1) Cel
- zbudować wejścia `kernel_derived` dla T2:
  - `m_up`, `m_down`, `m_strange` (QW-2088),
  - `m_h` (QW-2089),
- wykonać bramki strict bez retune i bez scan.

## 2) Zasady rygoru
- brak skanów i brak optymalizacji pod targety,
- formuły zamrożone w tym planie,
- źródła tylko z wcześniejszych raportów:
  - `report_qw2063_derivational_reconstruction_shared_flavor_basis.json`,
  - `report_qw2064_micro_derived_renormalization_constants_gate.json`,
  - plus stałe matematyczne/fizyczne.

## 3) Zmienne źródłowe
- `m_top_pred` z QW-2063,
- `omega`, `phi`, `beta_uv`, `z_beta`, `delta_eta` z QW-2064,
- `alpha_geo = 4 ln 2`.

## 4) Zamrożone formuły kernel-derived
### 4.1 Light-quark chain (QW-2088)
- `m_u = m_top * beta_uv^(2.5) * (1 + delta_eta/10)`
- `m_d = m_u * (1 + delta_eta)`
- `m_s = m_d * beta_uv^(-0.6) * (1 + omega)`

Niezależne cross-checki hadronowe (external constraints):
- `(m_d/m_u)_obs = 2.16 +/- 0.40`
- `(m_s/m_d)_obs = 19.9 +/- 3.0`

### 4.2 Higgs chain (QW-2089)
- `lambda_eff = (alpha_geo / (8*pi)) * (1 + delta_eta/10 + beta_uv + omega/2)`
- `v_eff = (m_top / (omega * sqrt(z_beta/100))) * beta_uv^(0.25) * (1 - delta_eta/10)`
- `m_h = sqrt(2*lambda_eff) * v_eff`

Niezależne punkty walidacji lambda (external constraints):
- `lambda(125 GeV) = 0.129 +/- 0.020`
- `lambda(500 GeV) = 0.120 +/- 0.030`

## 5) Kryteria wykonania
1. wygenerować wejścia:
   - `t2_nonanchor_light_quark_input_qw2088.json`
   - `t2_nonanchor_higgs_input_qw2089.json`
2. uruchomić:
   - `QW_2088_LIGHT_QUARK_MASS_NONANCHOR_GATE.py`
   - `QW_2089_HIGGS_SELFCOUPLING_STRICT_GATE.py`
   - `QW_2096_T2_NONANCHOR_STRICT_GATE.py`
3. przeliczyć integrację:
   - `QW_2069_FULL_SM_GR_DERIVATION_PACKAGE.py`
   - `QW_2071_SM_GR_FULL_PRECISION_CLOSURE_GATE.py`
   - `QW_2094_STRICT_RIGOR_DEFECT_SWEEP.py`

## 6) Oczekiwany typ wyniku
- wynik deterministyczny i audytowalny,
- PASS tylko przy spełnieniu hard checks (provenance + tolerance + consistency),
- FAIL jest wynikiem poprawnym naukowo, jeśli warunki nie są spełnione.

