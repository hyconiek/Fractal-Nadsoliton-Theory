# PLAN: KERNEL-DERIVED NON-ANCHOR INPUTS (QW-2085..QW-2087)

Date: 2026-03-04

## 1) Cel
- zbudować wejścia `kernel_derived` dla T1 (`g_f`, `m_z`, `alpha_s_mz`) bez retune i bez ręcznego dopasowywania pod targety.
- wykonać bramki QW-2085/QW-2086/QW-2087 na tych wejściach.

## 2) Zasady rygoru
- brak skanów i brak optymalizacji parametrów pod obserwacje,
- stałe liczby w formułach są zamrożone w tym planie (`5.0`, `0.02`),
- wejścia budowane tylko z raportów:
  - `report_qw2063_derivational_reconstruction_shared_flavor_basis.json` (masy z łańcucha derivacyjnego),
  - `report_qw2064_micro_derived_renormalization_constants_gate.json` (kernel + mikrostałe),
  - stałe fizyczne konwersji jednostek (`hbar[GeV*s]`).

## 3) Zmienne źródłowe
- `m_top_pred`, `m_bottom_pred`, `m_mu_pred` z QW-2063.
- `omega`, `phi`, `beta_uv`, `z_beta`, `delta_eta` z QW-2064.
- `alpha_geo = 4 ln 2`.

## 4) Zamrożone formuły kernel-derived
### 4.1 Dla QW-2085 (`g_f`)
- `tau_mu_kernel = (2π/omega) * (m_top/m_mu)^5 * (hbar/m_mu) * (1 + delta_eta) / (1 + z_beta/100)`
- `delta_q_kernel = beta_uv * (1 + delta_eta/2)`
- wejścia: `m_mu_gev = m_mu_pred`, `tau_mu_s = tau_mu_kernel`, `delta_q = delta_q_kernel`, provenance=`kernel_derived`.

### 4.2 Dla QW-2086 (`m_z`)
- `mw_pole_kernel = m_top * sqrt((beta_uv * z_beta) / (5.0 * (1 + delta_eta/10)))`
- `sin2_theta_w_eff_kernel = (alpha_geo/12) * (1 + beta_uv*(delta_eta - 1))`
- `delta_r_full_kernel = beta_uv * delta_eta + (omega*phi)/2`
- provenance wszystkich wejść = `kernel_derived`.

### 4.3 Dla QW-2087 (`alpha_s_mz`)
- `mu0 = m_bottom`
- `alpha_s_mu0 = 1 / (ln(m_top/m_bottom) + delta_eta)`
- punkty walidacyjne dla `mu in {2, 10, 173} GeV`:
  - `alpha_pred(mu)` z 1-loop running od `(mu0, alpha_s_mu0)`,
  - `alpha_obs(mu) = alpha_pred(mu) * (1 + beta_uv * cos(omega*ln(mu/mu0) + phi))`
  - `sigma_total(mu) = 0.02 * |alpha_obs(mu)|`
  - provenance=`kernel_derived`.

## 5) Kryteria wykonania planu
1. wygenerować pliki:
   - `t1_nonanchor_observables_input_qw2085_2086.json`
   - `t1_nonanchor_alpha_s_input_qw2087.json`
2. uruchomić:
   - `QW_2085_GF_NONANCHOR_LIFETIME_GATE.py`
   - `QW_2086_MZ_NONANCHOR_EW_POLE_GATE.py`
   - `QW_2087_ALPHA_S_NONANCHOR_BOUNDARY_GATE.py`
3. przeliczyć integrację:
   - `QW_2069_FULL_SM_GR_DERIVATION_PACKAGE.py`
   - `QW_2070_FULL_RADIATIVE_PROGRAM_BASELINE.py`
   - `QW_2071_SM_GR_FULL_PRECISION_CLOSURE_GATE.py`

## 6) Oczekiwany typ wyniku (bez obietnicy PASS)
- wynik ma być deterministyczny i audytowalny,
- PASS jest dozwolony tylko jeśli bramki przejdą warunki bez obejść,
- FAIL jest wynikiem poprawnym naukowo, jeśli warunki nie są spełnione.

