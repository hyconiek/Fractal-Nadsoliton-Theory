# RELEASE 5.0 TEXTBOOK EDITION (EN + PL)

**Version:** 5.0.0  
**Date:** 2026-03-04  
**Branch:** `main`

**Note:** filename remains `RELEASE_4_9_TEXTBOOK_EN_PL.md` for continuity, but content below is the current Release 5.0 textbook snapshot.

---

## ENGLISH VERSION

## 1) What Release 5.0 Is

Release 5.0 is the first branch state where the strict internal SM+GR package path is closed end-to-end under a locked no-scan/no-retune methodology.

In simple words:
- the theory is no longer only “strong internal partial,”
- in the current audited chain it reaches full internal strict closure.

## 2) One-Page Status (High-School Level)

### 2.1 What is now strong

1. One frozen kernel supports mass + flavor + GW sectors together.
2. Full package closure map passes (`QW-2069`).
3. Full radiative program gate passes (`QW-2070`).
4. Full precision closure gate passes (`QW-2071`).
5. Missing-14 strict frontier is fully closed (`QW-2081`).
6. CKM CP strict refinement now passes (`QW-2097`).
7. Defect sweep reports no critical failures (`QW-2094`, `130` checks, `0` failed).

### 2.2 What is still missing

1. Independent external multiteam replication.
2. Community-level acceptance based on independent reproduction.

### 2.3 Bottom line

- **Internal strict closure:** achieved in current audited chain.
- **Final community-confirmed ToE:** not yet; independent external replication is still required.

## 3) Core Idea in Textbook Language

The Nadsoliton program assumes that particles, interactions, and spacetime patterns are different manifestations of one deeper informational-dynamical substrate.

A common structural kernel is used across sectors:

$$
K(d)=\frac{\cos(\omega d+\phi)}{1+\beta d^{\eta}}
$$

For the working strict chain:
- $\omega = 0.18575$
- $\phi = 0.16250$
- $\beta = 1.00000$
- $\eta = 1.80000$

Interpretation:
- cosine term gives oscillatory structure,
- denominator damps amplitude with distance,
- same functional form is reused in all audited sectors.

## 4) Release 5.0 Closure Gates

### 4.1 QW-2069 (full SM+GR derivation package)

- Verdict: `FULL_SM_GR_DERIVATION_PACKAGE_PASS`
- Coverage:
  - `n_total_registry = 32`
  - `n_derived_strict_internal = 30`
  - `n_definition_constants = 2`
  - `n_missing = 0`
  - `n_strict_unresolved = 0`

### 4.2 QW-2070 (full radiative program)

- Verdict: `FULL_RADIATIVE_PROGRAM_PASS`
- Channels:
  - implemented: `7/7`
  - closure-ready: `7/7`
  - missing: `0`

### 4.3 QW-2071 (full precision closure gate)

- Verdict: `SM_GR_FULL_PRECISION_CLOSURE_PASS`
- Pass count: `6/6`

### 4.4 QW-2081 (missing-14 strict frontier)

- Verdict: `MISSING14_STRICT_RIGOR_FRONTIER_PASS_ALL_CLOSED`
- Unresolved in tracked scope: `0`

### 4.5 QW-2097 (CKM CP strict refinement)

- Verdict: `CKM_CP_TARGET_REFINEMENT_GATE_PASS_STRICT`
- Pass count: `6/6`
- Selected deterministic phase scheme: `kernel_cp_extension_v1`

### 4.6 QW-2094 (strict rigor defect sweep)

- Verdict: `STRICT_RIGOR_DEFECT_SWEEP_PASS_NO_CRITICAL_DEFECTS`
- Checks: `130`
- Failed: `0`

## 5) Empirical GW Status (Methodological)

- External blind chain remains strong:
  - `QW-2016`: `TRIAD_BLIND_EXTERNAL_VALIDATION_PASS_STRONG`
  - `QW-2017`: `BETA_OBSERVABLE_BLIND_EXTERNAL_INTERVENTION_PASS`
- Preregistered GW holdout branch:
  - `QW-2078` autocollector metrics satisfy hard thresholds,
  - `QW-2077`: GW branch supported, global verdict mixed because PMNS/cosmology still require complete external updates.

## 6) What Release 5.0 Means Scientifically

It means:
1. strict internal first-principles closure is achieved in the audited package chain,
2. no direct missing and no strict-unresolved entries remain in package scope.

It does not mean:
1. automatic community-level final ToE acceptance,
2. immunity from future falsification.

## 7) Exact Next Step

- Execute independent external multiteam replication on frozen manifests/artifacts,
- publish independent rerun reports.

## 8) Main Artifacts

- `report_qw2069_full_sm_gr_derivation_package.json`
- `report_qw2070_full_radiative_program_baseline.json`
- `report_qw2071_sm_gr_full_precision_closure_gate.json`
- `report_qw2081_missing14_strict_rigor_frontier.json`
- `report_qw2097_ckm_cp_target_refinement_gate.json`
- `report_qw2094_strict_rigor_defect_sweep.json`

---

## WERSJA POLSKA

## 1) Czym jest Release 5.0

Release 5.0 to pierwszy stan gałęzi, w którym ścieżka ścisłego domknięcia pakietu SM+GR przechodzi od początku do końca w metodologii locked no-scan/no-retune.

Po ludzku:
- teoria nie jest już tylko „mocnym częściowym domknięciem wewnętrznym”,
- w aktualnym audytowanym łańcuchu osiąga pełne ścisłe domknięcie wewnętrzne.

## 2) Status w skrócie (licealnie)

### 2.1 Co jest teraz mocne

1. Jedno zamrożone jądro wspiera razem sektor masy, flavor i GW.
2. Mapa domknięcia pełnego pakietu przechodzi (`QW-2069`).
3. Pełny program radiacyjny przechodzi (`QW-2070`).
4. Bramka pełnego domknięcia precyzyjnego przechodzi (`QW-2071`).
5. Frontier missing-14 jest domknięty (`QW-2081`).
6. Ścisłe doprecyzowanie CKM CP przechodzi (`QW-2097`).
7. Sweep usterek rygoru nie wykazuje krytycznych błędów (`QW-2094`, `130` kontroli, `0` failed).

### 2.2 Czego nadal brakuje

1. Niezależnej replikacji zewnętrznej multiteam.
2. Potwierdzenia społecznościowego opartego na niezależnym odtworzeniu wyników.

### 2.3 Wniosek

- **Ścisłe domknięcie wewnętrzne:** osiągnięte w aktualnym łańcuchu.
- **Finalne ToE potwierdzone przez społeczność:** jeszcze nie; potrzebna niezależna replikacja zewnętrzna.

## 3) Rdzeń teorii w języku podręcznikowym

Program nadsolitonowy zakłada, że cząstki, oddziaływania i wzory czasoprzestrzeni są różnymi przejawami jednej głębszej struktury informacyjno-dynamicznej.

Wspólne jądro strukturalne:

$$
K(d)=\frac{\cos(\omega d+\phi)}{1+\beta d^{\eta}}
$$

Dla aktualnego łańcucha strict:
- $\omega = 0.18575$
- $\phi = 0.16250$
- $\beta = 1.00000$
- $\eta = 1.80000$

Interpretacja:
- cosinus daje strukturę oscylacyjną,
- mianownik tłumi amplitudę z odległością,
- ta sama postać funkcji jest używana między sektorami.

## 4) Bramki domknięcia Release 5.0

### 4.1 QW-2069 (pełny pakiet derivacji SM+GR)

- Werdykt: `FULL_SM_GR_DERIVATION_PACKAGE_PASS`
- Pokrycie:
  - `n_total_registry = 32`
  - `n_derived_strict_internal = 30`
  - `n_definition_constants = 2`
  - `n_missing = 0`
  - `n_strict_unresolved = 0`

### 4.2 QW-2070 (pełny program radiacyjny)

- Werdykt: `FULL_RADIATIVE_PROGRAM_PASS`
- Kanały:
  - implemented: `7/7`
  - closure-ready: `7/7`
  - missing: `0`

### 4.3 QW-2071 (bramka pełnego domknięcia precyzyjnego)

- Werdykt: `SM_GR_FULL_PRECISION_CLOSURE_PASS`
- Wynik: `6/6`

### 4.4 QW-2081 (strict frontier missing-14)

- Werdykt: `MISSING14_STRICT_RIGOR_FRONTIER_PASS_ALL_CLOSED`
- Unresolved w śledzonym zakresie: `0`

### 4.5 QW-2097 (CKM CP strict refinement)

- Werdykt: `CKM_CP_TARGET_REFINEMENT_GATE_PASS_STRICT`
- Wynik: `6/6`
- Wybrany deterministyczny schemat fazy: `kernel_cp_extension_v1`

### 4.6 QW-2094 (sweep usterek rygoru)

- Werdykt: `STRICT_RIGOR_DEFECT_SWEEP_PASS_NO_CRITICAL_DEFECTS`
- Kontrole: `130`
- Failed: `0`

## 5) Status empiryczny GW (metodologicznie)

- Łańcuch blind external pozostaje silny:
  - `QW-2016`: `TRIAD_BLIND_EXTERNAL_VALIDATION_PASS_STRONG`
  - `QW-2017`: `BETA_OBSERVABLE_BLIND_EXTERNAL_INTERVENTION_PASS`
- Gałąź prerejestrowana GW holdout:
  - metryki z `QW-2078` spełniają twarde progi,
  - `QW-2077`: gałąź GW jest supported, a globalny wynik jest mixed, bo PMNS/kosmologia wymagają pełnych aktualizacji danych zewnętrznych.

## 6) Co naukowo oznacza Release 5.0

Znaczy:
1. ścisłe domknięcie first-principles jest osiągnięte w audytowanym łańcuchu pakietowym,
2. w zakresie pakietu nie ma już pozycji direct missing ani strict-unresolved.

Nie znaczy:
1. automatycznego finalnego uznania ToE przez społeczność,
2. odporności na przyszłą falsyfikację.

## 7) Najbliższy krok

- Wykonać niezależną replikację zewnętrzną multiteam na zamrożonych manifestach/artefaktach,
- opublikować niezależne raporty rerun.

## 8) Główne artefakty

- `report_qw2069_full_sm_gr_derivation_package.json`
- `report_qw2070_full_radiative_program_baseline.json`
- `report_qw2071_sm_gr_full_precision_closure_gate.json`
- `report_qw2081_missing14_strict_rigor_frontier.json`
- `report_qw2097_ckm_cp_target_refinement_gate.json`
- `report_qw2094_strict_rigor_defect_sweep.json`
