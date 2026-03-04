# Independent Check Guide (EN + PL)

## EN: How other people can verify this independently

## Goal
Run a truly independent confirmatory check with the locked protocol, without retuning.

## 1) Freeze reference
Use and publish:
- commit hash,
- lock file SHA256,
- exact Python version and OS.

Current lock:
- `external_confirmatory_v2/independent_multiteam_lock_qw2053/protocol_lock_qw2053.json`
- SHA256: `2b49385e38a9985b7cf031992da821c3da23a6c0390e4643d56b111adf07c880`

## 2) Independent team setup (must be separate)
- Separate machine/account (not your development environment).
- Fresh clone of repository.
- No local edits before run.

## 3) External data policy
- Do not use local private copies hidden from protocol.
- Download only from documented public sources:
  - `DATA_SOURCES_EXTERNAL_DOWNLOADS.md`
- Verify checksums before running gates.
- If NANOGrav archive is kept outside repository root, pass explicit path in QW-2014:
  - `python3 QW_2014_TRUE_EXTERNAL_BETA_CHANNEL_AUTOCOLLECTOR_V2.py --nanograv-archive <path_to_NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz>`

## 4) Required execution order (locked)
Run exactly:
1. `python3 QW_2048_SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION.py`
2. `python3 QW_2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.py`
3. `python3 QW_2051_INDEPENDENT_REHEARSAL_GATE.py`
4. `python3 QW_2052_EXTERNAL_SOURCE_ONLY_GOVERNANCE_GATE.py`
5. `python3 QW_2054_PROTOCOL_LOCK_INTEGRITY_GATE.py`
6. `python3 QW_2067_STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_STRENGTHENED_GATE.py`

## 5) Pass criteria
Minimum independent pass package:
- `QW-2049`: pass with all hard flags true,
- `QW-2051`: independent rehearsal pass,
- `QW-2052`: external-source-only governance pass,
- `QW-2054`: lock integrity pass,
- `QW-2067`: strengthened strict internal closure pass.

## 6) What each team should publish
Each team should publish:
- commit hash,
- lock SHA256,
- all generated report JSON files,
- exact command log,
- environment manifest (OS, Python, package versions),
- statement: "no retune/no threshold edits/no post-hoc selection".

## 7) What would count as a real independent confirmation
At least 2-3 independent teams should reproduce the same locked verdicts from fresh environments.

## 8) Extended full-precision package check (recommended)
After the locked confirmatory path, run:
1. `python3 QW_2068_SM_GR_PARAMETER_REGISTRY.py`
2. `python3 QW_2074_STRICT_NOFIT_MISSING_PARAMETER_DERIVATIONS.py`
3. `python3 QW_2075_STRICT_CP_PHASE_DERIVATION_GATE.py`
4. `python3 QW_2081_MISSING14_STRICT_RIGOR_FRONTIER.py`
5. `python3 QW_2083_MISSING14_EPISTEMIC_STATUS_GATE.py`
6. `python3 QW_2093_KERNEL_DERIVED_NONANCHOR_INPUTS_PLAN_EXECUTOR.py`
7. `python3 QW_2085_GF_NONANCHOR_LIFETIME_GATE.py`
8. `python3 QW_2086_MZ_NONANCHOR_EW_POLE_GATE.py`
9. `python3 QW_2087_ALPHA_S_NONANCHOR_BOUNDARY_GATE.py`
10. `python3 QW_2084_T1_NONANCHOR_STRICT_GATE.py`
11. `python3 QW_2069_FULL_SM_GR_DERIVATION_PACKAGE.py`
12. `python3 QW_2072_EW_YUKAWA_FLAVOR_RADIATIVE_BASELINES.py`
13. `python3 QW_2070_FULL_RADIATIVE_PROGRAM_BASELINE.py`
14. `python3 QW_2071_SM_GR_FULL_PRECISION_CLOSURE_GATE.py`
15. `python3 QW_2094_STRICT_RIGOR_DEFECT_SWEEP.py`

Expected current state (2026-03-04):
- `QW-2075`: strict CP phase round partial (`PMNS promoted`, `CKM still non-closing`),
- `QW-2081`: missing-14 strict frontier partial (`1` strict target-miss + `13` strict-unresolved),
- `QW-2083`: deterministic missing-14 epistemic integration pass-all-mapped (`0` direct-missing, `13` explicitly non-closing),
- `QW-2093`: frozen-plan kernel-derived input builder pass (`KERNEL_DERIVED_NONANCHOR_INPUTS_BUILT_FROZEN_PLAN`),
- `QW-2085`: G_F non-anchor lifetime gate pass (`5/6`) with strict non-anchor provenance,
- `QW-2086`: M_Z non-anchor EW-pole gate pass (`5/6`) with strict non-anchor provenance,
- `QW-2087`: alpha_s non-anchor boundary gate pass (`8/9`) with strict non-anchor provenance and validation consistency,
- `QW-2084`: T1 strict non-anchor aggregate gate pass (`6/6`) with no detected anchor leakage/circularity after upstream non-anchor gates,
- `QW-2069`: partial strong internal (`15/32` strict-derived, `0` direct missing, `16` strict-unresolved),
- `QW-2072`: EW/Yukawa/CKM-PMNS baselines created,
- `QW-2073`: channel upgrade checks executed,
- `QW-2070`: partial radiative baseline (`7/7` implemented and `7/7` closure-ready; radiative-sensitive strict-unresolved remains nonzero),
- `QW-2071`: partial strong internal (`3/6` closure flags; `0` direct missing parameters, strict-unresolved remains nonzero, `0` missing radiative channels).
- `QW-2094`: strict-rigor defect sweep pass (`27` checks, `0` failed).

This extension is the transparent path for community teams to track progress from strong internal closure toward full SM+GR precision closure.

## 9) Empirical falsification package (prospective)
Run:
1. `python3 QW_2076_EMPIRICAL_PREDICTION_PREREGISTRATION.py`
2. (GW block) `python3 QW_2078_GW_EXTERNAL_HOLDOUT_AUTOCOLLECTOR.py <external_gw_holdout_features.csv> empirical_observations_input_qw2077.gw_autocollected.json`
3. (PMNS block) `python3 QW_2079_PMNS_CP_EXTERNAL_AUTOCOLLECTOR.py --input-observation empirical_observations_input_qw2077.gw_autocollected.json --output-observation empirical_observations_input_qw2077.gw_pmns_autocollected.json --sin-delta-central <value> --sin-delta-ci95-low <value> --sin-delta-ci95-high <value> --source <citation>`
4. (Cosmology block) `python3 QW_2080_COSMO_WEFF_EXTERNAL_AUTOCOLLECTOR.py --nodes-csv <external_cosmo_nodes.csv> --input-observation empirical_observations_input_qw2077.gw_pmns_autocollected.json --output-observation empirical_observations_input_qw2077.full_autocollected.json --source <citation>`
5. `python3 QW_2077_EMPIRICAL_PREDICTION_VALIDATION_GATE.py empirical_observations_input_qw2077.full_autocollected.json`

What this gives:
- pre-registered predictions before new data arrival,
- reproducible GW metric extraction with locked weights (no manual transcription),
- reproducible, explicitly cited PMNS and cosmology block completion (no manual JSON editing),
- explicit support/falsification/pending verdicts,
- reduced risk of post-hoc interpretation.

---

## PL: Jak inni ludzie mają to niezależnie sprawdzić

## Cel
Wykonać prawdziwie niezależny check confirmatory na zablokowanym protokole, bez dostrajania.

## 1) Referencja freeze
Każdy zespół publikuje:
- hash commita,
- SHA256 lock-file,
- dokładną wersję Pythona i systemu.

Aktualny lock:
- `external_confirmatory_v2/independent_multiteam_lock_qw2053/protocol_lock_qw2053.json`
- SHA256: `2b49385e38a9985b7cf031992da821c3da23a6c0390e4643d56b111adf07c880`

## 2) Środowisko niezależnego zespołu
- Osobna maszyna/konto (nie Twoje środowisko deweloperskie).
- Świeży clone repo.
- Brak lokalnych zmian przed uruchomieniem.

## 3) Zasada danych zewnętrznych
- Nie używać ukrytych lokalnych kopii spoza protokołu.
- Pobierać tylko z udokumentowanych źródeł publicznych:
  - `DATA_SOURCES_EXTERNAL_DOWNLOADS.md`
- Sprawdzić checksumy przed uruchomieniem gate'ów.
- Jeśli archiwum NANOGrav nie leży w katalogu głównym repo, podać jawną ścieżkę w QW-2014:
  - `python3 QW_2014_TRUE_EXTERNAL_BETA_CHANNEL_AUTOCOLLECTOR_V2.py --nanograv-archive <sciezka_do_NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz>`

## 4) Wymagana kolejność uruchomienia (lock)
Uruchomić dokładnie:
1. `python3 QW_2048_SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION.py`
2. `python3 QW_2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.py`
3. `python3 QW_2051_INDEPENDENT_REHEARSAL_GATE.py`
4. `python3 QW_2052_EXTERNAL_SOURCE_ONLY_GOVERNANCE_GATE.py`
5. `python3 QW_2054_PROTOCOL_LOCK_INTEGRITY_GATE.py`
6. `python3 QW_2067_STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_STRENGTHENED_GATE.py`

## 5) Kryteria zaliczenia
Minimalny pakiet niezależnego PASS:
- `QW-2049`: PASS + wszystkie hard flagi true,
- `QW-2051`: PASS próby niezależnej,
- `QW-2052`: PASS governance external-source-only,
- `QW-2054`: PASS integralności locka,
- `QW-2067`: PASS strengthened strict internal closure.

## 6) Co każdy zespół powinien opublikować
Każdy zespół publikuje:
- hash commita,
- SHA256 locka,
- wszystkie wygenerowane raporty JSON,
- log komend,
- manifest środowiska (OS, Python, wersje pakietów),
- oświadczenie: "no retune/no threshold edits/no post-hoc selection".

## 7) Co uznajemy za realne niezależne potwierdzenie
Co najmniej 2-3 niezależne zespoły muszą odtworzyć te same zablokowane werdykty na świeżych środowiskach.

## 8) Rozszerzony check pełnego pakietu precyzyjnego (zalecane)
Po ścieżce locked confirmatory uruchomić:
1. `python3 QW_2068_SM_GR_PARAMETER_REGISTRY.py`
2. `python3 QW_2074_STRICT_NOFIT_MISSING_PARAMETER_DERIVATIONS.py`
3. `python3 QW_2075_STRICT_CP_PHASE_DERIVATION_GATE.py`
4. `python3 QW_2081_MISSING14_STRICT_RIGOR_FRONTIER.py`
5. `python3 QW_2083_MISSING14_EPISTEMIC_STATUS_GATE.py`
6. `python3 QW_2093_KERNEL_DERIVED_NONANCHOR_INPUTS_PLAN_EXECUTOR.py`
7. `python3 QW_2085_GF_NONANCHOR_LIFETIME_GATE.py`
8. `python3 QW_2086_MZ_NONANCHOR_EW_POLE_GATE.py`
9. `python3 QW_2087_ALPHA_S_NONANCHOR_BOUNDARY_GATE.py`
10. `python3 QW_2084_T1_NONANCHOR_STRICT_GATE.py`
11. `python3 QW_2069_FULL_SM_GR_DERIVATION_PACKAGE.py`
12. `python3 QW_2072_EW_YUKAWA_FLAVOR_RADIATIVE_BASELINES.py`
13. `python3 QW_2070_FULL_RADIATIVE_PROGRAM_BASELINE.py`
14. `python3 QW_2071_SM_GR_FULL_PRECISION_CLOSURE_GATE.py`
15. `python3 QW_2094_STRICT_RIGOR_DEFECT_SWEEP.py`

Oczekiwany obecny stan (2026-03-04):
- `QW-2075`: częściowy status ścisłej rundy faz CP (`PMNS promoted`, `CKM nadal non-closing`),
- `QW-2081`: częściowy frontier brakującej 14 (`1` strict target-miss + `13` strict-unresolved),
- `QW-2083`: deterministyczna integracja statusów epistemicznych brakującej 14 (`0` direct missing, `13` pozycji non-closing sklasyfikowanych),
- `QW-2093`: PASS budowy wejść kernel-derived z zamrożonego planu (`KERNEL_DERIVED_NONANCHOR_INPUTS_BUILT_FROZEN_PLAN`),
- `QW-2085`: PASS bramki G_F non-anchor lifetime (`5/6`) ze ścisłym pochodzeniem non-anchor,
- `QW-2086`: PASS bramki M_Z non-anchor EW-pole (`5/6`) ze ścisłym pochodzeniem non-anchor,
- `QW-2087`: PASS bramki alpha_s non-anchor boundary (`8/9`) ze ścisłym pochodzeniem non-anchor i spójnością walidacyjną,
- `QW-2084`: PASS bramki agregującej T1 strict non-anchor (`6/6`) bez wykrytego anchor leakage/circularity po bramkach upstream,
- `QW-2069`: status częściowy strong internal (`15/32` ścisłe derivacje, `0` direct missing, `16` strict-unresolved),
- `QW-2072`: utworzono baseline EW/Yukawa/CKM-PMNS,
- `QW-2073`: wykonano upgrade kanałów radiacyjnych,
- `QW-2070`: status częściowy radiative baseline (`7/7` zaimplementowanych i `7/7` closure-ready; radiative-sensitive strict-unresolved pozostaje > 0),
- `QW-2071`: status częściowy strong internal (`3/6` flag domknięcia; `0` direct missing, strict-unresolved pozostaje > 0, `0` brakujących kanałów radiacyjnych).
- `QW-2094`: PASS sweepu usterek rygoru ścisłego (`27` kontroli, `0` błędów).

To rozszerzenie daje przejrzystą ścieżkę dla zespołów zewnętrznych od mocnego domknięcia wewnętrznego do pełnego domknięcia precyzyjnego SM+GR.

## 9) Pakiet falsyfikacji empirycznej (prospektywny)
Uruchomić:
1. `python3 QW_2076_EMPIRICAL_PREDICTION_PREREGISTRATION.py`
2. (Blok GW) `python3 QW_2078_GW_EXTERNAL_HOLDOUT_AUTOCOLLECTOR.py <external_gw_holdout_features.csv> empirical_observations_input_qw2077.gw_autocollected.json`
3. (Blok PMNS) `python3 QW_2079_PMNS_CP_EXTERNAL_AUTOCOLLECTOR.py --input-observation empirical_observations_input_qw2077.gw_autocollected.json --output-observation empirical_observations_input_qw2077.gw_pmns_autocollected.json --sin-delta-central <wartość> --sin-delta-ci95-low <wartość> --sin-delta-ci95-high <wartość> --source <cytowanie>`
4. (Blok cosmology) `python3 QW_2080_COSMO_WEFF_EXTERNAL_AUTOCOLLECTOR.py --nodes-csv <external_cosmo_nodes.csv> --input-observation empirical_observations_input_qw2077.gw_pmns_autocollected.json --output-observation empirical_observations_input_qw2077.full_autocollected.json --source <cytowanie>`
5. `python3 QW_2077_EMPIRICAL_PREDICTION_VALIDATION_GATE.py empirical_observations_input_qw2077.full_autocollected.json`

Co to daje:
- pre-rejestrację predykcji przed pojawieniem się nowych danych,
- reprodukowalne wyliczenie metryk GW na zablokowanych wagach (bez ręcznego przepisywania),
- reprodukowalne, jawnie cytowane uzupełnienie bloków PMNS i cosmology bez ręcznej edycji JSON,
- jawne werdykty support/falsified/pending,
- ograniczenie ryzyka interpretacji post-hoc.
