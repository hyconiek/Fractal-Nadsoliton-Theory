# RELEASE 4.9 TEXTBOOK EDITION (EN + PL)

**Version:** 4.9  
**Date:** 2026-03-04  
**Branch:** `main`

---

## ENGLISH VERSION

## 1) What This Release Is

Release 4.9 is a **rigor release** focused on one central question:

> Can the same frozen Nadsoliton kernel be supported both by
> 1) macro-level sector gates (mass + flavor + GW), and
> 2) micro-level pointwise derivation,
> without sector-by-sector retuning?

In plain language:
- we test whether the theory is only "fitting data,"
- or whether it is also structurally consistent from the inside.

## 2) One-Page Status (For a High-School Reader)

### 2.1 What is already strong

1. One frozen kernel passes mass + flavor + GW gates together.
2. Blind external transfer tests pass.
3. Pointwise micro-identifiability was repaired in this release path.
4. A reproducibility rehearsal in isolated environment passes.
5. Internal first-principles closure has been strengthened by QW-2063..QW-2067.

### 2.2 What is still missing

1. A truly independent multiteam confirmatory run outside this local environment.
2. Full precision SM+GR package closure (QW-2069/2070/2071 shows partial status).

### 2.3 Bottom line

- **Internal scientific closure:** very strong.
- **Community-level final ToE claim:** not yet; external independent confirmation is still required.

## 3) The Core Idea in Textbook Language

### 3.1 Nadsoliton picture

The Nadsoliton framework assumes that what we call particles, forces, and spacetime patterns are different manifestations of one deeper informational-dynamical structure.

### 3.2 Frozen kernel (single structural function)

The common structural function is:

$$
K(d)=\frac{\cos(\omega d+\phi)}{1+\beta d^{\eta}}
$$

In Release 4.9 working gate (QW-2049), the selected kernel is:

- $\omega = 0.18575$
- $\phi = 0.16250$
- $\beta = 1.00000$
- $\eta = 1.80000$

Interpretation for students:
- cosine term gives oscillatory pattern,
- denominator damps amplitude with distance,
- same function is reused across sectors.

## 4) Main Breakthrough of Release 4.9

The key breakthrough is not "new retuning." It is a **methodological repair**:

- old phase estimation in micro-pointwise derivation tended to collapse into low-frequency proxy,
- new method locks phase from the signed-dynamic operator spectrum,
- this repaired the two missing strict micro flags:
  - enough pointwise bins,
  - strong phase conditioning.

This is why QW-2048 became full PASS.

## 5) Key New Gates and Results

### 5.1 QW-2048 (spectral phase-locked pointwise derivation)

- Verdict: `SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION_PASS`
- Pass count: `8/8`
- Important metrics:
  - `n_rows_total = 342`
  - `n_bins = 17`
  - `phase_min_median = 0.852` (above strict threshold 0.75)
  - joint target-in-CI95 bin coverage = `0.8235`

Meaning:
- pointwise micro-identifiability is no longer the blocker.

### 5.2 QW-2049 (micro -> Stage-C -> external intersection)

- Verdict: `SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE_PASS`
- Pass count: `7/7`
- All hard flags pass simultaneously.

Meaning:
- internal bridge between micro support and Stage-C macro closure is now strict-pass.

### 5.3 QW-2050 (freeze bundle)

- Verdict: `SPECTRAL_MICRO_BRIDGE_FREEZE_BUNDLE_READY`
- Package for independent execution was built:
  - manifest with SHA256,
  - runbook,
  - fixed kernel vector.

Meaning:
- handoff package exists in reproducible form.

### 5.4 QW-2051 (independent rehearsal in isolated `/tmp`)

- Final verdict after dependency fix: `INDEPENDENT_REHEARSAL_GATE_PASS`
- Pass count: `7/7`

Meaning:
- bundle was not only prepared, but also successfully dry-run under isolation rules.

### 5.5 QW-2063 (deterministic no-scan triad reconstruction)

- Verdict: `DERIVATIONAL_RECONSTRUCTION_TRIAD_PASS_PHYSICAL_PROVISIONAL_FIRST_PRINCIPLES`
- Pass count: `11/12`
- Physical triad: PASS

Meaning:
- mass + flavor + GW pass together under deterministic no-scan mapping,
- one formal first-principles flag still remained open before QW-2064.

### 5.6 QW-2064 (micro-derived renormalization constants gate)

- Verdict: `MICRO_DERIVED_RENORMALIZATION_CONSTANTS_GATE_PASS_WITH_WIDE_CI_WARNING`
- Pass count: `8/8`
- Key values:
  - target: `Z_beta=100.0`, `delta_eta=0.8`
  - micro median: `Z_beta=114.740`, `delta_eta=1.125`

Meaning:
- foundational renormalization constants are supported by micro-derivation,
- but `Z_beta` dispersion was still broad (warning).

### 5.7 QW-2065 (strict first-principles internal closure gate)

- Verdict: `STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_PASS`
- Pass count: `12/12`

Meaning:
- internal strict first-principles closure passes on the current gate chain.

### 5.8 QW-2066 (compatibility-filtered tightening)

- Verdict: `COMPATIBILITY_FILTERED_MICRO_CONSTANTS_TIGHTENING_PASS`
- Pass count: `6/6`
- Dispersion tightening:
  - `z_beta_log_iqr: 3.133 -> 2.124`

Meaning:
- the micro-constants warning is reduced by deterministic compatibility filtering.

### 5.9 QW-2067 (strengthened strict closure gate)

- Verdict: `STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_STRENGTHENED_PASS`
- Pass count: `3/3`

Meaning:
- strict internal closure remains PASS and is additionally strengthened.

### 5.10 QW-2068 (explicit SM+GR parameter registry)

- Artifact: `report_qw2068_sm_gr_parameter_registry.json`
- Target set size: `32` parameters

Meaning:
- the closure target is explicit, finite, and auditable.

### 5.11 QW-2069 (full SM+GR derivation package audit)

- Verdict: `FULL_SM_GR_DERIVATION_PACKAGE_PARTIAL_STRONG_INTERNAL`
- Coverage:
  - strict-derived: `28/32`
  - model-formula-only: `1/32`
  - anchor-dependent no-fit: `1/32`
  - coupled-anchor-dependent: `0/32`
  - model-assumption-nonclosing: `0/32`
  - SI-definition constants: `2/32`
  - missing direct derivation: `0/32`
  - strict-unresolved (closure criterion): `7/32`

Meaning:
- internal closure is strong in covered domains,
- full all-parameter precision closure is not yet reached.

### 5.12 QW-2070 (full radiative program baseline)

- Verdict: `FULL_RADIATIVE_PROGRAM_PARTIAL_BASELINE`
- Coverage:
  - implemented channels: `7/7`
  - closure-ready channels: `7/7` (after QW-2073 upgrade)
  - missing channels: `0/7`

Meaning:
- radiative program is now formalized,
- and channel-level closure is complete; remaining blocker is parameter-level direct derivation coverage.

### 5.13 QW-2071 (full precision closure gate)

- Verdict: `SM_GR_FULL_PRECISION_CLOSURE_PARTIAL_STRONG_INTERNAL`
- Gate flags pass: `3/6`
- Missing parameters: `0`
- Strict-unresolved parameters: `7`
- Missing radiative channels: `0`
- Implemented but non-closing radiative channels: `0`

Meaning:
- strict internal foundation exists,
- global full-precision SM+GR closure is still open.

### 5.14 QW-2074 (strict no-fit missing-parameter derivation round)

- Verdict: `STRICT_NOFIT_MISSING_PARAMETER_DERIVATION_ROUND1`
- Updates: `4`
- Epistemic impact on package map:
  - `2` entries moved to `physical_relation_anchor_dependent`,
  - `2` entries moved to `si_definition`,
  - missing direct derivations reduced from `19` to `15`.

Meaning:
- this is an anti-overclaim step: no retune, no scan, explicit epistemic labeling,
- these entries are **not** counted as strict first-principles derivations.

### 5.15 QW-2075 (strict CP phase derivation gate)

- Verdict: `STRICT_CP_PHASE_DERIVATION_PARTIAL_PMNS_ONLY`
- Pass count: `7/8`
- Updates promoted: `1` (`delta_cp_pmns`)

Meaning:
- PMNS CP phase is now derived in deterministic strict-internal mode,
- CKM CP phase remains non-closing in this round (outside registry tolerance branches).

### 5.16 QW-2078 (GW external holdout autocollector)

- Verdict: `GW_EXTERNAL_HOLDOUT_AUTOCOLLECTED`
- Function: auto-builds GW observation block for QW-2077 from locked metrics/weights

Meaning:
- GW branch can be evaluated reproducibly without manual metric transcription,
- the independence claim still depends on external data provenance and execution separation.

### 5.17 QW-2079 + QW-2080 (PMNS + cosmology external autocollectors)

- Verdicts:
  - `PMNS_CP_EXTERNAL_AUTOCOLLECTED` (QW-2079)
  - `COSMO_WEFF_EXTERNAL_AUTOCOLLECTED` (QW-2080)
- Function:
  - QW-2079 fills PMNS CI block,
  - QW-2080 fills preregistered cosmology z-nodes from validated CSV.

Meaning:
- QW-2077 no longer needs manual PMNS/cosmology value editing,
- pending status can now be removed strictly by supplying externally cited input files.

### 5.18 QW-2081 (missing-14 strict-rigor frontier audit)

- Verdict: `MISSING14_STRICT_RIGOR_FRONTIER_PARTIAL_ONLY`
- Classification counts:
  - strict candidate target-miss: `3` (`delta_cp_ckm`, `h0`, `lambda_cosmological`)
  - anchor-dependent baseline-only: `1` (`G_newton`)

Meaning:
- this gives an explicit, auditable boundary between what is currently derivable and what is not,
- it blocks false closure claims while preserving deterministic reproducibility.

### 5.19 QW-2082 (missing-14 strict closure roadmap)

- Verdict: `MISSING14_STRICT_CLOSURE_ROADMAP_READY`
- Tiered execution plan:
  - `T1`: `delta_cp_ckm`
  - `T3`: `h0`, `lambda_cosmological`
  - `T4`: `G_newton`

Meaning:
- closure work is now operationally ordered, not only diagnostically described,
- each unresolved parameter has explicit required observables and hard closure conditions.

### 5.20 QW-2083 (missing-14 epistemic status gate)

- Verdict: `MISSING14_EPISTEMIC_STATUS_GATE_PASS_WITH_TARGET_MISS`
- Deterministic status integration for all missing-14 IDs:
  - `3` strict target-miss (`delta_cp_ckm`, `h0`, `lambda_cosmological`)
  - `1` non-closing but explicitly classified entry (`G_newton`)
  - `0` still strictly missing

Meaning:
- all 14 missing IDs are now explicitly handled in the package map,
- but explicit handling is not equal to strict closure.

### 5.21 QW-2084 (T1 strict non-anchor gate)

- Verdict: `T1_NONANCHOR_STRICT_GATE_PASS`
- Pass count: `6/6`
- T1 strict non-anchor checks (aggregate from dedicated gates):
  - `alpha_s_mz`: PASS (from `QW-2087` strict non-anchor boundary gate),
  - `g_f`: PASS (from `QW-2085` strict non-anchor lifetime gate),
  - `m_z`: PASS (from `QW-2086` strict non-anchor EW-pole gate).

Meaning:
- T1 non-anchor closure is now achieved under deterministic, provenance-checked upstream gates,
- the aggregate gate remains strict because it propagates only non-anchor PASS artifacts.

### 5.22 QW-2085 (G_F non-anchor lifetime gate)

- Verdict: `GF_NONANCHOR_LIFETIME_GATE_PASS`
- Pass count: `5/6`
- Current result:
  - strict non-anchor PASS on kernel-derived muon-lifetime chain,
  - tolerance check remains satisfied (`rel_err_pct = 2.6524` vs tolerance `5.0%`).

Meaning:
- strict promotion is now achieved with explicit kernel-derived provenance,
- no retune/no-scan constraints remain active.

### 5.23 QW-2086 (M_Z non-anchor EW-pole gate)

- Verdict: `MZ_NONANCHOR_EW_POLE_GATE_PASS`
- Pass count: `5/6`
- Current result:
  - strict non-anchor PASS on kernel-derived EW-pole chain,
  - tolerance check remains satisfied (`rel_err_pct = 1.7202` vs tolerance `5.0%`).

Meaning:
- this T1 channel is now closure-ready under strict non-anchor criteria,
- provenance guard remains active to prevent future anchor leakage.

### 5.24 QW-2087 (alpha_s non-anchor boundary gate)

- Verdict: `ALPHA_S_NONANCHOR_BOUNDARY_GATE_PASS`
- Pass count: `8/9`
- Current result:
  - strict non-anchor PASS on kernel-derived boundary + validation points,
  - tolerance check remains satisfied (`rel_err_pct = 2.1338` vs tolerance `5.0%`),
  - validation consistency is strong (`z_mean = 0.4294`, `z_max = 0.4948`).

Meaning:
- strict promotion now includes both provenance and validation consistency,
- numeric agreement is accepted only together with those hard checks.

### 5.25 QW-2093 (kernel-derived non-anchor input plan executor)

- Verdict: `KERNEL_DERIVED_NONANCHOR_INPUTS_BUILT_FROZEN_PLAN`
- Output artifacts:
  - `t1_nonanchor_observables_input_qw2085_2086.json`
  - `t1_nonanchor_alpha_s_input_qw2087.json`
- Rule: frozen formula plan, no scan, no retune.

Meaning:
- T1 non-anchor input generation is now deterministic and auditable,
- QW-2085/2086/2087 PASS status is tied to a reproducible input-build step.

### 5.26 QW-2094 (strict-rigor defect sweep)

- Verdict: `STRICT_RIGOR_DEFECT_SWEEP_PASS_NO_CRITICAL_DEFECTS`
- Sweep scope:
  - internal consistency of QW-2084/2085/2086/2087/2090/2091/2092/2097/2098,
  - cross-report consistency between QW-2069 and QW-2071,
  - pre-gate consistency for QW-2102 and QW-2103,
  - merged preflight consistency for QW-2104,
  - intake/gap consistency for QW-2105 and QW-2106,
  - design-guidance consistency for QW-2107 and QW-2108,
  - evidence-manifest consistency for QW-2109,
  - strict-level/status integrity in package report.
- Result: `130` checks, `0` failed.

Meaning:
- no critical implementation mismatch was detected in the audited non-anchor T1 closure path,
- this reduces risk of false-success reporting caused by internal report inconsistency.

## 6) What "Internally Closed" Means (And What It Does Not Mean)

### 6.1 It means

1. The same frozen kernel survives multiple strict internal gates.
2. Micro-derivation and macro-gates are no longer disconnected.
3. Blind transfer metrics remain significant in primary and stress checks.

### 6.2 It does not mean

1. Final global ToE proof accepted by independent teams.
2. Immunity from future falsification.

Science rule: even very strong internal evidence is still internal until independent groups replicate it.

## 7) External Metrics Snapshot (From QW-2049)

Primary external holdout:
- corr = `0.054705` (q95 = `0.036683`, p = `0.006999`)
- gain = `0.000571` (q95 = `-0.001152`, p = `0.006999`)

Stress external holdout:
- corr = `0.330346` (q95 = `0.050699`, p = `0.000200`)
- gain = `0.055985` (q95 = `-0.032107`, p = `0.000200`)

Student-level interpretation:
- primary signal is modest but significant,
- stress signal is strong and highly significant,
- both satisfy gate conditions.

## 8) Current Scientific Position vs "Holy Grail" Claim

A fair statement in Release 4.9 is:

- This is now a **high-rigor unification candidate with strict internal closure path**.
- It is **not yet** the final community-confirmed ToE.
- The remaining step is inherently social-scientific, not just computational: independent external multiteam confirmation.

## 9) What Would Falsify It Now

The framework would be seriously weakened if independent multiteam runs show one of the following:

1. fixed-kernel gate fails under locked protocol,
2. external blind significance disappears under independent execution,
3. reproducibility breaks under manifest-verified artifacts.

## 10) Exact Next Step After Release 4.9

- Execute: `RUN_TRULY_INDEPENDENT_MULTITEAM_CONFIRMATORY_PACKAGE`
- With fixed bundle from QW-2050
- And rehearsal protocol validated by QW-2051

That is the shortest scientifically correct path from "strong internal closure" to "credible external confirmation".

## 11) External Data Rule (No Large Binaries in Git)

- Large raw payloads are treated as external sources, not frozen git artifacts.
- Canonical source list and download commands:
  - `DATA_SOURCES_EXTERNAL_DOWNLOADS.md`
- Practical meaning:
1. freeze scripts/manifests/reports in git,
2. fetch heavy public archives from official providers,
3. verify integrity locally, then run fixed gates.

## 12) Do We Already Derive *All* Known Physical Values?

Short answer: **No, not yet in full global strict form.**

What is already derived in strict internal gate scope:
1. mass-chain gate targets,
2. CKM/PMNS gate-level flavor targets,
3. PMNS CP phase from deterministic complex flavor operator,
4. GW discriminator targets,
5. micro-supported renormalization constants (`Z_beta`, `delta_eta`) with tightening.
6. explicit package audit map and closure gate logic (`QW-2069`, `QW-2071`).

What is still open for a full “all known values” claim:
1. package closure gap remains: `28/32` parameters are strict-derived, with `7/32` still strict-unresolved (even though direct `missing` is now `0/32`),
2. strict frontier over the current missing-14 remains partial (`QW-2081`: `3` strict target-miss + `1` non-closing),
3. closure roadmap exists (`QW-2082`) and is partially executed (T1/T2 + EW-secondary propagation), but T3/T4 and CKM CP tolerance closure are still open,
4. epistemic status integration (`QW-2083`) addresses all missing-14 IDs, but does not by itself convert non-closing classes into strict closure,
5. full exhaustive Standard Model + GR parameter derivation package at final precision,
6. independent external multiteam replication/audit.

Therefore:
- internal first-principles closure path is now very strong,
- but a final “all known physical values fully derived and externally confirmed” statement is still premature.

---

## WERSJA POLSKA

## 1) Czym jest ten release

Release 4.9 to **release rygorystyczny**. Odpowiada na jedno kluczowe pytanie:

> Czy to samo zamrożone jądro nadsolitonowe jest jednocześnie wspierane przez
> 1) bramki sektorowe (masa + flavor + GW), oraz
> 2) derywację punktową na poziomie mikro,
> bez dostrajania sektor po sektorze?

Po ludzku:
- sprawdzamy, czy teoria tylko "dopasowuje dane",
- czy też jest spójna konstrukcyjnie od środka.

## 2) Status w skrócie (dla laika)

### 2.1 Co jest już mocne

1. Jedno zamrożone jądro przechodzi wspólnie masę, flavor i GW.
2. Ślepe testy transferu zewnętrznego przechodzą.
3. Naprawiono historyczną lukę identyfikowalności mikro-punktowej.
4. Izolowana próba generalna reprodukcji przeszła.
5. Domknięcie first-principles zostało dodatkowo wzmocnione przez QW-2063..QW-2067.

### 2.2 Czego jeszcze brakuje

1. Prawdziwie niezależnej replikacji multiteam poza tym środowiskiem.
2. Pełnego domknięcia precyzyjnego pakietu SM+GR (QW-2069/2070/2071 pokazuje status częściowy).

### 2.3 Wniosek

- **Domknięcie wewnętrzne:** bardzo silne.
- **Finalny status ToE dla społeczności naukowej:** jeszcze nie; potrzebna niezależna walidacja zewnętrzna.

## 3) Rdzeń teorii w języku podręcznikowym

### 3.1 Obraz nadsolitonu

W tym podejściu cząstki, oddziaływania i wzory czasoprzestrzeni są różnymi przejawami jednej głębszej struktury informacyjno-dynamicznej.

### 3.2 Wspólne jądro (jedna funkcja strukturalna)

$$
K(d)=\frac{\cos(\omega d+\phi)}{1+\beta d^{\eta}}
$$

W bramce Release 4.9 (QW-2049) wybrane jądro to:

- $\omega = 0.18575$
- $\phi = 0.16250$
- $\beta = 1.00000$
- $\eta = 1.80000$

Interpretacja licealna:
- licznik (cosinus) daje rytm/oscylację,
- mianownik tłumi amplitudę wraz z odległością,
- ta sama funkcja jest używana we wszystkich sektorach.

## 4) Główny przełom Release 4.9

To nie jest przełom typu "kolejne strojenie".
To jest **naprawa metodologii**:

- stary odczyt fazy w mikro-derywacji miał tendencję do zapadania się do niskiej częstości,
- nowy odczyt fazy bierze nośną ze spektrum operatora signed-dynamic,
- to domknęło dwa brakujące warunki rygorystyczne:
  - liczbę binów punktowych,
  - siłę warunku fazowego.

Dlatego QW-2048 dał pełny PASS.

## 5) Najważniejsze nowe bramki i wyniki

### 5.1 QW-2048 (punktowa derywacja z blokadą fazy widmowej)

- Werdykt: `SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION_PASS`
- Wynik: `8/8`
- Kluczowe metryki:
  - `n_rows_total = 342`
  - `n_bins = 17`
  - `phase_min_median = 0.852` (powyżej progu 0.75)
  - pokrycie joint target-in-CI95 = `0.8235`

Znaczenie:
- identyfikowalność punktowa mikro przestała być blokerem.

### 5.2 QW-2049 (przecięcie mikro -> Stage-C -> external)

- Werdykt: `SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE_PASS`
- Wynik: `7/7`
- Wszystkie flagi twarde są TRUE.

Znaczenie:
- most między mikro i makro domknął się rygorystycznie.

### 5.3 QW-2050 (pakiet freeze)

- Werdykt: `SPECTRAL_MICRO_BRIDGE_FREEZE_BUNDLE_READY`
- Zbudowano pakiet do niezależnego uruchomienia:
  - manifest SHA256,
  - runbook,
  - zamrożony wektor jądra.

### 5.4 QW-2051 (izolowana próba generalna)

- Werdykt końcowy po poprawce zależności: `INDEPENDENT_REHEARSAL_GATE_PASS`
- Wynik: `7/7`

Znaczenie:
- pakiet nie tylko istnieje, ale jest też proceduralnie odtwarzalny 1:1 w izolacji.

### 5.5 QW-2063 (deterministyczna rekonstrukcja triady bez skanu)

- Werdykt: `DERIVATIONAL_RECONSTRUCTION_TRIAD_PASS_PHYSICAL_PROVISIONAL_FIRST_PRINCIPLES`
- Wynik: `11/12`
- Triada fizyczna: PASS

Znaczenie:
- masa + flavor + GW przechodzą razem w deterministycznym trybie no-scan,
- przed QW-2064 pozostawała jedna formalna flaga first-principles.

### 5.6 QW-2064 (bramka mikro-derywacji stałych renormalizacyjnych)

- Werdykt: `MICRO_DERIVED_RENORMALIZATION_CONSTANTS_GATE_PASS_WITH_WIDE_CI_WARNING`
- Wynik: `8/8`
- Kluczowe liczby:
  - target: `Z_beta=100.0`, `delta_eta=0.8`
  - mediana mikro: `Z_beta=114.740`, `delta_eta=1.125`

Znaczenie:
- stale renormalizacyjne mają wsparcie z mikro-derywacji,
- ale rozrzut `Z_beta` był jeszcze szeroki (warning).

### 5.7 QW-2065 (strict first-principles internal closure gate)

- Werdykt: `STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_PASS`
- Wynik: `12/12`

Znaczenie:
- wewnętrzne ścisłe domknięcie first-principles przechodzi w aktualnym łańcuchu bramek.

### 5.8 QW-2066 (zawężenie rozrzutu przez compatibility filtering)

- Werdykt: `COMPATIBILITY_FILTERED_MICRO_CONSTANTS_TIGHTENING_PASS`
- Wynik: `6/6`
- Zawężenie:
  - `z_beta_log_iqr: 3.133 -> 2.124`

Znaczenie:
- warning dla stałych mikro został istotnie zredukowany deterministycznym filtrem kompatybilności.

### 5.9 QW-2067 (wzmocniona bramka ścisłego domknięcia)

- Werdykt: `STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_STRENGTHENED_PASS`
- Wynik: `3/3`

Znaczenie:
- ścisłe domknięcie wewnętrzne pozostaje PASS i jest dodatkowo wzmocnione.

### 5.10 QW-2068 (jawny rejestr parametrów SM+GR)

- Artefakt: `report_qw2068_sm_gr_parameter_registry.json`
- Rozmiar zbioru docelowego: `32` parametry

Znaczenie:
- cel domknięcia jest jawny, skończony i audytowalny.

### 5.11 QW-2069 (audyt pełnego pakietu derivacji SM+GR)

- Werdykt: `FULL_SM_GR_DERIVATION_PACKAGE_PARTIAL_STRONG_INTERNAL`
- Pokrycie:
  - ścisłe derivacje: `28/32`
  - tylko model-formula: `1/32`
  - anchor-dependent no-fit: `1/32`
  - coupled-anchor-dependent: `0/32`
  - model-assumption-nonclosing: `0/32`
  - stałe SI-definition: `2/32`
  - brak bezpośredniej derivacji: `0/32`
  - strict-unresolved (kryterium domknięcia): `7/32`

Znaczenie:
- domknięcie wewnętrzne jest mocne w pokrytych domenach,
- pełne precyzyjne domknięcie wszystkich parametrów nie jest jeszcze osiągnięte.

### 5.12 QW-2070 (bazowy pełny program radiacyjny)

- Werdykt: `FULL_RADIATIVE_PROGRAM_PARTIAL_BASELINE`
- Pokrycie:
  - zaimplementowane kanały: `7/7`
  - kanały closure-ready: `7/7` (po upgrade QW-2073)
  - brakujące kanały: `0/7`

Znaczenie:
- program radiacyjny został formalnie zdefiniowany,
- a domknięcie na poziomie kanałów jest kompletne; blokerem pozostaje pokrycie bezpośrednich derivacji parametrów.

### 5.13 QW-2071 (bramka pełnego domknięcia precyzyjnego)

- Werdykt: `SM_GR_FULL_PRECISION_CLOSURE_PARTIAL_STRONG_INTERNAL`
- Przejście flag: `3/6`
- Brakujące parametry: `0`
- Parametry strict-unresolved: `7`
- Brakujące kanały radiacyjne: `0`
- Kanały radiacyjne zaimplementowane, ale non-closing: `0`

Znaczenie:
- ścisła baza wewnętrzna istnieje,
- globalne pełne domknięcie precyzyjne SM+GR pozostaje otwarte.

### 5.14 QW-2074 (ścisła runda no-fit dla brakujących parametrów)

- Werdykt: `STRICT_NOFIT_MISSING_PARAMETER_DERIVATION_ROUND1`
- Liczba aktualizacji: `4`
- Wpływ epistemiczny na mapę pakietu:
  - `2` pozycje przeniesiono do `physical_relation_anchor_dependent`,
  - `2` pozycje przeniesiono do `si_definition`,
  - liczba braków bezpośredniej derivacji spadła z `19` do `15`.

Znaczenie:
- to krok anty-nadinterpretacyjny: brak retune, brak skanu, jawne etykiety epistemiczne,
- te pozycje **nie** są liczone jako ścisłe derivacje first-principles.

### 5.15 QW-2075 (ścisła bramka derywacji faz CP)

- Werdykt: `STRICT_CP_PHASE_DERIVATION_PARTIAL_PMNS_ONLY`
- Wynik: `7/8`
- Promowane aktualizacje: `1` (`delta_cp_pmns`)

Znaczenie:
- faza CP PMNS została wyprowadzona w trybie deterministycznym strict-internal,
- faza CP CKM pozostaje niedomknięta w tej rundzie (poza gałęziami tolerancji rejestru).

### 5.16 QW-2078 (autocollector zewnętrznego holdoutu GW)

- Werdykt: `GW_EXTERNAL_HOLDOUT_AUTOCOLLECTED`
- Funkcja: automatycznie buduje blok obserwacyjny GW do QW-2077 z zablokowanych metryk/wag

Znaczenie:
- gałąź GW można oceniać odtwarzalnie bez ręcznego przepisywania metryk,
- teza o niezależności nadal zależy od pochodzenia danych zewnętrznych i separacji środowiska.

### 5.17 QW-2079 + QW-2080 (autocollectory PMNS + kosmologia)

- Werdykty:
  - `PMNS_CP_EXTERNAL_AUTOCOLLECTED` (QW-2079)
  - `COSMO_WEFF_EXTERNAL_AUTOCOLLECTED` (QW-2080)
- Funkcja:
  - QW-2079 wypełnia blok CI dla PMNS,
  - QW-2080 wypełnia prerejestrowane węzły cosmology z walidowanego CSV.

Znaczenie:
- QW-2077 nie wymaga już ręcznej edycji wartości PMNS/cosmology,
- status pending można zdjąć rygorystycznie po dostarczeniu zewnętrznych, cytowalnych danych.

### 5.18 QW-2081 (audyt missing-14 w rygorze ścisłym)

- Werdykt: `MISSING14_STRICT_RIGOR_FRONTIER_PARTIAL_ONLY`
- Liczności klas:
  - strict candidate target-miss: `3` (`delta_cp_ckm`, `h0`, `lambda_cosmological`)
  - anchor-dependent baseline-only: `1` (`G_newton`)

Znaczenie:
- to jawna, audytowalna granica między tym, co dziś jest wyprowadzalne, a tym, co jeszcze nie,
- zapobiega fałszywemu domknięciu, zachowując pełną odtwarzalność deterministyczną.

### 5.19 QW-2082 (roadmap ścisłego domykania missing-14)

- Werdykt: `MISSING14_STRICT_CLOSURE_ROADMAP_READY`
- Plan wykonania warstwowego:
  - `T1`: `delta_cp_ckm`
  - `T3`: `h0`, `lambda_cosmological`
  - `T4`: `G_newton`

Znaczenie:
- prace domykające są teraz uporządkowane operacyjnie, a nie tylko opisane diagnostycznie,
- każdy nierozwiązany parametr ma jawnie wskazane wymagane obserwable i twarde warunki domknięcia.

### 5.20 QW-2083 (bramka epistemicznego statusu missing-14)

- Werdykt: `MISSING14_EPISTEMIC_STATUS_GATE_PASS_WITH_TARGET_MISS`
- Deterministyczna integracja statusów dla całej brakującej 14:
  - `3` strict target-miss (`delta_cp_ckm`, `h0`, `lambda_cosmological`)
  - `1` pozycja non-closing z jawną klasyfikacją (`G_newton`)
  - `0` pozycji nadal ściśle brakujących

Znaczenie:
- cała brakująca 14 ma teraz jawne mapowanie statusu w pakiecie,
- ale jawna obsługa nie jest równoznaczna ze ścisłym domknięciem.

### 5.21 QW-2084 (bramka T1 strict non-anchor)

- Werdykt: `T1_NONANCHOR_STRICT_GATE_PASS`
- Wynik: `6/6`
- Kontrole T1 strict non-anchor (agregacja z dedykowanych bramek):
  - `alpha_s_mz`: PASS (z bramki `QW-2087` strict non-anchor boundary),
  - `g_f`: PASS (z bramki `QW-2085` strict non-anchor lifetime),
  - `m_z`: PASS (z bramki `QW-2086` strict non-anchor EW-pole).

Znaczenie:
- domknięcie T1 non-anchor jest teraz osiągnięte przy deterministycznych, provenance-checked bramkach upstream,
- bramka agregująca pozostaje ścisła, bo propaguje wyłącznie artefakty non-anchor PASS.

### 5.22 QW-2085 (bramka G_F non-anchor lifetime)

- Werdykt: `GF_NONANCHOR_LIFETIME_GATE_PASS`
- Wynik: `5/6`
- Wynik bieżący:
  - strict non-anchor PASS na kernel-derived łańcuchu czasu życia mionu,
  - test tolerancji pozostaje spełniony (`rel_err_pct = 2.6524` przy tolerancji `5.0%`).

Znaczenie:
- osiągnięto promocję strict przy jawnym pochodzeniu kernel-derived,
- utrzymano warunek no-retune/no-scan.

### 5.23 QW-2086 (bramka M_Z non-anchor EW-pole)

- Werdykt: `MZ_NONANCHOR_EW_POLE_GATE_PASS`
- Wynik: `5/6`
- Wynik bieżący:
  - strict non-anchor PASS na kernel-derived łańcuchu EW-pole,
  - test tolerancji pozostaje spełniony (`rel_err_pct = 1.7202` przy tolerancji `5.0%`).

Znaczenie:
- kanał T1 jest teraz closure-ready dla kryteriów strict non-anchor,
- osłona provenance nadal blokuje przyszły anchor leakage.

### 5.24 QW-2087 (bramka alpha_s non-anchor boundary)

- Werdykt: `ALPHA_S_NONANCHOR_BOUNDARY_GATE_PASS`
- Wynik: `8/9`
- Wynik bieżący:
  - strict non-anchor PASS na kernel-derived boundary + punktach walidacyjnych,
  - test tolerancji pozostaje spełniony (`rel_err_pct = 2.1338` przy tolerancji `5.0%`),
  - spójność walidacyjna jest silna (`z_mean = 0.4294`, `z_max = 0.4948`).

Znaczenie:
- promocja strict łączy pochodzenie non-anchor i spójność walidacyjną,
- sama zgodność numeryczna nadal nie wystarcza bez tych warunków.

### 5.25 QW-2093 (executor planu wejść kernel-derived non-anchor)

- Werdykt: `KERNEL_DERIVED_NONANCHOR_INPUTS_BUILT_FROZEN_PLAN`
- Artefakty wyjściowe:
  - `t1_nonanchor_observables_input_qw2085_2086.json`
  - `t1_nonanchor_alpha_s_input_qw2087.json`
- Reguła: zamrożony plan formuł, bez skanu i bez retune.

Znaczenie:
- generacja wejść T1 non-anchor jest teraz deterministyczna i audytowalna,
- status PASS w QW-2085/2086/2087 jest powiązany z reprodukowalnym krokiem budowy wejść.

### 5.26 QW-2094 (sweep usterek w rygorze ścisłym)

- Werdykt: `STRICT_RIGOR_DEFECT_SWEEP_PASS_NO_CRITICAL_DEFECTS`
- Zakres sweepu:
  - spójność wewnętrzna QW-2084/2085/2086/2087/2090/2091/2092/2097/2098,
  - spójność między raportami QW-2069 i QW-2071,
  - spójność pre-gate dla QW-2102 i QW-2103,
  - spójność scalonego preflightu QW-2104,
  - spójność intake/gap dla QW-2105 i QW-2106,
  - spójność guidance design dla QW-2107 i QW-2108,
  - spójność evidence-manifest dla QW-2109,
  - integralność strict-level/status w raporcie pakietowym.
- Wynik: `130` kontroli, `0` błędów.

Znaczenie:
- w audytowanej ścieżce domknięcia T1 non-anchor nie wykryto krytycznych niespójności implementacyjnych,
- maleje ryzyko fałszywych sukcesów wynikających z niespójności raportów.

## 6) Co znaczy "domknięte wewnętrznie" i czego nie znaczy

### 6.1 Znaczy

1. To samo zamrożone jądro przechodzi wiele bramek naraz.
2. Mikro-derywacja i makro-gate już się nie rozjeżdżają.
3. Ślepe metryki external pozostają istotne statystycznie.

### 6.2 Nie znaczy

1. Ostatecznego dowodu ToE uznanego przez niezależne zespoły.
2. Odporności na przyszłą falsyfikację.

Reguła nauki: bardzo mocne dowody wewnętrzne nadal są wewnętrzne, dopóki nie ma niezależnej replikacji.

## 7) Migawka metryk zewnętrznych (QW-2049)

Primary holdout:
- corr = `0.054705` (q95 = `0.036683`, p = `0.006999`)
- gain = `0.000571` (q95 = `-0.001152`, p = `0.006999`)

Stress holdout:
- corr = `0.330346` (q95 = `0.050699`, p = `0.000200`)
- gain = `0.055985` (q95 = `-0.032107`, p = `0.000200`)

Interpretacja licealna:
- sygnał primary jest umiarkowany, ale istotny,
- sygnał stress jest silny i bardzo istotny,
- oba spełniają kryteria bramki.

## 8) Obecna pozycja naukowa a "Święty Graal"

Uczciwe sformułowanie na Release 4.9:

- To jest teraz **kandydat unifikacji o bardzo wysokim rygorze wewnętrznym**.
- To **nie jest jeszcze** finalnie potwierdzone ToE przez społeczność.
- Brakujący krok ma naturę zewnętrznej, niezależnej replikacji multiteam.

## 9) Co mogłoby to teraz obalić

Teoria zostałaby istotnie osłabiona, gdyby niezależne zespoły wykazały:

1. padnięcie bramek przy stałym jądrze i zablokowanym protokole,
2. zanik istotności w blind external,
3. brak odtwarzalności przy zgodnych hashach manifestu.

## 10) Najbliższy krok po Release 4.9

- Wykonać: `RUN_TRULY_INDEPENDENT_MULTITEAM_CONFIRMATORY_PACKAGE`
- Na pakiecie freeze z QW-2050
- Z procedurą rehearsal potwierdzoną przez QW-2051

To jest najkrótsza poprawna naukowo droga od "bardzo mocnego domknięcia wewnętrznego" do "wiarygodnego potwierdzenia zewnętrznego".

## 11) Zasada danych zewnętrznych (bez dużych binarek w gicie)

- Duże surowe paczki danych traktujemy jako źródła zewnętrzne, nie jako artefakty freeze w repo.
- Kanoniczna lista źródeł i komendy pobrania:
  - `DATA_SOURCES_EXTERNAL_DOWNLOADS.md`
- Znaczenie praktyczne:
1. w gicie zamrażamy skrypty/manifesty/raporty,
2. ciężkie archiwa pobieramy z oficjalnych źródeł publicznych,
3. lokalnie sprawdzamy integralność i uruchamiamy niezmienione bramki.

## 12) Czy mamy już wyprowadzenie *wszystkich* znanych wartości fizycznych?

Krótka odpowiedź: **Nie, jeszcze nie w pełnej globalnej wersji ścisłej.**

Co jest już wyprowadzone w rygorze wewnętrznych bramek:
1. cele masowe (mass-chain),
2. cele flavor CKM/PMNS na poziomie gate,
3. faza CP PMNS z deterministycznego operatora flavor,
4. cele dyskryminatorów GW,
5. wsparte mikro-derywacją stałe renormalizacyjne (`Z_beta`, `delta_eta`) po zawężeniu.
6. jawna mapa pakietu i logika bramki domknięcia (`QW-2069`, `QW-2071`).

Czego wciąż brakuje do tezy „wszystkie znane wartości”:
1. luka domknięcia pakietu pozostaje: `28/32` parametrów jest ścisle wyprowadzonych, a `7/32` pozostaje strict-unresolved (mimo że bezpośrednio `missing` to już `0/32`),
2. rygorystyczny frontier dla brakującej 14 pozostaje częściowy (`QW-2081`: `3` kandydaty strict poza tolerancją + `1` non-closing),
3. roadmap domykania istnieje (`QW-2082`) i jest częściowo wykonana (T1/T2 + propagacja EW-secondary), ale T3/T4 oraz CKM CP tolerance nadal pozostają otwarte,
4. integracja statusów epistemicznych (`QW-2083`) domyka mapowanie, ale nie zamienia klas non-closing w strict closure,
5. kompletnego, wyczerpującego pakietu derivacji wszystkich parametrów SM+GR z finalną precyzją,
6. niezależnej replikacji/audytu multiteam poza tym środowiskiem.

Wniosek:
- ścieżka wewnętrznego domknięcia first-principles jest teraz bardzo mocna,
- ale finalne twierdzenie „wszystkie znane wartości fizyczne w pełni wyprowadzone i potwierdzone zewnętrznie” jest jeszcze przedwczesne.

## 13) Addendum Update (2026-03-04, QW-2090..QW-2098)

### EN

- T2 strict non-anchor branch is now closed internally:
  - `QW-2095` input builder PASS,
  - `QW-2088` light-quark gate PASS,
  - `QW-2089` Higgs self-coupling gate PASS,
  - `QW-2096` aggregate T2 gate PASS.
- CKM CP refinement gate (`QW-2097`) was executed in strict deterministic mode and remains a target-miss:
  - verdict: `CKM_CP_TARGET_REFINEMENT_GATE_TARGET_MISS`,
  - selected best strict candidate error: `54.9101%` (tolerance `20%`),
  - no false PASS was promoted.
- Remaining unresolved T3/T4 channels were executed under strict metadata-hardening:
  - `QW-2090`: `H0_LAMBDA_DECOUPLING_GATE_TARGET_MISS` (`7/9`) because strict decoupling candidate exists but H0/Lambda tolerances are not yet met on current H(z) snapshot,
    with new identifiability diagnostic showing weak two-parameter lever arm (`E` span `~0.487`); flatness-projection diagnostic is registry-compatible but non-closing for strict decoupling claim.
  - `QW-2091`: `NEUTRINO_ABSOLUTE_SCALE_GATE_PASS_STRICT` (`8/8`),
  - `QW-2092`: `GNEWTON_SI_BRIDGE_GATE_PENDING_NONCLOSING` (`6/8`) after tautology hardening blocks strict pass for backsolved bridge input.
- Input builders for these channels were additionally hardened:
  - `QW-2099` now reports strict identifiability readiness for H(z) input (`HZ_EXTERNAL_DECOUPLING_AUTOCOLLECTED_WEAK_LEVERARM` on current snapshot),
  - `QW-2101` now supports strict provenance preflight (`--strict-dimensionless-only --omit-g-si-optional --require-strict-ready`) and remains non-strict for backsolved bridge snapshots.
- New strict-input pre-gates (no-overclaim hardening):
  - `QW-2102`: `HZ_DECOUPLING_IDENTIFIABILITY_GATE_WEAK_LEVERARM_PENDING` (`3/7`), failed: `n_nodes_ge_5`, `z_span_ge_0p8`, `e_span_ge_1p0`, `design_condition_lt_8`.
  - `QW-2103`: `GNEWTON_DIMENSIONLESS_PROVENANCE_GATE_PENDING_NONCLOSING` (`5/8`), failed: `bridge_origin_external_dimensionless`, `provenance_anchor_free`, `g_si_not_primary_input`.
- New merged preflight gate:
  - `QW-2104`: `T3T4_STRICT_PREFLIGHT_GATE_PENDING` (`0/8`), combining intake readiness (`QW-2106`) with `QW-2099/2102/2090` and `QW-2101/2103/2092` plus logic-defect checks.
- New deterministic gap quantification:
  - `QW-2105`: `T3T4_STRICT_INPUT_GAPS_PRESENT`, explicitly listing what external inputs are still missing for strict-ready T3/T4.
- New deterministic H(z) acquisition design search:
  - `QW-2107`: `HZ_STRICT_DESIGN_FOUND`, showing strict identifiability can be reached with +2 external nodes and providing auditable redshift candidates (current top pair: `[0.10, 0.90]`).
- New deterministic G dimensionless acquisition specification:
  - `QW-2108`: `GNEWTON_DIMENSIONLESS_ACQUISITION_SPEC_READY`, defining strict external contract for `g_dimensionless_mu_ref` at `mu_ref=1 GeV`: target `6.708830750342e-39`, accepted range `[6.373389212825e-39, 7.044272287859e-39]`.
- New strict evidence-manifest gate for external artifacts:
  - `QW-2109`: `STRICT_EXTERNAL_EVIDENCE_MANIFEST_GATE_PASS` (`29/29`) after sidecar evidence completion (`acquired_utc` + `artifact_sha256` + acquisition protocol/command/analyst and schema/key consistency).
- New deterministic strict sidecar-template builder:
  - `QW-2110`: auto-builds strict sidecar templates with payload hashes for QW-2109:
    `external_hz_nodes_qw2099.metadata.strict.template.json`,
    `external_gnewton_bridge_qw2101.metadata.strict.template.json`.
- New deterministic external acquisition packet for blocked T3/T4 closure:
  - `QW-2111`: `T3T4_STRICT_EXTERNAL_ACQUISITION_PACKET_READY`, operationalizing exact H(z) + direct dimensionless G requirements and rerun protocol.
- New strict H(z) candidate-pack gate:
  - `QW-2112`: `HZ_STRICT_NODE_PACK_PENDING` (`2/12`), adding non-destructive merge + per-node provenance validation + strict threshold checks for candidate external H(z) nodes.
- New strict direct-dimensionless G candidate-pack gate:
  - `QW-2113`: `GNEWTON_DIRECT_DIMENSIONLESS_PACK_PENDING` (`1/16`), enforcing QW-2108 contract for direct external bridge candidate payloads and metadata.
- New raw-input intake gate:
  - `QW-2106`: `STRICT_EXTERNAL_INPUT_INTAKE_GATE_PENDING` (`10/18`), requiring metadata sidecars and strict-ready raw input structure before T3/T4 autocollectors.
- EW secondary propagation gate was executed:
  - `QW-2098`: `EW_SECONDARY_NONANCHOR_CLOSURE_GATE_TARGET_MISS` (`8/10`),
  - strict-derived promotion: `v_higgs`, `sin2_theta_w_mz`,
  - explicit strict target-miss remains: `m_w`, `alpha_em_inv_mz`.
- Current package closure snapshot after integration:
  - `QW-2069`: strict-derived `28/32`, direct missing `0/32`, strict-unresolved `7/32`.
  - `QW-2081`: strict-unresolved in missing-14 scope is `4/14`:
    `delta_cp_ckm`, `h0`, `lambda_cosmological`, `G_newton`.
  - `QW-2094` defect sweep remains clean: `130` checks, `0` failed (including `QW-2102/2103` pre-gates, `QW-2104` merged preflight consistency, `QW-2105/2106` intake-gap consistency, design-guidance consistency for `QW-2107/2108`, evidence-manifest consistency for `QW-2109`, and closure-tooling consistency for `QW-2111/2112/2113`).

### PL

- Gałąź T2 strict non-anchor jest już domknięta wewnętrznie:
  - PASS: `QW-2095` (budowa wejść),
  - PASS: `QW-2088` (lekkie kwarki),
  - PASS: `QW-2089` (samosprzężenie Higgsa),
  - PASS: `QW-2096` (agregat T2).
- Bramka doprecyzowania fazy CKM (`QW-2097`) została wykonana w trybie ścisłym deterministycznym i nadal daje target-miss:
  - werdykt: `CKM_CP_TARGET_REFINEMENT_GATE_TARGET_MISS`,
  - najlepszy kandydat strict ma błąd `54.9101%` (tolerancja `20%`),
  - nie promowano żadnego fałszywego PASS.
- Pozostałe kanały T3/T4 zostały uruchomione z utwardzeniem metadanych:
  - `QW-2090`: `H0_LAMBDA_DECOUPLING_GATE_TARGET_MISS` (`7/9`) bo kandydat strict istnieje, ale tolerancje H0/Lambda nie są jeszcze spełnione dla obecnej migawki H(z),
    z nową diagnostyką identyfikowalności: słaba dźwignia dwupunktowego rozdziału (`E` span `~0.487`); diagnostyka flatness-projection jest zgodna z rejestrem, ale nie domyka ścisłej tezy decoupling.
  - `QW-2091`: `NEUTRINO_ABSOLUTE_SCALE_GATE_PASS_STRICT` (`8/8`),
  - `QW-2092`: `GNEWTON_SI_BRIDGE_GATE_PENDING_NONCLOSING` (`6/8`) po utwardzeniu anty-tautologicznym dla wejścia backsolved.
- Budownicze wejść dla tych kanałów zostały dodatkowo uszczelnione:
  - `QW-2099` raportuje gotowość strict identyfikowalności dla wejścia H(z) (`HZ_EXTERNAL_DECOUPLING_AUTOCOLLECTED_WEAK_LEVERARM` dla obecnej migawki),
  - `QW-2101` wspiera ścisły preflight proweniencji (`--strict-dimensionless-only --omit-g-si-optional --require-strict-ready`) i dla migawki backsolved pozostaje non-strict.
- Nowe pre-gate’y jakości wejścia (utwardzenie anty-overclaim):
  - `QW-2102`: `HZ_DECOUPLING_IDENTIFIABILITY_GATE_WEAK_LEVERARM_PENDING` (`3/7`), fail dla: `n_nodes_ge_5`, `z_span_ge_0p8`, `e_span_ge_1p0`, `design_condition_lt_8`.
  - `QW-2103`: `GNEWTON_DIMENSIONLESS_PROVENANCE_GATE_PENDING_NONCLOSING` (`5/8`), fail dla: `bridge_origin_external_dimensionless`, `provenance_anchor_free`, `g_si_not_primary_input`.
- Nowa scalona bramka preflight:
  - `QW-2104`: `T3T4_STRICT_PREFLIGHT_GATE_PENDING` (`0/8`), łącząca gotowość intake (`QW-2106`) z `QW-2099/2102/2090` oraz `QW-2101/2103/2092` wraz z kontrolą defektów logicznych.
- Nowa deterministyczna kwantyfikacja luk:
  - `QW-2105`: `T3T4_STRICT_INPUT_GAPS_PRESENT`, jawnie wyliczająca, jakich danych zewnętrznych brakuje do strict-ready T3/T4.
- Nowy deterministyczny search projektu akwizycji H(z):
  - `QW-2107`: `HZ_STRICT_DESIGN_FOUND`, pokazujący, że strict identyfikowalność można osiągnąć przy +2 zewnętrznych węzłach i podający audytowalne kandydaty redshift (obecnie para top: `[0.10, 0.90]`).
- Nowa deterministyczna specyfikacja akwizycji bezwymiarowego mostu G:
  - `QW-2108`: `GNEWTON_DIMENSIONLESS_ACQUISITION_SPEC_READY`, definiująca ścisły kontrakt zewnętrzny dla `g_dimensionless_mu_ref` przy `mu_ref=1 GeV`: target `6.708830750342e-39`, zakres akceptacji `[6.373389212825e-39, 7.044272287859e-39]`.
- Nowa bramka strict evidence-manifest dla artefaktów zewnętrznych:
  - `QW-2109`: `STRICT_EXTERNAL_EVIDENCE_MANIFEST_GATE_PASS` (`29/29`) po uzupełnieniu pól dowodowych sidecar i weryfikacji zgodności hash/schema/key.
- Nowy deterministyczny generator strict sidecar-template:
  - `QW-2110`: automatycznie buduje template sidecar z hashami payloadów dla QW-2109:
    `external_hz_nodes_qw2099.metadata.strict.template.json`,
    `external_gnewton_bridge_qw2101.metadata.strict.template.json`.
- Nowy deterministyczny pakiet akwizycji danych dla blokowanych kanałów T3/T4:
  - `QW-2111`: `T3T4_STRICT_EXTERNAL_ACQUISITION_PACKET_READY`, zamienia luki diagnostyczne w operacyjną listę pomiarową i protokół rerun.
- Nowa bramka strict candidate-pack dla H(z):
  - `QW-2112`: `HZ_STRICT_NODE_PACK_PENDING` (`2/12`), dodająca nieinwazyjny merge + walidację proweniencji per-węzeł + kontrolę progów strict dla kandydatów H(z).
- Nowa bramka strict candidate-pack dla bezpośredniego mostu bezwymiarowego G:
  - `QW-2113`: `GNEWTON_DIRECT_DIMENSIONLESS_PACK_PENDING` (`1/16`), wymuszająca kontrakt QW-2108 dla kandydatów payload+metadata.
- Nowa bramka intake dla surowych wejść:
  - `QW-2106`: `STRICT_EXTERNAL_INPUT_INTAKE_GATE_PENDING` (`10/18`), wymagająca metadanych sidecar i struktury strict-ready przed autocollectorami T3/T4.
- Bramka propagacji wtórnych parametrów EW została wykonana:
  - `QW-2098`: `EW_SECONDARY_NONANCHOR_CLOSURE_GATE_TARGET_MISS` (`8/10`),
  - promocja do strict-derived: `v_higgs`, `sin2_theta_w_mz`,
  - jawny strict target-miss pozostaje dla: `m_w`, `alpha_em_inv_mz`.
- Aktualna migawka domknięcia pakietu po integracji:
  - `QW-2069`: strict-derived `28/32`, direct missing `0/32`, strict-unresolved `7/32`.
  - `QW-2081`: strict-unresolved w zakresie missing-14 to `4/14`:
    `delta_cp_ckm`, `h0`, `lambda_cosmological`, `G_newton`.
  - `QW-2094` pozostaje czysty: `130` kontroli, `0` błędów (w tym pre-gate `QW-2102/2103`, spójność scalonej bramki `QW-2104`, spójność intake/gap `QW-2105/2106`, spójność guidance dla `QW-2107/2108`, spójność evidence-manifest dla `QW-2109` oraz spójność narzędzi domykania `QW-2111/2112/2113`).

## 14) Web-Fetch Update (2026-03-04 UTC)

### EN

Post web-fetch strict-ingestion update:

1. `QW-2112` is now `HZ_STRICT_NODE_PACK_READY` (`12/12`).
2. `QW-2113` is now `GNEWTON_DIRECT_DIMENSIONLESS_PACK_READY` (`16/16`).
3. `QW-2099` + `QW-2102` are now strict-ready/pass for H(z) identifiability.
4. `QW-2101` + `QW-2103` + `QW-2092` are now strict-ready/pass for direct dimensionless `G_newton` bridge.
5. `QW-2106` (`18/18`) and `QW-2109` (`29/29`) both pass on strict evidence sidecars.
6. `QW-2105` is now gap-closed ready.
7. `QW-2104` remains pending only because `QW-2090` is still strict target-miss (`7/9`) for `h0/lambda`.
8. `QW-2094` remains clean (`130` checks, `0` failed).

### PL

Po web-fetch i wpięciu danych strict:

1. `QW-2112` ma teraz `HZ_STRICT_NODE_PACK_READY` (`12/12`).
2. `QW-2113` ma teraz `GNEWTON_DIRECT_DIMENSIONLESS_PACK_READY` (`16/16`).
3. `QW-2099` + `QW-2102` są już strict-ready/pass dla identyfikowalności H(z).
4. `QW-2101` + `QW-2103` + `QW-2092` są już strict-ready/pass dla bezpośredniego mostu bezwymiarowego `G_newton`.
5. `QW-2106` (`18/18`) i `QW-2109` (`29/29`) przechodzą z pełnym sidecar evidence.
6. `QW-2105` jest zamknięte jako gap-closed ready.
7. `QW-2104` pozostaje pending wyłącznie dlatego, że `QW-2090` nadal jest strict target-miss (`7/9`) dla `h0/lambda`.
8. `QW-2094` pozostaje czysty (`130` kontroli, `0` błędów).
