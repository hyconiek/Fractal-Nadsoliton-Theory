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
3. GW discriminator targets,
4. micro-supported renormalization constants (`Z_beta`, `delta_eta`) with tightening.

What is still open for a full “all known values” claim:
1. complete precision radiative program across all sectors,
2. full exhaustive Standard Model parameter derivation package at final precision,
3. independent external multiteam replication/audit.

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
3. cele dyskryminatorów GW,
4. wsparte mikro-derywacją stałe renormalizacyjne (`Z_beta`, `delta_eta`) po zawężeniu.

Czego wciąż brakuje do tezy „wszystkie znane wartości”:
1. pełnego precyzyjnego programu radiacyjnego we wszystkich sektorach,
2. kompletnego, wyczerpującego pakietu derivacji wszystkich parametrów SM z finalną precyzją,
3. niezależnej replikacji/audytu multiteam poza tym środowiskiem.

Wniosek:
- ścieżka wewnętrznego domknięcia first-principles jest teraz bardzo mocna,
- ale finalne twierdzenie „wszystkie znane wartości fizyczne w pełni wyprowadzone i potwierdzone zewnętrznie” jest jeszcze przedwczesne.
