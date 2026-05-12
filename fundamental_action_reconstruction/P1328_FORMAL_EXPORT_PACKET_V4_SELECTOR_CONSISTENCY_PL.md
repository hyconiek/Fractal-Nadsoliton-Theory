# P1328 — Formal export packet: `V4_selector_consistency` and residual obligations (PL)

Status: `FORMAL_EXPORT_PARTIAL_ONE_LEMMA_EXPORTED`
As of: `2026-05-12`
Depends on: `P1327`

## Cel
Zrealizować rekomendację z `P1327` i wyeksportować pierwszy theorem-level
blok dla O3-EXCLUSION: lemma `V4_selector_consistency`.

## Lemma L1 — `V4_selector_consistency`
**Statement (roboczy formalny):**
Dla klasy kandydatów opisywanych przez `R_common_v1` oraz selektora
`S_sel_strict_v4`, przy dopuszczalnych transformacjach `G_allowed_v1`, znak
kierunku jest zgodny z wynikiem selektora i stabilny lokalnie w granicy
perturbacyjnej użytej w `P1326`.

**Zakres eksportu:**
- obejmuje zgodność mapowania `R_common_v1 -> sign` dla `v4`,
- obejmuje lokalną stabilność (flip-rate bound) potwierdzoną w replay/adversarial,
- nie obejmuje jeszcze pełnej eliminacji residual loophole `open(Z2/eps)`.

**Status L1:** `EXPORTED_FORMAL_STATEMENT_WITH_PROOF_SKETCH`.

## Lemma L2 — `Residual_loophole_elimination`
**Statement target:**
Pokazać, że residual `open(Z2/eps)` nie tworzy legalnej alternatywnej gałęzi
po integracji z `v4` i constraints strict-core.

**Status L2:** `PENDING_NOT_EXPORTED`.

## Decyzja profesorska
- Mamy pierwszy formalny eksport (`L1`) i to jest realny postęp.
- Brak `L2` blokuje strict closure.
- `QW-2191` strict pozostaje `NOT_CLOSED`.

## Konsekwencja
Pipeline formalny przesuwa się z "zero lemmas" do "1/2 lemmas" wymaganych
przez finalną bramkę O3.

## Czego dokument nie twierdzi
- Nie twierdzi strict closure.
- Nie twierdzi ToE closure.
- Nie twierdzi, że lokalna stabilność numeryczna sama zamyka residual loophole.

## Rekomendowany następny uczciwy krok
Uruchomić **P1329 residual-loophole elimination proof attempt**
(z explicit case split dla `z2/eps`) i wyeksportować L2 albo uczciwie wykazać,
że L2 nadal nie jest gotowe.

## Dla laika
To tak, jakbyśmy dołożyli pierwszy brakujący dokument prawny: wiemy już
formalnie, że nowy kompas działa spójnie i stabilnie lokalnie. Ale drugi
kluczowy dokument (że nie ma ukrytej alternatywnej ścieżki) nadal czeka.
