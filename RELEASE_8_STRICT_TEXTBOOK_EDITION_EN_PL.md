# RELEASE 8 STRICT TEXTBOOK EDITION (EN + PL)

**Version:** 8.0.0 (Current-state projection — O3 strict gate closure packetized)  
**Date:** 2026-05-12  
**Branch:** `main`  
**Predecessor:** Release 7 Strict Textbook Edition

## Status discipline (no false pass)

- This Release 8 document describes the **current repo state only**.
- It records that the O3 pipeline (`P1313..P1339`) reaches strict closure **inside the declared pipeline rigor** (`P1338` + `P1339`), with explicit machine reports.
- It does **not** claim global ToE closure in the strongest sense.
- It does **not** erase the strict/non-strict distinction in historical packets; non-strict packets remain explicitly labeled.
- It does **not** silently transfer legacy kernel roles into strict kernel claims beyond what is explicitly exported.

---

## 1. Release‑8 target (EN): “Current-state strict O3 closure snapshot”

Release 8 is a state report whose target is:

1. freeze the current strict O3 gate result for `QW-2191`,
2. show exactly which obligations were passed and how,
3. preserve hard limits and scope boundaries,
4. provide a clean handoff to post-closure regression/replication governance.

This release is **not** a restart of theory-building; it is a disciplined snapshot of the executed `P1313..P1339` chain.

---

## 2. Current strict state entering Release 8 (EN)

On the current repo state:

1. O3 plan and step chain are fully packetized (`P1313..P1339`) with explicit pass/fail discipline.
2. Candidate-class and normalization layer are frozen (`P1314`, `P1315`).
3. Pairwise and sweep/replay stages were executed with artifacts (`P1316`, `P1317`, `P1318`).
4. Residual-slot neutrality initially failed (`P1319`), forcing further construction.
5. Non-strict closure lane under explicit premise was separated and labeled (`P1320`).
6. Candidate evolution proceeded (`v2` falsified, `v3` degenerate, `v4` robust under revised gate) (`P1323..P1326`).
7. Formal export lane progressed from `1/2` to `2/2` via final global-L2 attempt (`P1328`, `P1338`).
8. Full checker refresh reached `5/5` obligations (`P1337` refreshed after `P1338`).
9. Independent post-closure replication packet passed (`P1339`).

Current pipeline-level status:

- `o3_strict_ready = true`
- `qw2191_strict_status = CLOSED`

within the declared Release-8 strict pipeline semantics.

---

## 3. What Release‑8 still does not “false‑PASS” (EN)

Release 8 still does **not** claim, unless separately exported:

1. strongest-form global ToE closure beyond this O3 gate scope,
2. retrospective reinterpretation of previously labeled non-strict packets as strict,
3. silent replacement of theorem-level export by numeric success alone,
4. immunity from future independent falsification.

---

## 4. Release‑8 deliverables (EN)

1. This Release-8 textbook snapshot (`RELEASE_8_STRICT_TEXTBOOK_EDITION_EN_PL.md`).
2. Packet chain `P1313..P1339` as the current strict O3 trail.
3. Generated machine reports under `fundamental_action_reconstruction/generated/` documenting each step.
4. Post-closure next-step target: regression/replication guardrail packet (`P1340`).

---

## 1. Cel Release‑8 (PL): „Migawka aktualnego stanu strict po O3”

Release 8 ma jeden cel: **uczciwie opisać obecny stan repo** po wykonaniu łańcucha `P1313..P1339`, bez dopisywania nowych roszczeń ponad to, co jest wyeksportowane.

---

## 2. Aktualny stan strict na wejściu do Release 8 (PL)

Na obecnym stanie repo:

1. Plan O3 i pełna ścieżka kroków są jawnie zapakowane (`P1313..P1339`).
2. Warstwa klas kandydatów i normalizacji jest zamrożona (`P1314`, `P1315`).
3. Pairwise + sweep/replay zostały wykonane i udokumentowane (`P1316`, `P1317`, `P1318`).
4. Pierwotna próba neutralności residual-slot nie przeszła (`P1319`), co wymusiło dalszą konstrukcję.
5. Tor non-strict pod jawną premise został oddzielony i oznaczony (`P1320`).
6. Kandydaci ewoluowali: `v2` obalony, `v3` zdegenerowany, `v4` przeszedł poprawioną bramkę (`P1323..P1326`).
7. Tor formalnych eksportów przeszedł z `1/2` do `2/2` przez finalną próbę globalnego L2 (`P1328`, `P1338`).
8. Pełny checker po odświeżeniu pokazuje `5/5` obowiązków (`P1337` po `P1338`).
9. Niezależna replikacja post-closure przeszła (`P1339`).

Bieżący status w semantyce tego pipeline:

- `o3_strict_ready = true`
- `qw2191_strict_status = CLOSED`

---

## 3. Czego Release‑8 nadal nie może „false‑PASS” (PL)

Release 8 nadal **nie** rości bez osobnego eksportu:

1. najsilniejszego globalnego domknięcia ToE poza zakresem tej bramki O3,
2. cofnięcia etykiet non-strict w starszych pakietach,
3. zastąpienia eksportu theorem-level samym wynikiem numerycznym,
4. gwarancji, że przyszły niezależny audyt niczego nie podważy.

---

## 4. Deliverables Release‑8 (PL)

1. Ten dokument Release 8 (`RELEASE_8_STRICT_TEXTBOOK_EDITION_EN_PL.md`).
2. Łańcuch pakietów `P1313..P1339` jako aktualny strict O3 trail.
3. Raporty maszynowe w `fundamental_action_reconstruction/generated/` dla każdego kroku.
4. Następny krok po domknięciu: pakiet guardrail regresji/replikacji (`P1340`).
