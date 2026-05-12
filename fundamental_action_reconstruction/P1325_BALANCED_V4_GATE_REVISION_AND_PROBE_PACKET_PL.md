# P1325 — Balanced-v4 gate revision and probe packet (PL)

Status: `V4_PASSES_REVISED_NUMERICAL_GATE_NOT_YET_THEOREM_LEVEL`
As of: `2026-05-12`
Depends on: `P1324`

## Cel
Wykonać następny uczciwy krok po `P1324`: naprawić sprzeczną bramkę
(robust + informatywność) i sprawdzić nowego kandydata `S_sel_strict_v4` na
spójnym kryterium.

## Problem z P1324 (diagnoza)
W `P1324` użyto warunku `negative_count = 0` i jednocześnie
`sign_diversity = 2`, co jest wewnętrznie sprzeczne (brak znaków ujemnych
uniemożliwia dywersyfikację znaków).

## Rewizja bramki
Nowa bramka (`v2`) dla kandydata selektora:
1. **informatywność:** `sign_diversity = 2`,
2. **robustność lokalna:** `flip_rate <= 0.08` pod małą perturbacją.

## Artefakt wykonawczy
- skrypt: `p1325_balanced_v4_gate_revision_and_probe.py`
- raport: `generated/p1325_balanced_v4_gate_revision_and_probe_report_v1.json`

## Wynik
- `positive_count = 193`, `negative_count = 207`, `sign_diversity = 2`,
- `flip_rate = 0.03`,
- status: `CANDIDATE_V4_PASSES_REVISED_GATE`.

## Decyzja profesorska
To jest pierwszy kandydat (`v4`), który na obecnym teście jest jednocześnie:
- nietrywialny (rozróżnia znaki),
- i lokalnie stabilny.

Ale to **jeszcze nie jest** domknięcie `QW-2191` strict, bo brakuje:
1. niezależnego replay/adversarial dla `v4`,
2. theorem-level derivation strict-side,
3. integracji z pełnym O3-EXCLUSION bez residual loopholes.

## Konsekwencja
- Program przesuwa się jakościowo do etapu "kandydat mocny numerycznie".
- `QW-2191` strict pozostaje `NOT_CLOSED` do czasu domknięcia punktów 1–3.

## Czego dokument nie twierdzi
- Nie twierdzi strict closure.
- Nie twierdzi ToE closure.
- Nie twierdzi, że test numeryczny zastępuje dowód theorem-level.

## Rekomendowany następny uczciwy krok
Uruchomić **P1326 v4 independent replay + adversarial + O3 reintegration**:
- powtórzyć testy wieloseedowe,
- przeprowadzić ataki adversarial na granicy flip-rate,
- włączyć wynik do formalnego O3-EXCLUSION update.

## Dla laika
Poprzedni test miał wadliwą regułę. Naprawiliśmy ją i nowy kompas (`v4`) wypadł
pierwszy raz dobrze: potrafi odróżniać kierunki i jest stabilny na małych
zakłóceniach. To bardzo dobry sygnał, ale jeszcze trzeba go sprawdzić ostrzej i
wpiąć w pełny dowód.
