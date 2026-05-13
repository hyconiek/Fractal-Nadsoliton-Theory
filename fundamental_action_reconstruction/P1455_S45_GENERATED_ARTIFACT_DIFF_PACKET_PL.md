# P1455 — S4.5 Generated Artifact Diff Checkpoint (PL)

Status: `P1455_EXECUTED_LOCAL_ONLY_DIFF_NO_GLOBAL_CLAIM`
As of: `2026-05-13`

## Cel

Wykonać uczciwy, lokalny checkpoint różnic artefaktów `generated/*.json|*.csv`
między bazą `P1454` a aktualnym stanem po kolejnych krokach, bez jakichkolwiek
claimów global closure i bez legacy bridge.

## Wejście

- baseline manifest: `generated/p1454_generated_artifact_manifest.json`,
- skan bieżących artefaktów: `generated/*.json|*.csv`.

## Procedura (strict-only)

1. Odczytaj baseline `P1454`.
2. Zbuduj bieżący manifest (`P1455`) tylko dla `.json/.csv`.
3. Wylicz:
   - `added`,
   - `removed`,
   - `size_changed`.
4. Wyeksportuj:
   - `p1455_generated_artifact_manifest.json`,
   - `p1455_s45_generated_artifact_diff_summary.json`,
   - `p1455_s45_generated_artifact_diff_rows.csv`.

## Dyscyplina naukowa

- `scope = LOCAL_ONLY_NON_GLOBAL_CLAIM`,
- `strict_core_qw2191_closed = false`,
- `legacy_bridge_used = false`.

## Werdykt

- `PASS_DIFF_STABLE` — brak zmian,
- `PASS_DIFF_WITH_CHANGES` — zmiany wykryte i jawnie wylistowane.

W obu przypadkach to wyłącznie kontrola artefaktów; brak globalnych wniosków
fizycznych.
