# P1441 — Selector Source Bifurcation (Post-V3)

Status: `P1441_BIFURCATION_DECLARED_NO_FALSE_PASS`
As of: `2026-05-13`

## Professorial decision

After `P1440` fail, declare an explicit bifurcation:

1. **Branch A (strict)**: V4 strict internal selector-source attempt,
2. **Branch B (non-strict)**: explicit selector premise branch.

No branch may claim strict-core global closure at this stage.

## Why this is honest

- `QW-2191` remains open.
- We avoid mixing strict and non-strict evidence under one label.
- No legacy bridge import is admitted.

## Artifact

- `generated/p1441_selector_source_bifurcation_summary.json`

## PL (skrót)

Po niepowodzeniu v3 robimy jawne rozgałęzienie: strict-v4 vs gałąź non-strict.
Bez udawania, że problem QW-2191 jest domknięty.
