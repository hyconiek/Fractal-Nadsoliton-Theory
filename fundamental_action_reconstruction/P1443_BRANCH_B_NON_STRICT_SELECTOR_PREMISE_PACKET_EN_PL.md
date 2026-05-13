# P1443 — Branch B / Explicit Non-Strict Selector Premise

Status: `P1443_NON_STRICT_BRANCH_DECLARED`
As of: `2026-05-13`

## Professorial decision

After strict Branch A v4 failure (`P1442`), activate Branch B explicitly as non-strict fallback.

## Contract

- Premise is explicitly labeled non-strict.
- It may enable conditional projection checks only.
- It does **not** close strict-core `QW-2191`.
- It does **not** permit global strict closure claims.

## Artifact

- `generated/p1443_branch_b_non_strict_selector_premise_summary.json`

## PL (skrót)

Uruchomiono gałąź B jako jawnie non-strict: wolno robić tylko warunkowe testy projekcji,
bez twierdzenia, że strict-core jest domknięty.
