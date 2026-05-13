# P1444 — Branch B Conditional Projection Check (Non-Strict)

Status: `P1444_CONDITIONAL_CHECK_RECORDED`
As of: `2026-05-13`

## Professorial decision

Execute a conditional projection check under Branch B non-strict premise,
and keep strict/non-strict outputs explicitly separated.

## Result

`PASS_CONDITIONAL_NON_STRICT_PROJECTION_CHECK`

Interpretation:

- This PASS is local and conditional.
- It does not close strict-core `QW-2191`.
- It does not authorize global strict closure claim.

## Artifact

- `generated/p1444_branch_b_conditional_projection_check_summary.json`

## PL (skrót)

Warunkowy test projekcji przeszedł tylko w gałęzi non-strict.
To nie jest domknięcie strict-core.
