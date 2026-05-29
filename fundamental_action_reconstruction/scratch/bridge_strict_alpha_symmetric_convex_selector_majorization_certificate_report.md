# Scratch strict-alpha symmetric-convex selector majorization certificate probe

Status: majorization/convexity selector certificate for eta=9/5; no physical selector theorem.

- Conditional theorem: for fixed positive integer `m,n`, every canonical ledger majorizes the floor/ceil balanced ledger.
- Selector class: any symmetric strictly convex / Schur-convex branch action selects that balanced ledger.
- Target result: for `m=5, n=8`, checked powers `[2, 3, 4, 5, 6]` uniquely select `[2, 2, 2, 1, 1]` and `q^5=256/243`.
- Bounded scan: `180` `(m,n)` cases and `5583` canonical ledgers checked; all majorization/convex-power checks passed: `True`.
- Honest read: this reduces the selector source obligation to deriving a symmetric strictly convex branch action; it does not derive that action from strict geometry.
- No false pass: no strict convex-action source theorem, no QW-2191 discharge, no ToE closure.
