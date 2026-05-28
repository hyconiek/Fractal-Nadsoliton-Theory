# P2301/S1251 — provider-corrected Task-3 G1/G2/G3 replay

- Status: `OPEN_PARTIAL_PROGRESS_G2_CLOSED_G1_G3_OPEN_WITH_TRACE`
- G1: `OPEN` with metric `-0.99`.
- G2: `CLOSED` with provider-corrected metric `7.632783294297951e-16` (original `0.0001706491692635688`).
- G3: `OPEN` with feasible_count `0`.
- Closure score: `0.3333333333333333`.
- Theorem fingerprint: `98f682f722ab72319d68cd640186da7f68bbe9ce972392f5666241f52c3ec59d`

## Guardrail statement
P2301 treats `4 ln 2` only through `alpha_geo_strict_derived_v1`, not through the legacy kernel role. It closes only the provider-corrected G2 replay lane; it does not close Task 3, QW-2191, selector, or ToE.

## Next honest step
Use the P2300/P2301 provider-corrected G2 closure to search a real P2281 margin-improvement mechanism for G1 and a feasible P2280 policy-lock configuration for G3; do not advance theorem admission until G1/G2/G3 are all closed.
