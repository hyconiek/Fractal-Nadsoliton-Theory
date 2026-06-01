# P2385 S1335: exact Z12 support chamber theorem

Status: `OPEN_PROGRESS_EXACT_Z12_SUPPORT_CHAMBER_THEOREM_SOURCE_STILL_OPEN`

## Result

P2385/S1335 separates the finite support-selection theorem from the analytic cap proof.
For 5-subsets of `Z12` scored by `a*h1+b*h5`, the chamber

```text
b>0, a>=0, a/b<1/3
```

has exactly the 12 `(h1,h5)=(0,4)` supports as unique maximizers.
These supports are precisely the length-5 paths in the step-5 cycle on `Z12`.

## Finite certificate

- Supports checked: `792`.
- Full pair distribution: `{'0,0': 12, '0,2': 12, '0,4': 12, '1,1': 72, '1,2': 96, '1,3': 72, '2,0': 12, '2,1': 96, '2,2': 204, '2,3': 48, '3,1': 72, '3,2': 48, '3,3': 24, '4,0': 12}`.
- Target support count: `12`.
- Boundary tie pair at `a/b=1/3`: `[{'h1': 3, 'h5': 3, 'count': 24, 'target_pair_0_4': False, 'score_gap_from_0_4_over_b_at_r_equals_1_3_numerator': 0, 'strictly_beaten_for_r_lt_1_3': True, 'boundary_tie_pair_at_r_1_3': True}]`.

## Hard limits

- This is an exact finite chamber theorem for the support-scoring layer, not a strict source theorem for the cap or bang-bang profile.
- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.
