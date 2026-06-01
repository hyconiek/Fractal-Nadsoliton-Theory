# P2369 S1319: self-recorded ledger closed-form uniqueness theorem

Status: `OPEN_PROGRESS_CLOSED_FORM_LEDGER_UNIQUENESS_NO_QW2191_DISCHARGE`

## Result

P2369/S1319 replaces the P2368 ledger enumeration with a closed-form integer proof: convex smoothing gives the ripple lower bound, and the arrow penalty uniquely selects the nonincreasing ordered minimizer.

## Certificate

- Ripple lower bound: `14`.
- Ripple minimizer permutations: `10`.
- Brute-force cross-check over 35 compositions: `True`.
- Arrow tiebreak unique winner: `True` with `[[2, 2, 2, 1, 1]]`.
- Endpoint values distinct: `True`.
- Exact eta identity holds: `True`.

## Hard limits

- This proves finite ledger uniqueness, not strict derivation of the ordered d5 support or arrow action.
- No `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.
