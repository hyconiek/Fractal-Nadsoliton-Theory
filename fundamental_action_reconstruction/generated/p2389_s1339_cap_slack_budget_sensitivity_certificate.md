# P2389 S1339: cap slack budget sensitivity certificate

Status: `OPEN_PROGRESS_CAP_SLACK_BUDGET_SENSITIVITY_CERTIFICATE_SOURCE_STILL_OPEN`

## Result

P2389/S1339 quantifies how much room the accepted cap `M=1.6` has above the unique P2388 cap threshold.
It records the cap surplus, scalar margin `F(1.6)`, derivative/sensitivity band, and source-geometry surplus against the just-threshold profile.

## Checks

- Cap surplus `1.6-root`: `0.02517864256463742`.
- Scalar margin `F(1.6)`: `0.02569532538409236`.
- Mean-value slope: `1.0205206781155312`.
- Root sensitivity interval: `{'dM_dT_lower_1_over_max_F_prime': 0.976632666780895, 'dM_dT_upper_1_over_min_F_prime': 0.9831667514907347, 'meaning': "for a small additive increase in the scalar target, the unique cap root moves by approximately dT/F'(M*)"}`.
- Early interval shortening vs threshold: `0.009992658233011209`.
- Early-half mass surplus vs threshold: `0.01258932128231871`.

## Hard limits

- This is a slack/sensitivity acceptance budget, not a strict source theorem deriving `M=1.6` or the density.
- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.
