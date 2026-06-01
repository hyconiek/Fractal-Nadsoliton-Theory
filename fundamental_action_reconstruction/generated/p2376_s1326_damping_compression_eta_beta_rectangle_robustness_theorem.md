# P2376 S1326: damping-compression eta-beta rectangle robustness theorem

Status: `OPEN_PROGRESS_DAMPING_COMPRESSION_ETA_BETA_RECTANGLE_ROBUST_NO_QW2191_DISCHARGE`

## Result

P2376/S1326 proves that the P2374/P2375 compression-polarity candidate is robust on `eta in [9/5,2]` and `beta_tors in [0,0.1]`.

## Certificate

- Margin: `F(eta,x)=(1+5^eta)*(1+x)^3 - 8*(1+5*x) > 0`.
- Minimum corner margin: `11.119491591942388`.
- Minimum corner dF/dx: `17.358474775827162`.
- Minimum corner dF/deta: `29.16219672210299`.
- Grid support scans select: `[{'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}]`.

## Hard limits

- This is a bivariate robustness theorem for a candidate feature, not a strict dynamical source theorem.
- No `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, L_total promotion, or ToE closure is claimed.
