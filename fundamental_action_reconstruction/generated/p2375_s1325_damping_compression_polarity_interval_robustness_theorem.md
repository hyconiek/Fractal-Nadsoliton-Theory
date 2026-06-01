# P2375 S1325: damping-compression polarity interval robustness theorem

Status: `OPEN_PROGRESS_DAMPING_COMPRESSION_POLARITY_INTERVAL_ROBUST_NO_QW2191_DISCHARGE`

## Result

P2375/S1325 proves that the P2374 compression-polarity candidate is robust over `beta_tors in [0,0.1]`.

## Certificate

- Equivalent margin: `F(x)=(1+5^(9/5))*(1+x)^3 - 8*(1+5*x) > 0`.
- Margin minimum at x=0: `11.119491591942388`.
- Derivative minimum at x=0: `17.358474775827162`.
- Canonical C1/C5: `0.2354294000560429`.
- Interval endpoint beta_tors=0.1 C1/C5: `0.2348840372868392`.
- All sampled support scans select: `[{'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}]`.

## Hard limits

- This is an interval robustness theorem for a candidate feature, not a strict dynamical source theorem.
- No `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, L_total promotion, or ToE closure is claimed.
