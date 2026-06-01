# P2387 S1337: bathtub exact KKT branch certificate

Status: `OPEN_PROGRESS_EXACT_KKT_BRANCH_CERTIFICATE_SOURCE_STILL_OPEN`

## Result

P2387/S1337 upgrades the P2386 LP-dual audit from sampled KKT evidence to a branchwise KKT certificate.
Given the already-certified monotonicity `q'(s)<0`, the cut `t=1/M=0.625` gives `q(s)>lambda` for `s<t`, `q(s)=lambda` at the cut, and `q(s)<lambda` for `s>t`.
Therefore `mu=max(q-lambda,0)` and `rho*=M*1_[0,t)` make dual feasibility and complementarity algebraic branch identities.

## Checks

- `q(0)-lambda`: `9.507655486694514`.
- `lambda-q(1)`: `0.2125807556876289`.
- Closed value gap: `0.0`.
- Audit-grid max `q_prime`: `-0.24174708731451278`.

## Hard limits

- This is an exact KKT acceptance certificate, not a strict source theorem deriving `rho*` or `M`.
- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.
