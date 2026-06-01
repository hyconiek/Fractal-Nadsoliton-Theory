# P2386 S1336: bathtub LP dual certificate

Status: `OPEN_PROGRESS_PROOF_SIDE_LP_DUAL_CERTIFICATE_SOURCE_STILL_OPEN`

## Result

P2386/S1336 rewrites the bounded-density bathtub step as a primal/dual LP certificate.
The primal maximizes `int q(s)*rho(s) ds` under `0<=rho<=M` and `int rho=1`; the dual minimizes `lambda + M*int mu(s) ds` under `lambda+mu(s)>=q(s)` and `mu>=0`.
For `M=1.6`, the certificate uses `lambda=q(1/M)` and `mu(s)=max(q(s)-lambda,0)`.

## Numerical/closed-form checks

- Early cut `1/M`: `0.625`.
- Primal value: `1.4115211199555233`.
- Dual value: `1.4115211199555233`.
- Closed-form `W_M(5)-3*W_M(1)`: `1.411521119955523`.
- Absolute primal-dual gap: `0.0`.
- Sampled KKT complementarity max error: `7.105427357601002e-16`.

## Hard limits

- This is a proof-side LP/KKT certificate for the bounded-density acceptance criterion, not a strict source theorem deriving the cap or density.
- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.
