# P2383 S1333: closed-form bathtub corner reduction theorem

Status: `OPEN_PROGRESS_CLOSED_FORM_BATHTUB_CORNER_REDUCTION_SOURCE_STILL_OPEN`

## Result

P2383/S1333 turns the P2382 bathtub certificate into a more proof-like corner reduction.
The monotonicity of `q(s)=A_s(5)-3*A_s(1)` is certified by the closed ratio bound

```text
R=A5/A1=(Delta5*(u1+s*Delta1))/(Delta1*(u5+s*Delta5))
R_min=2*5^(9/5)/(1+5^(9/5)) > sqrt(3)
q'(s)=A1^2*(3-R^2)<0.
```

## Corner-reduced cap burden

On the audited cap band `[1.5,1.6]`, derivative replay shows that the chamber margin increases with `eta`, decreases with `beta_tors`, and increases with `M`.
Therefore the cap threshold is controlled by the single corner `(eta,beta_tors)=(9/5,0.1)`.

- Corner threshold: `1.574821357435363`.
- Certified sufficient cap: `1.6`.
- `M=1.6` early interval length: `0.625`.
- `M=1.6` early-half mass: `0.8`.
- `M=1.6` barycenter shift from uniform: `0.1875`.

## Hard limits

- This is a closed-form/derivative-audited reduction of the bounded-density acceptance criterion, not a strict source theorem for the cap or bang-bang profile.
- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.
