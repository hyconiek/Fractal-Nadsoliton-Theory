# P2388 S1338: cap-threshold root uniqueness certificate

Status: `OPEN_PROGRESS_UNIQUE_CAP_THRESHOLD_ROOT_CERTIFICATE_SOURCE_STILL_OPEN`

## Result

P2388/S1338 turns the worst-corner cap threshold into a unique-root certificate for the scalar equation

```text
F(M)=W_M(5)-3*W_M(1)-(3*K_strict(1)-K_strict(5))=0.
```

The sign-changing bracket plus the P2384 positive cap-derivative proof gives uniqueness; bisection and Newton are retained as reproducible computations rather than as the theorem core.

## Checks

- Bracket: `[1.5748213564353628, 1.574821358435363]`.
- Root midpoint: `1.5748213574353627`.
- Bisection width: `2.220446049250313e-16`.
- Newton final residual: `1.1102230246251565e-15`.
- `F(1.6)`: `0.02569532538409236`.

## Hard limits

- This is a unique threshold acceptance certificate, not a strict source theorem deriving `M` or the density.
- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.
