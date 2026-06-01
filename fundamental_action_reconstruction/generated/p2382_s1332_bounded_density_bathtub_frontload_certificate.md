# P2382 S1332: bounded-density bathtub frontload certificate

Status: `OPEN_PROGRESS_BOUNDED_DENSITY_BATHTUB_FRONTLOAD_CERTIFICATE_SOURCE_STILL_OPEN`

## Result

P2382/S1332 generalizes the P2379-P2381 affine profile work to all normalized densities with a pointwise cap `0<=rho<=M`.
Because the audited contrast `q(s)=A_s(5)-3*A_s(1)` is strictly decreasing, the bathtub maximizer is the earliest bang-bang profile.

## Certificate

- Worst-corner cap threshold: `1.574821357435363`.
- Certified sufficient cap: `1.6`.
- Worst-corner early interval length for `M=1.6`: `0.625`.
- Below-threshold negative control selects: `{'3,3': 24}`.
- Cap `M=1.6` replay selects: `{'0,4': 12}`.

## Hard limits

- This is a bounded-density variational reduction/acceptance certificate, not a strict source theorem deriving the density cap or bang-bang profile.
- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.
