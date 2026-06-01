# P2384 S1334: symbolic bathtub inequality proof packet

Status: `OPEN_PROGRESS_SYMBOLIC_BATHTUB_INEQUALITY_PROOF_SOURCE_STILL_OPEN`

## Result

P2384/S1334 extracts an inequality-level proof core from P2383.
The ratio part uses only the coarse bound `5^eta >= 5^(3/2)=sqrt(125)` and the algebraic threshold `3+2*sqrt(3)`.
The derivative part records closed identities for the eta, beta_tors, and M directions and supplies positive coarse endpoint gaps on the cap band `[1.5,1.6]`.

## Certified coarse gaps

- `sqrt(125)-(3+2*sqrt(3)) = 4.7162382723611938949909756606446364433174812904368673652428722681487385943714240`.
- beta-direction gap `N5-3*N1 >= 0.7039814949405757301853178046162302479753196963436893275333213168479736341442150`.
- cap-direction gap `h(x5)-3*h(x1) >= 0.54664619303656935396159122398125481702355274557740194895287403066171062071284188`.
- `M=1.6` early interval length: `0.625`.
- `M=1.6` early-half mass: `0.8`.
- `M=1.6` barycenter shift from uniform: `0.1875`.

## Hard limits

- This is an inequality proof packet for the bounded-density acceptance criterion, not a strict source theorem for the cap or bang-bang profile.
- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.
