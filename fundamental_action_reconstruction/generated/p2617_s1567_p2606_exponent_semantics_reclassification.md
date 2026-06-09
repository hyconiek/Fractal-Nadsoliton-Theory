# P2617/S1567 P2606 exponent-semantics reclassification

Status: `P2617_P2606_EXPONENT_SEMANTICS_RECLASSIFIED_DELTA_PROBE_RETAINED_STRICT_RESIDUAL_EXPORT_QUARANTINED`

## Theorem

The P2603 value 4/5 is the codimension/log-slope delta, not the strict denominator exponent eta. In the strict denominator layer the audited exponent is eta=1+delta=9/5, so a P2606 residual computed with denominator power 4/5 is not the strict-side nonlinear compression residual for K_strict_gate.

## Proof

- P2603 supplies the logarithmic damping slope/codimension delta = D_f - 1 = 4/5.
- The strict denominator target audited by P2414 is S(d)=1+d^(9/5), so the denominator power is eta_strict=9/5.
- The relation between the linear backbone and the codimension slope is eta_strict = 1 + delta = 1 + 4/5 = 9/5.
- P2606 used the label eta=4/5 inside the denominator power, thereby using delta where the strict kernel notation reserves eta for the full denominator exponent.
- Therefore P2606 may be retained only as a numerical codimension-slope perturbation probe, not as an exported strict-side nonlinear residual addition for the strict denominator.

## Computed checks

- `eta_strict = 1 + delta`: `True`.
- P2606 used delta as denominator eta: `True`.
- Max strict eta=9/5 vs P2606 eta=4/5 kernel difference: `0.3520850545472334`.
- P2606 strict residual interpretation quarantined: `True`.
- P2606 codimension probe retained: `True`.

## Scope guards

P2617 reclassifies only the P2606 exponent semantics. It does not revalidate P2606 as a strict-side residual addition, does not revalidate P2607/P2608, and does not export QW-2191 discharge, APD sourcehood, legacy role transfer, or ToE closure.

## Fingerprint

`519305cbc95f9e5427497c7c6aadccbb42c66546ed33436796e9dc122f707027`
