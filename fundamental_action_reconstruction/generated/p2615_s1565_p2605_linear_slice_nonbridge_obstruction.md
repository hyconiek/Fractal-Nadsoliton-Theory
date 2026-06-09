# P2615/S1565 P2605 linear-slice non-bridge obstruction

Status: `P2615_P2605_LINEAR_SLICE_RETAINED_ONLY_AS_BOUNDARY_NEGATIVE_CONTROL_NONBRIDGE_OBSTRUCTION_RECORDED`

## Theorem

For a node-preserving constant-beta denominator map 1+beta_tors*d = 1+beta*d^eta, equality at any two distinct positive nodes a,b forces eta=1 and beta=beta_tors.

## Exact proof

- Cancel the shared unit term to obtain beta_tors*d = beta*d^eta at each audited positive node.
- At nodes a and b, divide the two equations: a/b = (a/b)^eta because beta and beta_tors are nonzero.
- Since a/b is positive and not equal to 1, taking logs gives (eta-1)*log(a/b)=0, hence eta=1.
- Substituting eta=1 back into beta_tors*d = beta*d gives beta=beta_tors.
- Therefore the P2605 eta=1 slice is an exact boundary equality only; it cannot by itself be a bridge to any nonlinear strict-side compression exponent eta != 1.

## Computed checks

- P2605 quarantine retained: `True`.
- P2606 nonlinear residual component retained: `True`.
- Non-bridge damping source revalidation inherited from P2613/P2614: `True`.
- Remaining quarantines after P2615: `['P2605', 'P2607', 'P2608']`.

## Scope guards

P2615 does not revalidate P2605 as a full bridge, does not revalidate P2607, does not re-enable P2608 role-bearing L_total, and does not export QW-2191 discharge, APD sourcehood, legacy physical-role transfer, or ToE closure.

## Fingerprint

`16a0f97a65393f80599388dd74fce83a9e4f9fb49e1dd5f3c5ecf7f7bd1a9afb`
