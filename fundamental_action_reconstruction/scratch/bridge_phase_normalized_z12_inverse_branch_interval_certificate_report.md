# Z12 inverse-branch interval certificate

Status: `finite-z12-monotone-inverse-branch-certified-by-interval-brackets`

## Interval certificate

- All targets inside legacy branch range: `True`
- Legacy derivative upper bound on [0,2]: `-0.447676338367`
- Derivative bound strictly negative: `True`
- All roots bracketed/endpoints: `True`
- Max root interval width: `2.220e-16`
- Max midpoint residual: `3.548e-16`
- Strict targets decreasing: `True`
- Roots monotone nondecreasing: `True`
- Compression x(11)/11: `0.122016558479`

## Blocker refinement

- Previous monotone status: `OPEN_MONOTONE_INVERSE_BRANCH_EXISTENCE_WITNESS_NO_BRIDGE_THEOREM`
- Previous truncated-formal-inverse admissible: `False`
- Refined status: `global_z12_output_matching_map_certified; strict_transport_derivation remains open`

## Proof certificate

- `existence`: Every strict normalized target S(d), d=0..11, lies between L(0) and L(2), so IVT gives at least one root in [0,2].
- `uniqueness`: The analytic derivative upper bound for L on [0,2] is strictly negative, so each root is unique on the first legacy branch.
- `monotonicity`: The finite strict targets are strictly decreasing and L is strictly decreasing, hence the inverse roots are monotone nondecreasing; bisection intervals numerically certify the finite roots.
- `blocker_refinement`: This certifies a global finite Z12 output-matching map, but it does not derive that map from strict nadsoliton dynamics.

## Hard limits

- No unqualified identity K_legacy_ont == K_strict_gate is asserted.
- No strict-side derivation of the inverse branch or transport ODE is claimed.
- No beta_tors -> chi_11 theorem is asserted.
- No legacy physical-role transfer onto K_strict_gate is used without an explicit bridge theorem.
- No QW-2191 discharge is claimed.
- No ToE closure is claimed.
