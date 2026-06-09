# P2619/S1569 P2618 selector-source obligation lattice

Status: `P2619_SELECTOR_SOURCE_OBLIGATION_LATTICE_EXACT_C2_ENUMERATION_NO_SELECTOR_EXPORT_NO_LTOTAL_NO_QW2191_NO_TOE`

## Theorem

For a C2-odd strict phase sign, no selector map from orientation-invariant legacy scalar data to {+1,-1} is C2-equivariant. A C2-equivariant selector exists only after the input already contains a freely transforming orientation/sign torsor or an equivalent source premise.

## Proof

- Let C2={1,r} act on the strict sign set Sigma={+1,-1} by r·sigma=-sigma.
- If an input x is legacy-scalar/invariant, then r·x=x.
- Equivariance of a selector f requires f(r·x)=r·f(x), hence f(x)=-f(x).
- No element of Sigma satisfies sigma=-sigma, so no such selector exists from invariant input data.
- If the input is an orientation torsor X={x_+,x_-} with r·x_+=x_-, then f(x_+)=+1, f(x_-)=-1 and its negation are exactly the two equivariant maps; the missing object is therefore the source torsor, not a numerical fit.

## Computed C2 enumeration

- `1_orientation_invariant_legacy_scalar_classes`: candidates `2`, equivariant maps `0`, selector available `False`.
- `2_orientation_invariant_legacy_scalar_classes`: candidates `4`, equivariant maps `0`, selector available `False`.
- `3_orientation_invariant_legacy_scalar_classes`: candidates `8`, equivariant maps `0`, selector available `False`.
- `4_orientation_invariant_legacy_scalar_classes`: candidates `16`, equivariant maps `0`, selector available `False`.
- `1_orientation_odd_torsor_source_pairs`: candidates `4`, equivariant maps `2`, selector available `True`.
- `2_orientation_odd_torsor_source_pairs`: candidates `16`, equivariant maps `4`, selector available `True`.
- `3_orientation_odd_torsor_source_pairs`: candidates `64`, equivariant maps `8`, selector available `True`.

## Minimal source obligations

- `['orientation_odd_source']`
- `['symmetry_breaking_boundary']`
- `['spin_pin_orientation_source']`

## Bridge policy

- beta_tors: A beta_tors magnitude may remain legacy damping input, but as an orientation-invariant scalar it cannot by itself be the strict phase-sign source.
- axis-only data: Axis-only or quotient selector data can reduce continuous degeneracy but still leaves the Z2 sign unresolved.
- role transfer: Still blocked; after any real selector premise, rerun bridge validity and role-transfer audit separately.

## Scope guards

No strict selector source, no `beta_tors -> chi11` route reopening, no GF(2) bridge revalidation, no role-transfer revalidation, no role-bearing `L_total`, no `QW-2191` discharge, and no ToE closure are exported.

## Fingerprint

`c41aa8769c9bb0b819ed8701c16449817988b4bf055f5b17a3d3ab982cc851a6`
