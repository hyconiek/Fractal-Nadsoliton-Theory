# P2571/S1521 APD Sobolev measure/boundary dependence audit

Status: `P2571_APD_SOBOLEV_MEASURE_BOUNDARY_DEPENDENCE_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE`

## Result

- Frontier atom under attack: `strict_dynamical_source_for_A_P_D`.
- Candidate principle: `fix Sobolev derivative order k=2, then vary positive measure weights and endpoint slope penalties on the P2569 interpolation family`.
- Fixed derivative order: `2`.
- Variants audited: `7`.
- Distinct minimizers after rounding: `7`.
- All variants preserve APD nodes: `True`.
- Measure or boundary changes selector: `True`.
- Finite APD values select measure/boundary class: `False`.

## Interpretation

Even after fixing the Sobolev derivative order to `k=2`, changing the positive measure weight or endpoint slope-penalty class changes the minimizing `lambda` on the same APD interpolation family.  The finite APD nodes remain exact in every audited variant, so measure and boundary data are additional strict source obligations.

## Recommended next honest step

Do not treat k=2 hydrodynamic/Sobolev intuition as an APD source theorem. The next honest step is to derive the measure and boundary class of the APD action from strict nadsoliton dynamics; without that derivation, even a fixed second-order roughness law is a family of conditional selectors rather than strict A/P/D dynamics.

## Negative controls

No APD measure source, boundary-class source, weighted roughness selector source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.

## Fingerprint

`f3d0999ada41b472dab05ca0010a2f4ad87a9e1f7df4a9957c6337c514a1f3df`
