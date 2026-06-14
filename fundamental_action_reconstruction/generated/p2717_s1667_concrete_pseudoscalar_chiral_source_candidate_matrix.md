# P2717/S1667 concrete pseudoscalar/chiral source candidate matrix

Status: `P2717_CONCRETE_PSEUDOSCALAR_CHIRAL_SOURCE_MATRIX_NO_ACCEPTED_SOURCE`

## Candidate matrix
- `levi_civita_volume_orientation_density`: accepted=False, missing=['strict_artifact_exports_source_law', 'nonzero_signed_value_exported', 'sign_not_orientation_convention', 'coupling_to_p2708_p2714_orientation_torsor_exported']. A Levi-Civita density is sign-odd, but its sign is an orientation convention unless a strict law chooses the orientation.
- `pontryagin_or_chiral_anomaly_density`: accepted=False, missing=['strict_artifact_exports_source_law', 'nonzero_signed_value_exported', 'coupling_to_p2708_p2714_orientation_torsor_exported']. A Pontryagin/anomaly-type density would be the right parity class, but no strict artifact exports such a nonzero signed density or its coupling to the Z12 torsor.
- `eta_or_spectral_asymmetry_index`: accepted=False, missing=['strict_artifact_exports_source_law', 'nonzero_signed_value_exported', 'coupling_to_p2708_p2714_orientation_torsor_exported']. A spectral asymmetry could source chirality, but no current strict artifact exports an eta/index sign coupled to the boundary-cocycle orientation torsor.
- `oriented_z12_cycle_cup_product`: accepted=False, missing=['nonzero_signed_value_exported', 'sign_not_orientation_convention']. The oriented Z12 cycle is the existing P2708/P2714 torsor itself; it supplies the two signs but not a non-premise rule selecting one.

## Decision
P2717 audits concrete pseudoscalar/chiral source classes after P2716.  Levi-Civita orientation density, Pontryagin/anomaly density, eta/spectral asymmetry, and the oriented Z12 cup-product candidate all have the right sign intuition or parity, but none satisfies all strict criteria: exported source law, nonzero signed value, nonconventional sign, and coupling to the P2708/P2714 orientation torsor.  Therefore no candidate fixes lambda or discharges QW-2191.

## Next honest step
Do not keep enumerating generic pseudoscalar names.  A next admissible move must provide one explicit formula/artifact computing a nonzero signed pseudoscalar/chiral value and its coupling to the orientation torsor, or pivot to a different genuinely new typed object outside the closed lanes.  Otherwise preserve the P2697-P2717 no-new-live-frontier certificate.
