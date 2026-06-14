# P2726/S1676 chiral-bispectrum affine orientation-flow transition matrix

Status: `P2726_AFFINE_ORIENTATION_FLOW_NONZERO_BUT_UNSOURCED_NO_STRICT_DYNAMIC_CHIRAL_SOURCE`

## Finite transition result
- checked_transition_rows=576
- orientation_preserving_delta_values=[0.0]
- orientation_flipping_delta_values=[-4.0, 4.0]
- orientation_flipping_delta_histogram={-4.0: 144, 4.0: 144}
The only nonzero signed jumps in the complete affine transition matrix occur when the transition flips the orientation torsor.  Since current artifacts do not export a non-premise law selecting such a flip direction or its P2721 coupling polarity, the nonzero +/-4 jumps are conditional evidence, not a strict source.

## Acceptance
- accepted_as_strict_dynamic_chiral_source=False
- missing=['strict_orientation_flip_source_law_exported', 'nonpremise_orientation_flip_direction_selected', 'p2721_coupling_polarity_selected']

## Recommendation
The next admissible step is no longer to search translation-flow variants.  It must either export a strict orientation-flip/chiral-flow source law that selects one branch and one P2721 polarity, then rerun this matrix as an acceptance test, or preserve the P2697-P2726 no-new-live-frontier certificate.
