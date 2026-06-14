# P2718/S1668 chiral-bispectrum explicit signed formula torsor-coupling audit

Status: `P2718_CHIRAL_BISPECTRUM_SIGNED_FORMULA_POSITIVE_BUT_NO_STRICT_TORSOR_SOURCE`

## Formula result
- checked_rows=24
- orientation_separating=True
- nonzero_signed_value_on_all_rows=True

## Acceptance
- accepted_as_strict_pseudoscalar_source=False
- missing=['translation_or_source_localizer_strictly_exported', 'phase_origin_reference_nonpremise', 'coupling_to_p2708_p2714_torsor_exported', 'qw2191_safe_nonpremise_selector_source']
The chiral-bispectrum imaginary marker is a real nonzero signed formula and separates orientation, but it is translation/source-blind without a non-premise phase-origin/source localizer and has no exported coupling theorem to the P2708/P2714 torsor.

## Decision
P2718 tests the concrete chiral-bispectrum formula demanded after P2717.  Its imaginary marker is nonzero on all 24 source/orientation rows and separates orientation with values +2 and -2.  However, the marker is still not a strict torsor-breaking source: it lacks a non-premise phase-origin/source localizer and no theorem couples it to the P2708/P2714 orientation torsor as a QW-2191-safe selector source.

## Next honest step
Do not promote the chiral-bispectrum marker by itself.  The next admissible move is a narrow phase-origin/source-localizer theorem audit for this exact formula, proving a non-premise origin reference and exported torsor coupling; if that cannot be supplied, preserve the P2697-P2718 no-new-live-frontier certificate.
