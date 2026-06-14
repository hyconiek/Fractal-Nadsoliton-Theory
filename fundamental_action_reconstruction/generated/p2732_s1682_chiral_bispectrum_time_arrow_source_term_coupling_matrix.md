# P2732/S1682 chiral-bispectrum time-arrow source-term coupling matrix

Status: `P2732_CHIRAL_BISPECTRUM_TIME_ARROW_SOURCE_TERM_CONDITIONAL_NO_GO`

## Coupling matrix
- coupling_row_count=4
- field_count_per_row=4096
- all_rows_have_unique_constant_tau_ground_state=True
- selected_tau_sign_histogram={-1: 2, 1: 2}
- all_lambda_flip_pairs_reverse_tau=True
- all_orientation_flip_pairs_reverse_tau=True
The chiral-bispectrum source term is strong enough to choose a tau ground state only after lambda and orientation/polarity are fixed.  The finite matrix is exactly paired: flipping lambda or the orientation torsor reverses the selected tau sign, and current artifacts do not export a non-premise law selecting one lambda/P2721 polarity.

## Recommendation
Do not replay direct Im(B_{1,5}) tau-coupling as closure.  A next admissible proof-grade move must export a strict law fixing the coupling sign lambda or the P2721 polarity from inside the theory; otherwise preserve the P2697-P2732 no-new-live-frontier certificate.
