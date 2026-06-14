# P2743/S1693 affine-frame transition character no-go

Status: `P2743_AFFINE_FRAME_TRANSITION_CHARACTER_NO_GO`

## Finite transition-character audit
- frame_count=48
- transition_count=2304
- unit_counts={1: 576, 5: 576, 7: 576, 11: 576}
- transition_unit_orbit_count_under_simultaneous_affine_action=4
- inversion_odd_character_count=2
- characters_with_nonzero_global_signed_sum=0
- characters_with_balanced_positive_negative_transitions=2
- all_unit_11_sign_flip_witnesses_pass=True

## Theorem statement
The two inversion-odd U(12) characters are valid pointwise signs on affine-frame transition units, but the 48x48 transition ensemble contains each unit exactly 576 times.  Hence every inversion-odd character has 1152 positive and 1152 negative transitions and global signed sum zero; multiplication by unit 11 pairs character signs.  Without a strict source selecting one transition unit or one character polarity, this frame-character observable does not export a non-premise orbit-safe signed value or fix P2721/lambda.

## Recommendation
Do not promote affine-frame transition characters: P2743 shows the two inversion-odd U(12) character signs are pointwise real but globally balanced across all 2304 frame transitions, and current artifacts export no strict transition-unit source or P2721 polarity coupling.  The next proof-grade move must supply a strict source law selecting a nonzero transition-unit/character polarity with P2721 coupling, or pivot outside finite Z12 sign-character/frame observables; otherwise preserve the P2697-P2743 no-new-live-frontier certificate.
