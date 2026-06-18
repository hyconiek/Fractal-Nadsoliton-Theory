# P2851/S1801 amplitude-normalization kernel bridge-atom audit

Status: `P2851_AMPLITUDE_NORMALIZATION_KERNEL_BRIDGE_ATOM_AUDIT_NO_CLOSURE`

## Ratio stats
- min_ratio=-1237.9850503876855
- max_ratio=793.8604020936507
- max_minus_min=2031.8454524813362
- sign_changes_present=True
- constant_ratio=False

## Least-squares constant-amplitude fit
- best_fit_amplitude=1.721420697189017
- sse=42.34890837776295
- max_abs_residual=2.7574988103791687

## Two-point constant-amplitude matrix
A target-independent amplitude-only bridge A*K_strict(d)=K_legacy(d) requires identical pointwise ratios K_legacy(d)/K_strict(d) at every nonzero strict point; two distinct ratios refute an amplitude-only bridge.
- pair=[0, 1]: ratio_a=2.4331873129223913; ratio_b=1.5117350775413128; same=False
- pair=[0, 2]: ratio_a=2.4331873129223913; ratio_b=-7.077103640210554; same=False
- pair=[1, 2]: ratio_a=1.5117350775413128; ratio_b=-7.077103640210554; same=False
- pair=[2, 4]: ratio_a=-7.077103640210554; ratio_b=-49.09253317078477; same=False
- pair=[4, 8]: ratio_a=-49.09253317078477; ratio_b=-1237.9850503876855; same=False
- pair=[8, 12]: ratio_a=-1237.9850503876855; ratio_b=259.635881756419; same=False

## Premise matrix
- accepted_as_amplitude_normalization_bridge_atom=False
- missing_premises=['target_independent_constant_amplitude_map', 'zero_residual_after_constant_amplitude_fit', 'phase_frequency_compatibility', 'damping_compression_compatibility', 'alpha_geo_source_law_safe_for_strict_kernel', 'role_transfer_theorem_available']

## Boundary
P2851 tests exactly one bridge atom: amplitude/normalization passage from legacy alpha_geo=4 ln 2 to the unit-amplitude strict gate kernel.  Pointwise legacy/strict amplitude ratios are not constant, two-point ratio checks reject an amplitude-only bridge, and the best constant-amplitude least-squares fit leaves a nonzero residual.  Thus alpha_geo cannot be silently absorbed into K_strict_gate without the missing phase/frequency and damping/compression bridge atoms plus a role-safe alpha_geo source law.

## Recommendation
Do not claim amplitude passage, full bridge, role transfer, L_total, EOM, Hamiltonian, or ToE closure from alpha_geo rescaling.  The next admissible move is a combined bridge-obligation reconciliation matrix over the now-tested atoms (damping/compression, amplitude, EML syntax) to identify whether any single remaining typed source atom is still live; otherwise preserve no-new-live-frontier.
