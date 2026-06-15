# P2761/S1711 kernel-moment to physical-coupling provenance obstruction

Status: `P2761_KERNEL_MOMENT_COUPLING_PROVENANCE_OBSTRUCTION_NO_CLOSURE`

## Moment/coupling audit
- moments={'M0': 1.2300294722058696, 'M1': 0.40861872472814803, 'M2': -14.5188745761074, 'M3': -285.6563371577367}
- derived_lagrangian_coefficients={'lambda_sm_eff': 0.33220238535858243, 'kappa_gr_eff': 11.91403866517398, 'epsilon_mix_eff': 18.13817242656519}
- row_count=3
- accepted_physical_coupling_count=0
- missing_global_provenance_atom_count=5
- stale_flags_quarantined=True

## Coupling rows
- lambda_sm_eff: accepted=False; expected_mass_dimension_4d=0; blocker=P1562 supplies a number, but current artifacts do not export the reference unit, sign convention, and variational normalization needed to interpret it as this physical coupling.
- kappa_gr_eff: accepted=False; expected_mass_dimension_4d=2; blocker=P1562 supplies a number, but current artifacts do not export the reference unit, sign convention, and variational normalization needed to interpret it as this physical coupling.
- epsilon_mix_eff: accepted=False; expected_mass_dimension_4d=1; blocker=P1562 supplies a number, but current artifacts do not export the reference unit, sign convention, and variational normalization needed to interpret it as this physical coupling.

## Missing global provenance atoms
- canonical physical length/reference cell mapping for moment powers
- action-density normalization and hbar/c unit convention
- field normalization Z_phi/Z_A/Z_psi and curvature normalization
- sign convention theorem for negative moments and positive couplings
- variational insertion theorem tying coefficients to nonproxy EOM residuals

## Decision
The moment-to-coupling numerics are present, but every target coupling lacks exported unit/reference, sign, and variational-normalization provenance; stale P1562 closure flags remain quarantined by later artifacts.

## Recommendation
Do not promote P1562 moment-derived numbers into physical couplings yet.  The next proof-grade move should target exactly one missing global provenance atom: construct a canonical physical length/reference-cell and action-density normalization theorem for the strict moment map; only after that rerun the sign and variational-normalization rows.  Without that theorem, preserve the P2697-P2761 no-closure certificate.
