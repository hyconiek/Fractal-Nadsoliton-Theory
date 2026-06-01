# P2381 S1331: affine frontload source obligation theorem

Status: `OPEN_PROGRESS_AFFINE_FRONTLOAD_SOURCE_OBLIGATION_QUANTIFIED_SOURCE_STILL_OPEN`

## Result

P2381/S1331 turns the P2380 sufficient affine front-loaded profile into a quantified source-obligation statement.
Any rectangle-uniform affine source must derive `lambda` above the worst-corner threshold, not merely name a normalized profile.

## Certificate

- Worst-corner lambda threshold: `0.7916644842269429`.
- Necessary obligations at that threshold: `{'lambda': 0.7916644842269429, 'rho_s0': 1.791664484226943, 'rho_s1': 0.2083355157730571, 'endpoint_contrast_rho0_over_rho1': 8.599899434231027, 'early_half_mass_int_0_to_1_2': 0.6979161210567357, 'late_half_mass_int_1_2_to_1': 0.3020838789432643, 'transport_time_barycenter_int_s_rho': 0.3680559192955095, 'barycenter_shift_from_uniform': 0.13194408070449049, 'l1_distance_from_uniform': 0.39583224211347146, 'l2_squared_distance_from_uniform': 0.20891088519543718, 'profile_variance_int_rho_minus_1_sq': 0.20891088519543718}`.
- Lambda=0.8 obligations: `{'lambda': 0.8, 'rho_s0': 1.8, 'rho_s1': 0.19999999999999996, 'endpoint_contrast_rho0_over_rho1': 9.000000000000002, 'early_half_mass_int_0_to_1_2': 0.7, 'late_half_mass_int_1_2_to_1': 0.3, 'transport_time_barycenter_int_s_rho': 0.3666666666666667, 'barycenter_shift_from_uniform': 0.13333333333333333, 'l1_distance_from_uniform': 0.4, 'l2_squared_distance_from_uniform': 0.21333333333333337, 'profile_variance_int_rho_minus_1_sq': 0.21333333333333337}`.
- Below-threshold negative control selects: `{'3,3': 24}`.
- Lambda=0.8 corner replay selects: `{'0,4': 12}`.

## Hard limits

- This is a necessary source-obligation theorem for the affine family, not a derivation of that source from strict dynamics.
- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.
