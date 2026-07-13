# P3143/S2093 V_lift weight/scale strict-source audit

Status: `P3143_VLIFT_WEIGHT_SCALE_SOURCE_AUDIT_BOUNDED_NO_GO`

## Constructed object
- `Omega_VL source-obligation matrix for V_lift weights/scale`
- Targets: `mu, w_theta, w_s, kappa`
- Classification: `strict_source_audit_matrix_for_axiom_branch_parameters`

## Finite theorem
`P3143_T1_vlift_weight_scale_source_obstruction_matrix`: The audited repo-backed scalar candidates provide positive receiver values in 16/20 target rows, but none supplies the full strict-source package for any V_lift coefficient: strict source law, unit-bearing normalization, global Z12 compatibility, noncircularity, and no closed-lane replay.  Therefore P3142's weights/scale remain axiom-branch parameters on current artifacts.

## Finite counts
- `candidates_tested`: `5`
- `targets_tested`: `4`
- `candidate_target_rows`: `20`
- `positive_value_rows`: `16`
- `accepted_strict_source_rows`: `0`

## Source-defect table
- `positive scalar value`: `16/20` rows pass
- `strict source law`: `0/20` rows pass
- `unit-bearing normalization`: `0/20` rows pass
- `global Z12 compatibility`: `12/20` rows pass
- `noncircularity`: `16/20` rows pass
- `no closed-lane replay`: `8/20` rows pass

## Candidate summary
- `C1_entropy_bit_cell` (P2689): value `1.0`, strict source `False`, unit-bearing `False`
- `C2_alpha_geo_shape_norm` (P2691): value `2.772588722239781`, strict source `False`, unit-bearing `False`
- `C3_beta_z_beta_positive_orbit` (P2692): value `1.0`, strict source `False`, unit-bearing `False`
- `C4_vlift_hessian_self_normalization` (P3142): value `2.0`, strict source `False`, unit-bearing `False`
- `C5_dhl_receiver_obligation_gap` (P3139): value `0.0`, strict source `False`, unit-bearing `False`

## Decision
P3143 performs the narrow source audit requested by P3142.  Existing positive scalar lanes can provide useful receiver magnitudes, but every row fails at least one strict obligation, and no row exports a strict unit-bearing source for mu, w_theta, w_s, or kappa.

## Recommendation
Do not keep recycling scalar magnitude sources for V_lift.  The next honest proof-grade step is exactly one new typed unit-measure object Upsilon_sel that couples the selector local chart to an action measure without importing A_origin/A_lambda; otherwise preserve the P3140-P3143 non-strict axiom-branch/no-strict-source boundary.
