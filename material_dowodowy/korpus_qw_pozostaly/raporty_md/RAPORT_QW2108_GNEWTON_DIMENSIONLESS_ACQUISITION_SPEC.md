# RAPORT QW-2108: GNEWTON DIMENSIONLESS ACQUISITION SPEC

- Date UTC: 2026-03-04T12:35:49.486183+00:00
- Verdict: **GNEWTON_DIMENSIONLESS_ACQUISITION_SPEC_READY**

## Registry Target
- G_ref (SI): `6.674300000000e-11`
- tolerance_rel_pct: `5`

## Dimensionless Bridge Spec
- mu_ref_gev: `1`
- g_dimensionless_target: `6.708830750342e-39`
- acceptance_range: `[6.373389212825e-39, 7.044272287859e-39]`

## Strict Contract
- bridge_observable_origin must be `external_dimensionless_observable`
- value must NOT be backsolved from G_SI
- provenance_anchor_free must be true
- seeded_from_registry must be false
- g_si_input_optional must be null for strict path

## Artifact
- JSON: `report_qw2108_gnewton_dimensionless_acquisition_spec.json`
