# RAPORT QW-2105: T3T4 STRICT INPUT GAP REPORT

- Date UTC: 2026-03-04T13:28:22.202181+00:00
- Verdict: **T3T4_STRICT_INPUT_GAPS_CLOSED_READY_FOR_STRICT_RERUN**

## H(z) Path
- strict_ready: `True`
- current n_nodes: `5`
- current z_span: `1.15`
- current e_span: `4.276715605330501`
- current cond([E,1]): `5.7956053033256225`
- gaps:
  - none

## H(z) Design Guidance
- QW-2107 verdict: `HZ_STRICT_DESIGN_FOUND`
- suggested_added_z_top5:
  - [0.1, 0.9]
  - [0.1, 0.92]
  - [0.12, 0.92]
  - [0.1, 0.94]
  - [0.12, 0.94]

## G_newton Path
- strict_ready: `True`
- bridge_observable_origin: `external_dimensionless_observable`
- gaps:
  - none

## G_newton Design Guidance
- QW-2108 verdict: `GNEWTON_DIMENSIONLESS_ACQUISITION_SPEC_READY`
- mu_ref_gev: `1.0`
- g_dimensionless_target: `6.708830750342087e-39`
- g_dimensionless_acceptance_range: `[6.373389212824983e-39, 7.044272287859191e-39]`

## Meta
- QW-2104 verdict: `T3T4_STRICT_PREFLIGHT_GATE_PASS` (8/8)
- QW-2106 verdict: `STRICT_EXTERNAL_INPUT_INTAKE_GATE_PASS` (18/18)

## Artifact
- JSON: `report_qw2105_t3t4_strict_input_gap_report.json`
- required_next_step: `RERUN_T3T4_STRICT_CHAIN_AND_PROMOTE_CLOSURE`
