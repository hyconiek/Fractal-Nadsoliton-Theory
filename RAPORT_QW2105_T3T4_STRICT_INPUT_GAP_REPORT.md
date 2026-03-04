# RAPORT QW-2105: T3T4 STRICT INPUT GAP REPORT

- Date UTC: 2026-03-04T11:59:42.198881+00:00
- Verdict: **T3T4_STRICT_INPUT_GAPS_PRESENT**

## H(z) Path
- strict_ready: `False`
- current n_nodes: `3`
- current z_span: `0.22999999999999998`
- current e_span: `0.48701913687450005`
- current cond([E,1]): `10.956177494448182`
- gaps:
  - add_nodes: at least 2 more valid H(z) nodes
  - extend_redshift_span: include at least one validated node at z >= 1.180
  - increase_E_span: add validated nodes at sufficiently separated z so E(z) span >= 1.0
  - improve_design_condition: reduce cond([E,1]) below 8 (current 10.956) via wider z coverage

## G_newton Path
- strict_ready: `False`
- bridge_observable_origin: `backsolved_from_g_si`
- gaps:
  - provide_direct_dimensionless_bridge: g_dimensionless_mu_ref must come from external observable, not backsolve
  - set_bridge_observable_origin_external_dimensionless
  - set_provenance_anchor_free_true
  - set_g_si_input_optional_null

## Meta
- QW-2104 verdict: `T3T4_STRICT_PREFLIGHT_GATE_PENDING` (0/6)

## Artifact
- JSON: `report_qw2105_t3t4_strict_input_gap_report.json`
- required_next_step: `COLLECT_STRICT_READY_HZ_AND_GNEWTON_EXTERNAL_INPUTS_THEN_RERUN_QW2099_QW2101_QW2102_QW2103_QW2090_QW2092_QW2104_QW2094`
