# RAPORT QW-2120: SCALAR SCALE VACUUM CLOSURE STRICT GATE

- Date UTC: 2026-03-04T18:38:59.526970+00:00
- Verdict: **SCALAR_SCALE_VACUUM_CLOSURE_STRICT_FAIL_INSUFFICIENT_FLOOR**
- pass_count: `7/8`

## Inputs
- v_higgs [GeV]: `243.017802492`
- m_h [GeV]: `126.363590676`
- required shift >= `0.681874763`

## Strict floor
- lambda_h = m_h^2/(2 v^2): `0.135187875`
- higgs_curvature_floor = m_h^2/v^2: `0.270375750`
- max g_i = max(m_i^2/v^2): `0.506775986`
- strict_floor_used_for_verdict: `0.506775986`

## Exploratory nonclosing channel (not used in verdict)
- micro_renorm_floor_candidate: `1.088153846`

## Artifact
- JSON: `report_qw2120_scalar_scale_vacuum_closure_strict_gate.json`
