# RAPORT QW-1970: STRUCTURAL GW CONTROL TERM GATE

- Data UTC: 2026-03-03T08:37:41.159985+00:00
- Verdict: **STRUCTURAL_CONTROL_TERM_INSUFFICIENT**

## Baseline vs Structural Term
- baseline bootstrap pass (5000): 71.94%
- best structural-term bootstrap pass (5000): 79.94%
- delta: 8.00 pp

## Best Structural Candidate
- xi: 0.000400
- GW auc/adv/sep/gap: 0.8298/0.3221/0.002839/0.002206

## Xi Local Robustness
- xi=-0.000600: det_pass=False, boot3000=27.73%
- xi=0.000400: det_pass=True, boot3000=79.47%
- xi=0.001400: det_pass=False, boot3000=35.43%

## Required Next Step
- NEED_DEEPER_SHARED_DYNAMICS_FOR_CONTROL_CHANNEL

## Artifacts
- JSON: `report_qw1970_structural_gw_control_term_gate.json`
