# RAPORT QW-1971: TWO-TERM STRUCTURAL CONTROL DYNAMICS GATE

- Data UTC: 2026-03-03T08:42:03.760663+00:00
- Verdict: **TWO_TERM_STRUCTURAL_LOCKABLE**

## Baseline vs Best Two-Term
- baseline (QW-1970) bootstrap pass 5000: 79.94%
- best two-term bootstrap pass 5000: 100.00%
- delta: 20.06 pp

## Best Two-Term Candidate
- xi1: 0.000675
- xi2: -0.001000
- GW auc/adv/sep/gap: 0.8986/0.5316/0.003120/0.001357

## Local Robustness (xi1, xi2 neighborhood)
- (xi1,xi2)=(0.000375,-0.001000): det_pass=True, boot3000=100.00%
- (xi1,xi2)=(0.000675,-0.001300): det_pass=True, boot3000=100.00%
- (xi1,xi2)=(0.000675,-0.001000): det_pass=True, boot3000=100.00%
- (xi1,xi2)=(0.000675,-0.000700): det_pass=True, boot3000=100.00%
- (xi1,xi2)=(0.000975,-0.001000): det_pass=True, boot3000=100.00%

## Required Next Step
- FREEZE_TWO_TERM_STRUCTURE_AND_RUN_TRUE_EXTERNAL_CONFIRMATORY

## Artifacts
- JSON: `report_qw1971_two_term_structural_control_dynamics_gate.json`
