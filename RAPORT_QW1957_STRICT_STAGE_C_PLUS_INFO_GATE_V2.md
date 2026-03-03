# RAPORT QW-1957: STRICT STAGE-C+INFO GATE V2

- Data UTC: 2026-03-03T07:29:26.078407+00:00
- Verdict: **STRICT_STAGE_C_PLUS_INFO_GATE_V2_FAIL**
- Readiness: **TOE_STAGE_C_PLUS_INFO_BLOCKED_V2**

## Flags
- stage_b_closed: True
- triad_strict_pass_q1946: False
- minimal_repair_pass_q1955: False
- two_state_repair_pass_q1956: False

## Top Blockers
- triad::mass_mean_rel_pct: 4.2162
- triad::ckm_mean_rel_pct: 1.5996
- info1956::closed_info_gain_vs_control: 1.1778
- info1955::info_gain_vs_control: 1.1361
- triad::q_assignment_joint_feasibility: 1.0000
- info1955::acc_gain_vs_control: 0.8611
- info1956::channel_complementarity: 0.4735
- info1955::channel_complementarity: 0.4648

## Required Next Step
- REWORK_MINIMAL_REPAIR_OPERATOR_AND_MASS_FLAVOR_LINK

## Artifacts
- JSON: `report_qw1957_strict_stage_c_plus_info_gate_v2.json`
