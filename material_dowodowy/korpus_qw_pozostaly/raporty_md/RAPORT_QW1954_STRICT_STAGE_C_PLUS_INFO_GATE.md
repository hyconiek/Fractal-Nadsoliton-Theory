# RAPORT QW-1954: STRICT STAGE-C+INFO GATE

- Data UTC: 2026-03-03T07:17:58.240996+00:00
- Verdict: **STRICT_STAGE_C_PLUS_INFO_GATE_FAIL**
- Readiness: **TOE_STAGE_C_PLUS_INFO_BLOCKED**

## Flags
- stage_b_closed: True
- triad_strict_pass_q1946: False
- info_dedeg_pass_q1952: False
- info_two_state_pass_q1953: False

## Top Blockers
- triad::mass_mean_rel_pct: 4.2162
- triad::ckm_mean_rel_pct: 1.5996
- info1952::acc_gain_vs_control: 1.2778
- info1953::closed_info_gain_vs_control: 1.1832
- info1952::info_gain_vs_control: 1.1780
- triad::q_assignment_joint_feasibility: 1.0000
- info1953::channel_complementarity: 0.9901
- info1952::channel_complementarity: 0.9900

## Required Next Step
- REWORK_INFORMATION_CHANNEL_DEGENERACY_AND_MASS_LINK_UNDER_FROZEN_KERNEL

## Artifacts
- JSON: `report_qw1954_strict_stage_c_plus_info_gate.json`
