# RAPORT QW-1789: COHORT CRITERIA SWEEP

- Data UTC: 2026-03-02T19:12:36.536246+00:00
- Precomputed pairs: 376
- Fixed protocol: frac=0.95, q_width=0.20
- Recommendation strength: STRONG
- Verdict: **COHORT_CRITERIA_SELECTION_SUPPORTED**

## Criteria Results
- K0_base: n_pairs=91 | full_reparam=0.186 | full_delta=0.644 | P(reparam>0)=1.000 | P(delta>0)=1.000 | std_rep=0.070 | std_delta=0.132 | score=0.791 | pass_basic=True
- K1_relax_n: n_pairs=91 | full_reparam=0.193 | full_delta=0.610 | P(reparam>0)=0.938 | P(delta>0)=1.000 | std_rep=0.070 | std_delta=0.116 | score=0.778 | pass_basic=True
- K2_relax_stab: n_pairs=94 | full_reparam=0.080 | full_delta=0.672 | P(reparam>0)=0.938 | P(delta>0)=1.000 | std_rep=0.075 | std_delta=0.121 | score=0.695 | pass_basic=True
- K3_balanced: n_pairs=94 | full_reparam=0.100 | full_delta=0.690 | P(reparam>0)=0.938 | P(delta>0)=1.000 | std_rep=0.049 | std_delta=0.126 | score=0.732 | pass_basic=True
- K4_wider: n_pairs=96 | full_reparam=0.019 | full_delta=0.738 | P(reparam>0)=0.750 | P(delta>0)=1.000 | std_rep=0.064 | std_delta=0.083 | score=0.634 | pass_basic=False
- K5_widest: n_pairs=96 | full_reparam=0.139 | full_delta=0.851 | P(reparam>0)=0.625 | P(delta>0)=1.000 | std_rep=0.044 | std_delta=0.113 | score=0.709 | pass_basic=False

## Recommended
- name=K0_base, n_match_min=120, stability_max=0.65, n_pairs=91, score=0.791

## Artifacts
- JSON: `report_qw1789_cohort_criteria_sweep.json`
