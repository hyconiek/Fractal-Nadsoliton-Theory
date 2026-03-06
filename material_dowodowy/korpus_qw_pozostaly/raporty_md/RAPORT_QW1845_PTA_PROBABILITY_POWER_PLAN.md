# RAPORT QW-1845: PTA PROBABILITY POWER PLAN

- Data UTC: 2026-03-02T22:53:42.822590+00:00
- Test: H0 p<=0.900 vs H1 p>0.900
- alpha=0.050, power_target=0.80
- current n=14

## Grid (p_true assumptions)
- p_true=0.93: n=549, k_crit=506, alpha_act=0.0488, power=0.8038, add_vs_current=535
- p_true=0.95: n=179, k_crit=168, alpha_act=0.0486, power=0.8129, add_vs_current=165
- p_true=0.97: n=76, k_crit=73, alpha_act=0.0470, power=0.8057, add_vs_current=62
- p_true=0.99: n=46, k_crit=45, alpha_act=0.0480, power=0.9225, add_vs_current=32
- p_true=1.00: n=29, k_crit=29, alpha_act=0.0471, power=1.0000, add_vs_current=15

## Recommendation
- recommended_total_n_for_p_true_up_to_0.97: 549
- additional_needed_vs_current: 535

## Artifacts
- JSON: `report_qw1845_pta_probability_power_plan.json`
