import json
import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2025_s975_strict_cutkosky_same_scheme_cohomology_amplitude_bridge_seed.py"
OUT = ROOT / "generated" / "p2025_s975_strict_cutkosky_same_scheme_cohomology_amplitude_bridge_seed.json"
OUT_CSV = ROOT / "generated" / "p2025_s975_strict_cutkosky_same_scheme_cohomology_amplitude_bridge_seed_per_channel_power_aware_verdicts.csv"
OUT_QUALITY_CSV = ROOT / "generated" / "p2025_s975_strict_cutkosky_same_scheme_cohomology_amplitude_bridge_seed_per_channel_wilcoxon_quality.csv"


def test_p2025_exports_same_scheme_bridge_seed_without_false_closure():
    subprocess.run([sys.executable, str(SCRIPT)], check=True)
    data = json.loads(OUT.read_text(encoding="utf-8"))
    assert data["schema_version"] == "p2025_s975_v202"
    assert data["status"] == "OPEN_OBSTRUCTION_WITH_TRACE"
    assert all(data["gatekeeper_checks"].values())
    assert len(data["toe_closure_gaps_7tasks"]) == 7
    assert len(data["task_numeric_evidence_7"]) == 7
    for row in data["task_numeric_evidence_7"]:
        assert row["honest_verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE"
        assert len(row["method_stack"]) >= 2
        assert 0.0 <= row["local_readiness_score_0_1"] <= 1.0
    task2 = next(r for r in data["task_numeric_evidence_7"] if r["task_id"] == 2)
    assert task2["metrics"]["q95_margin_before_abs"] >= 0.0
    assert task2["metrics"]["q95_margin_after_abs"] >= 0.0
    assert task2["metrics"]["q95_margin_improvement_abs"] >= 0.0
    assert 0.0 <= task2["metrics"]["q95_progress_score_0_1"] <= 1.0
    q95_proj = data["task2_strict_unitarity_witness"]["q95_blocker_direct_relief_projection"]
    assert q95_proj["scope"] == "STRICT_TASK2_Q95_BLOCKER_DIRECT_RELIEF_PROJECTION"
    assert q95_proj["theorem_target"] == "ONE_MORE_LOCAL_STEP_CAN_CLOSE_Q95_MARGIN"
    assert "projected_margin_after_next_step_le_zero" in q95_proj["pass_fail_criteria"]
    assert q95_proj["verdict"] in {"PREDICTED_CLOSE_IN_NEXT_STEP", "OPEN_OBSTRUCTION_WITH_TRACE"}
    if q95_proj["projected_margin_after_next_step"] is not None:
        assert q95_proj["estimated_steps_to_closure"] >= 0.0
        if q95_proj["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
            assert "projected_margin_after_next_step=" in q95_proj["fail_trace"]
            assert " > 0" in q95_proj["fail_trace"]
    ls_exec = data["task2_strict_unitarity_witness"]["q95_blocker_local_line_search_execution"]
    assert ls_exec["scope"] == "STRICT_TASK2_Q95_BLOCKER_LOCAL_LINE_SEARCH_EXECUTION"
    assert ls_exec["theorem_target"] == "EXISTS_LOCAL_S_STEP_WITH_Q95_GAP_NOT_WORSE_THAN_BASELINE"
    assert isinstance(ls_exec["computed_rows"], list)
    assert ls_exec["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    if ls_exec["computed_rows"]:
        best_gap = ls_exec["aggregate_metrics"]["best_gap_abs_quad"]
        baseline_gap = ls_exec["aggregate_metrics"]["baseline_gap_abs_quad"]
        assert best_gap is not None and baseline_gap is not None
        assert best_gap <= baseline_gap + 1e-18
    if ls_exec["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "=" in ls_exec["fail_trace"]
        assert ">" in ls_exec["fail_trace"]
    opt_exec = data["task2_strict_unitarity_witness"]["q95_blocker_continuous_optimization_execution"]
    assert opt_exec["scope"] == "STRICT_TASK2_Q95_BLOCKER_CONTINUOUS_OPTIMIZATION_EXECUTION"
    assert opt_exec["theorem_target"] == "EXISTS_S_IN_RANGE_WITH_Q95_GAP_BELOW_THRESHOLD_UNDER_CROSSCHECK"
    assert opt_exec["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    if opt_exec["computed_rows"]:
        row = opt_exec["computed_rows"][0]
        assert 0.05 <= row["s_opt"] <= 3.5
        assert row["gap_abs_quad_opt"] >= 0.0
        assert row["cross_integrator_gap_abs_n400_vs_n800"] >= 0.0
    if opt_exec["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "=" in opt_exec["fail_trace"]
        assert ">" in opt_exec["fail_trace"]
    refined_exec = data["task2_strict_unitarity_witness"]["q95_blocker_refined_window_execution"]
    assert refined_exec["scope"] == "STRICT_TASK2_Q95_BLOCKER_REFINED_WINDOW_EXECUTION"
    assert refined_exec["theorem_target"] == "EXISTS_S_IN_LOCAL_WINDOW_WITH_Q95_GAP_AND_CROSSCHECK_BELOW_THRESHOLDS"
    assert refined_exec["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    if refined_exec["computed_rows"]:
        br = refined_exec["computed_rows"][0]
        assert br["gap_abs_quad"] >= 0.0
        assert br["cross_integrator_gap_abs_n600_vs_n1200"] >= 0.0
        assert br["gap_uncertainty_std_across_3_integrators"] >= 0.0
    if refined_exec["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "=" in refined_exec["fail_trace"]
        assert ">" in refined_exec["fail_trace"]
    sub_attempt = data["task2_strict_unitarity_witness"]["q95_blocker_single_row_substitution_attempt"]
    assert sub_attempt["scope"] == "STRICT_TASK2_Q95_BLOCKER_SINGLE_ROW_SUBSTITUTION_ATTEMPT"
    assert sub_attempt["theorem_target"] == "SINGLE_DOMINANT_ROW_REPLACEMENT_CAN_CLOSE_GLOBAL_Q95"
    assert sub_attempt["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    if sub_attempt["computed_rows"]:
        row = sub_attempt["computed_rows"][0]
        assert row["q95_gap_abs_baseline"] >= 0.0
        assert row["q95_gap_abs_after_single_substitution"] >= 0.0
    if sub_attempt["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "=" in sub_attempt["fail_trace"]
        assert ">" in sub_attempt["fail_trace"]
    sub2_attempt = data["task2_strict_unitarity_witness"]["q95_blocker_two_row_substitution_attempt"]
    assert sub2_attempt["scope"] == "STRICT_TASK2_Q95_BLOCKER_TWO_ROW_SUBSTITUTION_ATTEMPT"
    assert sub2_attempt["theorem_target"] == "TOP2_DOMINANT_ROW_REPLACEMENT_CAN_CLOSE_GLOBAL_Q95"
    assert sub2_attempt["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    if sub2_attempt["computed_rows"]:
        assert 1 <= len(sub2_attempt["computed_rows"]) <= 2
        for rr in sub2_attempt["computed_rows"]:
            assert rr["baseline_gap_abs"] >= 0.0
            assert rr["replacement_gap_abs"] >= 0.0
            assert rr["replacement_cross_integrator_gap_abs_n600_vs_n1200"] >= 0.0
        am = sub2_attempt["aggregate_metrics"]
        assert am["q95_gap_abs_baseline"] >= 0.0
        assert am["q95_gap_abs_after_two_row_substitution"] >= 0.0
        assert 0.0 <= am["improvement_fraction_of_baseline"] <= 1.0
    if sub2_attempt["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "=" in sub2_attempt["fail_trace"]
        assert ">" in sub2_attempt["fail_trace"]
    mk_scan = data["task2_strict_unitarity_witness"]["q95_blocker_min_k_substitution_scan"]
    assert mk_scan["scope"] == "STRICT_TASK2_Q95_BLOCKER_MIN_K_SUBSTITUTION_SCAN"
    assert mk_scan["theorem_target"] == "FIND_MIN_K_DOMINANT_ROW_SUBSTITUTIONS_TO_CLOSE_GLOBAL_Q95"
    assert mk_scan["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    if mk_scan["computed_rows"]:
        assert mk_scan["domain"]["k_scan_range_inclusive"][0] == 1
        for rr in mk_scan["computed_rows"]:
            assert rr["k"] >= 1
            assert rr["q95_gap_abs_after_k_row_substitution"] >= 0.0
            assert rr["improvement_abs"] >= 0.0
            assert isinstance(rr["closes_threshold"], bool)
            assert len(rr["replacement_pairs"]) == rr["k"]
        am = mk_scan["aggregate_metrics"]
        assert am["q95_gap_abs_baseline"] >= 0.0
        assert am["best_q95_after_k_substitution"] >= 0.0
        assert am["best_k_by_q95"] >= 1
    if mk_scan["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "=" in mk_scan["fail_trace"]
        assert ">" in mk_scan["fail_trace"]
    joint_exec = data["task2_strict_unitarity_witness"]["q95_blocker_joint_topk_continuous_execution"]
    assert joint_exec["scope"] == "STRICT_TASK2_Q95_BLOCKER_JOINT_TOPK_CONTINUOUS_EXECUTION"
    assert joint_exec["theorem_target"] == "JOINT_TOPK_LOCAL_CONTINUOUS_OPTIMIZATION_CAN_CLOSE_GLOBAL_Q95"
    assert joint_exec["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    if joint_exec["computed_rows"]:
        assert 1 <= len(joint_exec["computed_rows"]) <= 3
        for rr in joint_exec["computed_rows"]:
            assert rr["gap_abs_baseline"] >= 0.0
            assert rr["gap_abs_optimized"] >= 0.0
            assert rr["improvement_abs"] >= 0.0
            assert rr["s_window"][0] <= rr["s_optimized"] <= rr["s_window"][1]
            assert isinstance(rr["optimizer_success"], bool)
            assert rr["optimizer_nfev"] >= 1
        am = joint_exec["aggregate_metrics"]
        assert am["q95_gap_abs_baseline"] >= 0.0
        assert am["q95_gap_abs_after_joint_topk_optimization"] >= 0.0
        assert 0.0 <= am["improvement_fraction_of_baseline"] <= 1.0
    if joint_exec["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "=" in joint_exec["fail_trace"]
        assert ">" in joint_exec["fail_trace"]
    joint4_exec = data["task2_strict_unitarity_witness"]["q95_blocker_joint_top4_vector_optimization_execution"]
    assert joint4_exec["scope"] == "STRICT_TASK2_Q95_BLOCKER_JOINT_TOP4_VECTOR_OPTIMIZATION_EXECUTION"
    assert joint4_exec["theorem_target"] == "JOINT_TOP4_VECTOR_OPTIMIZATION_CAN_CLOSE_GLOBAL_Q95"
    assert joint4_exec["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    if joint4_exec["computed_rows"]:
        assert 1 <= len(joint4_exec["computed_rows"]) <= 4
        for rr in joint4_exec["computed_rows"]:
            assert rr["gap_abs_baseline"] >= 0.0
            assert rr["gap_abs_optimized"] >= 0.0
            assert rr["improvement_abs"] >= 0.0
            assert rr["s_window"][0] <= rr["s_optimized"] <= rr["s_window"][1]
        am = joint4_exec["aggregate_metrics"]
        assert am["q95_gap_abs_baseline"] >= 0.0
        assert am["q95_gap_abs_after_joint_top4_vector_optimization"] >= 0.0
        assert 0.0 <= am["improvement_fraction_of_baseline"] <= 1.0
        assert isinstance(am["optimizer_success"], bool)
        assert am["optimizer_nfev"] >= 0
    if joint4_exec["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "=" in joint4_exec["fail_trace"]
        assert ">" in joint4_exec["fail_trace"]
    joint_profile = data["task2_strict_unitarity_witness"]["q95_blocker_adaptive_joint_topk_profile_execution"]
    assert joint_profile["scope"] == "STRICT_TASK2_Q95_BLOCKER_ADAPTIVE_JOINT_TOPK_PROFILE_EXECUTION"
    assert joint_profile["theorem_target"] == "EXISTS_K_IN_{2,4,6}_JOINT_VECTOR_OPTIMIZATION_THAT_CLOSES_GLOBAL_Q95"
    assert joint_profile["domain"]["k_grid"] == [2, 4, 6]
    assert joint_profile["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    if joint_profile["computed_rows"]:
        for rr in joint_profile["computed_rows"]:
            assert rr["k"] in {2, 4, 6}
            assert rr["q95_gap_abs_after_joint_vector_optimization"] >= 0.0
            assert rr["improvement_abs"] >= 0.0
            assert 0.0 <= rr["improvement_fraction_of_baseline"] <= 1.0
            assert isinstance(rr["closes_threshold"], bool)
            assert isinstance(rr["optimizer_success"], bool)
            assert rr["optimizer_nfev"] >= 0
        am = joint_profile["aggregate_metrics"]
        assert am["q95_gap_abs_baseline"] >= 0.0
        assert am["best_q95_after_joint_vector_optimization"] >= 0.0
        assert am["best_k_by_q95"] in {2, 4, 6}
        assert (am["first_k_that_closes_threshold"] in {2, 4, 6}) or (am["first_k_that_closes_threshold"] is None)
    if joint_profile["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "=" in joint_profile["fail_trace"]
        assert ">" in joint_profile["fail_trace"]
    bestk_cert = data["task2_strict_unitarity_witness"]["q95_blocker_joint_bestk_exact_recompute_certificate"]
    assert bestk_cert["scope"] == "STRICT_TASK2_Q95_BLOCKER_JOINT_BESTK_EXACT_RECOMPUTE_CERTIFICATE"
    assert bestk_cert["theorem_target"] == "BEST_K_JOINT_OPTIMIZATION_RECOMPUTE_CLOSES_Q95_UNDER_INDEPENDENT_INTEGRATOR_REPLAY"
    assert bestk_cert["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    if bestk_cert["computed_rows"]:
        for rr in bestk_cert["computed_rows"]:
            assert rr["gap_abs_quad_optimized"] >= 0.0
            assert rr["gap_abs_fixed_quad_n1200_optimized"] >= 0.0
            assert rr["cross_integrator_gap_abs"] >= 0.0
        am = bestk_cert["aggregate_metrics"]
        assert am["q95_gap_abs_after_joint_bestk_quad"] >= 0.0
        assert am["q95_gap_abs_after_joint_bestk_fixed_quad_n1200"] >= 0.0
        assert am["q95_cross_integrator_gap_abs"] >= 0.0
    if bestk_cert["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "=" in bestk_cert["fail_trace"]
        assert ">" in bestk_cert["fail_trace"]
    knee_cert = data["task2_strict_unitarity_witness"]["q95_blocker_adaptive_profile_knee_certificate"]
    assert knee_cert["scope"] == "STRICT_TASK2_Q95_BLOCKER_ADAPTIVE_PROFILE_KNEE_CERTIFICATE"
    assert knee_cert["theorem_target"] == "ADAPTIVE_PROFILE_EXHIBITS_NONZERO_MARGINAL_Q95_GAIN_AND_IDENTIFIES_BEST_K"
    assert knee_cert["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    if knee_cert["computed_rows"]:
        for rr in knee_cert["computed_rows"]:
            assert rr["k_from"] in {2, 4, 6}
            assert rr["k_to"] in {2, 4, 6}
            assert rr["k_from"] < rr["k_to"]
            assert rr["delta_k"] > 0
            assert rr["marginal_q95_gain_per_k"] >= 0.0
        assert knee_cert["aggregate_metrics"]["best_q95_after_joint_vector_optimization"] >= 0.0
    if knee_cert["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "=" in knee_cert["fail_trace"]
        assert ">" in knee_cert["fail_trace"]
    blocker_rows = data["task2_strict_unitarity_witness"]["q95_blocker_choice_panel"]["rows"]
    assert any(r["criterion"] == "q95_refined_window_gap_abs" for r in blocker_rows)
    assert data["verdict"] == data["task2_strict_unitarity_witness"]["verdict"]
    assert data["fail_trace"] == data["task2_strict_unitarity_witness"]["fail_trace"]
    assert "task_priority_decision_panel" in data
    tpp = data["task_priority_decision_panel"]
    assert len(tpp["rows"]) == 7
    assert tpp["recommended_next_task_id"] in {1, 2, 3, 4, 5, 6, 7}
    assert tpp["recommended_lane"] in {"kernel_split_robust_discm_integration", "non_selector_fallback_due_to_qw2191_guardrail"}
    assert tpp["symbolic_normalization_certificate"]["exactly_one"] is True
    assert tpp["normalized_weight_entropy_nats"] >= 0.0
    assert tpp["score_dispersion_l2"] >= 0.0
    assert tpp["score_mad"] >= 0.0
    assert tpp["score_std"] >= 0.0
    assert tpp["score_cv"] >= 0.0
    assert tpp["score_mean_ci95_t_interval"]["lower"] <= tpp["score_mean_ci95_t_interval"]["upper"]
    assert -1.0 <= tpp["score_spearman_rank_stability"] <= 1.0
    assert tpp["score_covariance_scalar"] >= 0.0
    assert tpp["score_pca_effective_rank"] >= 0.0
    assert len(tpp["score_pca_variance_ratio"]) == 7
    assert tpp["score_centering_symbolic_certificate"]["exactly_zero"] is True
    assert tpp["bootstrap_readiness_summary"]["bootstrap_size"] == 512
    assert len(tpp["bootstrap_readiness_summary"]["mean_q05_q50_q95"]) == 3
    assert len(tpp["bootstrap_readiness_summary"]["std_q05_q50_q95"]) == 3
    assert len(tpp["bootstrap_readiness_summary"]["top_index_frequency_over_resamples"]) == 7
    assert abs(sum(tpp["bootstrap_readiness_summary"]["top_index_frequency_over_resamples"]) - 1.0) < 1e-12
    assert tpp["robust_spread"]["iqr"] >= 0.0
    assert tpp["robust_spread"]["mad_scaled"] >= 0.0
    assert "branch_integrator_replay_task7_panel" in tpp
    rpt = tpp["branch_integrator_replay_task7_panel"]
    assert len(rpt["seeds"]) == 4
    assert len(rpt["stress_abs_gap_q05_q50_q95"]) == 3
    assert len(rpt["rows"]) == 4
    for rr in rpt["rows"]:
        assert rr["bootstrap_size"] == 512
        assert len(rr["top_index_frequency_over_resamples"]) == 7
        assert abs(sum(rr["top_index_frequency_over_resamples"]) - 1.0) < 1e-12
        assert 0.0 <= rr["task7_frequency"] <= 1.0
    assert rpt["task7_frequency_span_over_seeds"] >= 0.0
    assert len(rpt["leader_task_id_per_seed"]) == 4
    assert "controlled_substitution_replay" in rpt
    csr = rpt["controlled_substitution_replay"]
    assert csr["guard_threshold_task7_frequency_span"] >= 0.0
    assert "preconditions" in csr
    assert "leader_stable_over_seeds" in csr["preconditions"]
    assert "task7_frequency_span_bounded" in csr["preconditions"]
    assert csr["status"] in {"SKIPPED_DUE_TO_STABILITY_GUARD", "EXECUTED_LOCAL_ONLY"}
    if csr["executed"]:
        assert isinstance(csr["leader_changed"], bool)
    assert "controlled_substitution_guard_sensitivity" in rpt
    csg = rpt["controlled_substitution_guard_sensitivity"]
    assert len(csg["rows"]) == 4
    for gr in csg["rows"]:
        assert gr["threshold"] > 0.0
        assert isinstance(gr["allow_controlled_replay"], bool)
    assert 0.0 <= csg["allow_frequency_over_threshold_grid"] <= 1.0
    assert "substitution_replay_governance" in rpt
    srg = rpt["substitution_replay_governance"]
    assert isinstance(srg["go_for_actual_substitution_replay"], bool)
    assert 0.0 <= srg["allow_frequency_observed"] <= 1.0
    assert 0.0 <= srg["allow_frequency_threshold"] <= 1.0
    assert srg["reason"] in {"GO", "HOLD_AND_RECALIBRATE"}
    assert "actual_substitution_replay" in rpt
    asr = rpt["actual_substitution_replay"]
    assert isinstance(asr["executed"], bool)
    assert isinstance(asr["governance_go"], bool)
    assert 0.0 <= asr["allow_frequency_observed"] <= 1.0
    assert 0.0 <= asr["allow_frequency_threshold"] <= 1.0
    assert asr["status"] in {"SKIPPED_DUE_TO_GOVERNANCE_HOLD", "EXECUTED_SINGLE_BRANCH_ROBUST_SUBSTITUTION_REPLAY_LOCAL_ONLY"}
    if asr["executed"]:
        assert asr["status"] == "EXECUTED_SINGLE_BRANCH_ROBUST_SUBSTITUTION_REPLAY_LOCAL_ONLY"
        assert isinstance(asr["leader_changed"], bool)
    else:
        assert asr["status"] == "SKIPPED_DUE_TO_GOVERNANCE_HOLD"
    assert "actual_substitution_replay_comparative_report" in rpt
    asrcr = rpt["actual_substitution_replay_comparative_report"]
    assert isinstance(asrcr["executed"], bool)
    assert asrcr["report_scope"] == "LOCAL_SEQUENCING_DIAGNOSTIC_ONLY"
    assert asrcr["status"] in {"SKIPPED_DUE_TO_GOVERNANCE_HOLD", "EXECUTED_COMPARATIVE_REPORT_LOCAL_ONLY"}
    if asrcr["executed"]:
        assert asrcr["status"] == "EXECUTED_COMPARATIVE_REPORT_LOCAL_ONLY"
        assert isinstance(asrcr["leader_changed"], bool)
        assert asrcr["task7_rank_delta_abs"] >= 0.0
        assert asrcr["task7_score_delta_abs"] >= 0.0
        assert asrcr["stability_verdict"] in {"LEADER_STABLE_UNDER_SINGLE_REPLAY", "LEADER_SHIFTED_UNDER_SINGLE_REPLAY"}
    else:
        assert asrcr["status"] == "SKIPPED_DUE_TO_GOVERNANCE_HOLD"
    assert "cross_seed_actual_substitution_replay_panel" in rpt
    csrp = rpt["cross_seed_actual_substitution_replay_panel"]
    assert csrp["report_scope"] == "LOCAL_SEQUENCING_DIAGNOSTIC_ONLY"
    assert csrp["status"] in {"SKIPPED_DUE_TO_GOVERNANCE_HOLD", "EXECUTED_CROSS_SEED_COMPARATIVE_REPORT_LOCAL_ONLY"}
    assert len(csrp["seeds"]) == 4
    if csrp["status"] == "EXECUTED_CROSS_SEED_COMPARATIVE_REPORT_LOCAL_ONLY":
        assert len(csrp["rows"]) == 4
        assert 0.0 <= csrp["leader_changed_frequency_over_seeds"] <= 1.0
        assert len(csrp["task7_rank_delta_abs_q05_q50_q95"]) == 3
        assert len(csrp["task7_score_delta_abs_q05_q50_q95"]) == 3
        assert csrp["stability_verdict"] in {"SEED_ROBUST_LEADER_STABILITY", "SEED_SENSITIVE_LEADER_SHIFT_DETECTED"}
    assert "cross_seed_substitution_governance" in rpt
    cssg = rpt["cross_seed_substitution_governance"]
    assert isinstance(cssg["ready_for_costlier_next_replay_step"], bool)
    assert cssg["reason"] in {"GO_CROSS_SEED_STABLE", "HOLD_AND_RECALIBRATE"}
    assert "criteria" in cssg
    assert "cross_seed_panel_executed" in cssg["criteria"]
    assert "leader_change_frequency_zero" in cssg["criteria"]
    assert "rank_delta_q95_bounded" in cssg["criteria"]
    assert "score_delta_q95_bounded" in cssg["criteria"]
    assert "thresholds" in cssg
    assert cssg["thresholds"]["leader_changed_frequency_max"] >= 0.0
    assert cssg["thresholds"]["task7_rank_delta_abs_q95_max"] >= 0.0
    assert cssg["thresholds"]["task7_score_delta_abs_q95_max"] >= 0.0
    if cssg["criteria"]["cross_seed_panel_executed"]:
        assert "observed" in cssg
        assert 0.0 <= cssg["observed"]["leader_changed_frequency_over_seeds"] <= 1.0
        assert cssg["observed"]["task7_rank_delta_abs_q95"] >= 0.0
        assert cssg["observed"]["task7_score_delta_abs_q95"] >= 0.0
    assert "nonclosure_lock_after_governance" in rpt
    nlg = rpt["nonclosure_lock_after_governance"]
    assert nlg["global_status_must_remain_open_obstruction_with_trace"] is True
    assert nlg["actual_substitution_replay_is_local_only"] is True
    assert nlg["cross_seed_replay_is_local_only"] is True
    assert nlg["costlier_step_readiness_is_not_closure_claim"] is True
    assert "task7_attack_and_task4_verification_packet" in rpt
    t74 = rpt["task7_attack_and_task4_verification_packet"]
    assert t74["status"] in {"HOLD_DUE_TO_GOVERNANCE", "EXECUTED_LOCAL_STRICT_GOVERNANCE_STEP"}
    assert t74["scope"] == "SEQUENCING_EXECUTION_ONLY_NOT_CLOSURE"
    assert "task7_discm_common_basis_attack" in t74
    assert "task4_po3_nonempty_verification" in t74
    if t74["status"] == "EXECUTED_LOCAL_STRICT_GOVERNANCE_STEP":
        t7 = t74["task7_discm_common_basis_attack"]
        t4 = t74["task4_po3_nonempty_verification"]
        assert t7["executed"] is True and t7["result_kind"] == "OPEN_PRECURSOR_NOT_CLOSURE"
        assert t4["executed"] is True and t4["result_kind"] == "OPEN_PRECURSOR_NOT_CLOSURE"
        assert t7["basis_condition_number"] >= 0.0
        assert t7["max_bootstrap_coef_std"] >= 0.0
        assert t7["max_channel_residual_l2"] >= 0.0
        assert isinstance(t4["solver_success"], bool)
        assert t4["objective_value"] >= 0.0
        assert t4["covariant_proxy_d1"] > 0.0
    assert "governance_result_discussion" in rpt
    grd = rpt["governance_result_discussion"]
    assert grd["status"] == "SEQUENCING_DISCUSSION_ONLY_NOT_CLOSURE"
    assert grd["cross_seed_governance_reason"] in {"GO_CROSS_SEED_STABLE", "HOLD_AND_RECALIBRATE"}
    assert isinstance(grd["task7_discm_attack_executed"], bool)
    assert isinstance(grd["task4_po3_verification_executed"], bool)
    assert grd["task7_result_snapshot"]["basis_condition_number"] >= 0.0
    assert grd["task7_result_snapshot"]["max_bootstrap_coef_std"] >= 0.0
    assert grd["task7_result_snapshot"]["max_channel_residual_l2"] >= 0.0
    assert isinstance(grd["task4_result_snapshot"]["solver_success"], bool)
    assert grd["task4_result_snapshot"]["objective_value"] >= 0.0
    assert grd["task4_result_snapshot"]["covariant_proxy_d1"] > 0.0
    assert grd["task4_result_snapshot"]["constraints_hold"] is True
    assert "task7_task4_trend_panel" in rpt
    ttp = rpt["task7_task4_trend_panel"]
    assert ttp["status"] == "LOCAL_TREND_ESTIMATE_NOT_CLOSURE"
    assert ttp["num_runs"] == 3
    assert len(ttp["rows"]) == 3
    assert ttp["task7_residual_l2_span"] >= 0.0
    assert ttp["task7_uncertainty_span"] >= 0.0
    assert ttp["task4_objective_span"] >= 0.0
    assert ttp["task4_covariant_proxy_span"] >= 0.0
    assert ttp["stability_snapshot"] in {"STABLE_LOCAL_TREND", "DRIFT_REVIEW_NEEDED"}
    assert "trend_gate_for_costlier_step" in rpt
    tg = rpt["trend_gate_for_costlier_step"]
    assert tg["scope"] == "SEQUENCING_GOVERNANCE_ONLY_NOT_CLOSURE"
    assert tg["status"] in {"GO_COMPOSITE_GOVERNANCE_STABLE", "HOLD_DUE_TO_COMPOSITE_GOVERNANCE"}
    assert isinstance(tg["ready_for_costlier_step"], bool)
    assert isinstance(tg["criteria"]["cross_seed_governance_go"], bool)
    assert isinstance(tg["criteria"]["nonclosure_lock_active"], bool)
    assert isinstance(tg["criteria"]["trend_stable"], bool)
    assert "composite_nonclosure_enforcement" in rpt
    cne = rpt["composite_nonclosure_enforcement"]
    assert cne["status"] == "ENFORCED"
    assert cne["scope"] == "STRICT_NONCLOSURE_GUARD"
    assert cne["checks"]["global_payload_status_open"] is True
    assert cne["checks"]["all_7_tasks_open"] is True
    assert cne["checks"]["composite_governance_not_interpreted_as_closure"] is True
    assert "nonclosure_status_history_audit" in rpt
    nsha = rpt["nonclosure_status_history_audit"]
    assert nsha["status"] == "AUDIT_TRAIL_LOCAL_PACKET"
    assert nsha["scope"] == "SEQUENCING_AUDIT_ONLY_NOT_CLOSURE"
    assert len(nsha["rows"]) == 5
    assert nsha["all_rows_global_open"] is True
    assert nsha["all_rows_all7_open"] is True
    for hr in nsha["rows"]:
        assert hr["global_status"] == "OPEN_OBSTRUCTION_WITH_TRACE"
        assert hr["all_7_tasks_status"] == "OPEN_OBSTRUCTION_WITH_TRACE"
        assert hr["nonclosure_guard_active"] is True
    assert "governance_nonclosure_consistency_gate" in rpt
    gncg = rpt["governance_nonclosure_consistency_gate"]
    assert gncg["status"] in {"CONSISTENT", "INCONSISTENT"}
    assert gncg["scope"] == "SEQUENCING_GOVERNANCE_ONLY_NOT_CLOSURE"
    assert isinstance(gncg["checks"]["if_go_then_all7_open"], bool)
    assert isinstance(gncg["checks"]["global_payload_open"], bool)
    assert isinstance(gncg["checks"]["history_all_rows_global_open"], bool)
    assert isinstance(gncg["checks"]["history_all_rows_all7_open"], bool)
    assert gncg["status"] == "CONSISTENT"
    assert "governance_nonclosure_failure_simulation" in rpt
    gnfs = rpt["governance_nonclosure_failure_simulation"]
    assert gnfs["status"] == "SIMULATED_FAILURE_DETECTED"
    assert gnfs["scope"] == "TEST_ONLY_DIAGNOSTIC_NOT_RUNTIME_CLAIM"
    assert gnfs["simulated_case"]["ready_for_costlier_step"] is True
    assert gnfs["simulated_case"]["all_7_tasks_open"] is False
    assert gnfs["would_be_consistent_under_simulation"] is False
    assert "governance_history_nonclosure_failure_simulation" in rpt
    ghnfs = rpt["governance_history_nonclosure_failure_simulation"]
    assert ghnfs["status"] == "SIMULATED_FAILURE_DETECTED"
    assert ghnfs["scope"] == "TEST_ONLY_DIAGNOSTIC_NOT_RUNTIME_CLAIM"
    assert ghnfs["simulated_case"]["global_payload_open"] is True
    assert ghnfs["simulated_case"]["history_all_rows_all7_open"] is False
    assert ghnfs["would_be_consistent_under_simulation"] is False
    assert "governance_history_global_nonclosure_failure_simulation" in rpt
    ghgnfs = rpt["governance_history_global_nonclosure_failure_simulation"]
    assert ghgnfs["status"] == "SIMULATED_FAILURE_DETECTED"
    assert ghgnfs["scope"] == "TEST_ONLY_DIAGNOSTIC_NOT_RUNTIME_CLAIM"
    assert ghgnfs["simulated_case"]["if_go_then_all7_open"] is True
    assert ghgnfs["simulated_case"]["global_payload_open"] is True
    assert ghgnfs["simulated_case"]["history_all_rows_global_open"] is False
    assert ghgnfs["would_be_consistent_under_simulation"] is False
    assert "governance_nonclosure_single_flip_matrix" in rpt
    gnsm = rpt["governance_nonclosure_single_flip_matrix"]
    assert gnsm["status"] == "SIMULATED_FAILURE_MATRIX_EXPORTED"
    assert gnsm["scope"] == "TEST_ONLY_DIAGNOSTIC_NOT_RUNTIME_CLAIM"
    assert gnsm["baseline_checks"] == {
        "if_go_then_all7_open": True,
        "global_payload_open": True,
        "history_all_rows_global_open": True,
        "history_all_rows_all7_open": True,
    }
    assert len(gnsm["rows"]) == 4
    seen = set()
    for rr in gnsm["rows"]:
        assert rr["status"] == "SIMULATED_FAILURE_DETECTED"
        assert rr["would_be_consistent_under_simulation"] is False
        assert rr["flipped_check"] in gnsm["baseline_checks"]
        seen.add(rr["flipped_check"])
        false_count = sum(1 for _, v in rr["simulated_checks"].items() if not bool(v))
        assert false_count == 1
    assert seen == set(gnsm["baseline_checks"].keys())
    assert "governance_nonclosure_two_flip_matrix" in rpt
    gntm = rpt["governance_nonclosure_two_flip_matrix"]
    assert gntm["status"] == "SIMULATED_FAILURE_MATRIX_EXPORTED"
    assert gntm["scope"] == "TEST_ONLY_DIAGNOSTIC_NOT_RUNTIME_CLAIM"
    assert gntm["baseline_checks"] == gnsm["baseline_checks"]
    assert gntm["coverage_summary"]["num_checks"] == 4
    assert gntm["coverage_summary"]["expected_rows_n_choose_2"] == 6
    assert gntm["coverage_summary"]["exported_rows"] == 6
    assert gntm["coverage_summary"]["all_rows_have_exactly_two_false"] is True
    assert len(gntm["rows"]) == 6
    seen_pairs = set()
    for rr in gntm["rows"]:
        assert rr["status"] == "SIMULATED_FAILURE_DETECTED"
        assert rr["would_be_consistent_under_simulation"] is False
        assert len(rr["flipped_checks"]) == 2
        pair = tuple(sorted(rr["flipped_checks"]))
        seen_pairs.add(pair)
        false_count = sum(1 for _, v in rr["simulated_checks"].items() if not bool(v))
        assert false_count == 2
    assert len(seen_pairs) == 6
    for row in tpp["rows"]:
        assert "zscore_vs_task_mean" in row
    assert data["backend_loop_fit_precursor"]["loss_l2"] > 0.0
    assert data["backend_loop_fit_precursor"]["loss_gap"] < 1.0
    assert data["backend_loop_fit_precursor"]["multistart_loss_span"] < 1.0
    assert len(data["backend_loop_fit_precursor"]["multistart_rows"]) == 4
    assert len(data["backend_multichannel_discm_loop_precursor"]["channels"]) == 3
    assert data["backend_multichannel_discm_loop_precursor"]["max_solver_loss_gap"] < 1.0
    assert data["backend_multichannel_discm_loop_precursor"]["global_loss_spread"] > 0.0
    assert data["backend_renormalization_b1_precursor"]["max_abs_solver_gap"] < 1e-9
    assert data["backend_renormalization_b1_precursor"]["fit_residual_l2"] > 0.0
    assert len(data["backend_renormalization_b1_precursor"]["rows"]) == 5
    assert data["po3_nonempty_certifier_precursor"]["solver_success"] is True
    assert all(data["po3_nonempty_certifier_precursor"]["constraints"].values())
    assert data["po3_nonempty_certifier_precursor"]["covariant_consistency_proxy_value_d1"] > 0.0
    assert len(data["background_transport_nu_precursor"]["rows"]) == 4
    assert data["background_transport_nu_precursor"]["max_closure_fro_error"] < 1e-10
    assert data["background_transport_nu_precursor"]["max_symmetry_commutator_fro_error"] < 1.0
    assert data["po2_sufficiency_trace_precursor"]["hessian_rank_symbolic"] == 4
    assert data["po2_sufficiency_trace_precursor"]["max_abs_delta_bg_yf_under_constraints"] < 1e-12
    assert len(data["po2_sufficiency_trace_precursor"]["numeric_rows"]) == 8
    assert data["qw2191_selector_premise_precursor"]["status"] == "NON_STRICT_PREMISE_PACKET"
    assert data["qw2191_selector_premise_precursor"]["strict_core_closure_claimed"] is False
    assert data["qw2191_selector_premise_precursor"]["max_ratio_w1_w0"] < 1.0
    assert data["discm_common_basis_integration_precursor"]["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert len(data["discm_common_basis_integration_precursor"]["rows"]) == 3
    assert data["discm_common_basis_integration_precursor"]["common_basis_condition_number"] < 1e4
    assert len(data["channel_phase_space_cutkosky_precursor"]["rows"]) == 3
    assert data["channel_phase_space_cutkosky_precursor"]["global_min_integral"] > 0.0
    assert len(data["channel_phase_space_cutkosky_precursor"]["tolerance_sweep_rows"]) == 3
    assert data["channel_phase_space_cutkosky_precursor"]["tolerance_span_max"] < 1e-10
    assert "ur_uncertainty_transport_bridge_precursor" in data
    uutbp = data["ur_uncertainty_transport_bridge_precursor"]
    assert uutbp["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert uutbp["scope"] == "STRICT_UR_LINK_UNCERTAINTY_SYNTHESIS"
    assert uutbp["rows_count"] == 10
    assert uutbp["median_abs_delta_center"] >= 0.0
    assert uutbp["p95_abs_delta_center"] >= 0.0
    assert uutbp["median_delta_std"] >= 0.0
    assert uutbp["p95_delta_std"] >= 0.0
    assert 0.0 <= uutbp["residue_positive_rate"] <= 1.0
    assert uutbp["bounded_p95_abs_delta_center"] is True
    assert uutbp["bounded_p95_delta_std"] is True
    assert "ur_transport_cross_source_agreement_precursor" in data
    utcsa = data["ur_transport_cross_source_agreement_precursor"]
    assert utcsa["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert utcsa["scope"] == "STRICT_UR_LINK_CROSS_SOURCE_AGREEMENT"
    assert utcsa["common_s_count"] == 5
    assert len(utcsa["rows"]) == 5
    assert utcsa["max_delta_center_abs_gap"] >= 0.0
    assert utcsa["p95_delta_center_abs_gap"] >= 0.0
    assert utcsa["median_delta_std_ratio_p2016_over_p2015"] >= 0.0
    assert utcsa["max_delta_std_ratio_p2016_over_p2015"] >= 0.0
    assert utcsa["center_gap_bounded_p95"] is True
    assert utcsa["std_ratio_bounded_max"] is True
    assert "ur_channel_trace_budget_precursor" in data
    uctbp = data["ur_channel_trace_budget_precursor"]
    assert uctbp["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert uctbp["scope"] == "STRICT_UR_CHANNEL_TRACE_BUDGET"
    assert uctbp["num_channels"] == 4
    assert len(uctbp["rows"]) == 4
    assert uctbp["total_trace_all_channels"] > 0.0
    assert abs(uctbp["trace_share_sum"] - 1.0) < 1e-12
    assert uctbp["all_channels_monotone_nonincreasing"] is True
    assert 0.0 <= uctbp["max_channel_trace_share"] <= 1.0
    assert "ur_channel_class_mapping_precursor" in data
    uccmp = data["ur_channel_class_mapping_precursor"]
    assert uccmp["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert uccmp["scope"] == "STRICT_CHANNEL_MAP_TO_TASK2_CLASSES"
    assert uccmp["mapping_kind"] == "EXPLICIT_WEIGHTED_PRECURSOR_NOT_UNIQUENESS_THEOREM"
    assert len(uccmp["rows"]) == 4
    assert set(uccmp["class_trace_budget"].keys()) == {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
    assert abs(sum(v["trace_share"] for v in uccmp["class_trace_budget"].values()) - 1.0) < 1e-12
    assert uccmp["trace_conservation_gap_abs"] < 1e-12
    assert "ur_class_bounded_uncertainty_residual_budget_precursor" in data
    ucburbp = data["ur_class_bounded_uncertainty_residual_budget_precursor"]
    assert ucburbp["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert ucburbp["scope"] == "STRICT_TASK2_CLASS_BOUNDED_UNCERTAINTY_RESIDUAL_BUDGET"
    assert len(ucburbp["rows"]) == 3
    for rr in ucburbp["rows"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert rr["trace_sum"] >= 0.0
        assert 0.0 <= rr["trace_share"] <= 1.0
        assert rr["residual_l2_backend_sub"] >= 0.0
        assert rr["uncertainty_p95_delta_std"] >= 0.0
        assert rr["risk_proxy_residual_uncertainty_trace"] >= 0.0
    assert ucburbp["risk_proxy_min"] >= 0.0
    assert ucburbp["risk_proxy_max"] >= ucburbp["risk_proxy_min"]
    assert ucburbp["risk_proxy_span"] >= 0.0
    assert ucburbp["all_rows_bounded_uncertainty"] is True
    assert "ur_class_readiness_gate_precursor" in data
    ucrgp = data["ur_class_readiness_gate_precursor"]
    assert ucrgp["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert ucrgp["scope"] == "STRICT_TASK2_CLASS_SEQUENCING_GATE"
    assert len(ucrgp["rows"]) == 3
    assert ucrgp["risk_threshold_go"] >= 0.0
    assert ucrgp["uncertainty_threshold_go"] > 0.0
    assert 0 <= ucrgp["go_count"] <= 3
    assert set(ucrgp["priority_order_low_risk_to_high_risk"]) == {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
    assert "ur_class_first_exact_integration_replay_precursor" in data
    ucfirp = data["ur_class_first_exact_integration_replay_precursor"]
    assert ucfirp["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert ucfirp["scope"] == "STRICT_TASK2_CLASS_FIRST_REPLAY_PACKET"
    assert ucfirp["selected_class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
    assert ucfirp["selection_basis"] == "min_risk_proxy_from_ur_class_readiness_gate_precursor"
    assert ucfirp["baseline_mean_risk_proxy"] >= 0.0
    assert ucfirp["selected_class_risk_proxy"] >= 0.0
    assert ucfirp["selected_class_uncertainty_p95_delta_std"] >= 0.0
    assert 0.0 <= ucfirp["selected_class_trace_share"] <= 1.0
    assert isinstance(ucfirp["ready_for_costlier_exact_integration_replay"], bool)
    assert "ur_class_first_replay_delta_precursor" in data
    ucfrdp = data["ur_class_first_replay_delta_precursor"]
    assert ucfrdp["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert ucfrdp["scope"] == "STRICT_TASK2_CLASS_FIRST_REPLAY_DELTA_PANEL"
    assert ucfrdp["selected_class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
    assert len(ucfrdp["baseline_integrals_over_s_grid"]) == 5
    assert len(ucfrdp["replay_integrals_over_s_grid"]) == 5
    assert ucfrdp["delta_l2_replay_minus_baseline"] >= 0.0
    assert ucfrdp["delta_linf_replay_minus_baseline"] >= 0.0
    assert ucfrdp["replay_settings"]["epsabs"] < 1e-12
    assert "ur_class_first_vs_all_class_replay_comparison_precursor" in data
    ucfvac = data["ur_class_first_vs_all_class_replay_comparison_precursor"]
    assert ucfvac["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert ucfvac["scope"] == "STRICT_TASK2_CLASS_FIRST_VS_ALL_CLASS_REPLAY_COMPARISON"
    assert ucfvac["selected_class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
    assert len(ucfvac["rows"]) == 3
    assert ucfvac["selected_class_delta_l2"] >= 0.0
    assert ucfvac["mean_other_classes_delta_l2"] >= 0.0
    assert isinstance(ucfvac["selected_is_min_delta_l2"], bool)
    assert ucfvac["replay_settings"]["epsabs"] < 1e-12
    assert "ur_cost_vs_gain_precursor" in data
    ucvg = data["ur_cost_vs_gain_precursor"]
    assert ucvg["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert ucvg["scope"] == "STRICT_TASK2_COST_VS_GAIN_PANEL"
    assert ucvg["class_first"]["estimated_cost_units"] > 0.0
    assert ucvg["all_class"]["estimated_cost_units"] > ucvg["class_first"]["estimated_cost_units"]
    assert ucvg["class_first"]["gain_proxy_delta_l2"] >= 0.0
    assert ucvg["all_class"]["gain_proxy_delta_l2_mean"] >= 0.0
    assert ucvg["class_first"]["gain_per_cost"] >= 0.0
    assert ucvg["all_class"]["gain_per_cost"] >= 0.0
    assert isinstance(ucvg["class_first_more_cost_efficient"], bool)
    assert ucvg["cost_ratio_all_over_class_first"] >= 1.0
    assert "ur_runtime_tolerance_benchmark_precursor" in data
    urtbp = data["ur_runtime_tolerance_benchmark_precursor"]
    assert urtbp["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert urtbp["scope"] == "STRICT_TASK2_RUNTIME_TOLERANCE_BENCHMARK"
    assert urtbp["selected_class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
    assert len(urtbp["rows"]) == 3
    for rr in urtbp["rows"]:
        assert rr["tol"] > 0.0
        assert rr["runtime_seconds"] > 0.0
        assert rr["delta_l2_vs_baseline"] >= 0.0
        assert rr["gain_per_second"] >= 0.0
    assert urtbp["slowest_runtime_seconds"] >= urtbp["fastest_runtime_seconds"]
    assert "ur_all_class_exact_integration_sweep_precursor" in data
    uaeisp = data["ur_all_class_exact_integration_sweep_precursor"]
    assert uaeisp["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert uaeisp["scope"] == "STRICT_TASK2_ALL_CLASS_EXACT_INTEGRATION_SWEEP"
    assert uaeisp["num_rows"] == 15
    assert len(uaeisp["rows"]) == 15
    for rr in uaeisp["rows"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert rr["epsabs"] > 0.0
        assert rr["epsrel"] > 0.0
        assert rr["limit"] in {600, 1200, 2000}
        assert rr["delta_l2_vs_baseline"] >= 0.0
        assert rr["delta_linf_vs_baseline"] >= 0.0
        assert rr["integration_warning_count"] >= 0
        assert isinstance(rr["numerical_stress_flag"], bool)
    assert uaeisp["delta_l2_min"] >= 0.0
    assert uaeisp["delta_l2_max"] >= uaeisp["delta_l2_min"]
    assert uaeisp["integration_warning_total"] >= 0
    assert isinstance(uaeisp["any_numerical_stress_flag"], bool)
    assert "ur_numerical_stress_ranking_precursor" in data
    unsrp = data["ur_numerical_stress_ranking_precursor"]
    assert unsrp["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsrp["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_RANKING"
    assert unsrp["ranking_key"] == "integration_warning_count_desc_then_delta_l2_desc"
    assert 1 <= unsrp["top_k"] <= 5
    assert len(unsrp["rows_top_k"]) == unsrp["top_k"]
    for rr in unsrp["rows_top_k"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert rr["integration_warning_count"] >= 0
        assert rr["delta_l2_vs_baseline"] >= 0.0
    assert "ur_numerical_stress_alt_parameterization_precursor" in data
    unsapp = data["ur_numerical_stress_alt_parameterization_precursor"]
    assert unsapp["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsapp["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_ALT_PARAMETERIZATION"
    assert unsapp["transform"] == "x_equals_u_squared"
    assert unsapp["num_rows"] == len(unsapp["rows"])
    assert len(unsapp["rows"]) == uaeisp["num_rows"]
    for rr in unsapp["rows"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert rr["original_integration_warning_count"] >= 0
        assert rr["alt_integration_warning_count"] >= 0
        assert rr["original_delta_l2_vs_baseline"] >= 0.0
        assert rr["alt_delta_l2_vs_baseline"] >= 0.0
    assert "ur_numerical_stress_alt_transform_comparison_precursor" in data
    unsatc = data["ur_numerical_stress_alt_transform_comparison_precursor"]
    assert unsatc["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsatc["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_ALT_TRANSFORM_COMPARISON"
    assert 1 <= unsatc["top_k"] <= 5
    assert len(unsatc["rows_u1"]) == unsatc["top_k"]
    assert len(unsatc["rows_u2"]) == unsatc["top_k"]
    assert len(unsatc["rows_u4"]) == unsatc["top_k"]
    assert unsatc["ranking_key"] == "warning_count_then_delta_l2_then_runtime"
    assert len(unsatc["ranking_rows"]) == unsatc["top_k"]
    for rr in unsatc["rows_u1"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert rr["u1_integration_warning_count"] >= 0
        assert rr["u1_delta_l2_vs_baseline"] >= 0.0
        assert rr["u1_runtime_seconds"] > 0.0
    for rr in unsatc["rows_u4"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert rr["u4_integration_warning_count"] >= 0
        assert rr["u4_delta_l2_vs_baseline"] >= 0.0
    for rr in unsatc["ranking_rows"]:
        assert rr["winner"] in {"u1", "u2", "u4"}
        assert len(rr["ranked_transforms"]) == 3
        assert set(rr["ranked_transforms"]) == {"u1", "u2", "u4"}
    assert "ur_numerical_stress_alt_fullgrid_tritransform_precursor" in data
    unsaft = data["ur_numerical_stress_alt_fullgrid_tritransform_precursor"]
    assert unsaft["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsaft["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_ALT_FULLGRID_TRITRANSFORM"
    assert unsaft["ranking_key"] == "warning_count_then_delta_l2_then_runtime"
    assert unsaft["num_rows"] == 15
    assert len(unsaft["rows"]) == 15
    assert len(unsaft["by_class"]) == 3
    for rr in unsaft["rows"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert rr["winner"] in {"u1", "u2", "u4"}
        assert len(rr["ranked_transforms"]) == 3
        assert set(rr["ranked_transforms"]) == {"u1", "u2", "u4"}
        assert rr["u1_runtime_seconds"] > 0.0 and rr["u2_runtime_seconds"] > 0.0 and rr["u4_runtime_seconds"] > 0.0
    for rr in unsaft["by_class"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert rr["num_rows"] == 5
        assert set(rr["winner_counts"].keys()) == {"u1", "u2", "u4"}
        assert 0.0 <= rr["recommended_transform_frequency"] <= 1.0
    assert "ur_numerical_stress_class_conditional_replay_precursor" in data
    unsccr = data["ur_numerical_stress_class_conditional_replay_precursor"]
    assert unsccr["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsccr["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_CLASS_CONDITIONAL_REPLAY"
    assert unsccr["num_rows"] == 15
    assert len(unsccr["rows"]) == 15
    assert set(unsccr["class_transform_policy"].keys()) == {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
    assert set(unsccr["class_transform_policy"].values()).issubset({"u1", "u2", "u4"})
    for rr in unsccr["rows"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert rr["chosen_transform"] in {"u1", "u2", "u4"}
        assert rr["chosen_warning_count"] >= 0
        assert rr["baseline_warning_count"] >= 0
        assert rr["chosen_delta_l2_vs_baseline"] >= 0.0
        assert rr["chosen_runtime_seconds"] > 0.0
    assert "ur_numerical_stress_policy_counterfactual_precursor" in data
    unspc = data["ur_numerical_stress_policy_counterfactual_precursor"]
    assert unspc["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unspc["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_COUNTERFACTUAL"
    assert unspc["ranking_key"] == "warning_total_then_delta_l2_span_then_runtime_total"
    assert len(unspc["rows"]) == 4
    assert unspc["best_policy"] in {"always_u1", "always_u2", "always_u4", "class_conditional"}
    names = [r["policy"] for r in unspc["rows"]]
    assert set(names) == {"always_u1", "always_u2", "always_u4", "class_conditional"}
    for rr in unspc["rows"]:
        assert rr["warning_total"] >= 0
        assert rr["delta_l2_span"] >= 0.0
        assert rr["runtime_total_seconds"] > 0.0
    assert "ur_numerical_stress_policy_pareto_front_precursor" in data
    unsp_pf = data["ur_numerical_stress_policy_pareto_front_precursor"]
    assert unsp_pf["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_pf["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_PARETO_FRONT"
    assert unsp_pf["axes"] == ["warning_total", "delta_l2_span", "runtime_total_seconds"]
    assert len(unsp_pf["rows"]) == 4
    assert 1 <= unsp_pf["pareto_frontier_count"] <= 4
    assert len(unsp_pf["pareto_frontier_policies"]) == unsp_pf["pareto_frontier_count"]
    for rr in unsp_pf["rows"]:
        assert rr["policy"] in {"always_u1", "always_u2", "always_u4", "class_conditional"}
        assert rr["warning_total"] >= 0
        assert rr["delta_l2_span"] >= 0.0
        assert rr["runtime_total_seconds"] > 0.0
        assert isinstance(rr["pareto_frontier"], bool)
        assert rr["dominated_by_policy"] in {"none", "always_u1", "always_u2", "always_u4", "class_conditional"}
        assert rr["dominance_margin"]["warning_total"] >= 0.0
        assert rr["dominance_margin"]["delta_l2_span"] >= 0.0
        assert rr["dominance_margin"]["runtime_total_seconds"] >= 0.0
    assert "ur_numerical_stress_policy_pareto_stability_precursor" in data
    unsp_ps = data["ur_numerical_stress_policy_pareto_stability_precursor"]
    assert unsp_ps["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_ps["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_PARETO_STABILITY"
    assert unsp_ps["bootstrap_size"] == 512
    assert unsp_ps["resampling_rule"] == "iid_row_resample_with_replacement_over_fullgrid_tri_rows"
    assert len(unsp_ps["rows"]) == 4
    assert unsp_ps["most_stable_frontier_policy"] in {"always_u1", "always_u2", "always_u4", "class_conditional"}
    for rr in unsp_ps["rows"]:
        assert rr["policy"] in {"always_u1", "always_u2", "always_u4", "class_conditional"}
        assert 0.0 <= rr["pareto_front_frequency"] <= 1.0
        assert 0.0 <= rr["pareto_front_frequency_ci95_jeffreys"]["lower"] <= rr["pareto_front_frequency_ci95_jeffreys"]["upper"] <= 1.0
        assert 0 <= rr["bootstrap_successes"] <= rr["bootstrap_trials"]
        assert rr["bootstrap_trials"] == 512
    assert "ur_numerical_stress_policy_budgeted_selection_precursor" in data
    unsp_bs = data["ur_numerical_stress_policy_budgeted_selection_precursor"]
    assert unsp_bs["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_bs["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_BUDGETED_SELECTION"
    assert unsp_bs["selection_rule"] == "argmin_warning_then_delta_l2_then_runtime_under_runtime_cap"
    assert unsp_bs["runtime_budget_quantiles"] == [0.25, 0.5, 0.75, 1.0]
    assert len(unsp_bs["rows"]) == 4
    assert len(unsp_bs["budget_vote_rows"]) == 4
    assert unsp_bs["budget_recommended_policy"] in {"always_u1", "always_u2", "always_u4", "class_conditional"}
    for rr in unsp_bs["rows"]:
        assert rr["runtime_budget_quantile"] in {0.25, 0.5, 0.75, 1.0}
        assert rr["runtime_budget_cap_seconds"] >= 0.0
        assert rr["eligible_policy_count"] >= 0
        assert rr["selected_policy"] in {"none", "always_u1", "always_u2", "always_u4", "class_conditional"}
        assert rr["selected_warning_total"] >= 0
        assert rr["selected_delta_l2_span"] >= 0.0
        assert rr["selected_runtime_total_seconds"] >= 0.0
    for rr in unsp_bs["budget_vote_rows"]:
        assert rr["policy"] in {"always_u1", "always_u2", "always_u4", "class_conditional"}
        assert rr["wins_across_budgets"] >= 0
    assert "ur_numerical_stress_policy_budget_fragility_precursor" in data
    unsp_bf = data["ur_numerical_stress_policy_budget_fragility_precursor"]
    assert unsp_bf["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_bf["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_BUDGET_FRAGILITY"
    assert unsp_bf["base_policy"] in {"always_u1", "always_u2", "always_u4", "class_conditional"}
    assert unsp_bf["cap_scale_grid"] == [0.85, 1.0, 1.15]
    assert unsp_bf["objective_jitter_grid"] == [0.0, 0.01, 0.02]
    assert unsp_bf["bootstrap_size_per_cell"] == 256
    assert len(unsp_bf["rows"]) == 9
    assert 0 <= unsp_bf["global_flip_count_vs_base_policy"] <= unsp_bf["global_trials"]
    assert 0.0 <= unsp_bf["global_flip_frequency_vs_base_policy"] <= 1.0
    assert 0.0 <= unsp_bf["global_flip_frequency_ci95_jeffreys"]["lower"] <= unsp_bf["global_flip_frequency_ci95_jeffreys"]["upper"] <= 1.0
    assert set(unsp_bf["policy_win_counter_global"].keys()) == {"always_u1", "always_u2", "always_u4", "class_conditional"}
    for rr in unsp_bf["rows"]:
        assert rr["runtime_cap_scale"] in {0.85, 1.0, 1.15}
        assert rr["objective_jitter_frac"] in {0.0, 0.01, 0.02}
        assert rr["bootstrap_size"] == 256
        assert 0 <= rr["flip_count_vs_base_policy"] <= rr["local_trials"] <= 256
        assert 0.0 <= rr["flip_frequency_vs_base_policy"] <= 1.0
        assert set(rr["winner_counts"].keys()) == {"always_u1", "always_u2", "always_u4", "class_conditional"}
    assert "ur_numerical_stress_policy_budget_fragility_by_class_precursor" in data
    unsp_bfc = data["ur_numerical_stress_policy_budget_fragility_by_class_precursor"]
    assert unsp_bfc["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_bfc["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_BUDGET_FRAGILITY_BY_CLASS"
    assert len(unsp_bfc["rows"]) == 3
    for rr in unsp_bfc["rows"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert rr["base_policy"] in {"always_u1", "always_u2", "always_u4", "class_conditional", "none"}
        assert 0 <= rr["global_flip_count_vs_base_policy"] <= rr["global_trials"]
        assert 0.0 <= rr["global_flip_frequency_vs_base_policy"] <= 1.0
        assert 0.0 <= rr["global_flip_frequency_ci95_jeffreys"]["lower"] <= rr["global_flip_frequency_ci95_jeffreys"]["upper"] <= 1.0
        assert set(rr["policy_win_counter_global"].keys()) == {"always_u1", "always_u2", "always_u4", "class_conditional"}
    assert "ur_numerical_stress_policy_class_adaptive_fallback_precursor" in data
    unsp_caf = data["ur_numerical_stress_policy_class_adaptive_fallback_precursor"]
    assert unsp_caf["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_caf["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_CLASS_ADAPTIVE_FALLBACK"
    assert 0.0 <= unsp_caf["class_fragility_lb_threshold"] <= 1.0
    assert set(unsp_caf["class_policy_map"].keys()) == {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
    assert len(unsp_caf["rows"]) == 3
    for rr in unsp_caf["rows"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert 0.0 <= rr["fragility_flip_frequency_wilson_lb95_proxy"] <= 1.0
        assert rr["base_policy"] in {"always_u1", "always_u2", "always_u4", "class_conditional", "none"}
        assert rr["robust_winner_policy"] in {"always_u1", "always_u2", "always_u4", "class_conditional"}
        assert isinstance(rr["use_fallback"], bool)
        assert rr["selected_policy"] in {"always_u1", "always_u2", "always_u4", "class_conditional"}
    assert unsp_caf["base_budget_policy"] in {"always_u1", "always_u2", "always_u4", "class_conditional"}
    assert unsp_caf["adaptive_warning_total"] >= 0
    assert unsp_caf["adaptive_delta_l2_span"] >= 0.0
    assert unsp_caf["adaptive_runtime_total_seconds"] >= 0.0
    assert "ur_numerical_stress_policy_ablation_precursor" in data
    unsp_ab = data["ur_numerical_stress_policy_ablation_precursor"]
    assert unsp_ab["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_ab["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_ABLATION"
    assert len(unsp_ab["rows"]) == 3
    assert unsp_ab["best_regime_lexicographic"] in {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only"}
    for rr in unsp_ab["rows"]:
        assert rr["regime"] in {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only"}
        assert set(rr["class_policy_map"].keys()) == {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert rr["warning_total"] >= 0
        assert rr["delta_l2_span"] >= 0.0
        assert rr["runtime_total_seconds"] >= 0.0
    assert len(unsp_ab["pairwise_dominance_rows"]) == 3
    for rr in unsp_ab["pairwise_dominance_rows"]:
        assert rr["regime_a"] in {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only"}
        assert rr["regime_b"] in {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only"}
        assert rr["regime_a"] != rr["regime_b"]
        assert rr["bootstrap_size"] == 256
        assert 0.0 <= rr["a_dominates_b_frequency"] <= 1.0
        assert 0.0 <= rr["b_dominates_a_frequency"] <= 1.0
        assert 0.0 <= rr["tie_or_incomparable_frequency"] <= 1.0
        assert abs(rr["a_dominates_b_frequency"] + rr["b_dominates_a_frequency"] + rr["tie_or_incomparable_frequency"] - 1.0) < 1e-12
    assert "ur_numerical_stress_policy_cross_class_constrained_ablation_precursor" in data
    unsp_cc = data["ur_numerical_stress_policy_cross_class_constrained_ablation_precursor"]
    assert unsp_cc["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_cc["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_CROSS_CLASS_CONSTRAINED_ABLATION"
    assert unsp_cc["constraint_rule"] == "row_runtime_le_class_q75_u2_and_row_warning_le_class_q75_u2"
    assert unsp_cc["feasibility_pass_rate_threshold"] == 0.80
    assert set(unsp_cc["class_runtime_caps_seconds"].keys()) == {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
    assert set(unsp_cc["class_warning_caps"].keys()) == {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
    assert len(unsp_cc["rows"]) == 3
    assert unsp_cc["best_feasible_regime_lexicographic"] in {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only", "none"}
    for rr in unsp_cc["rows"]:
        assert rr["regime"] in {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only"}
        assert 0 <= rr["constraint_pass_rows"] <= rr["constraint_total_rows"]
        assert rr["constraint_total_rows"] > 0
        assert 0.0 <= rr["constraint_pass_rate"] <= 1.0
        assert rr["warning_total"] >= 0
        assert rr["delta_l2_span"] >= 0.0
        assert rr["runtime_total_seconds"] >= 0.0
    assert "ur_numerical_stress_policy_cross_class_constrained_bootstrap_dominance_precursor" in data
    unsp_ccd = data["ur_numerical_stress_policy_cross_class_constrained_bootstrap_dominance_precursor"]
    assert unsp_ccd["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_ccd["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_CROSS_CLASS_CONSTRAINED_BOOTSTRAP_DOMINANCE"
    assert len(unsp_ccd["rows"]) == 3
    for rr in unsp_ccd["rows"]:
        assert rr["regime_a"] in {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only"}
        assert rr["regime_b"] in {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only"}
        assert rr["regime_a"] != rr["regime_b"]
        assert 0 <= rr["usable_bootstrap_trials"] <= rr["bootstrap_size_requested"] == 256
        assert 0.0 <= rr["a_dominates_b_frequency"] <= 1.0
        assert 0.0 <= rr["b_dominates_a_frequency"] <= 1.0
        assert 0.0 <= rr["tie_or_incomparable_frequency"] <= 1.0
        assert abs(rr["a_dominates_b_frequency"] + rr["b_dominates_a_frequency"] + rr["tie_or_incomparable_frequency"] - 1.0) < 1e-12
    assert "ur_numerical_stress_policy_cross_class_threshold_sweep_precursor" in data
    unsp_th = data["ur_numerical_stress_policy_cross_class_threshold_sweep_precursor"]
    assert unsp_th["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_th["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_CROSS_CLASS_THRESHOLD_SWEEP"
    assert unsp_th["threshold_grid"] == [0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
    assert len(unsp_th["rows"]) == len(unsp_th["threshold_grid"]) == 6
    assert len(unsp_th["best_regime_sequence"]) == 6
    assert isinstance(unsp_th["best_regime_stable_over_nonempty_thresholds"], bool)
    for rr in unsp_th["rows"]:
        assert rr["feasibility_pass_rate_threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
        assert 0 <= rr["num_feasible_regimes"] <= 3
        assert rr["best_feasible_regime_lexicographic"] in {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only", "none"}
    assert "ur_numerical_stress_policy_joint_stress_map_precursor" in data
    unsp_js = data["ur_numerical_stress_policy_joint_stress_map_precursor"]
    assert unsp_js["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_js["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_JOINT_STRESS_MAP"
    assert unsp_js["threshold_grid"] == [0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
    assert unsp_js["cap_scale_grid"] == [0.85, 1.0, 1.15]
    assert unsp_js["jitter_grid"] == [0.0, 0.01, 0.02]
    assert len(unsp_js["rows"]) == unsp_js["num_cells"] == 54
    assert 0 <= unsp_js["num_stable_cells_winner_freq_ge_070"] <= unsp_js["num_cells"]
    assert 0.0 <= unsp_js["stable_cell_frequency"] <= 1.0
    for rr in unsp_js["rows"]:
        assert rr["threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
        assert rr["cap_scale"] in {0.85, 1.0, 1.15}
        assert rr["jitter"] in {0.0, 0.01, 0.02}
        assert rr["bootstrap_size_requested"] == 128
        assert 0 <= rr["usable_trials"] <= rr["bootstrap_size_requested"]
        assert rr["winner"] in {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only", "none"}
        assert 0.0 <= rr["winner_frequency"] <= 1.0
        assert 0.0 <= rr["winner_frequency_ci95_jeffreys"]["lower"] <= rr["winner_frequency_ci95_jeffreys"]["upper"] <= 1.0
        assert 0.0 <= rr["winner_entropy_norm"] <= 1.0
        assert set(rr["winner_counts"].keys()) == {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only"}
    assert "ur_numerical_stress_policy_stability_topology_precursor" in data
    unsp_top = data["ur_numerical_stress_policy_stability_topology_precursor"]
    assert unsp_top["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_top["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_STABILITY_TOPOLOGY"
    assert unsp_top["stable_cell_definition"] == "winner_frequency_ge_0p70_and_winner_not_none"
    assert unsp_top["num_components"] >= 0
    assert unsp_top["largest_component_size"] >= 0
    assert len(unsp_top["components"]) == unsp_top["num_components"]
    for cc in unsp_top["components"]:
        assert cc["size"] == len(cc["rows"])
        assert cc["dominant_winner"] in {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only"}
        assert set(cc["winner_counts"].keys()) == {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only"}
        for rr in cc["rows"]:
            assert rr["threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
            assert rr["cap_scale"] in {0.85, 1.0, 1.15}
            assert rr["jitter"] in {0.0, 0.01, 0.02}
            assert rr["winner"] in {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only"}
    assert "ur_numerical_stress_policy_stability_boundary_margin_precursor" in data
    unsp_bm = data["ur_numerical_stress_policy_stability_boundary_margin_precursor"]
    assert unsp_bm["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_bm["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_STABILITY_BOUNDARY_MARGIN"
    assert unsp_bm["distance_metric"] == "manhattan_on_threshold_cap_jitter_grid"
    assert unsp_bm["num_stable_rows"] == len(unsp_bm["rows"])
    for rr in unsp_bm["rows"]:
        assert rr["threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
        assert rr["cap_scale"] in {0.85, 1.0, 1.15}
        assert rr["jitter"] in {0.0, 0.01, 0.02}
        assert rr["winner"] in {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only"}
        assert rr["boundary_manhattan_margin"] >= 0
    for rr in unsp_bm["component_margin_rows"]:
        assert rr["dominant_winner"] in {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only"}
        assert rr["size"] >= 0
        q = rr["boundary_margin_q05_q50_q95"]
        assert len(q) == 3 and q[0] <= q[1] <= q[2]
        assert rr["boundary_margin_min"] <= rr["boundary_margin_max"]
    assert "ur_numerical_stress_policy_weighted_boundary_risk_precursor" in data
    unsp_wr = data["ur_numerical_stress_policy_weighted_boundary_risk_precursor"]
    assert unsp_wr["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_wr["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_WEIGHTED_BOUNDARY_RISK"
    assert "winner_frequency" in unsp_wr["risk_definition"]
    assert unsp_wr["margin_normalization_divisor"] >= 1.0
    assert len(unsp_wr["rows"]) == 54
    assert len(unsp_wr["best_corridor_rows"]) <= 8
    for rr in unsp_wr["rows"]:
        assert rr["threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
        assert rr["cap_scale"] in {0.85, 1.0, 1.15}
        assert rr["jitter"] in {0.0, 0.01, 0.02}
        assert rr["winner"] in {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only", "none"}
        assert 0.0 <= rr["winner_frequency"] <= 1.0
        assert 0.0 <= rr["winner_entropy_norm"] <= 1.0
        assert rr["boundary_manhattan_margin"] >= 0
        assert 0.0 <= rr["boundary_margin_norm"] <= 1.0
        assert 0.0 <= rr["weighted_boundary_risk_score"] <= 1.0
    assert "ur_numerical_stress_policy_weighted_boundary_risk_sensitivity_precursor" in data
    unsp_ws = data["ur_numerical_stress_policy_weighted_boundary_risk_sensitivity_precursor"]
    assert unsp_ws["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_ws["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_WEIGHTED_BOUNDARY_RISK_SENSITIVITY"
    assert len(unsp_ws["rows"]) == 5
    assert isinstance(unsp_ws["winner_set_stable_over_weight_grid"], bool)
    for rr in unsp_ws["rows"]:
        w = rr["weights"]
        assert abs((w["wf"] + w["ent"] + w["margin"]) - 1.0) < 1e-12
        assert 0.0 <= w["wf"] <= 1.0 and 0.0 <= w["ent"] <= 1.0 and 0.0 <= w["margin"] <= 1.0
        assert len(rr["best_corridor_rows"]) <= 8
        for cr in rr["best_corridor_rows"]:
            assert cr["winner"] in {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only"}
            assert 0.0 <= cr["weighted_boundary_risk_score"] <= 1.0
    assert "ur_numerical_stress_policy_weighted_boundary_risk_bayesian_precursor" in data
    unsp_wb = data["ur_numerical_stress_policy_weighted_boundary_risk_bayesian_precursor"]
    assert unsp_wb["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_wb["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_WEIGHTED_BOUNDARY_RISK_BAYESIAN"
    assert len(unsp_wb["dirichlet_alpha"]) == 3
    assert unsp_wb["posterior_sample_size"] == 512
    assert len(unsp_wb["winner_set_posterior_rows"]) >= 1
    p_sum = 0.0
    for rr in unsp_wb["winner_set_posterior_rows"]:
        assert rr["count"] >= 0
        assert 0.0 <= rr["posterior_probability"] <= 1.0
        assert 0.0 <= rr["posterior_probability_ci95_jeffreys"]["lower"] <= rr["posterior_probability_ci95_jeffreys"]["upper"] <= 1.0
        p_sum += rr["posterior_probability"]
    assert abs(p_sum - 1.0) < 1e-12
    assert 0.0 <= unsp_wb["most_probable_winner_set_probability"] <= 1.0
    assert set(unsp_wb["best_cell_winner_posterior_counts"].keys()) == {"base_budget_policy", "class_adaptive_fallback", "robust_winner_only", "none"}
    assert "task2_strict_unitarity_witness" in data
    t2w = data["task2_strict_unitarity_witness"]
    assert t2w["theorem_target"].startswith("bounded DiscM-CutSum gap")
    assert t2w["channel"] == "graviton_to_gauge_gauge"
    assert len(t2w["s_grid"]) >= 10
    assert len(t2w["rows"]) == len(t2w["s_grid"])
    assert len(t2w["residue_or_weight_sign_proxy_rows"]) >= 5
    for rr in t2w["rows"]:
        assert rr["s"] in set(t2w["s_grid"])
        assert rr["disc_value"] >= 0.0 and rr["cutsum_value"] >= 0.0
        assert rr["gap_abs"] >= 0.0 and rr["gap_rel"] >= 0.0 and rr["uncertainty_estimate"] >= 0.0
        assert rr["cutsum_native"] >= 0.0 and rr["cutsum_warped"] >= 0.0
        assert rr["cutsum_scheme_gap_abs"] >= 0.0
    for rr in t2w["residue_or_weight_sign_proxy_rows"]:
        assert rr["s"] in set(t2w["s_grid"])
        assert rr["min_effective_weight"] <= rr["max_effective_weight"]
        assert isinstance(rr["all_nonnegative"], bool)
    assert t2w["aggregate_metrics"]["max_gap_abs"] >= 0.0
    assert t2w["aggregate_metrics"]["q95_gap_abs"] >= 0.0
    assert t2w["aggregate_metrics"]["max_gap_rel"] >= 0.0
    ci = t2w["aggregate_metrics"]["consistency_ci95"]
    assert 0.0 <= ci["lower"] <= ci["upper"]
    assert isinstance(t2w["aggregate_metrics"]["all_nonnegative_weights"], bool)
    assert isinstance(t2w["aggregate_metrics"]["min_effective_weight_global"], float)
    assert isinstance(t2w["aggregate_metrics"]["obstruction_is_numerically_sensitive"], bool)
    assert t2w["aggregate_metrics"]["max_bin_disc_proxy_gap_abs"] >= 0.0
    assert t2w["aggregate_metrics"]["max_bin_scheme_gap_abs"] >= 0.0
    assert t2w["aggregate_metrics"]["q95_cross_integrator_gap_abs"] >= 0.0
    assert t2w["aggregate_metrics"]["q95_convergence_delta_n400_to_n800_abs"] >= 0.0
    assert len(t2w["phase_space_bin_contribution_rows"]) == len(t2w["rows"])
    assert "bin_obstruction_ranking" in t2w
    bor = t2w["bin_obstruction_ranking"]
    assert len(bor["rows"]) == 8
    assert bor["total_disc_proxy_gap_sum"] >= 0.0
    assert 0.0 <= bor["top2_disc_proxy_gap_share"] <= 1.0
    share_sum = 0.0
    for rr in bor["rows"]:
        assert 0 <= rr["bin_index"] <= 7
        assert 0.0 <= rr["x_left"] < rr["x_right"] <= 1.0
        assert rr["disc_proxy_gap_sum"] >= 0.0
        assert rr["scheme_gap_sum"] >= 0.0
        assert rr["mean_disc_proxy_gap"] >= 0.0
        assert 0.0 <= rr["disc_proxy_gap_share"] <= 1.0
        share_sum += rr["disc_proxy_gap_share"]
    assert abs(share_sum - 1.0) < 1e-9 or bor["total_disc_proxy_gap_sum"] == 0.0
    assert "endpoint_refinement" in t2w
    erf = t2w["endpoint_refinement"]
    assert erf["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert erf["scope"] == "STRICT_TASK2_ENDPOINT_REFINEMENT_TOP2_BINS"
    assert len(erf["top2_bins"]) == 2
    assert all(0 <= int(b) <= 7 for b in erf["top2_bins"])
    assert len(erf["rows"]) == len(t2w["s_grid"])
    assert erf["q95_gap_abs_baseline"] >= 0.0
    assert erf["q95_gap_abs_refined"] >= 0.0
    for rr in erf["rows"]:
        assert rr["s"] in set(t2w["s_grid"])
        assert rr["disc_value"] >= 0.0
        assert rr["cutsum_value_refined"] >= 0.0
        assert rr["gap_abs_refined"] >= 0.0
        assert len(rr["rows"]) == 8
        for br in rr["rows"]:
            assert 0 <= br["bin_index"] <= 7
            assert isinstance(br["is_top2"], bool)
            assert br["bin_integral_refined"] >= 0.0
    assert "endpoint_split_domain" in t2w
    esd = t2w["endpoint_split_domain"]
    assert esd["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert esd["scope"] == "STRICT_TASK2_TOP2_BIN_ENDPOINT_SPLIT_DOMAIN"
    assert len(esd["top2_bins"]) == 2
    assert esd["left_total_abs"] >= 0.0 and esd["right_total_abs"] >= 0.0
    assert esd["left_to_right_abs_ratio"] >= 0.0
    assert esd["dominant_endpoint_half"] in {"LEFT", "RIGHT", "BALANCED"}
    assert len(esd["rows"]) == len(t2w["s_grid"]) * 2
    assert "endpoint_adaptive_transform" in t2w
    eat = t2w["endpoint_adaptive_transform"]
    assert eat["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert eat["scope"] == "STRICT_TASK2_ENDPOINT_ADAPTIVE_TRANSFORM_SELECTION"
    assert eat["dominant_endpoint_half"] in {"LEFT", "RIGHT", "BALANCED"}
    assert eat["transform_grid"] == ["native", "left_focus", "right_focus"]
    assert len(eat["rows"]) == 3
    assert eat["recommended_transform"] in {"native", "left_focus", "right_focus"}
    assert eat["q95_gap_abs_baseline"] >= 0.0
    for rr in eat["rows"]:
        assert rr["transform"] in {"native", "left_focus", "right_focus"}
        assert rr["q95_gap_abs"] >= 0.0
    assert eat["bootstrap_size"] == 256
    assert len(eat["bootstrap_rows"]) == 3
    assert eat["most_frequent_bootstrap_transform"] in {"native", "left_focus", "right_focus"}
    assert "conditional_recompute_gate" in eat
    crg = eat["conditional_recompute_gate"]
    assert crg["status"] in {"EXECUTED_LOCAL_RECOMPUTE", "SKIPPED_DUE_TO_STABILITY_CI95_GATE"}
    assert crg["scope"] == "STRICT_TASK2_ENDPOINT_ADAPTIVE_RECOMPUTE_GATE"
    assert crg["ci95_lower_threshold"] == 0.6
    assert crg["dominant_transform"] in {"native", "left_focus", "right_focus"}
    assert 0.0 <= crg["dominant_transform_ci95_lower"] <= 1.0
    assert isinstance(crg["gate_passed"], bool)
    if crg["gate_passed"]:
        assert crg["status"] == "EXECUTED_LOCAL_RECOMPUTE"
        assert crg["q95_gap_abs_recompute"] >= 0.0
        assert crg["max_gap_rel_recompute"] >= 0.0
        assert isinstance(crg["delta_q95_gap_abs_recompute_minus_baseline"], float)
        assert isinstance(crg["delta_max_gap_rel_recompute_minus_baseline"], float)
    else:
        assert crg["status"] == "SKIPPED_DUE_TO_STABILITY_CI95_GATE"
        assert crg["q95_gap_abs_recompute"] is None
        assert crg["max_gap_rel_recompute"] is None
        assert crg["delta_q95_gap_abs_recompute_minus_baseline"] is None
        assert crg["delta_max_gap_rel_recompute_minus_baseline"] is None
    total_freq = 0.0
    for rr in eat["bootstrap_rows"]:
        assert rr["transform"] in {"native", "left_focus", "right_focus"}
        assert 0 <= rr["selection_count"] <= 256
        assert 0.0 <= rr["selection_frequency"] <= 1.0
        assert 0.0 <= rr["selection_frequency_ci95_jeffreys"]["lower"] <= rr["selection_frequency_ci95_jeffreys"]["upper"] <= 1.0
        assert rr["q95_gap_abs_point_estimate"] >= 0.0
        total_freq += rr["selection_frequency"]
    assert abs(total_freq - 1.0) < 1e-12
    for rr in esd["rows"]:
        assert rr["s"] in set(t2w["s_grid"])
        assert rr["bin_index"] in set(esd["top2_bins"])
        assert rr["x_left"] < rr["x_mid"] < rr["x_right"]
        assert rr["left_half_integral_abs"] >= 0.0
        assert rr["right_half_integral_abs"] >= 0.0
        assert rr["left_minus_right_abs"] >= 0.0
        assert isinstance(rr["left_dominates"], bool)
    for sr in t2w["phase_space_bin_contribution_rows"]:
        assert sr["s"] in set(t2w["s_grid"])
        assert len(sr["rows"]) == 8
        assert sr["max_bin_disc_proxy_gap_abs"] >= 0.0
        assert sr["max_bin_scheme_gap_abs"] >= 0.0
        for br in sr["rows"]:
            assert 0 <= br["bin_index"] <= 7
            assert 0.0 <= br["x_left"] < br["x_right"] <= 1.0
            assert br["native_bin_integral"] >= 0.0
            assert br["warped_bin_integral"] >= 0.0
            assert br["bin_scheme_gap_abs"] >= 0.0
            assert br["bin_disc_proxy_gap_abs"] >= 0.0
    assert t2w["verdict"] in {"OPEN_OBSTRUCTION_WITH_TRACE", "CLOSED_NUMERICAL_WITNESS_TASK2"}
    if t2w["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert ">" in t2w["fail_trace"]
        assert any(
            tag in t2w["fail_trace"]
            for tag in {"q95_gap_abs=", "max_gap_rel=", "q95_cross_integrator_gap_abs=", "q95_convergence_delta_n400_to_n800_abs=", "q95_refined_window_best_gap_abs=", "q95_gap_abs_quad_high_precision_top3=", "q95_gap_abs_quad_high_precision_upper_envelope_top3=", "q95_gap_abs_upper_tail_envelope_top3=", "q95_gap_abs_n2400_top3=", "q95_delta_gap_abs_n2400_minus_n6400_abs_top3=", "q95_delta_gap_abs_n3200_minus_n6400_abs_top3=", "q95_delta_gap_abs_n1600_minus_n6400_abs_top3=", "q95_delta_gap_abs_n800_minus_n6400_abs_top3=", "q95_delta_gap_abs_n12800_minus_n6400_abs_top3=", "q95_total_monotone_violation_top3=", "min_effective_weight_global="}
        )
    assert "closure_consistency" in t2w
    cc = t2w["closure_consistency"]
    assert set(cc["criteria_evaluation"].keys()) == {
        "q95_gap_abs_le_threshold",
        "q95_gap_abs_n2400_le_threshold",
        "q95_gap_abs_n2400_vs_n6400_delta_le_threshold",
        "q95_gap_abs_n3200_vs_n6400_delta_le_threshold",
        "q95_gap_abs_n1600_vs_n6400_delta_le_threshold",
        "q95_gap_abs_n800_vs_n6400_delta_le_threshold",
        "q95_gap_abs_n12800_vs_n6400_delta_le_threshold",
        "q95_convergence_tail_ratio_n12800_le_one",
        "q95_gap_abs_n25600_vs_n12800_delta_le_threshold",
        "q95_monotone_violation_n25600_le_zero",
        "max_gap_rel_le_threshold",
        "all_nonnegative_weights",
        "q95_cross_integrator_gap_le_threshold",
        "q95_convergence_delta_le_threshold",
    }
    assert isinstance(cc["all_criteria_satisfied"], bool)
    assert cc["dominant_blocker"] in {"q95_gap_abs", "max_gap_rel", "weight_sign_nonnegativity", "q95_cross_integrator_gap", "q95_convergence_delta_n400_to_n800_abs", "none", "q95_refined_window_gap_abs", "q95_tail_budget_upper_envelope_gap_abs", "q95_n2400_gap_abs", "q95_n2400_vs_n6400_delta_abs", "q95_n3200_vs_n6400_delta_abs", "q95_n1600_vs_n6400_delta_abs", "q95_n800_vs_n6400_delta_abs", "q95_n12800_vs_n6400_delta_abs", "q95_n12800_tail_ratio", "q95_n25600_vs_n12800_delta_abs", "q95_n25600_monotone_violation"}
    assert isinstance(cc["dominant_inequality"], str) and len(cc["dominant_inequality"]) > 0
    assert "falsifier_trace_consistency" in t2w
    ftc = t2w["falsifier_trace_consistency"]
    assert ftc["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert ftc["scope"] == "STRICT_TASK2_FALSIFIER_TRACE_CONSISTENCY"
    assert set(ftc["checks"].keys()) == {
        "if_open_then_fail_trace_equals_dominant_inequality",
        "if_closed_then_fail_trace_empty",
        "if_open_then_dominant_not_none",
    }
    assert all(isinstance(v, bool) for v in ftc["checks"].values())
    assert ftc["all_checks_pass"] is True
    assert "q95_dominant_s_attribution" in t2w
    qsa = t2w["q95_dominant_s_attribution"]
    assert qsa["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qsa["scope"] == "STRICT_TASK2_Q95_BLOCKER_S_POINT_ATTRIBUTION"
    assert qsa["q95_gap_abs_reference"] >= 0.0
    assert len(qsa["rows_sorted_by_gap_abs_desc"]) == len(t2w["s_grid"])
    assert len(qsa["top3_rows"]) <= 3
    assert 0.0 <= qsa["top3_gap_abs_share_of_total"] <= 1.0
    assert 0 <= qsa["q95_tail_count"] <= len(t2w["s_grid"])
    assert qsa["q95_tail_abs_mean"] >= 0.0
    assert qsa["q95_tail_abs_std"] >= 0.0
    if qsa["dominant_s_value"] is not None:
        assert qsa["dominant_s_value"] in set(t2w["s_grid"])
    prev_gap = None
    for rr in qsa["rows_sorted_by_gap_abs_desc"]:
        assert rr["s"] in set(t2w["s_grid"])
        assert rr["gap_abs"] >= 0.0
        assert rr["gap_rel"] >= 0.0
        assert isinstance(rr["is_q95_tail"], bool)
        if prev_gap is not None:
            assert prev_gap + 1e-18 >= rr["gap_abs"]
        prev_gap = rr["gap_abs"]
    assert "q95_dominant_s_crosscheck" in t2w
    qsc = t2w["q95_dominant_s_crosscheck"]
    assert qsc["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qsc["scope"] == "STRICT_TASK2_Q95_DOMINANT_S_DUAL_INTEGRATOR_CROSSCHECK"
    assert len(qsc["rows"]) <= 3
    assert qsc["crosscheck_gap_abs_max"] >= 0.0
    assert qsc["crosscheck_gap_abs_q95"] >= 0.0
    assert qsc["crosscheck_gap_abs_max"] >= qsc["crosscheck_gap_abs_q95"]
    assert qsc["q95_gap_abs_quad_top3"] >= 0.0
    assert qsc["q95_gap_abs_fixed_quad_top3"] >= 0.0
    assert isinstance(qsc["delta_q95_fixed_minus_quad_top3"], float)
    for rr in qsc["rows"]:
        assert rr["s"] in set(t2w["s_grid"])
        assert rr["disc_quad"] >= 0.0
        assert rr["cutsum_quad"] >= 0.0
        assert rr["cutsum_fixed_quad_n400"] >= 0.0
        assert rr["gap_abs_quad"] >= 0.0
        assert rr["gap_abs_fixed_quad"] >= 0.0
        assert rr["crosscheck_gap_abs"] >= 0.0
    assert "q95_dominant_s_sign_check" in t2w
    qss = t2w["q95_dominant_s_sign_check"]
    assert qss["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qss["scope"] == "STRICT_TASK2_Q95_DOMINANT_S_SIGN_RESIDUE_CHECK"
    assert isinstance(qss["symbolic_integrand"], str) and len(qss["symbolic_integrand"]) > 0
    assert isinstance(qss["symbolic_nonnegativity_certificate"], bool)
    assert isinstance(qss["all_rows_nonnegative_numeric"], bool)
    assert len(qss["rows"]) <= 3
    for rr in qss["rows"]:
        assert rr["s"] in set(t2w["s_grid"])
        assert rr["integrand_min_over_x_grid"] >= -1e-12
        assert rr["integrand_max_over_x_grid"] >= rr["integrand_min_over_x_grid"]
        assert isinstance(rr["all_nonnegative_over_x_grid"], bool)
    assert "q95_dominant_s_convergence" in t2w
    qcv = t2w["q95_dominant_s_convergence"]
    assert qcv["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qcv["scope"] == "STRICT_TASK2_Q95_DOMINANT_S_FIXED_QUAD_CONVERGENCE"
    assert len(qcv["rows"]) <= 3
    assert qcv["max_delta_n400_to_n800_abs"] >= 0.0
    assert qcv["median_convergence_ratio_400_800_over_200_400"] >= 0.0
    assert qcv["q95_gap_abs_n800_top3"] >= 0.0
    for rr in qcv["rows"]:
        assert rr["s"] in set(t2w["s_grid"])
        assert rr["disc_reference"] >= 0.0
        assert rr["cutsum_fixed_quad_n200"] >= 0.0
        assert rr["cutsum_fixed_quad_n400"] >= 0.0
        assert rr["cutsum_fixed_quad_n800"] >= 0.0
        assert rr["gap_abs_n200"] >= 0.0
        assert rr["gap_abs_n400"] >= 0.0
        assert rr["gap_abs_n800"] >= 0.0
        assert rr["delta_n200_to_n400_abs"] >= 0.0
        assert rr["delta_n400_to_n800_abs"] >= 0.0
        assert rr["convergence_ratio_400_800_over_200_400"] >= 0.0
    assert "q95_blocker_n1600_recompute_certificate" in t2w
    qn16 = t2w["q95_blocker_n1600_recompute_certificate"]
    assert qn16["scope"] == "STRICT_TASK2_Q95_BLOCKER_N1600_RECOMPUTE_CERTIFICATE"
    assert qn16["theorem_target"] == "Q95_GAP_ABS_N1600_TOP3_LE_THRESHOLD"
    assert qn16["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    for rr in qn16["computed_rows"]:
        assert rr["cutsum_fixed_quad_n1600"] >= 0.0
        assert rr["gap_abs_n1600"] >= 0.0
    if qn16["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "q95_gap_abs_n1600_top3=" in qn16["fail_trace"] and ">" in qn16["fail_trace"]
    assert "q95_blocker_n3200_recompute_certificate" in t2w
    qn32 = t2w["q95_blocker_n3200_recompute_certificate"]
    assert qn32["scope"] == "STRICT_TASK2_Q95_BLOCKER_N3200_RECOMPUTE_CERTIFICATE"
    assert qn32["theorem_target"] == "Q95_GAP_ABS_N3200_TOP3_LE_THRESHOLD"
    assert qn32["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    for rr in qn32["computed_rows"]:
        assert rr["cutsum_fixed_quad_n3200"] >= 0.0
        assert rr["gap_abs_n3200"] >= 0.0
    if qn32["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "q95_gap_abs_n3200_top3=" in qn32["fail_trace"] and ">" in qn32["fail_trace"]
    assert "q95_blocker_n6400_recompute_certificate" in t2w
    qn64 = t2w["q95_blocker_n6400_recompute_certificate"]
    assert qn64["scope"] == "STRICT_TASK2_Q95_BLOCKER_N6400_RECOMPUTE_CERTIFICATE"
    assert qn64["theorem_target"] == "Q95_GAP_ABS_N6400_TOP3_LE_THRESHOLD"
    assert qn64["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    for rr in qn64["computed_rows"]:
        assert rr["cutsum_fixed_quad_n6400"] >= 0.0
        assert rr["gap_abs_n6400"] >= 0.0
    if qn64["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "q95_gap_abs_n6400_top3=" in qn64["fail_trace"] and ">" in qn64["fail_trace"]
    assert "q95_blocker_ninf_extrapolation_certificate" in t2w
    qninf = t2w["q95_blocker_ninf_extrapolation_certificate"]
    assert qninf["scope"] == "STRICT_TASK2_Q95_BLOCKER_NINF_EXTRAPOLATION_CERTIFICATE"
    assert qninf["theorem_target"] == "Q95_GAP_ABS_EXTRAPOLATED_NINF_UPPER_TOP3_LE_THRESHOLD"
    assert qninf["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    assert qninf["domain"]["fixed_quad_n_levels"] == [1600, 3200, 6400]
    for rr in qninf["computed_rows"]:
        assert rr["gap_abs_n1600"] >= 0.0
        assert rr["gap_abs_n3200"] >= 0.0
        assert rr["gap_abs_n6400"] >= 0.0
        assert rr["delta_1600_3200_abs"] >= 0.0
        assert rr["delta_3200_6400_abs"] >= 0.0
        assert rr["effective_order_p"] >= 1.0
        assert rr["gap_abs_extrapolated_ninf"] >= 0.0
        assert rr["gap_abs_extrapolation_error_abs"] >= 0.0
    assert qninf["aggregate_metrics"]["q95_gap_abs_extrapolated_ninf_upper_top3"] >= qninf["aggregate_metrics"]["q95_gap_abs_extrapolated_ninf_top3"]
    if qninf["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "q95_gap_abs_extrapolated_ninf_upper_top3=" in qninf["fail_trace"] and ">" in qninf["fail_trace"]
    assert "q95_blocker_uniform_top3_obstruction_certificate" in t2w
    qunif = t2w["q95_blocker_uniform_top3_obstruction_certificate"]
    assert qunif["scope"] == "STRICT_TASK2_Q95_BLOCKER_UNIFORM_TOP3_OBSTRUCTION_CERTIFICATE"
    assert qunif["theorem_target"] == "MIN_TOP3_EXTRAPOLATED_NINF_UPPER_MARGIN_GT_ZERO_IMPLIES_OPEN_OBSTRUCTION"
    assert qunif["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    for rr in qunif["computed_rows"]:
        assert rr["gap_abs_extrapolated_ninf_upper"] >= 0.0
        assert isinstance(rr["signed_margin_upper_minus_threshold"], float)
    assert isinstance(qunif["aggregate_metrics"]["all_top3_upper_bounds_above_threshold"], bool)
    if qunif["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "min_signed_margin_upper_minus_threshold=" in qunif["fail_trace"] and "> 0" in qunif["fail_trace"]
    assert "q95_blocker_quad_hp_top3_certificate" in t2w
    qhp = t2w["q95_blocker_quad_hp_top3_certificate"]
    assert qhp["scope"] == "STRICT_TASK2_Q95_BLOCKER_QUAD_HP_TOP3_CERTIFICATE"
    assert qhp["theorem_target"] == "Q95_GAP_ABS_QUAD_HP_TOP3_AND_CROSSCHECK_LE_THRESHOLDS"
    assert qhp["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    assert qhp["domain"]["quad_limit"] == 400
    for rr in qhp["computed_rows"]:
        assert rr["cutsum_quad_high_precision"] >= 0.0
        assert rr["cutsum_quad_high_precision_abs_error_estimate"] >= 0.0
        assert rr["gap_abs_quad_high_precision"] >= 0.0
        assert rr["gap_abs_fixed_quad_n6400"] >= 0.0
        assert rr["cross_integrator_gap_abs_quad_hp_vs_n6400"] >= 0.0
    if qhp["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "=" in qhp["fail_trace"] and ">" in qhp["fail_trace"]
    assert "q95_blocker_quad_hp_error_envelope_certificate" in t2w
    qhpe = t2w["q95_blocker_quad_hp_error_envelope_certificate"]
    assert qhpe["scope"] == "STRICT_TASK2_Q95_BLOCKER_QUAD_HP_ERROR_ENVELOPE_CERTIFICATE"
    assert qhpe["theorem_target"] == "Q95_GAP_ABS_QUAD_HP_UPPER_ENVELOPE_TOP3_LE_THRESHOLD"
    assert qhpe["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    for rr in qhpe["computed_rows"]:
        assert rr["gap_abs_quad_high_precision"] >= 0.0
        assert rr["quad_abs_error_estimate"] >= 0.0
        assert rr["gap_abs_quad_high_precision_upper_envelope"] >= rr["gap_abs_quad_high_precision"]
    if qhpe["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "q95_gap_abs_quad_high_precision_upper_envelope_top3=" in qhpe["fail_trace"] and ">" in qhpe["fail_trace"]
    assert "q95_blocker_tail_budget_certificate" in t2w
    qtb = t2w["q95_blocker_tail_budget_certificate"]
    assert qtb["scope"] == "STRICT_TASK2_Q95_BLOCKER_TAIL_BUDGET_CERTIFICATE"
    assert qtb["theorem_target"] == "Q95_GAP_ABS_UPPER_TAIL_ENVELOPE_TOP3_LE_THRESHOLD"
    assert qtb["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    assert qtb["domain"]["n_levels"] == [1600, 3200, 6400]
    for rr in qtb["computed_rows"]:
        assert rr["gap_abs_n1600"] >= 0.0
        assert rr["gap_abs_n3200"] >= 0.0
        assert rr["gap_abs_n6400"] >= 0.0
        assert rr["delta_n1600_n3200_abs"] >= 0.0
        assert rr["delta_n3200_n6400_abs"] >= 0.0
        assert rr["delta_ratio_32_64_over_16_32"] >= 0.0
        assert rr["tail_budget_beyond_n6400"] >= 0.0
        assert rr["gap_abs_upper_tail_envelope"] >= rr["gap_abs_n6400"]
    if qtb["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "q95_gap_abs_upper_tail_envelope_top3=" in qtb["fail_trace"] and ">" in qtb["fail_trace"]
    assert "q95_blocker_n2400_recompute_certificate" in t2w
    qn96 = t2w["q95_blocker_n2400_recompute_certificate"]
    assert qn96["scope"] == "STRICT_TASK2_Q95_BLOCKER_N2400_RECOMPUTE_CERTIFICATE"
    assert qn96["theorem_target"] == "Q95_GAP_ABS_N2400_TOP3_LE_THRESHOLD_WITH_SMALL_DELTA_VS_N6400"
    assert qn96["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    assert qn96["domain"]["n_levels"] == [2400, 6400]
    for rr in qn96["computed_rows"]:
        assert rr["cutsum_fixed_quad_n2400"] >= 0.0
        assert rr["gap_abs_n6400"] >= 0.0
        assert rr["gap_abs_n2400"] >= 0.0
        assert rr["delta_gap_abs_n2400_minus_n6400_abs"] >= 0.0
        assert isinstance(rr["signed_margin_n2400_minus_threshold"], float)
    if qn96["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "=" in qn96["fail_trace"] and ">" in qn96["fail_trace"]
    assert "q95_blocker_margin" in t2w
    qbm = t2w["q95_blocker_margin"]
    assert qbm["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qbm["scope"] == "STRICT_TASK2_Q95_BLOCKER_MARGIN_WITH_UNCERTAINTY"
    assert qbm["q95_gap_abs_observed"] >= 0.0
    assert qbm["q95_gap_abs_threshold"] > 0.0
    assert isinstance(qbm["q95_margin_observed_minus_threshold"], float)
    assert qbm["uncertainty_abs_from_fixed_quad_convergence"] >= 0.0
    iv = qbm["q95_gap_abs_interval_from_uncertainty"]
    assert 0.0 <= iv["lower"] <= iv["upper"]
    assert qbm["margin_robust_sign"] in {"ABOVE_THRESHOLD_ROBUST", "BELOW_THRESHOLD_ROBUST", "AMBIGUOUS_WITHIN_UNCERTAINTY"}
    assert "q95_blocker_interval_separation_certificate" in t2w
    qisc = t2w["q95_blocker_interval_separation_certificate"]
    assert qisc["scope"] == "STRICT_TASK2_Q95_BLOCKER_INTERVAL_SEPARATION_CERTIFICATE"
    assert qisc["theorem_target"] == "Q95_LOWER_BOUND_EXCEEDS_THRESHOLD_IMPLIES_OPEN_OBSTRUCTION"
    assert qisc["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    assert len(qisc["computed_rows"]) == 1
    qisc_row = qisc["computed_rows"][0]
    assert qisc_row["q95_gap_abs_lower_bound"] <= qisc_row["q95_gap_abs_upper_bound"]
    if qisc["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "q95_gap_abs_lower_bound=" in qisc["fail_trace"] and ">" in qisc["fail_trace"]
    assert "q95_blocker_counterfactual" in t2w
    qbc = t2w["q95_blocker_counterfactual"]
    assert qbc["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qbc["scope"] == "STRICT_TASK2_Q95_BLOCKER_COUNTERFACTUAL_THRESHOLD_PRESSURE"
    assert qbc["threshold_scale_required_for_observed_crossing"] >= 0.0
    assert qbc["threshold_scale_required_for_upper_bound_crossing"] >= 0.0
    assert isinstance(qbc["observed_would_close_under_current_threshold"], bool)
    assert isinstance(qbc["upper_bound_would_close_under_current_threshold"], bool)
    assert qbc["pressure_interpretation"] in {
        "THRESHOLD_MUCH_LOOSER_NEEDED_FOR_CLOSURE",
        "CLOSURE_WITHIN_CURRENT_OR_STRICTER_THRESHOLD",
    }
    assert "q95_blocker_sensitivity" in t2w
    qbs = t2w["q95_blocker_sensitivity"]
    assert qbs["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qbs["scope"] == "STRICT_TASK2_Q95_BLOCKER_LOCAL_S_SENSITIVITY"
    assert qbs["s_perturbation_eps"] > 0.0
    assert len(qbs["rows"]) <= 3
    assert qbs["max_local_slope_abs_dgap_ds"] >= 0.0
    assert qbs["median_local_slope_abs_dgap_ds"] >= 0.0
    for rr in qbs["rows"]:
        assert rr["s_center"] in set(t2w["s_grid"])
        assert rr["s_minus"] <= rr["s_center"] <= rr["s_plus"]
        assert rr["gap_abs_s_minus"] >= 0.0
        assert rr["gap_abs_s_center"] >= 0.0
        assert rr["gap_abs_s_plus"] >= 0.0
        assert rr["local_slope_abs_dgap_ds"] >= 0.0
    assert "q95_blocker_directionality" in t2w
    qbd = t2w["q95_blocker_directionality"]
    assert qbd["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qbd["scope"] == "STRICT_TASK2_Q95_BLOCKER_LOCAL_DIRECTIONALITY"
    assert 0.0 <= qbd["upward_step_reduces_gap_frequency"] <= 1.0
    assert 0.0 <= qbd["downward_step_reduces_gap_frequency"] <= 1.0
    assert len(qbd["rows"]) <= 3
    for rr in qbd["rows"]:
        assert rr["s_center"] in set(t2w["s_grid"])
        assert isinstance(rr["delta_plus"], float)
        assert isinstance(rr["delta_minus"], float)
        assert isinstance(rr["upward_step_reduces_gap"], bool)
        assert isinstance(rr["downward_step_reduces_gap"], bool)
    assert "q95_blocker_actionability" in t2w
    qba = t2w["q95_blocker_actionability"]
    assert qba["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qba["scope"] == "STRICT_TASK2_Q95_BLOCKER_LOCAL_ACTIONABILITY"
    assert len(qba["rows"]) <= 3
    if qba["best_action_s"] is not None:
        assert qba["best_action_s"] in set(t2w["s_grid"])
    assert qba["best_action_move"] in {"INCREASE_S", "DECREASE_S", "NO_LOCAL_RELIEF"}
    assert qba["best_action_estimated_local_gap_reduction"] >= 0.0
    prev_est = None
    for rr in qba["rows"]:
        assert rr["s_center"] in set(t2w["s_grid"])
        assert rr["recommended_move"] in {"INCREASE_S", "DECREASE_S", "NO_LOCAL_RELIEF"}
        assert rr["estimated_local_gap_reduction"] >= 0.0
        if prev_est is not None:
            assert prev_est + 1e-18 >= rr["estimated_local_gap_reduction"]
        prev_est = rr["estimated_local_gap_reduction"]
    assert "q95_blocker_action_execution" in t2w
    qex = t2w["q95_blocker_action_execution"]
    assert qex["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qex["scope"] == "STRICT_TASK2_Q95_BLOCKER_ONE_STEP_ACTION_EXECUTION"
    assert qex["move"] in {"INCREASE_S", "DECREASE_S", "NO_LOCAL_RELIEF"}
    assert isinstance(qex["improves_locally"], bool)
    if qex["move"] == "NO_LOCAL_RELIEF":
        assert qex["s_before"] is None
        assert qex["s_after"] is None
        assert qex["gap_abs_before"] is None
        assert qex["gap_abs_after"] is None
        assert qex["delta_gap_abs_after_minus_before"] is None
    else:
        assert qex["s_before"] in set(t2w["s_grid"])
        assert 0.0 <= qex["s_after"] <= 3.5
        assert qex["gap_abs_before"] >= 0.0
        assert qex["gap_abs_after"] >= 0.0
        assert isinstance(qex["delta_gap_abs_after_minus_before"], float)
    assert "q95_action_effect_crosscheck" in t2w
    qec = t2w["q95_action_effect_crosscheck"]
    assert qec["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qec["scope"] == "STRICT_TASK2_Q95_ACTION_EFFECT_FIXED_QUAD_CROSSCHECK"
    assert isinstance(qec["effect_sign_consistent_across_orders"], bool)
    assert isinstance(qec["both_orders_improve"], bool)
    if qex["move"] == "NO_LOCAL_RELIEF":
        assert qec["delta_gap_abs_n400_after_minus_before"] is None
        assert qec["delta_gap_abs_n800_after_minus_before"] is None
    else:
        assert isinstance(qec["delta_gap_abs_n400_after_minus_before"], float)
        assert isinstance(qec["delta_gap_abs_n800_after_minus_before"], float)
    assert "q95_action_effect_bootstrap" in t2w
    qeb = t2w["q95_action_effect_bootstrap"]
    assert qeb["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qeb["scope"] == "STRICT_TASK2_Q95_ACTION_EFFECT_BOOTSTRAP_ROBUSTNESS"
    assert qeb["n_range_inclusive"] == [300, 900]
    assert 0.0 <= qeb["p_improve"] <= 1.0
    assert 0.0 <= qeb["p_improve_ci95_jeffreys"]["lower"] <= qeb["p_improve_ci95_jeffreys"]["upper"] <= 1.0
    assert len(qeb["delta_q05_q50_q95"]) == 3
    if qex["move"] == "NO_LOCAL_RELIEF":
        assert qeb["bootstrap_size"] == 0
    else:
        assert qeb["bootstrap_size"] == 128
    assert "q95_action_gate" in t2w
    qag = t2w["q95_action_gate"]
    assert qag["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qag["scope"] == "STRICT_TASK2_Q95_ACTION_BOOTSTRAP_DECISION_GATE"
    assert qag["p_improve_lb95_threshold"] == 0.60
    assert isinstance(qag["criterion_p_improve_lb95_ge_threshold"], bool)
    assert isinstance(qag["criterion_fixed_quad_sign_consistent"], bool)
    assert isinstance(qag["go_for_next_local_step"], bool)
    assert qag["reason"] in {"GO", "HOLD_AND_RECALIBRATE"}
    assert "q95_action_gate_consistency" in t2w
    qgc = t2w["q95_action_gate_consistency"]
    assert qgc["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qgc["scope"] == "STRICT_TASK2_Q95_ACTION_GATE_CONSISTENCY"
    assert set(qgc["checks"].keys()) == {"if_go_then_action_exists", "if_go_then_bootstrap_executed", "if_go_then_crosscheck_not_none"}
    assert all(isinstance(v, bool) for v in qgc["checks"].values())
    assert qgc["all_checks_pass"] is True
    assert "q95_blocker_step_efficiency" in t2w
    qse = t2w["q95_blocker_step_efficiency"]
    assert qse["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qse["scope"] == "STRICT_TASK2_Q95_BLOCKER_ONE_STEP_EFFICIENCY"
    assert qse["step_abs"] >= 0.0
    assert qse["gap_abs_before"] >= 0.0
    assert qse["gap_abs_after"] >= 0.0
    assert qse["gap_abs_reduction"] >= 0.0
    assert qse["reduction_per_unit_s_step"] >= 0.0
    assert 0.0 <= qse["relative_reduction_fraction"] <= 1.0
    assert "q95_after_one_step_local_margin" in t2w
    qlm = t2w["q95_after_one_step_local_margin"]
    assert qlm["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qlm["scope"] == "STRICT_TASK2_Q95_AFTER_ONE_STEP_LOCAL_MARGIN"
    assert qlm["q95_gap_abs_threshold"] > 0.0
    assert isinstance(qlm["moves_toward_closure"], bool)
    if qex["move"] == "NO_LOCAL_RELIEF":
        assert qlm["gap_abs_before"] is None
        assert qlm["gap_abs_after"] is None
        assert qlm["margin_before"] is None
        assert qlm["margin_after"] is None
        assert qlm["margin_delta_after_minus_before"] is None
    else:
        assert qlm["gap_abs_before"] >= 0.0
        assert qlm["gap_abs_after"] >= 0.0
        assert isinstance(qlm["margin_before"], float)
        assert isinstance(qlm["margin_after"], float)
        assert isinstance(qlm["margin_delta_after_minus_before"], float)
    assert "q95_after_one_step_progress_score" in t2w
    qps = t2w["q95_after_one_step_progress_score"]
    assert qps["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qps["scope"] == "STRICT_TASK2_Q95_AFTER_ONE_STEP_PROGRESS_SCORE"
    assert qps["margin_before_abs"] >= 0.0
    assert qps["margin_after_abs"] >= 0.0
    assert qps["absolute_margin_improvement"] >= 0.0
    assert 0.0 <= qps["normalized_progress_score_0_1"] <= 1.0
    assert "q95_blocker_choice_panel" in t2w
    qcp = t2w["q95_blocker_choice_panel"]
    assert qcp["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qcp["scope"] == "STRICT_TASK2_BLOCKER_CHOICE_NORMALIZED_OVERSHOOT"
    assert len(qcp["rows"]) == 18
    for rr in qcp["rows"]:
        assert rr["criterion"] in {"q95_gap_abs", "q95_refined_window_gap_abs", "q95_quad_hp_top3_gap_abs", "q95_quad_hp_upper_envelope_gap_abs", "q95_tail_budget_upper_envelope_gap_abs", "q95_n2400_gap_abs", "q95_n2400_vs_n6400_delta_abs", "q95_n3200_vs_n6400_delta_abs", "q95_n1600_vs_n6400_delta_abs", "q95_n800_vs_n6400_delta_abs", "q95_n12800_vs_n6400_delta_abs", "q95_n12800_tail_ratio", "q95_n25600_vs_n12800_delta_abs", "q95_n25600_monotone_violation", "max_gap_rel", "weight_sign_nonnegativity", "q95_cross_integrator_gap", "q95_convergence_delta_n400_to_n800_abs"}
        assert rr["normalized_overshoot"] >= 0.0
        assert isinstance(rr["is_satisfied"], bool)
    assert qcp["easiest_unresolved_blocker"] in {"q95_gap_abs", "max_gap_rel", "weight_sign_nonnegativity", "q95_cross_integrator_gap", "q95_convergence_delta_n400_to_n800_abs", "none", "q95_refined_window_gap_abs", "q95_quad_hp_top3_gap_abs", "q95_quad_hp_upper_envelope_gap_abs", "q95_tail_budget_upper_envelope_gap_abs", "q95_n2400_gap_abs", "q95_n2400_vs_n6400_delta_abs", "q95_n3200_vs_n6400_delta_abs", "q95_n1600_vs_n6400_delta_abs", "q95_n800_vs_n6400_delta_abs", "q95_n12800_vs_n6400_delta_abs", "q95_n12800_tail_ratio", "q95_n25600_vs_n12800_delta_abs", "q95_n25600_monotone_violation"}
    assert qcp["dominant_unresolved_blocker"] in {"q95_gap_abs", "max_gap_rel", "weight_sign_nonnegativity", "q95_cross_integrator_gap", "q95_convergence_delta_n400_to_n800_abs", "none", "q95_refined_window_gap_abs", "q95_quad_hp_top3_gap_abs", "q95_quad_hp_upper_envelope_gap_abs", "q95_tail_budget_upper_envelope_gap_abs", "q95_n2400_gap_abs", "q95_n2400_vs_n6400_delta_abs", "q95_n3200_vs_n6400_delta_abs", "q95_n1600_vs_n6400_delta_abs", "q95_n800_vs_n6400_delta_abs", "q95_n12800_vs_n6400_delta_abs", "q95_n12800_tail_ratio", "q95_n25600_vs_n12800_delta_abs", "q95_n25600_monotone_violation"}
    assert "q95_blocker_choice_consistency" in t2w
    qcc = t2w["q95_blocker_choice_consistency"]
    assert qcc["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert qcc["scope"] == "STRICT_TASK2_BLOCKER_CHOICE_CONSISTENCY"
    assert qcc["dominant_blocker"] in {"q95_gap_abs", "max_gap_rel", "weight_sign_nonnegativity", "q95_cross_integrator_gap", "q95_convergence_delta_n400_to_n800_abs", "none", "pending", "q95_refined_window_gap_abs", "q95_quad_hp_top3_gap_abs", "q95_quad_hp_upper_envelope_gap_abs", "q95_tail_budget_upper_envelope_gap_abs", "q95_n2400_gap_abs", "q95_n2400_vs_n6400_delta_abs", "q95_n3200_vs_n6400_delta_abs", "q95_n1600_vs_n6400_delta_abs", "q95_n800_vs_n6400_delta_abs", "q95_n12800_vs_n6400_delta_abs", "q95_n12800_tail_ratio", "q95_n25600_vs_n12800_delta_abs", "q95_n25600_monotone_violation"}
    assert qcc["easiest_unresolved_blocker"] in {"q95_gap_abs", "max_gap_rel", "weight_sign_nonnegativity", "q95_cross_integrator_gap", "q95_convergence_delta_n400_to_n800_abs", "none", "q95_refined_window_gap_abs", "q95_quad_hp_top3_gap_abs", "q95_quad_hp_upper_envelope_gap_abs", "q95_tail_budget_upper_envelope_gap_abs", "q95_n2400_gap_abs", "q95_n2400_vs_n6400_delta_abs", "q95_n3200_vs_n6400_delta_abs", "q95_n1600_vs_n6400_delta_abs", "q95_n800_vs_n6400_delta_abs", "q95_n12800_vs_n6400_delta_abs", "q95_n12800_tail_ratio", "q95_n25600_vs_n12800_delta_abs", "q95_n25600_monotone_violation"}
    assert isinstance(qcc["is_consistent_when_q95_dominates"], bool)
    if qcc["dominant_blocker"] != "pending":
        assert qcc["dominant_blocker"] == qcp["dominant_unresolved_blocker"]
    assert "dominant_blocker_numeric_margin" in t2w
    dbm = t2w["dominant_blocker_numeric_margin"]
    assert dbm["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert dbm["scope"] == "STRICT_TASK2_DOMINANT_BLOCKER_NUMERIC_MARGIN"
    assert dbm["dominant_blocker"] in {"q95_gap_abs", "max_gap_rel", "weight_sign_nonnegativity", "q95_cross_integrator_gap", "q95_convergence_delta_n400_to_n800_abs", "none", "q95_refined_window_gap_abs", "q95_quad_hp_top3_gap_abs", "q95_quad_hp_upper_envelope_gap_abs", "q95_tail_budget_upper_envelope_gap_abs", "q95_n2400_gap_abs", "q95_n2400_vs_n6400_delta_abs", "q95_n3200_vs_n6400_delta_abs", "q95_n1600_vs_n6400_delta_abs", "q95_n800_vs_n6400_delta_abs", "q95_n12800_vs_n6400_delta_abs", "q95_n12800_tail_ratio", "q95_n25600_vs_n12800_delta_abs", "q95_n25600_monotone_violation"}
    assert isinstance(dbm["observed_value"], float)
    assert isinstance(dbm["threshold_value"], float)
    assert isinstance(dbm["signed_margin_observed_minus_threshold"], float)
    assert dbm["normalized_overshoot_vs_threshold"] >= 0.0
    assert "dominant_blocker_selection_consistency" in t2w
    dsc = t2w["dominant_blocker_selection_consistency"]
    assert dsc["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert dsc["scope"] == "STRICT_TASK2_DOMINANT_BLOCKER_SELECTION_CONSISTENCY"
    assert dsc["dominant_blocker"] in {"q95_gap_abs", "max_gap_rel", "weight_sign_nonnegativity", "q95_cross_integrator_gap", "q95_convergence_delta_n400_to_n800_abs", "none", "q95_refined_window_gap_abs", "q95_quad_hp_top3_gap_abs", "q95_quad_hp_upper_envelope_gap_abs", "q95_tail_budget_upper_envelope_gap_abs", "q95_n2400_gap_abs", "q95_n2400_vs_n6400_delta_abs", "q95_n3200_vs_n6400_delta_abs", "q95_n1600_vs_n6400_delta_abs", "q95_n800_vs_n6400_delta_abs", "q95_n12800_vs_n6400_delta_abs", "q95_n12800_tail_ratio", "q95_n25600_vs_n12800_delta_abs", "q95_n25600_monotone_violation"}
    assert dsc["dominant_unresolved_expected"] in {"q95_gap_abs", "max_gap_rel", "weight_sign_nonnegativity", "q95_cross_integrator_gap", "q95_convergence_delta_n400_to_n800_abs", "none", "q95_refined_window_gap_abs", "q95_quad_hp_top3_gap_abs", "q95_quad_hp_upper_envelope_gap_abs", "q95_tail_budget_upper_envelope_gap_abs", "q95_n2400_gap_abs", "q95_n2400_vs_n6400_delta_abs", "q95_n3200_vs_n6400_delta_abs", "q95_n1600_vs_n6400_delta_abs", "q95_n800_vs_n6400_delta_abs", "q95_n12800_vs_n6400_delta_abs", "q95_n12800_tail_ratio", "q95_n25600_vs_n12800_delta_abs", "q95_n25600_monotone_violation"}
    assert dsc["is_argmax_overshoot"] is True
    assert "dominant_blocker_robustness_certificate" in t2w
    drc = t2w["dominant_blocker_robustness_certificate"]
    assert drc["scope"] == "STRICT_TASK2_DOMINANT_BLOCKER_ROBUSTNESS_CERTIFICATE"
    assert drc["theorem_target"] == "DOMINANT_BLOCKER_NORMALIZED_OVERSHOOT_GAP_GT_ZERO"
    assert drc["verdict"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    assert len(drc["computed_rows"]) == 1
    rr = drc["computed_rows"][0]
    assert rr["dominant_normalized_overshoot"] >= 0.0
    assert rr["second_largest_normalized_overshoot"] >= 0.0
    assert isinstance(rr["dominance_gap"], float)
    if drc["verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE":
        assert "dominance_gap=" in drc["fail_trace"] and "<= 0" in drc["fail_trace"]
    assert "criterion_coherence_sign" in t2w
    ccs = t2w["criterion_coherence_sign"]
    assert ccs["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert ccs["scope"] == "STRICT_TASK2_SIGN_CRITERION_COHERENCE"
    assert isinstance(ccs["min_effective_weight_global"], float)
    assert isinstance(ccs["min_effective_weight_global_min"], float)
    assert isinstance(ccs["all_nonnegative_weights_flag"], bool)
    assert ccs["flag_equals_numeric_inequality"] is True
    assert "criterion_coherence_convergence" in t2w
    ccc = t2w["criterion_coherence_convergence"]
    assert ccc["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert ccc["scope"] == "STRICT_TASK2_CONVERGENCE_CRITERION_COHERENCE"
    assert ccc["q95_convergence_delta_n400_to_n800_abs"] >= 0.0
    assert ccc["q95_convergence_delta_n400_to_n800_abs_max"] >= 0.0
    assert isinstance(ccc["q95_convergence_delta_le_threshold_flag"], bool)
    assert ccc["flag_equals_numeric_inequality"] is True
    assert "criterion_coherence_cross_integrator" in t2w
    cci = t2w["criterion_coherence_cross_integrator"]
    assert cci["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert cci["scope"] == "STRICT_TASK2_CROSS_INTEGRATOR_CRITERION_COHERENCE"
    assert cci["q95_cross_integrator_gap_abs"] >= 0.0
    assert cci["q95_cross_integrator_gap_abs_max"] >= 0.0
    assert isinstance(cci["q95_cross_integrator_gap_le_threshold_flag"], bool)
    assert cci["flag_equals_numeric_inequality"] is True
    assert "criterion_coherence_q95_gap" in t2w
    ccq = t2w["criterion_coherence_q95_gap"]
    assert ccq["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert ccq["scope"] == "STRICT_TASK2_Q95_GAP_CRITERION_COHERENCE"
    assert ccq["q95_gap_abs"] >= 0.0
    assert ccq["q95_gap_abs_max"] >= 0.0
    assert isinstance(ccq["q95_gap_abs_le_threshold_flag"], bool)
    assert ccq["flag_equals_numeric_inequality"] is True
    assert "criterion_coherence_max_gap_rel" in t2w
    ccm = t2w["criterion_coherence_max_gap_rel"]
    assert ccm["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert ccm["scope"] == "STRICT_TASK2_MAX_GAP_REL_CRITERION_COHERENCE"
    assert ccm["max_gap_rel"] >= 0.0
    assert ccm["max_gap_rel_max"] >= 0.0
    assert isinstance(ccm["max_gap_rel_le_threshold_flag"], bool)
    assert ccm["flag_equals_numeric_inequality"] is True
    assert "criterion_coherence_global_closure" in t2w
    ccg = t2w["criterion_coherence_global_closure"]
    assert ccg["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert ccg["scope"] == "STRICT_TASK2_GLOBAL_CLOSURE_CRITERION_COHERENCE"
    assert isinstance(ccg["all_criteria_satisfied_flag"], bool)
    assert isinstance(ccg["numeric_conjunction_recomputed"], bool)
    assert ccg["flag_equals_numeric_conjunction"] is True
    assert "criterion_coherence_weight_sign" in t2w
    cws = t2w["criterion_coherence_weight_sign"]
    assert cws["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert cws["scope"] == "STRICT_TASK2_WEIGHT_SIGN_CRITERION_COHERENCE"
    assert isinstance(cws["min_effective_weight_global"], float)
    assert isinstance(cws["min_effective_weight_global_min"], float)
    assert isinstance(cws["weight_sign_nonnegative_flag"], bool)
    assert cws["flag_equals_numeric_inequality"] is True
    assert "criterion_coherence_verdict_gate" in t2w
    cvg = t2w["criterion_coherence_verdict_gate"]
    assert cvg["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert cvg["scope"] == "STRICT_TASK2_VERDICT_GATE_COHERENCE"
    assert isinstance(cvg["all_criteria_satisfied_flag"], bool)
    assert cvg["verdict_task2"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE"}
    assert cvg["closed_label"] == "CLOSED_NUMERICAL_WITNESS_TASK2"
    assert cvg["open_label"] == "OPEN_OBSTRUCTION_WITH_TRACE"
    assert cvg["flag_matches_verdict_label"] is True
    assert "criterion_coherence_fail_trace_numeric" in t2w
    cfn = t2w["criterion_coherence_fail_trace_numeric"]
    assert cfn["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert cfn["scope"] == "STRICT_TASK2_FAIL_TRACE_NUMERIC_COHERENCE"
    assert cfn["verdict_task2"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE", "PENDING"}
    assert cfn["dominant_blocker"] in {"q95_gap_abs", "max_gap_rel", "q95_cross_integrator_gap", "q95_convergence_delta_n400_to_n800_abs", "weight_sign_nonnegativity", "none", "pending", "q95_refined_window_gap_abs", "q95_tail_budget_upper_envelope_gap_abs", "q95_n2400_gap_abs", "q95_n2400_vs_n6400_delta_abs", "q95_n3200_vs_n6400_delta_abs", "q95_n1600_vs_n6400_delta_abs", "q95_n800_vs_n6400_delta_abs", "q95_n12800_vs_n6400_delta_abs", "q95_n12800_tail_ratio", "q95_n25600_vs_n12800_delta_abs", "q95_n25600_monotone_violation"}
    assert isinstance(cfn["fail_trace"], str)
    assert cfn["trace_prefix_matches_dominant_blocker"] is True
    assert "criterion_coherence_dominant_margin_sign" in t2w
    cds = t2w["criterion_coherence_dominant_margin_sign"]
    assert cds["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert cds["scope"] == "STRICT_TASK2_DOMINANT_MARGIN_SIGN_COHERENCE"
    assert cds["verdict_task2"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE", "PENDING"}
    assert cds["dominant_blocker"] in {"q95_gap_abs", "max_gap_rel", "q95_cross_integrator_gap", "q95_convergence_delta_n400_to_n800_abs", "weight_sign_nonnegativity", "none", "pending", "q95_refined_window_gap_abs", "q95_tail_budget_upper_envelope_gap_abs", "q95_n2400_gap_abs", "q95_n2400_vs_n6400_delta_abs", "q95_n3200_vs_n6400_delta_abs", "q95_n1600_vs_n6400_delta_abs", "q95_n800_vs_n6400_delta_abs", "q95_n12800_vs_n6400_delta_abs", "q95_n12800_tail_ratio", "q95_n25600_vs_n12800_delta_abs", "q95_n25600_monotone_violation"}
    assert isinstance(cds["signed_margin_observed_minus_threshold"], float)
    assert cds["open_requires_positive_margin"] is True
    assert "criterion_coherence_open_trace_inequality" in t2w
    cot = t2w["criterion_coherence_open_trace_inequality"]
    assert cot["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert cot["scope"] == "STRICT_TASK2_OPEN_TRACE_INEQUALITY_COHERENCE"
    assert cot["verdict_task2"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE", "PENDING"}
    assert isinstance(cot["fail_trace"], str)
    assert cot["open_trace_contains_numeric_inequality_token"] is True
    assert "criterion_coherence_dominant_inequality_prefix" in t2w
    cdp = t2w["criterion_coherence_dominant_inequality_prefix"]
    assert cdp["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert cdp["scope"] == "STRICT_TASK2_DOMINANT_INEQUALITY_PREFIX_COHERENCE"
    assert cdp["dominant_blocker"] in {"q95_gap_abs", "max_gap_rel", "q95_cross_integrator_gap", "q95_convergence_delta_n400_to_n800_abs", "weight_sign_nonnegativity", "none", "pending", "q95_refined_window_gap_abs", "q95_tail_budget_upper_envelope_gap_abs", "q95_n2400_gap_abs", "q95_n2400_vs_n6400_delta_abs", "q95_n3200_vs_n6400_delta_abs", "q95_n1600_vs_n6400_delta_abs", "q95_n800_vs_n6400_delta_abs", "q95_n12800_vs_n6400_delta_abs", "q95_n12800_tail_ratio", "q95_n25600_vs_n12800_delta_abs", "q95_n25600_monotone_violation"}
    assert isinstance(cdp["dominant_inequality"], str)
    assert cdp["prefix_matches_dominant_blocker"] is True
    assert "criterion_coherence_fail_trace_equals_dominant" in t2w
    cfe = t2w["criterion_coherence_fail_trace_equals_dominant"]
    assert cfe["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert cfe["scope"] == "STRICT_TASK2_FAIL_TRACE_EQUALS_DOMINANT_COHERENCE"
    assert cfe["verdict_task2"] in {"CLOSED_NUMERICAL_WITNESS_TASK2", "OPEN_OBSTRUCTION_WITH_TRACE", "PENDING"}
    assert isinstance(cfe["fail_trace"], str)
    assert isinstance(cfe["dominant_inequality"], str)
    assert cfe["open_requires_exact_equality"] is True
    assert cc["verdict_matches_criteria"] is True

    assert "ur_numerical_stress_policy_weighted_boundary_risk_posterior_predictive_precursor" in data
    unsp_pp = data["ur_numerical_stress_policy_weighted_boundary_risk_posterior_predictive_precursor"]
    assert unsp_pp["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_pp["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_WEIGHTED_BOUNDARY_RISK_POSTERIOR_PREDICTIVE"
    assert len(unsp_pp["rows"]) == 9
    assert 0.0 <= unsp_pp["global_p_stays_optimal_mean"] <= 1.0
    for rr in unsp_pp["rows"]:
        assert rr["cap_scale"] in {0.85, 1.0, 1.15}
        assert rr["jitter"] in {0.0, 0.01, 0.02}
        assert rr["posterior_samples"] == 256
        assert 0 <= rr["valid_samples"] <= rr["posterior_samples"]
        assert 0 <= rr["stay_count"] <= rr["valid_samples"]
        assert 0.0 <= rr["p_winner_set_stays_optimal"] <= 1.0
        assert 0.0 <= rr["p_winner_set_stays_optimal_ci95_jeffreys"]["lower"] <= rr["p_winner_set_stays_optimal_ci95_jeffreys"]["upper"] <= 1.0
    assert "ur_numerical_stress_policy_posterior_predictive_decision_gate_precursor" in data
    unsp_pg = data["ur_numerical_stress_policy_posterior_predictive_decision_gate_precursor"]
    assert unsp_pg["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_pg["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_POSTERIOR_PREDICTIVE_DECISION_GATE"
    assert unsp_pg["rule"] == "GO_if_all_critical_cells_lb95_ge_threshold_else_HOLD"
    assert 0.0 <= unsp_pg["lb95_threshold"] <= 1.0
    assert len(unsp_pg["rows"]) == 9
    assert unsp_pg["decision"] in {"GO", "HOLD_AND_RECALIBRATE"}
    assert isinstance(unsp_pg["ready_for_next_costlier_policy_step"], bool)
    for rr in unsp_pg["rows"]:
        assert rr["cap_scale"] in {0.85, 1.0, 1.15}
        assert rr["jitter"] in {0.0, 0.01, 0.02}
        assert 0.0 <= rr["p_stay"] <= 1.0
        assert 0.0 <= rr["p_stay_lb95"] <= 1.0
        assert isinstance(rr["criterion_lb95_ge_threshold"], bool)
    assert "ur_numerical_stress_policy_posterior_predictive_gate_calibration_precursor" in data
    unsp_gc = data["ur_numerical_stress_policy_posterior_predictive_gate_calibration_precursor"]
    assert unsp_gc["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_gc["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_POSTERIOR_PREDICTIVE_GATE_CALIBRATION"
    assert unsp_gc["threshold_grid"] == [0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
    assert len(unsp_gc["rows"]) == 6
    assert 0.0 <= unsp_gc["go_rate_over_threshold_grid"] <= 1.0
    assert 0.70 <= unsp_gc["recommended_threshold_max_go_with_min_pass_rate_0p90"] <= 0.95
    for rr in unsp_gc["rows"]:
        assert rr["lb95_threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
        assert 0 <= rr["critical_cells_pass_count"] <= rr["critical_cells_total"] == 9
        assert 0.0 <= rr["critical_cells_pass_rate"] <= 1.0
        assert isinstance(rr["decision_go"], bool)
    assert "ur_numerical_stress_policy_posterior_predictive_gate_cost_calibration_precursor" in data
    unsp_gcc = data["ur_numerical_stress_policy_posterior_predictive_gate_cost_calibration_precursor"]
    assert unsp_gcc["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_gcc["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_POSTERIOR_PREDICTIVE_GATE_COST_CALIBRATION"
    assert "go_rate_local" in unsp_gcc["utility_definition"]
    assert unsp_gcc["base_runtime_reference_seconds"] > 0.0
    assert unsp_gcc["lambda_grid"] == [0.0, 0.1, 0.2, 0.4]
    assert len(unsp_gcc["rows"]) == 24
    assert len(unsp_gcc["selected_rows"]) == 4
    for rr in unsp_gcc["rows"]:
        assert rr["lambda_cost"] in {0.0, 0.1, 0.2, 0.4}
        assert rr["lb95_threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
        assert rr["go_rate_local"] in {0.0, 1.0}
        assert 0.0 <= rr["runtime_uplift_proxy"] <= 1.0
    for rr in unsp_gcc["selected_rows"]:
        assert rr["lambda_cost"] in {0.0, 0.1, 0.2, 0.4}
        assert rr["lb95_threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
    assert "ur_numerical_stress_policy_gate_frontier_utility_precursor" in data
    unsp_fu = data["ur_numerical_stress_policy_gate_frontier_utility_precursor"]
    assert unsp_fu["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_fu["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_GATE_FRONTIER_UTILITY"
    assert unsp_fu["axes"] == ["false_hold_risk", "false_go_risk", "runtime_uplift_proxy"]
    assert len(unsp_fu["rows"]) == 6
    assert 1 <= unsp_fu["pareto_frontier_count"] <= 6
    assert len(unsp_fu["pareto_frontier_thresholds"]) == unsp_fu["pareto_frontier_count"]
    for rr in unsp_fu["rows"]:
        assert rr["lb95_threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
        assert 0.0 <= rr["false_hold_risk"] <= 1.0
        assert 0.0 <= rr["false_go_risk"] <= 1.0
        assert 0.0 <= rr["runtime_uplift_proxy"] <= 1.0
        assert isinstance(rr["pareto_frontier"], bool)
    assert "ur_numerical_stress_policy_gate_frontier_knee_precursor" in data
    unsp_fk = data["ur_numerical_stress_policy_gate_frontier_knee_precursor"]
    assert unsp_fk["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_fk["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_GATE_FRONTIER_KNEE"
    assert unsp_fk["ideal_point"] == {"false_hold_risk": 0.0, "false_go_risk": 0.0, "runtime_uplift_proxy": 0.0}
    assert unsp_fk["recommended_knee_threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
    for rr in unsp_fk["rows"]:
        assert rr["lb95_threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
        assert isinstance(rr["pareto_frontier"], bool) and rr["pareto_frontier"] is True
        assert rr["ideal_point_distance_l2"] >= 0.0
    assert "ur_numerical_stress_policy_gate_frontier_knee_stability_precursor" in data
    unsp_fks = data["ur_numerical_stress_policy_gate_frontier_knee_stability_precursor"]
    assert unsp_fks["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_fks["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_GATE_FRONTIER_KNEE_STABILITY"
    assert unsp_fks["bootstrap_size"] == 512
    assert len(unsp_fks["rows"]) >= 1
    assert unsp_fks["most_stable_knee_threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
    for rr in unsp_fks["rows"]:
        assert rr["lb95_threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
        assert 0 <= rr["knee_selection_count"] <= 512
        assert 0.0 <= rr["knee_selection_frequency"] <= 1.0
        assert 0.0 <= rr["knee_selection_frequency_ci95_jeffreys"]["lower"] <= rr["knee_selection_frequency_ci95_jeffreys"]["upper"] <= 1.0
    assert "ur_numerical_stress_policy_gate_frontier_knee_cross_seed_stability_precursor" in data
    unsp_fkcs = data["ur_numerical_stress_policy_gate_frontier_knee_cross_seed_stability_precursor"]
    assert unsp_fkcs["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_fkcs["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_GATE_FRONTIER_KNEE_CROSS_SEED_STABILITY"
    assert unsp_fkcs["seed_grid"] == [975166, 975167, 975168, 975169]
    assert len(unsp_fkcs["rows"]) == 4
    assert len(unsp_fkcs["span_rows"]) == 6
    assert 0.0 <= unsp_fkcs["max_span_over_thresholds"] <= 1.0
    for sr in unsp_fkcs["rows"]:
        assert sr["seed"] in {975166, 975167, 975168, 975169}
        assert sr["most_stable_knee_threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
        assert len(sr["rows"]) == 6
        for rr in sr["rows"]:
            assert rr["lb95_threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
            assert 0.0 <= rr["knee_selection_frequency"] <= 1.0
    for rr in unsp_fkcs["span_rows"]:
        assert rr["lb95_threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
        assert 0.0 <= rr["knee_selection_frequency_span_over_seeds"] <= 1.0
        assert 0.0 <= rr["knee_selection_frequency_mean_over_seeds"] <= 1.0
    assert "ur_numerical_stress_policy_gate_frontier_knee_cross_seed_consensus_precursor" in data
    unsp_fkcc = data["ur_numerical_stress_policy_gate_frontier_knee_cross_seed_consensus_precursor"]
    assert unsp_fkcc["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_fkcc["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_GATE_FRONTIER_KNEE_CROSS_SEED_CONSENSUS"
    assert unsp_fkcc["seed_grid"] == [975166, 975167, 975168, 975169]
    assert unsp_fkcc["consensus_min_frequency_for_go"] == 0.75
    assert len(unsp_fkcc["rows"]) == 6
    assert 0.0 <= unsp_fkcc["consensus_strength"] <= 1.0
    assert unsp_fkcc["consensus_knee_threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
    assert isinstance(unsp_fkcc["consensus_go"], bool)
    for rr in unsp_fkcc["rows"]:
        assert rr["lb95_threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
        assert 0 <= rr["most_stable_vote_count"] <= 4
        assert 0.0 <= rr["most_stable_vote_frequency"] <= 1.0
    assert "ur_numerical_stress_policy_gate_frontier_knee_consensus_threshold_sweep_precursor" in data
    unsp_fkcts = data["ur_numerical_stress_policy_gate_frontier_knee_consensus_threshold_sweep_precursor"]
    assert unsp_fkcts["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_fkcts["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_GATE_FRONTIER_KNEE_CONSENSUS_THRESHOLD_SWEEP"
    assert unsp_fkcts["seed_grid"] == [975166, 975167, 975168, 975169]
    assert 0.0 <= unsp_fkcts["consensus_strength"] <= 1.0
    assert len(unsp_fkcts["rows"]) == 4
    assert isinstance(unsp_fkcts["go_is_robust_across_threshold_grid"], bool)
    assert isinstance(unsp_fkcts["hold_is_robust_across_threshold_grid"], bool)
    for rr in unsp_fkcts["rows"]:
        assert rr["consensus_min_frequency_for_go"] in {0.50, 0.60, 0.75, 0.90}
        assert isinstance(rr["consensus_go"], bool)
        assert -1.0 <= rr["consensus_margin"] <= 1.0
    assert "ur_numerical_stress_policy_gate_frontier_knee_weighted_cross_seed_consensus_precursor" in data
    unsp_fkw = data["ur_numerical_stress_policy_gate_frontier_knee_weighted_cross_seed_consensus_precursor"]
    assert unsp_fkw["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsp_fkw["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_POLICY_GATE_FRONTIER_KNEE_WEIGHTED_CROSS_SEED_CONSENSUS"
    assert unsp_fkw["seed_grid"] == [975166, 975167, 975168, 975169]
    assert unsp_fkw["weight_rule"] == "stability_weight = 1/(1+local_frequency_span)"
    assert unsp_fkw["consensus_min_frequency_for_go"] == 0.75
    assert len(unsp_fkw["seed_weight_rows"]) == 4
    assert len(unsp_fkw["rows"]) == 6
    assert 0.0 <= unsp_fkw["weighted_consensus_strength"] <= 1.0
    assert unsp_fkw["weighted_consensus_knee_threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
    assert isinstance(unsp_fkw["weighted_consensus_go"], bool)
    assert isinstance(unsp_fkw["unweighted_weighted_threshold_agreement"], bool)
    assert unsp_fkw["weighted_consensus_bootstrap_size"] == 512
    assert 0.0 <= unsp_fkw["weighted_consensus_strength_ci95_bootstrap"]["lower"] <= unsp_fkw["weighted_consensus_strength_ci95_bootstrap"]["upper"] <= 1.0
    assert isinstance(unsp_fkw["weighted_consensus_go_ci95_lb"], bool)
    for rr in unsp_fkw["seed_weight_rows"]:
        assert rr["seed"] in {975166, 975167, 975168, 975169}
        assert 0.0 <= rr["local_frequency_span"] <= 1.0
        assert 0.0 < rr["stability_weight"] <= 1.0
        assert rr["most_stable_knee_threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
    for rr in unsp_fkw["rows"]:
        assert rr["lb95_threshold"] in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95}
        assert rr["weighted_vote_score"] >= 0.0
        assert 0.0 <= rr["weighted_vote_frequency"] <= 1.0
    assert "ur_numerical_stress_alt_replay_trend_precursor" in data
    unsart = data["ur_numerical_stress_alt_replay_trend_precursor"]
    assert unsart["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsart["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_ALT_REPLAY_TREND"
    assert unsart["selection_rule"] == "per_improved_class_min_alt_warning_then_min_abs_delta_shift"
    assert unsart["num_improved_classes_replayed"] == len(unsart["rows"])
    for rr in unsart["rows"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert rr["selected_alt_integration_warning_count"] >= 0
        assert rr["selected_original_integration_warning_count"] >= rr["selected_alt_integration_warning_count"]
        assert rr["selected_alt_delta_l2_vs_baseline"] >= 0.0
        assert rr["selected_original_delta_l2_vs_baseline"] >= 0.0
    assert "ur_numerical_stress_alt_dominance_map_precursor" in data
    unsad = data["ur_numerical_stress_alt_dominance_map_precursor"]
    assert unsad["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsad["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_ALT_DOMINANCE_MAP"
    assert unsad["num_rows"] == len(unsad["rows"]) == uaeisp["num_rows"]
    assert len(unsad["by_class"]) == 3
    for rr in unsad["rows"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert isinstance(rr["alt_nonworse_both_axes"], bool)
        assert isinstance(rr["alt_strictly_better_on_at_least_one_axis"], bool)
        assert isinstance(rr["alt_pareto_dominates_original"], bool)
    for rr in unsad["by_class"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert rr["num_cases"] > 0
        assert 0.0 <= rr["pareto_dominance_frequency"] <= 1.0
        assert 0.0 <= rr["pareto_dominance_frequency_wilson_interval95"]["lower"] <= rr["pareto_dominance_frequency_wilson_interval95"]["upper"] <= 1.0
        assert 0.0 <= rr["nonworse_both_axes_frequency"] <= 1.0
    assert "ur_numerical_stress_alt_decision_gate_precursor" in data
    unsadg = data["ur_numerical_stress_alt_decision_gate_precursor"]
    assert unsadg["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsadg["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_ALT_DECISION_GATE"
    assert len(unsadg["rows"]) == 3
    assert 0 <= unsadg["num_recommended_classes"] <= 3
    for rr in unsadg["rows"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert 0.0 <= rr["dominance_wilson_lb95"] <= 1.0
        assert rr["dominance_lb_threshold"] == 0.5
        assert rr["span_worsening_tolerance"] >= 0.0
        assert isinstance(rr["criteria"]["dominance_wilson_lb95_ge_threshold"], bool)
        assert isinstance(rr["criteria"]["selected_delta_l2_alt_minus_original_le_tolerance"], bool)
        assert isinstance(rr["recommend_alt_parameterization_for_class"], bool)
    assert "ur_numerical_stress_alt_hysteresis_gate_precursor" in data
    unsah = data["ur_numerical_stress_alt_hysteresis_gate_precursor"]
    assert unsah["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsah["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_ALT_HYSTERESIS_GATE"
    assert len(unsah["rows"]) == 3
    assert unsah["num_force_on_classes"] + unsah["num_hold_classes"] + unsah["num_force_off_classes"] == 3
    for rr in unsah["rows"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert rr["hysteresis_threshold_on"] > rr["hysteresis_threshold_off"]
        assert isinstance(rr["span_ok"], bool)
        assert isinstance(rr["states"]["force_on"], bool)
        assert isinstance(rr["states"]["hold_previous_state"], bool)
        assert isinstance(rr["states"]["force_off"], bool)
        assert rr["state_partition_valid"] is True
    assert "ur_numerical_stress_alt_hysteresis_time_stability_precursor" in data
    unsaht = data["ur_numerical_stress_alt_hysteresis_time_stability_precursor"]
    assert unsaht["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsaht["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_ALT_HYSTERESIS_TIME_STABILITY"
    assert len(unsaht["seeds"]) == 6
    assert len(unsaht["rows"]) == 18
    assert len(unsaht["by_class"]) == 3
    assert len(unsaht["transition_matrix_by_class"]) == 3
    assert len(unsaht["entropy_rate_by_class"]) == 3
    for rr in unsaht["rows"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert rr["state"] in {"ON", "HOLD", "OFF"}
        assert rr["n_cases"] > 0
        assert 0 <= rr["k_sim"] <= rr["n_cases"]
        assert 0.0 <= rr["p_sim"] <= 1.0
        assert 0.0 <= rr["lb95_sim"] <= 1.0
    for rr in unsaht["by_class"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert rr["num_replays"] == 6
        assert rr["state_counts"]["ON"] + rr["state_counts"]["HOLD"] + rr["state_counts"]["OFF"] == 6
        assert 0.0 <= rr["transition_frequency"] <= 1.0
    for rr in unsaht["transition_matrix_by_class"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert rr["states_order"] == ["ON", "HOLD", "OFF"]
        for s in rr["states_order"]:
            assert rr["row_totals"][s] >= 0
            assert 0.0 <= rr["self_transition_wilson_lb95"][s] <= 1.0
            rowsum = sum(rr["counts"][s][t] for t in rr["states_order"])
            assert rowsum == rr["row_totals"][s]
            if rr["row_totals"][s] > 0:
                psum = sum(rr["transition_probabilities"][s][t] for t in rr["states_order"])
                assert abs(psum - 1.0) < 1e-12
            else:
                psum = sum(rr["transition_probabilities"][s][t] for t in rr["states_order"])
                assert psum == 0.0
    for rr in unsaht["entropy_rate_by_class"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert 0.0 <= rr["entropy_rate_bits_per_step"] <= rr["max_entropy_bits_per_step_for_3states"] + 1e-12
        assert 0.0 <= rr["normalized_entropy_rate_0_1"] <= 1.0 + 1e-12
        assert set(rr["state_occupancy_pi"].keys()) == {"ON", "HOLD", "OFF"}
        assert abs(sum(rr["state_occupancy_pi"].values()) - 1.0) < 1e-12
    assert 0.0 <= unsaht["entropy_rate_global_max_bits_per_step"] <= 2.0
    assert 0.0 <= unsaht["global_transition_frequency_max"] <= 1.0
    assert "ur_numerical_stress_alt_entropy_gate_precursor" in data
    unsaeg = data["ur_numerical_stress_alt_entropy_gate_precursor"]
    assert unsaeg["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsaeg["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_ALT_ENTROPY_GATE"
    assert len(unsaeg["rows"]) == 3
    assert 0 <= unsaeg["num_recommended_classes"] <= 3
    for rr in unsaeg["rows"]:
        assert rr["class"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
        assert 0.0 <= rr["dominance_wilson_lb95"] <= 1.0
        assert rr["dominance_lb_threshold"] == 0.5
        assert rr["span_worsening_tolerance"] >= 0.0
        assert 0.0 <= rr["normalized_entropy_rate_0_1"] <= 1.0 + 1e-12
        assert rr["entropy_threshold_norm"] == 0.6
        assert isinstance(rr["criteria"]["dominance_wilson_lb95_ge_threshold"], bool)
        assert isinstance(rr["criteria"]["selected_delta_l2_alt_minus_original_le_tolerance"], bool)
        assert isinstance(rr["criteria"]["normalized_entropy_rate_le_threshold"], bool)
        assert isinstance(rr["recommend_alt_parameterization_entropy_gated"], bool)
    assert "ur_numerical_stress_alt_entropy_threshold_calibration_precursor" in data
    unsetc = data["ur_numerical_stress_alt_entropy_threshold_calibration_precursor"]
    assert unsetc["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert unsetc["scope"] == "STRICT_TASK2_NUMERICAL_STRESS_ALT_ENTROPY_THRESHOLD_CALIBRATION"
    assert len(unsetc["rows"]) == 5
    assert unsetc["selected_entropy_threshold_norm"] == 0.6
    for rr in unsetc["rows"]:
        assert 0.0 <= rr["entropy_threshold_norm"] <= 1.0
        assert 0 <= rr["num_recommended_classes"] <= 3
    assert 0.0 <= unsetc["recommendation_count_span"] <= 3.0
    assert data["phase_common_basis_link_precursor"]["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert data["phase_common_basis_link_precursor"]["condition_number"] < 1e6
    assert data["phase_common_basis_link_precursor"]["residual_l2"] > 0.0
    assert len(data["phase_common_basis_link_precursor"]["loocv_rows"]) == 3
    assert data["phase_common_basis_link_precursor"]["loocv_residual_l2_max"] < 1.0
    assert len(data["phase_common_basis_link_precursor"]["bootstrap_rows_preview"]) == 8
    assert data["phase_common_basis_link_precursor"]["bootstrap_residual_l2_p95"] < 1.0
    assert data["phase_common_basis_link_precursor"]["stability_envelope_max"] < 1.0
    assert data["phase_common_basis_link_precursor"]["cross_channel_coef_spread_l2"] < 10.0
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["status"] == "OPEN_PRECURSOR_NOT_CLOSURE"
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["residual_l2"] > 0.0
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["solver_crosscheck_objective_gap"] < 1.0
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["lambda_sweep_residual_span"] < 1.0
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["holdout_rotation_residual_l2_max"] < 1.0
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["multistart_residual_l2_span"] < 1.0
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["perturbation_residual_l2_span"] < 1.0
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["combined_stress_panel"]["worst_case_residual_envelope"] < 1.0
    assert len(data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["cross_background_stress_panel_rows"]) == 2
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["cross_background_envelope_span"] < 0.2
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["cross_background_scale_source"]["method"] == "mean_det_frw_to_bianchi_over_nu_grid"
    assert len(data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["operator_transport_replay"]["rows"]) == 2
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["operator_transport_replay"]["residual_l2_span"] < 0.2
    assert len(data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["operator_transport_nu_sweep"]["rows"]) == 4
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["operator_transport_nu_sweep"]["residual_l2_span"] < 0.2
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["operator_transport_nu_sweep"]["solver_objective_gap_max"] < 1.0
    assert len(data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["operator_transport_nu_lambda_panel"]["rows"]) == 12
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["operator_transport_nu_lambda_panel"]["residual_l2_span"] < 0.3
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["operator_transport_nu_lambda_panel"]["solver_objective_gap_max"] < 1.0
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["operator_transport_nu_lambda_panel"]["weighted_envelope_method"] == "residual_l2 * abs(det(T_frw_to_bianchi(nu)))"
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["operator_transport_nu_lambda_panel"]["weighted_residual_l2_max"] < 1.0
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["operator_transport_nu_lambda_panel"]["condition_weighted_envelope_method"] == "residual_l2 * cond(T_frw_to_bianchi(nu))"
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["operator_transport_nu_lambda_panel"]["condition_weighted_residual_l2_max"] < 2.0
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["operator_transport_dual_criterion_frontier"]["pareto_frontier_count"] >= 1
    assert len(data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["operator_transport_dual_criterion_frontier"]["frontier_continuity_rows"]) == 3
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_gauge_gauge"]["residual_l2_after_substitution"] >= 0.0
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_fermion_fermion"]["residual_l2_after_substitution"] >= 0.0
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_scalar_scalar"]["residual_l2_after_substitution"] >= 0.0
    assert "delta_residual_l2_backend_minus_baseline" in data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_delta_report"]
    assert len(data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_delta_rows"]) == 3
    assert len(data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_transport_channel_delta_rows"]) == 12
    assert len(data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_rows"]) == 3
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_best"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
    assert len(data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_rank_robustness_rows"]) == 7
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_winner_count"] >= 1
    assert len(data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_bootstrap_rows_preview"]) == 8
    assert len(data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_bootstrap_winner_frequency_rows"]) == 3
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_bootstrap_winner_frequency_max_wilson_interval95"]["upper"] <= 1.0
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_bootstrap_winner_frequency_entropy_norm"] <= 1.0
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_bootstrap_winner_frequency_top2_margin"] >= 0.0
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_dirichlet_posterior_p_best_gt_050"] <= 1.0
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_dirichlet_posterior_best_quantiles"]["q05"] <= data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_dirichlet_posterior_best_quantiles"]["q95"]
    assert len(data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_bootstrap_size_rows"]) == 3
    assert len(data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_bootstrap_size_loo_rows"]) == 3
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_bootstrap_size_trend"]["r2"] <= 1.0
    assert "quadratic_coef" in data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_bootstrap_size_curvature"]
    assert len(data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_bootstrap_size_extrapolation_rows"]) == 2
    assert len(data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_priority_bootstrap_seed_rows"]) == 3
    assert data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_first_simulation_panel"]["selected_channel"] in {"gauge_gauge", "fermion_fermion", "scalar_scalar"}
    assert len(data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_first_simulation_panel"]["transport_rows"]) == 4
    assert "replay_metrics" in data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_first_simulation_panel"]
    rm = data["phase_common_basis_link_precursor"]["joint_coupled_fit"]["backend_substitution_channel_first_simulation_panel"]["replay_metrics"]
    assert len(rm["channel_first_bootstrap_winner_frequency_rows"]) == 3
    assert len(rm["channel_first_bootstrap_size_rows"]) == 3
    assert len(rm["channel_first_bootstrap_size_extrapolation_rows"]) == 2
    assert len(rm["channel_first_seed_rows"]) == 3
    assert "paired_delta_panel" in rm
    assert rm["paired_delta_panel"]["method"] == "common_index_paired_bootstrap_baseline_vs_channel_first"
    assert rm["paired_delta_panel"]["num_pairs"] == 128
    assert len(rm["paired_delta_panel"]["rows_preview"]) == 8
    assert len(rm["paired_delta_panel"]["per_channel_rows"]) == 3
    assert "wilcoxon_quality" in rm["paired_delta_panel"]
    assert set(rm["paired_delta_panel"]["wilcoxon_quality"].keys()) == {"seed", "extrap", "q05"}
    assert "low_power_summary" in rm["paired_delta_panel"]
    assert "power_aware_verdict" in rm["paired_delta_panel"]
    assert "nonworse_probability_ci95_conditions_met" in rm["paired_delta_panel"]["power_aware_verdict"]
    assert "per_channel_power_aware_verdicts" in rm["paired_delta_panel"]
    assert len(rm["paired_delta_panel"]["per_channel_power_aware_verdicts"]) == 3
    assert "mixed_verdict_regime" in rm["paired_delta_panel"]
    assert rm["paired_delta_panel"]["mixed_verdict_regime"]["num_channels"] == 3
    assert "per_channel_risk_ranking" in rm["paired_delta_panel"]
    assert len(rm["paired_delta_panel"]["per_channel_risk_ranking"]["rows"]) == 3
    assert "time_stability_risk_panel" in rm["paired_delta_panel"]
    assert rm["paired_delta_panel"]["time_stability_risk_panel"]["num_seeds"] == 6
    assert len(rm["paired_delta_panel"]["time_stability_risk_panel"]["rows"]) == 6
    assert "time_stability_seed_robust_gate" in rm["paired_delta_panel"]
    assert "risk_signal_stable" in rm["paired_delta_panel"]["time_stability_seed_robust_gate"]
    assert "branch_cut_sensitivity_panel" in rm["paired_delta_panel"]
    assert len(rm["paired_delta_panel"]["branch_cut_sensitivity_panel"]["rows"]) == 9
    assert "loglog_slope_rows" in rm["paired_delta_panel"]["branch_cut_sensitivity_panel"]
    assert len(rm["paired_delta_panel"]["branch_cut_sensitivity_panel"]["loglog_slope_rows"]) == 5
    assert "branch_cross_integrator_panel" in rm["paired_delta_panel"]
    assert len(rm["paired_delta_panel"]["branch_cross_integrator_panel"]["rows"]) == 5
    assert rm["paired_delta_panel"]["branch_cross_integrator_panel"]["methods"] == ["scipy.integrate.quad", "numpy.trapezoid"]
    assert rm["paired_delta_panel"]["branch_cross_integrator_panel"]["trapz_grid_points"] == 4001
    assert "branch_integrator_stress_matrix" in rm["paired_delta_panel"]
    assert len(rm["paired_delta_panel"]["branch_integrator_stress_matrix"]["rows"]) == 45
    assert "branch_integrator_cross_seed_envelope" in rm["paired_delta_panel"]
    assert len(rm["paired_delta_panel"]["branch_integrator_cross_seed_envelope"]["rows"]) == 3
    assert "branch_integrator_threshold_calibration_panel" in rm["paired_delta_panel"]
    assert rm["paired_delta_panel"]["branch_integrator_threshold_calibration_panel"]["bootstrap_size"] == 256
    assert "branch_robust_substitution_decision" in rm["paired_delta_panel"]
    assert "ready_for_branch_robust_substitution" in rm["paired_delta_panel"]["branch_robust_substitution_decision"]
    for vrow in rm["paired_delta_panel"]["per_channel_power_aware_verdicts"]:
        assert "nonworse_probability_conditions_met" in vrow
        assert "nonworse_probability_ci95_conditions_met" in vrow
        assert 0.0 <= vrow["prob_seed_nonworse"] <= 1.0
        assert 0.0 <= vrow["prob_extrap_nonworse"] <= 1.0
        assert 0.0 <= vrow["prob_q05_nonworse"] <= 1.0
        assert 0.0 <= vrow["prob_seed_nonworse_ci95"]["lower"] <= vrow["prob_seed_nonworse_ci95"]["upper"] <= 1.0
        assert 0.0 <= vrow["prob_extrap_nonworse_ci95"]["lower"] <= vrow["prob_extrap_nonworse_ci95"]["upper"] <= 1.0
        assert 0.0 <= vrow["prob_q05_nonworse_ci95"]["lower"] <= vrow["prob_q05_nonworse_ci95"]["upper"] <= 1.0
    for qrow in rm["paired_delta_panel"]["per_channel_wilcoxon_quality"]:
        assert qrow["seed"]["n_pairs"] == 64
        assert qrow["extrap"]["n_pairs"] == 64
        assert qrow["q05"]["n_pairs"] == 64
    assert data["gatekeeper_checks"]["phase_joint_stress_panel_envelope_bounded"] is True
    assert data["gatekeeper_checks"]["phase_joint_cross_background_envelope_span_bounded"] is True
    assert data["gatekeeper_checks"]["phase_joint_operator_transport_replay_bounded"] is True
    assert data["gatekeeper_checks"]["phase_joint_operator_transport_nu_sweep_bounded"] is True
    assert data["gatekeeper_checks"]["phase_joint_operator_transport_nu_sweep_solver_agreement"] is True
    assert data["gatekeeper_checks"]["phase_joint_operator_transport_nu_lambda_panel_bounded"] is True
    assert data["gatekeeper_checks"]["phase_joint_operator_transport_nu_lambda_solver_agreement"] is True
    assert data["gatekeeper_checks"]["phase_joint_operator_transport_nu_lambda_weighted_envelope_bounded"] is True
    assert data["gatekeeper_checks"]["phase_joint_operator_transport_nu_lambda_condition_weighted_envelope_bounded"] is True
    assert data["gatekeeper_checks"]["phase_joint_operator_transport_dual_frontier_nonempty"] is True
    assert data["gatekeeper_checks"]["phase_joint_operator_transport_dual_frontier_continuity_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_gauge_gauge_finite"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_fermion_fermion_finite"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_scalar_scalar_finite"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_delta_report_finite"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_delta_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_transport_channel_delta_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_panel_nonempty"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_span_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_rank_robustness_rows_nonempty"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_winner_set_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_rows_nonempty"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_winner_freq_max_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_winner_freq_wilson_lb_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_entropy_norm_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_top2_margin_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_dirichlet_p_best_gt_050_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_dirichlet_q05_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_size_rows_nonempty"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_size_span_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_size_monotone_guard"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_size_loo_rows_nonempty"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_size_loo_span_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_size_slope_finite"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_size_r2_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_size_curvature_finite"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_size_aic_delta_finite"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_size_extrap_rows_nonempty"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_size_extrap_gap_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_seed_rows_nonempty"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_priority_bootstrap_seed_span_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_selected_valid"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_transport_rows_nonempty"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_cond_weighted_median_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_replay_seed_span_nonworse"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_replay_extrap_gap_nonworse"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_replay_dirichlet_q05_nonworse"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_paired_seed_nonworse_prob_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_paired_extrap_nonworse_prob_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_paired_q05_nonworse_prob_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_paired_seed_wilcoxon_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_paired_extrap_wilcoxon_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_paired_q05_wilcoxon_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_paired_per_channel_rows_nonempty"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_paired_holm_pvalues_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_paired_wilcoxon_quality_exported"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_paired_wilcoxon_effective_pairs_nonzero"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_paired_wilcoxon_effective_pairs_threshold"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_paired_low_power_summary_exported"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_power_aware_verdict_exported"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_power_aware_ready_flag_consistent"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_per_channel_power_aware_verdicts_exported"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_per_channel_power_aware_ready_flag_consistent"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_per_channel_quality_csv_json_consistent"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_per_channel_verdict_csv_json_consistent"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_mixed_verdict_regime_exported"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_per_channel_risk_ranking_exported"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_time_stability_risk_panel_exported"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_time_stability_seed_robust_gate_exported"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_branch_cut_sensitivity_exported"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_branch_loglog_slope_span_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_branch_cross_integrator_exported"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_branch_cross_integrator_agreement_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_branch_integrator_stress_matrix_exported"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_branch_integrator_worst_case_gap_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_branch_integrator_cross_seed_envelope_exported"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_branch_integrator_cross_seed_envelope_bounded"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_branch_integrator_threshold_calibration_exported"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_branch_integrator_threshold_calibration_consistent"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_branch_robust_substitution_decision_exported"] is True
    assert data["gatekeeper_checks"]["phase_backend_substitution_channel_first_branch_robust_substitution_decision_consistent"] is True
    assert all(x["status"] == "OPEN_OBSTRUCTION_WITH_TRACE" for x in data["toe_closure_gaps_7tasks"])
    assert data["depends_on"]["same_scheme_tag"] == "STRICT_P2020_PHASESPACE_SCHEME_V1"
    assert data["upstream_manifest"]["same_scheme_tag"] == "STRICT_P2020_PHASESPACE_SCHEME_V1"
    assert isinstance(data["upstream_manifest_digest_sha256"], str) and len(data["upstream_manifest_digest_sha256"]) == 64
    assert data["phase_space_grid_refinement"]["max_abs_gap"] < 1e-12
    assert data["quadrature_tolerance_robustness"]["max_span"] < 1e-10
    assert data["finite_difference_slope_consistency"]["max_abs_gap"] < 1e-5
    assert len(data["scipy_numpy_sympy_calibration"]["feature_matrix"]) == 10
    assert data["scipy_numpy_sympy_calibration"]["weighted_design_rank_numeric"] == 3
    assert data["scipy_numpy_sympy_calibration"]["weighted_design_rank_symbolic"] == 3

    assert data["scipy_numpy_sympy_calibration"]["condition_robustness"]["p95"] < 1e8
    assert data["scipy_numpy_sympy_calibration"]["bootstrap_seed_robustness"]["max_span"] < 10.0
    assert data["gatekeeper_checks"]["theorem_digest_self_consistent"] is True
    assert data["gatekeeper_checks"]["reproducibility_digest_self_consistent"] is True
    assert data["environment_lock"]["python_major"] == 3
    assert data["reproducibility_probe"]["digest_1"] == data["reproducibility_probe"]["digest_2"]
    assert data["theorem_core_digest_sha256"] == data["theorem_core_digest_recomputed_sha256"]
    assert OUT_CSV.exists()
    with OUT_CSV.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 3
    assert OUT_QUALITY_CSV.exists()
    with OUT_QUALITY_CSV.open("r", encoding="utf-8", newline="") as f:
        qrows = list(csv.DictReader(f))
    assert len(qrows) == 3
