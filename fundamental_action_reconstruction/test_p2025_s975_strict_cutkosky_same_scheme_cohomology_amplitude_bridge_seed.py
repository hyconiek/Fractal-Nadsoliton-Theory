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
    assert data["schema_version"] == "p2025_s975_v144"
    assert data["status"] == "OPEN_OBSTRUCTION_WITH_TRACE"
    assert all(data["gatekeeper_checks"].values())
    assert len(data["toe_closure_gaps_7tasks"]) == 7
    assert len(data["task_numeric_evidence_7"]) == 7
    for row in data["task_numeric_evidence_7"]:
        assert row["honest_verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE"
        assert len(row["method_stack"]) >= 2
        assert 0.0 <= row["local_readiness_score_0_1"] <= 1.0
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
