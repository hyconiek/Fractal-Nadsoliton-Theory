import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2025_s975_strict_cutkosky_same_scheme_cohomology_amplitude_bridge_seed.py"
OUT = ROOT / "generated" / "p2025_s975_strict_cutkosky_same_scheme_cohomology_amplitude_bridge_seed.json"


def test_p2025_exports_same_scheme_bridge_seed_without_false_closure():
    subprocess.run([sys.executable, str(SCRIPT)], check=True)
    data = json.loads(OUT.read_text(encoding="utf-8"))
    assert data["schema_version"] == "p2025_s975_v67"
    assert data["status"] == "OPEN_OBSTRUCTION_WITH_TRACE"
    assert all(data["gatekeeper_checks"].values())
    assert len(data["toe_closure_gaps_7tasks"]) == 7
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
