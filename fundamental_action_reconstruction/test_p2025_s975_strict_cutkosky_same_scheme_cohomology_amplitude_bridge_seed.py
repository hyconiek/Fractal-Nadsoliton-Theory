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
    assert data["schema_version"] == "p2025_s975_v37"
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
