import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2024_s974_strict_cutkosky_minimal_channel_brst_data_object.py"
OUT = ROOT / "generated" / "p2024_s974_strict_cutkosky_minimal_channel_brst_data_object.json"


def test_p2024_exports_minimal_finite_brst_quartet_v11_without_cutkosky_promotion():
    subprocess.run([sys.executable, str(SCRIPT)], check=True)
    data = json.loads(OUT.read_text(encoding="utf-8"))

    assert data["result_kind"] == "PARTIAL_FINITE_BRST_QUARTET_MODEL_V11_EXPORTED__NOT_CHANNEL_PROJECTOR__DISCM_GHOST_CUT_STILL_OPEN_WITH_TRACE"
    assert data["status"] == "OPEN_OBSTRUCTION_WITH_TRACE"
    assert data["legacy_bridge_used"] is False
    assert data["model_scope"] == "finite_algebraic_quartet_not_full_channel_projector_v11"
    assert all(data["gatekeeper_checks"].values())
    assert data["symbolic_lock_guard"]["present_in_formal_required_loop_discontinuity"] is True
    assert data["symbolic_lock_guard"]["same_scheme_tag"] == "STRICT_P2020_PHASESPACE_SCHEME_V1"
    assert data["symbolic_lock_guard"]["same_scheme_tag_lock"] is True
    assert data["symbolic_lock_guard"]["python_major_version_required"] == 3
    assert data["gatekeeper_checks"]["python_major_version_lock"] is True
    assert isinstance(data["theorem_core_digest_sha256"], str) and len(data["theorem_core_digest_sha256"]) == 64
    assert data["gatekeeper_checks"]["theorem_core_digest_self_consistent"] is True
    assert data["theorem_core_digest_sha256"] == data["theorem_core_digest_recomputed_sha256"]
    nonuniq = data["transverse_label_nonuniqueness_witness"]
    assert nonuniq["invariant_under_T1_T2_swap"] is True

    quartet = data["FiniteBRSTQuartetModel_strict_gauge_polarization_v11"]
    assert quartet["basis_order"] == ["gauge_T1", "gauge_T2", "gauge_L", "ghost", "antighost", "B"]
    assert quartet["q_action"]["Q(gauge_L)"] == "ghost"
    assert quartet["q_action"]["Q(antighost)"] == "B"
    assert quartet["exact_checks"]["Q_squared_zero"] is True
    assert quartet["exact_checks"]["Q2_characteristic_polynomial"] == "lambda_q2**6"
    assert quartet["exact_checks"]["rank_Q"] == 2
    assert quartet["exact_checks"]["kernel_dimension"] == 4
    assert quartet["exact_checks"]["image_dimension"] == 2
    assert quartet["exact_checks"]["minimal_cohomology_dimension_kernel_mod_image"] == 2
    assert quartet["ghost_sector_trace_convention_candidate"]["quartet_supertrace"] == "0"
    assert quartet["ghost_sector_trace_convention_candidate"]["physical_transverse_supertrace"] == "2"

    numeric = quartet["numeric_cross_checks"]
    assert numeric["rank_numpy"] == 2
    assert numeric["kernel_dimension_scipy"] == 4
    assert numeric["cohomology_dimension_numeric"] == 2
    assert numeric["q2_frobenius_norm"] == 0.0
    assert numeric["q2_max_abs_eigenvalue"] == 0.0

    gap = data["symbolic_cutkosky_gap_after_minimal_BRST_witness"]
    assert gap["optimistic_unproved_assignment_defect"] == "0"
    assert gap["ghost_shifted_assignment_defect"] == "-1/11"
    assert gap["bridge_shifted_assignment_defect"] == "-1/17"
    assert gap["cohomology_trace_not_identified_with_p2020_amplitude_trace"] is True
    assert "undertermined" not in gap["meaning"]
    assert "underdetermined" in gap["meaning"]

    assert "does not identify gauge_T1/gauge_T2 with P2020 graviton plus/cross labels" in data["false_pass_guard"]
    assert gap["p2020_tree_cutsum_symbol"] == "CutSum_tree"
    assert gap["amplitude_normalization_symbol"] == "AmpNorm"
    assert "Build P2025" in data["next_honest_step"]
