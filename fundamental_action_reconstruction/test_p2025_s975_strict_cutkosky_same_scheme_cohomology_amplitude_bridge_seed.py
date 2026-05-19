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
    assert data["schema_version"] == "p2025_s975_v16"
    assert data["status"] == "OPEN_OBSTRUCTION_WITH_TRACE"
    assert all(data["gatekeeper_checks"].values())
    assert len(data["toe_closure_gaps_7tasks"]) == 7
    assert data["backend_loop_fit_precursor"]["loss_l2"] > 0.0
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
