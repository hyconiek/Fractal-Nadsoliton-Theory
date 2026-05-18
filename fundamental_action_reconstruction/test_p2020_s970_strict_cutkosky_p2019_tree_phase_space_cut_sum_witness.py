import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2020_s970_strict_cutkosky_p2019_tree_phase_space_cut_sum_witness.py"
OUT = ROOT / "generated" / "p2020_s970_strict_cutkosky_p2019_tree_phase_space_cut_sum_witness.json"


def test_p2020_integrates_p2019_tree_component_without_cutkosky_closure() -> None:
    subprocess.run([sys.executable, str(ROOT / "p2019_s969_strict_cutkosky_first_transverse_tree_amplitude_component.py")], check=True)
    subprocess.run([sys.executable, str(SCRIPT)], check=True)
    data = json.loads(OUT.read_text(encoding="utf-8"))

    assert data["packet_id"] == "P2020"
    assert data["stage_id"] == "S970"
    assert data["result_kind"] == "PASS_TREE_PHASE_SPACE_LINEAR_POLARIZATION_CUT_SUM_COMPONENT_WITNESS"
    assert data["status"] == "OPEN_OBSTRUCTION_WITH_TRACE"
    assert data["amplitude_input"]["polarization_sum_over_kappa2_Zgauge2"] == "16"
    assert data["exact_phase_space_cut_sum"]["CutSum_tree_identical_final_state_over_kappa2_Zgauge2"] == "1/pi"
    assert data["exact_phase_space_cut_sum"]["CutSum_tree_no_identical_symmetry_over_kappa2_Zgauge2"] == "2/pi"
    assert data["gatekeeper_checks"]["scipy_quadrature_matches_exact_integral"] is True
    assert data["gatekeeper_checks"]["generic_angular_transport_symbolically_verified"] is True
    assert data["angular_transport_certificate"]["transported_polarization_sum_over_kappa2_Zgauge2"] == "16"
    assert data["angular_transport_certificate"]["matches_p2019_axis_frame"] is True
    assert "does not claim circular helicity" in data["angular_transport_certificate"]["nomenclature_guard"]
    assert data["gatekeeper_checks"]["linear_polarization_resolved_cut_sum_exported"] is True
    assert data["gatekeeper_checks"]["linear_polarization_resolved_cut_sum_psd"] is True
    assert data["amplitude_input"]["final_state_summed_graviton_linear_polarization_gram_over_kappa2_Zgauge2"] == [["8", "0"], ["0", "8"]]
    assert data["exact_phase_space_cut_sum"]["linear_polarization_resolved_CutSum_no_identical_symmetry_over_kappa2_Zgauge2"] == [["1/pi", "0"], ["0", "1/pi"]]
    assert data["exact_phase_space_cut_sum"]["linear_polarization_resolved_CutSum_identical_final_state_over_kappa2_Zgauge2"] == [["1/(2*pi)", "0"], ["0", "1/(2*pi)"]]
    assert data["p1953_contract_update"]["CutSum_common_basis"] == "PARTIAL_TREE_LINEAR_POLARIZATION_RESOLVED_CUTSUM_CONVENTION_WINDOW_NOT_DRESSED_NOT_DISC_MATCHED"
    assert "not DiscM=CutSum" in data["false_pass_guard"]
