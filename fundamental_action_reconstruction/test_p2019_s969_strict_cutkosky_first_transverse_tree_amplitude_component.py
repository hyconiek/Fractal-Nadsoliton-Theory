import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2019_s969_strict_cutkosky_first_transverse_tree_amplitude_component.py"
OUT = ROOT / "generated" / "p2019_s969_strict_cutkosky_first_transverse_tree_amplitude_component.json"


def test_p2019_exports_first_local_transverse_tree_component_without_closure() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True)
    data = json.loads(OUT.read_text(encoding="utf-8"))

    assert data["packet_id"] == "P2019"
    assert data["stage_id"] == "S969"
    assert data["result_kind"] == "PASS_FIRST_TREE_TRANSVERSE_COMPONENT_WITNESS"
    assert data["status"] == "OPEN_OBSTRUCTION_WITH_TRACE"
    assert data["gatekeeper_checks"]["p1955_minimal_vertex_available"] is True
    assert data["gatekeeper_checks"]["p1956_local_transverse_projector_available"] is True
    assert data["M_tree_transverse_common_basis_component"]["amplitude_matrices_over_kappa_Zgauge"]["plus"] == [["-2", "0"], ["0", "2"]]
    assert data["M_tree_transverse_common_basis_component"]["amplitude_matrices_over_kappa_Zgauge"]["cross"] == [["0", "-2"], ["-2", "0"]]
    assert data["M_tree_transverse_common_basis_component"]["AbsM_tree_transverse_sum_over_local_polarizations_over_kappa2_Zgauge2"] == "16"
    assert data["p1953_contract_update"]["M_dressed_common_basis"] == "PARTIAL_TREE_TRANSVERSE_COMPONENT_ONLY_NOT_DRESSED"
