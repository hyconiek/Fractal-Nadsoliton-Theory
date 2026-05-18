import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.py"
OUT = ROOT / "generated" / "p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.json"


def test_p2017_packet_shape_and_provenance_gatekeeper_checks() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True)
    data = json.loads(OUT.read_text(encoding="utf-8"))

    assert data["packet_id"] == "P2017"
    assert data["stage_id"] == "S967"
    assert data["route"] == "strict_only"
    assert len(data["tensor_candidate_table"]) == 5
    assert data["coupled_quadrature_tensor_candidate_scan"]["count"] == 21
    assert data["diagnostic_result_kind"] in {"PASS_STRICT_KERNEL_QUADRATURE_NUMERICS", "OPEN_NUMERICAL_GAP"}
    assert data["result_kind"] == "OPEN_PROVENANCE_GAP_WITH_STRICT_QUADRATURE_TRACE"
    assert data["status"] == "OPEN_OBSTRUCTION_WITH_TRACE"
    assert data["gatekeeper_checks"]["provenance_gap_declared"] is True
    assert data["gatekeeper_checks"]["false_pass_blocked_by_provenance_gate"] is True
    assert data["provenance_gatekeeper_checks"]["strict_feynman_rule_integrand_exported"] is False
    assert "quadrature_error_bounded" in data["numerical_gatekeeper_checks"]
