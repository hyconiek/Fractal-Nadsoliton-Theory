import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
P2017_SCRIPT = ROOT / "p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.py"
SCRIPT = ROOT / "p2018_s968_strict_cutkosky_p2017_provenance_nonpromotion_audit.py"
OUT = ROOT / "generated" / "p2018_s968_strict_cutkosky_p2017_provenance_nonpromotion_audit.json"


def test_p2018_blocks_p2017_promotion_when_p1953_contract_is_unsatisfied() -> None:
    subprocess.run([sys.executable, str(P2017_SCRIPT)], check=True)
    subprocess.run([sys.executable, str(SCRIPT)], check=True)
    data = json.loads(OUT.read_text(encoding="utf-8"))

    assert data["packet_id"] == "P2018"
    assert data["stage_id"] == "S968"
    assert data["result_kind"] == "PASS_PROVENANCE_NONPROMOTION_AUDIT"
    assert data["status"] == "OPEN_OBSTRUCTION_WITH_TRACE"
    assert data["promotion_verdict"] == "BLOCK_PROMOTION_TO_BACKEND_CUTKOSKY_CLOSURE"
    assert data["provenance_gap_count"] >= 5
    assert data["gatekeeper_checks"]["p1953_contract_present"] is True
    assert data["gatekeeper_checks"]["p2017_fails_p1953_dressed_amplitude_contract"] is True
    assert data["p1953_contract_audit"]["contract_complete"] is False
    assert "M_dressed_common_basis" in data["p1953_contract_audit"]["missing_contract_fields"]
    assert data["gatekeeper_checks"]["numerical_diagnostics_do_not_override_provenance"] is True
    assert data["symbolic_nonpromotion_rule"]["current_promote_allowed"] is False
