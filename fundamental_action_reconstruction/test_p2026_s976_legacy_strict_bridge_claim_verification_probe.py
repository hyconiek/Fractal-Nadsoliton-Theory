import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2026_s976_legacy_strict_bridge_claim_verification_probe.py"
OUT = ROOT / "generated" / "p2026_s976_legacy_strict_bridge_claim_verification_probe.json"


def test_p2026_bridge_claim_probe_keeps_open_verdict():
    subprocess.run([sys.executable, str(SCRIPT)], check=True)
    data = json.loads(OUT.read_text(encoding="utf-8"))
    assert data["schema_version"] == "p2026_s976_v1"
    assert data["status"] == "OPEN_OBSTRUCTION_WITH_TRACE"
    assert data["result_kind"] in {"OPEN_NONBRIDGE_EVIDENCE", "CANDIDATE_BRIDGE_NEEDS_THEOREM"}
    assert data["fit_summary"]["r2"] <= 1.0
    assert data["fit_summary"]["max_abs_gap"] >= 0.0
    assert "identity_possible_without_new_axiom" in data["symbolic_check"]
    assert len(data["hard_limits"]) == 3

