import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2013_s963_strict_cutkosky_pipeline_rigor_audit_witness.py"
OUT = ROOT / "generated" / "p2013_s963_strict_cutkosky_pipeline_rigor_audit_witness.json"


def test_p2013_rigor_audit_packet_shape_and_gatekeeper_checks() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True)
    data = json.loads(OUT.read_text(encoding="utf-8"))

    assert data["packet_id"] == "P2013"
    assert data["stage_id"] == "S963"
    assert data["result_kind"] in {"PASS_RIGOR_METADATA_AUDIT_WITNESS", "OPEN_RIGOR_METADATA_GAPS"}
    assert "audit_table" in data and len(data["audit_table"]) == 15
    assert "gatekeeper_checks" in data
    assert "all_artifacts_present" in data["gatekeeper_checks"]
    assert "all_required_fields_present" in data["gatekeeper_checks"]
