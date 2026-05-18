import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2015_s965_strict_cutkosky_ur_link_uncertainty_witness.py"
OUT = ROOT / "generated" / "p2015_s965_strict_cutkosky_ur_link_uncertainty_witness.json"


def test_p2015_packet_shape_and_gatekeeper_checks() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True)
    data = json.loads(OUT.read_text(encoding="utf-8"))

    assert data["packet_id"] == "P2015"
    assert data["stage_id"] == "S965"
    assert len(data["uncertainty_table"]) == 5
    assert data["scan_parameters"]["samples"] == 77
    assert "residue_positive_under_scan" in data["gatekeeper_checks"]
    assert data["result_kind"] in {"PASS_UR_LINK_UNCERTAINTY_WITNESS", "OPEN_OBSTRUCTION_WITH_TRACE"}
