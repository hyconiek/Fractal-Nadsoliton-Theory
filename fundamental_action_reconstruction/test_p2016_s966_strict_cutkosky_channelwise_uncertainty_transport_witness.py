import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2016_s966_strict_cutkosky_channelwise_uncertainty_transport_witness.py"
OUT = ROOT / "generated" / "p2016_s966_strict_cutkosky_channelwise_uncertainty_transport_witness.json"


def test_p2016_packet_shape_and_gatekeeper_checks() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True)
    data = json.loads(OUT.read_text(encoding="utf-8"))

    assert data["packet_id"] == "P2016"
    assert data["stage_id"] == "S966"
    assert len(data["uncertainty_table"]) == 5
    assert data["scan_parameters"]["samples"] == 343
    assert "residue_positive_under_channelwise_scan" in data["gatekeeper_checks"]
    assert data["result_kind"] in {"PASS_CHANNELWISE_UNCERTAINTY_TRANSPORT_WITNESS", "OPEN_OBSTRUCTION_WITH_TRACE"}
