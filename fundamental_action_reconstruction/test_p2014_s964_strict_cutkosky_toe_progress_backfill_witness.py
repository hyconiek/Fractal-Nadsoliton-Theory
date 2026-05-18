import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT_2014 = ROOT / "p2014_s964_strict_cutkosky_toe_progress_backfill_witness.py"
SCRIPT_2013 = ROOT / "p2013_s963_strict_cutkosky_pipeline_rigor_audit_witness.py"
OUT_2014 = ROOT / "generated" / "p2014_s964_strict_cutkosky_toe_progress_backfill_witness.json"
OUT_2013 = ROOT / "generated" / "p2013_s963_strict_cutkosky_pipeline_rigor_audit_witness.json"


def test_p2014_backfill_then_p2013_rigor_audit_passes_required_fields() -> None:
    subprocess.run([sys.executable, str(SCRIPT_2014)], check=True)
    subprocess.run([sys.executable, str(SCRIPT_2013)], check=True)

    data_2014 = json.loads(OUT_2014.read_text(encoding="utf-8"))
    assert data_2014["packet_id"] == "P2014"
    assert data_2014["gatekeeper_checks"]["all_targets_present"] is True
    assert data_2014["gatekeeper_checks"]["toe_progress_present_everywhere"] is True

    data_2013 = json.loads(OUT_2013.read_text(encoding="utf-8"))
    assert data_2013["packet_id"] == "P2013"
    assert data_2013["gatekeeper_checks"]["all_required_fields_present"] is True
    assert data_2013["summary"]["passed"] == data_2013["summary"]["total"]
