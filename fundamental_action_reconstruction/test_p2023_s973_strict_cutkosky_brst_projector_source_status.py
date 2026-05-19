import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2023_s973_strict_cutkosky_brst_projector_source_status.py"
OUT = ROOT / "generated" / "p2023_s973_strict_cutkosky_brst_projector_source_status.json"


def test_p2023_separates_local_transverse_projector_from_brst_projector():
    subprocess.run([sys.executable, str(SCRIPT)], check=True)
    data = json.loads(OUT.read_text(encoding="utf-8"))

    assert data["result_kind"] == "OPEN_BRST_PROJECTOR_SOURCE_STATUS_LOCAL_TRANSVERSE_ONLY_WITH_TRACE"
    assert data["status"] == "OPEN_OBSTRUCTION_WITH_TRACE"
    assert data["legacy_bridge_used"] is False
    assert all(data["gatekeeper_checks"].values())

    assert data["local_projector_summary"]["all_local_projector_checks_pass"] is True
    assert data["local_pass_count"] == 1
    assert data["hard_failure_count"] == 4

    requirements = data["brst_projector_requirements"]
    assert {row["req_id"] for row in requirements} == {"B1", "B2", "B3", "B4", "B5"}
    b1 = next(row for row in requirements if row["req_id"] == "B1")
    assert b1["verdict"] == "PASS_LOCAL_TRANSVERSE_PROJECTOR_AVAILABLE_NOT_BRST"
    hard = [row for row in requirements if row["verdict"].startswith("FAIL")]
    assert {row["req_id"] for row in hard} == {"B2", "B3", "B4", "B5"}

    witness = data["symbolic_projector_gap_witness"]
    assert witness["local_transverse_trace_from_P2020_no_symmetry"] == "2/pi"
    assert witness["optimistic_unproved_assignment_defect"] == "0"
    assert witness["ghost_shifted_assignment_defect"] == "1/5"
    assert witness["q_shifted_assignment_defect"] == "1/7"
    assert "underdetermined" in witness["meaning"]

    assert "does not promote it to a BRST physical-state projector" in data["false_pass_guard"]
    assert "Q_BRST action" in data["next_honest_step"]
