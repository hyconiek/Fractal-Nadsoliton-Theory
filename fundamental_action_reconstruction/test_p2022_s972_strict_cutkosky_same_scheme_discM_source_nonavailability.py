import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2022_s972_strict_cutkosky_same_scheme_discM_source_nonavailability.py"
OUT = ROOT / "generated" / "p2022_s972_strict_cutkosky_same_scheme_discM_source_nonavailability.json"


def test_p2022_exports_same_scheme_discm_source_nonavailability_without_promotion():
    subprocess.run([sys.executable, str(SCRIPT)], check=True)
    data = json.loads(OUT.read_text(encoding="utf-8"))

    assert data["result_kind"] == "OPEN_SAME_SCHEME_DISCM_SOURCE_PARTIAL_TREE_VERTEX_AVAILABLE_WITH_TRACE"
    assert data["status"] == "OPEN_OBSTRUCTION_WITH_TRACE"
    assert data["legacy_bridge_used"] is False
    assert all(data["gatekeeper_checks"].values())

    requirements = data["discM_source_requirements"]
    assert len(requirements) == 5
    assert data["hard_failure_count"] == 4
    assert data["partial_source_count"] == 1
    assert {row["req_id"] for row in requirements} == {"D1", "D2", "D3", "D4", "D5"}

    d1 = next(row for row in requirements if row["req_id"] == "D1")
    assert d1["verdict"] == "PARTIAL_TREE_LEVEL_SOURCE_AVAILABLE_NOT_LOOP_DISCM_READY"
    assert any("P1955 minimal hAA vertex" in item for item in d1["evidence"])
    assert any("P2019 P1955/P1956 transverse tree component" in item for item in d1["evidence"])

    hard = [row for row in requirements if row["verdict"].startswith("FAIL")]
    assert {row["req_id"] for row in hard} == {"D2", "D3", "D4", "D5"}

    d3 = next(row for row in requirements if row["req_id"] == "D3")
    assert d3["verdict"] == "FAIL_PROXY_OR_SEED_ONLY"
    assert any("P2021 rejects" in item for item in d3["evidence"])

    witness = data["symbolic_underdetermination_witness"]
    assert witness["cutSum_P2020_no_symmetry_matrix"] == [["1/pi", "0"], ["0", "1/pi"]]
    assert witness["assignment_zero_defect_trace"] == "-2/pi"
    assert witness["assignment_cut_like_defect_trace"] == "0"
    assert "cannot decide Cutkosky equality" in witness["meaning"]

    assert "not DiscM=CutSum" in data["false_pass_guard"]
    assert "BRST physical-state projector" in data["next_honest_step"]
    assert "existing tree hAA chain" in data["next_honest_step"]
