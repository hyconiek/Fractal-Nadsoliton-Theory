#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"
SRC = {
    "p1779": G / "p1779_s729_strict_current_priority_success_condition_tracker_checkpoint.json",
    "p1782": G / "p1782_s732_strict_priority_closure_gap_matrix_checkpoint.json",
    "p1780": G / "p1780_s730_strict_theorem_gate_freeze_until_joint_residual_witness_checkpoint.json",
}
OUT = G / "p1808_s758_strict_w1_full_export_gate_reconciliation_checkpoint.json"


def detect(path: Path) -> str:
    txt = path.read_text(encoding="utf-8")
    if "W1 not FULL_EXPORT" in txt:
        return "NOT_FULL_EXPORT"
    if "W1 accepted as FULL_EXPORT" in txt:
        return "FULL_EXPORT_ACCEPTED"
    return "UNKNOWN"


def main() -> None:
    states = {k: detect(v) for k, v in SRC.items()}
    uniq = sorted(set(states.values()))
    conflict = len(uniq) > 1
    out = {
        "packet_id": "P1808",
        "stage_id": "S758",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE" if conflict else "PASS_ZERO_RECONCILED",
        "w1_source_states": states,
        "w1_consistency": "CONFLICT" if conflict else "CONSISTENT",
        "gate_vector": {
            "TG1_BW": "OPEN_LOCKED_BY_W1_CONSISTENCY_CONFLICT" if conflict else "OPEN_LOCKED_BY_UNIFIED_NONPROXY_RESIDUAL"
        }
    }
    OUT.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(str(OUT))

if __name__ == "__main__":
    main()
