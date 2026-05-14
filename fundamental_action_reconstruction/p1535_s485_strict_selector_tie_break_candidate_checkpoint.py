from __future__ import annotations

import json
from pathlib import Path

from p1527_s477_physics_level_justification_of_branch_reduction_checkpoint import branch_physics_score


def tie_break_by_provenance_depth(selected_branch: str, alternatives: list[dict[str, object]]) -> dict[str, object]:
    selected_score = branch_physics_score(selected_branch)
    tied = [a for a in alternatives if abs(branch_physics_score(str(a["branch_id"])) - selected_score) < 1e-9]

    if not tied:
        return {
            "tie_detected": False,
            "tie_break_status": "no_tie",
            "selected_branch_after_tie_break": selected_branch,
        }

    max_depth = max(int(a["provenance_depth"]) for a in tied)
    winners = [a for a in tied if int(a["provenance_depth"]) == max_depth]

    if len(winners) == 1:
        return {
            "tie_detected": True,
            "tie_break_status": "resolved_by_provenance_depth",
            "selected_branch_after_tie_break": winners[0]["branch_id"],
        }

    return {
        "tie_detected": True,
        "tie_break_status": "unresolved_tie",
        "selected_branch_after_tie_break": None,
    }


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    selected_branch = "strict_anchor_branch_C_trace"
    alternatives = [
        {"branch_id": "strict_anchor_branch_candidate_trace", "provenance_depth": 3},
        {"branch_id": "strict_anchor_branch_candidate_trace_alt", "provenance_depth": 2},
        {"branch_id": "branch_B_noncyclic", "provenance_depth": 5},
    ]

    tie_break = tie_break_by_provenance_depth(selected_branch, alternatives)

    lem4_status_update = "partial" if tie_break["tie_break_status"] == "resolved_by_provenance_depth" else "open"

    summary = {
        "checkpoint": "P1535_S485",
        "status": "PASS_STRICT_SELECTOR_TIE_BREAK_CANDIDATE",
        "date_utc": "2026-05-14",
        "route": "F_nadsoliton_to_LSM_plus_LGR",
        "strict_only": True,
        "legacy_bridge_used": False,
        "selected_branch_input": selected_branch,
        "alternatives": alternatives,
        "tie_break_result": tie_break,
        "lem4_status_update": lem4_status_update,
        "qw2191_closed": False,
        "next_required_objects": [
            "formal_tie_break_theorem_proof",
            "lem4_upgrade_proof_with_tie_break_soundness",
        ],
    }

    out_path = out_dir / "p1535_s485_strict_selector_tie_break_candidate_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1535] wrote {out_path}")


if __name__ == "__main__":
    main()
