from __future__ import annotations

import json
from pathlib import Path

from p1527_s477_physics_level_justification_of_branch_reduction_checkpoint import branch_physics_score

EPS = 1e-9


def find_equivalent_alternative(selected_branch: str, alternatives: list[str]) -> dict[str, object]:
    selected_score = branch_physics_score(selected_branch)
    for alt in alternatives:
        alt_score = branch_physics_score(alt)
        if abs(selected_score - alt_score) <= EPS:
            return {
                "equivalent_alternative_found": True,
                "equivalent_branch": alt,
                "selected_score": selected_score,
                "equivalent_score": alt_score,
            }

    return {
        "equivalent_alternative_found": False,
        "equivalent_branch": None,
        "selected_score": selected_score,
        "equivalent_score": None,
    }


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    selected_branch = "strict_anchor_branch_C_trace"
    alternatives = [
        "branch_A",
        "branch_B_noncyclic",
        "strict_anchor_branch_candidate_trace",  # equal score to selected
    ]

    scan = find_equivalent_alternative(selected_branch, alternatives)
    lem4_status_update = "open" if scan["equivalent_alternative_found"] else "partial"

    summary = {
        "checkpoint": "P1534_S484",
        "status": "PASS_LEM4_COUNTEREXAMPLE_SCANNER",
        "date_utc": "2026-05-14",
        "route": "F_nadsoliton_to_LSM_plus_LGR",
        "strict_only": True,
        "legacy_bridge_used": False,
        "selected_branch": selected_branch,
        "alternative_branches": alternatives,
        "counterexample_scan": scan,
        "lem4_status_update": lem4_status_update,
        "qw2191_closed": False,
        "next_required_objects": [
            "strict_selector_tie_break_theorem_or_internal_source",
            "lem4_upgrade_proof_after_tie_break_resolution",
        ],
    }

    out_path = out_dir / "p1534_s484_lem4_counterexample_scanner_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1534] wrote {out_path}")


if __name__ == "__main__":
    main()
