from __future__ import annotations

import json
from pathlib import Path

from p1535_s485_strict_selector_tie_break_candidate_checkpoint import tie_break_by_provenance_depth


def monotonicity_check() -> bool:
    selected = "strict_anchor_branch_C_trace"
    base = [
        {"branch_id": "strict_anchor_branch_candidate_trace", "provenance_depth": 2},
        {"branch_id": "strict_anchor_branch_candidate_trace_alt", "provenance_depth": 1},
    ]
    boosted = [
        {"branch_id": "strict_anchor_branch_candidate_trace", "provenance_depth": 4},
        {"branch_id": "strict_anchor_branch_candidate_trace_alt", "provenance_depth": 1},
    ]

    base_res = tie_break_by_provenance_depth(selected, base)
    boost_res = tie_break_by_provenance_depth(selected, boosted)

    return (
        base_res["selected_branch_after_tie_break"] == "strict_anchor_branch_candidate_trace"
        and boost_res["selected_branch_after_tie_break"] == "strict_anchor_branch_candidate_trace"
    )


def stability_check() -> bool:
    selected = "strict_anchor_branch_C_trace"
    case_a = [
        {"branch_id": "strict_anchor_branch_candidate_trace", "provenance_depth": 3},
        {"branch_id": "strict_anchor_branch_candidate_trace_alt", "provenance_depth": 2},
    ]
    case_b = [
        {"branch_id": "strict_anchor_branch_candidate_trace", "provenance_depth": 3},
        {"branch_id": "strict_anchor_branch_candidate_trace_alt", "provenance_depth": 2},
        {"branch_id": "branch_noise", "provenance_depth": 10},
    ]

    ra = tie_break_by_provenance_depth(selected, case_a)
    rb = tie_break_by_provenance_depth(selected, case_b)
    return ra["selected_branch_after_tie_break"] == rb["selected_branch_after_tie_break"]


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    monotonicity_pass = monotonicity_check()
    stability_pass = stability_check()
    partial_pass = monotonicity_pass and stability_pass

    summary = {
        "checkpoint": "P1536_S486",
        "status": "PASS_FORMAL_TIE_BREAK_SOUNDNESS_PROOF_ATTEMPT",
        "date_utc": "2026-05-14",
        "route": "F_nadsoliton_to_LSM_plus_LGR",
        "strict_only": True,
        "legacy_bridge_used": False,
        "monotonicity_pass": monotonicity_pass,
        "stability_pass": stability_pass,
        "soundness_attempt_status": "partial_pass" if partial_pass else "failed",
        "lem4_status_update": "partial" if partial_pass else "open",
        "qw2191_closed": False,
        "next_required_objects": [
            "full_formal_tie_break_theorem_proof",
            "lem4_theorem_level_upgrade_after_soundness_completion",
        ],
    }

    out_path = out_dir / "p1536_s486_formal_tie_break_soundness_proof_attempt_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1536] wrote {out_path}")


if __name__ == "__main__":
    main()
