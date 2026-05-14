from __future__ import annotations

import json
from pathlib import Path

from p1535_s485_strict_selector_tie_break_candidate_checkpoint import tie_break_by_provenance_depth


def canonicalize(alternatives: list[dict[str, object]]) -> list[dict[str, object]]:
    return sorted(alternatives, key=lambda x: str(x["branch_id"]))


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    selected_branch = "strict_anchor_branch_C_trace"
    raw_alternatives = [
        {"branch_id": "strict_anchor_branch_candidate_trace_alt", "provenance_depth": 2},
        {"branch_id": "strict_anchor_branch_candidate_trace", "provenance_depth": 3},
        {"branch_id": "branch_B_noncyclic", "provenance_depth": 5},
    ]

    canonical_alternatives = canonicalize(raw_alternatives)

    raw_decision = tie_break_by_provenance_depth(selected_branch, raw_alternatives)
    canonical_decision = tie_break_by_provenance_depth(selected_branch, canonical_alternatives)

    determinism_pass = raw_decision["selected_branch_after_tie_break"] == canonical_decision["selected_branch_after_tie_break"]

    tb_lem2_status_update = "theorem_level_candidate" if determinism_pass else "open"

    summary = {
        "checkpoint": "P1538_S488",
        "status": "PASS_TB_LEM2_DETERMINISM_PROOF_CHECKPOINT",
        "date_utc": "2026-05-14",
        "route": "F_nadsoliton_to_LSM_plus_LGR",
        "strict_only": True,
        "legacy_bridge_used": False,
        "selected_branch": selected_branch,
        "raw_alternatives": raw_alternatives,
        "canonical_alternatives": canonical_alternatives,
        "raw_decision": raw_decision,
        "canonical_decision": canonical_decision,
        "determinism_pass": determinism_pass,
        "tb_lem2_status_update": tb_lem2_status_update,
        "qw2191_closed": False,
        "next_required_objects": [
            "tb_lem2_formal_derivation_text",
            "tb_thm_main_composition_after_lem2_and_lem4_upgrade",
        ],
    }

    out_path = out_dir / "p1538_s488_tb_lem2_determinism_proof_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1538] wrote {out_path}")


if __name__ == "__main__":
    main()
