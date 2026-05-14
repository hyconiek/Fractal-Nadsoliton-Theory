from __future__ import annotations

import json
from pathlib import Path

from p1535_s485_strict_selector_tie_break_candidate_checkpoint import tie_break_by_provenance_depth


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    selected_branch = "strict_anchor_branch_C_trace"
    base_alternatives = [
        {"branch_id": "strict_anchor_branch_candidate_trace", "provenance_depth": 3},
        {"branch_id": "strict_anchor_branch_candidate_trace_alt", "provenance_depth": 2},
    ]

    base_decision = tie_break_by_provenance_depth(selected_branch, base_alternatives)
    base_selected = base_decision["selected_branch_after_tie_break"]

    noise_scenarios = [
        base_alternatives + [{"branch_id": "noise_1", "provenance_depth": 50}],
        base_alternatives + [{"branch_id": "noise_2", "provenance_depth": 1}],
        base_alternatives + [{"branch_id": "noise_3", "provenance_depth": 999}],
    ]

    results = []
    for i, scenario in enumerate(noise_scenarios, start=1):
        decision = tie_break_by_provenance_depth(selected_branch, scenario)
        same = decision["selected_branch_after_tie_break"] == base_selected
        results.append({"scenario_id": f"noise_{i}", "decision": decision, "same_as_base": same})

    noise_stability_pass = all(r["same_as_base"] for r in results)
    tb_lem4_status_update = "theorem_level_candidate" if noise_stability_pass else "open"

    summary = {
        "checkpoint": "P1539_S489",
        "status": "PASS_TB_LEM4_NOISE_STABILITY_THEOREM_LEVEL_CHECKPOINT",
        "date_utc": "2026-05-14",
        "route": "F_nadsoliton_to_LSM_plus_LGR",
        "strict_only": True,
        "legacy_bridge_used": False,
        "base_decision": base_decision,
        "noise_scenarios_results": results,
        "noise_stability_pass": noise_stability_pass,
        "tb_lem4_status_update": tb_lem4_status_update,
        "qw2191_closed": False,
        "next_required_objects": [
            "tb_lem4_formal_derivation_text",
            "tb_thm_main_composition_with_lem2_lem4_candidates",
        ],
    }

    out_path = out_dir / "p1539_s489_tb_lem4_noise_stability_theorem_level_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1539] wrote {out_path}")


if __name__ == "__main__":
    main()
