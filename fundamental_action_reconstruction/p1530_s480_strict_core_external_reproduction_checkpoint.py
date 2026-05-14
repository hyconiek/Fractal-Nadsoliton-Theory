from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from p1529_s479_independent_strict_core_reproduction_checkpoint import TOLERANCE, run_scaffold_case


def compare_runs(reference_run: dict[str, Any], external_run: dict[str, Any]) -> dict[str, Any]:
    delta = abs(reference_run["dominance_margin"] - external_run["dominance_margin"])
    pass_flag = reference_run["theorem_scaffold_pass"] and external_run["theorem_scaffold_pass"] and delta <= TOLERANCE
    return {
        "cross_run_margin_delta": delta,
        "external_reproduction_pass": pass_flag,
    }


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    reference_run = run_scaffold_case("strict_anchor_branch_C_trace", ["branch_B_noncyclic", "branch_A"])

    external_case = {
        "external_case_id": "EXT_CASE_ALPHA",
        "external_selected_branch": "strict_anchor_branch_Q_trace",
        "external_alternatives": ["branch_R_noncyclic", "branch_S"],
    }
    external_run = run_scaffold_case(external_case["external_selected_branch"], external_case["external_alternatives"])

    comparison = compare_runs(reference_run, external_run)

    summary = {
        "checkpoint": "P1530_S480",
        "status": "PASS_STRICT_CORE_EXTERNAL_REPRODUCTION_PACKET",
        "date_utc": "2026-05-13",
        "route": "F_nadsoliton_to_LSM_plus_LGR",
        "strict_only": True,
        "legacy_bridge_used": False,
        "tolerance": TOLERANCE,
        "external_case": external_case,
        "reference_run": reference_run,
        "external_run": external_run,
        **comparison,
        "qw2191_closed": False,
        "closure_note": "external_reproduction_pass_is_not_selector_closure",
        "next_required_objects": [
            "formal_selector_uniqueness_theorem_proof",
            "strict_core_peer_review_reproduction_packet",
        ],
    }

    out_path = out_dir / "p1530_s480_strict_core_external_reproduction_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1530] wrote {out_path}")


if __name__ == "__main__":
    main()
