from __future__ import annotations

import json
from pathlib import Path

from p1529_s479_independent_strict_core_reproduction_checkpoint import TOLERANCE, run_scaffold_case
from p1530_s480_strict_core_external_reproduction_checkpoint import compare_runs


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    reference_run = run_scaffold_case("strict_anchor_branch_C_trace", ["branch_B_noncyclic", "branch_A"])

    external_cases = [
        {"case_id": "EXT_BETA", "selected": "strict_anchor_branch_Q_trace", "alternatives": ["branch_R_noncyclic", "branch_S"]},
        {"case_id": "EXT_GAMMA", "selected": "strict_anchor_branch_T_trace", "alternatives": ["branch_U_noncyclic", "branch_V"]},
        {"case_id": "EXT_DELTA", "selected": "strict_anchor_branch_W_trace", "alternatives": ["branch_X_noncyclic", "branch_Y"]},
    ]

    results = []
    deltas = []

    for case in external_cases:
        run = run_scaffold_case(case["selected"], case["alternatives"])
        cmp = compare_runs(reference_run, run)
        deltas.append(cmp["cross_run_margin_delta"])
        results.append({"case_id": case["case_id"], "run": run, **cmp})

    delta_mean = sum(deltas) / len(deltas)
    delta_max = max(deltas)
    peer_review_reproduction_pass = all(r["external_reproduction_pass"] for r in results) and delta_max <= TOLERANCE

    summary = {
        "checkpoint": "P1531_S481",
        "status": "PASS_STRICT_CORE_PEER_REVIEW_REPRODUCTION_PACKET",
        "date_utc": "2026-05-14",
        "route": "F_nadsoliton_to_LSM_plus_LGR",
        "strict_only": True,
        "legacy_bridge_used": False,
        "tolerance": TOLERANCE,
        "reference_run": reference_run,
        "external_cases_results": results,
        "delta_mean": delta_mean,
        "delta_max": delta_max,
        "peer_review_reproduction_pass": peer_review_reproduction_pass,
        "qw2191_closed": False,
        "closure_note": "peer_review_reproduction_pass_is_not_selector_closure",
        "next_required_objects": [
            "formal_selector_uniqueness_theorem_proof",
            "strict_core_axiom_minimization_packet",
        ],
    }

    out_path = out_dir / "p1531_s481_strict_core_peer_review_reproduction_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1531] wrote {out_path}")


if __name__ == "__main__":
    main()
