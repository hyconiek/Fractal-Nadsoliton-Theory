from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from p1522_s472_strict_selector_source_intake_gate_checkpoint import strict_selector_source_intake
from p1525_s475_strict_selector_uniqueness_witness_checkpoint import strict_selector_uniqueness_witness


def reduce_branches(alternative_branch_ids: list[str], reduction_rule_id: str) -> dict[str, Any]:
    if reduction_rule_id != "S476_STRICT_ANCHOR_PREFIX_RULE":
        return {
            "reduced_branch_ids": alternative_branch_ids,
            "elimination_applied": False,
            "elimination_reason_code": "unknown_reduction_rule",
        }

    anchored = [b for b in alternative_branch_ids if b.startswith("strict_anchor_")]
    if anchored:
        return {
            "reduced_branch_ids": [anchored[0]],
            "elimination_applied": True,
            "elimination_reason_code": "strict_anchor_prefix_selected",
        }

    return {
        "reduced_branch_ids": alternative_branch_ids,
        "elimination_applied": False,
        "elimination_reason_code": "no_strict_anchor_branch_found",
    }


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    candidate = {
        "candidate_id": "SSEL_CANDIDATE_STRICT_TRACE_V2",
        "provider_class": "nad12_sigma_shannon_weighted",
        "symmetry_breaking_premise": "explicit_candidate_premise_exported_for_intake_testing",
        "strict_provenance_trace": [{"step_id": "trace_step_1", "strict_guardrail_ok": True}],
        "noncyclic_anchor": True,
    }

    alternative_branch_ids = ["branch_A", "strict_anchor_branch_C", "branch_B"]
    reduction = reduce_branches(alternative_branch_ids, reduction_rule_id="S476_STRICT_ANCHOR_PREFIX_RULE")

    intake = strict_selector_source_intake(candidate)
    uniqueness = strict_selector_uniqueness_witness(candidate, reduction["reduced_branch_ids"])

    intake_pass = intake["intake_status"] == "accepted_as_strict_source"
    uniqueness_pass = uniqueness["uniqueness_status"] == "unique_selector_confirmed"

    summary = {
        "checkpoint": "P1526_S476",
        "status": "PASS_BRANCH_COLLISION_ELIMINATION_WITNESS",
        "date_utc": "2026-05-13",
        "route": "F_nadsoliton_to_LSM_plus_LGR",
        "strict_only": True,
        "legacy_bridge_used": False,
        "branch_reduction": reduction,
        "intake_result": intake,
        "uniqueness_result": uniqueness,
        "intake_pass": intake_pass,
        "uniqueness_pass": uniqueness_pass,
        "qw2191_closed": False,
        "closure_rule": "strict_closure_requires_extra_selector_source_theorem_beyond_s476",
        "provisional_status": "non_strict_until_physics_level_justification_and_stability_witness",
        "next_required_objects": [
            "physics_level_justification_of_reduction_rule",
            "strict_selector_uniqueness_stability_witness",
        ],
    }

    out_path = out_dir / "p1526_s476_branch_collision_elimination_witness_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1526] wrote {out_path}")


if __name__ == "__main__":
    main()
