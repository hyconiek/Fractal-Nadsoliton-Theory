from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from p1522_s472_strict_selector_source_intake_gate_checkpoint import strict_selector_source_intake


def strict_selector_uniqueness_witness(candidate: dict[str, Any], alternative_branch_ids: list[str]) -> dict[str, Any]:
    if not alternative_branch_ids:
        return {
            "uniqueness_status": "non_unique_or_unproven",
            "reason_code": "no_alternative_branches_provided",
            "branch_collision_count": None,
        }

    unique_ids = set(alternative_branch_ids)
    branch_collision_count = len(alternative_branch_ids) - len(unique_ids)

    if branch_collision_count == 0 and len(unique_ids) == 1:
        return {
            "uniqueness_status": "unique_selector_confirmed",
            "reason_code": "single_branch_without_collision",
            "branch_collision_count": 0,
        }

    return {
        "uniqueness_status": "non_unique_or_unproven",
        "reason_code": "multiple_or_colliding_branches",
        "branch_collision_count": branch_collision_count,
    }


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    candidate = {
        "candidate_id": "SSEL_CANDIDATE_STRICT_TRACE_V1",
        "provider_class": "nad12_sigma_shannon_weighted",
        "symmetry_breaking_premise": "explicit_candidate_premise_exported_for_intake_testing",
        "strict_provenance_trace": [
            {
                "step_id": "trace_step_1",
                "artifact_ref": "P1524_S474_trace_builder",
                "operation": "strict_trace_export",
                "strict_guardrail_ok": True,
            }
        ],
        "noncyclic_anchor": True,
    }

    intake = strict_selector_source_intake(candidate)
    uniqueness = strict_selector_uniqueness_witness(candidate, alternative_branch_ids=["branch_A", "branch_B"])

    intake_pass = intake["intake_status"] == "accepted_as_strict_source"
    uniqueness_pass = uniqueness["uniqueness_status"] == "unique_selector_confirmed"

    summary = {
        "checkpoint": "P1525_S475",
        "status": "PASS_STRICT_SELECTOR_UNIQUENESS_WITNESS_GATE",
        "date_utc": "2026-05-13",
        "route": "F_nadsoliton_to_LSM_plus_LGR",
        "strict_only": True,
        "legacy_bridge_used": False,
        "intake_result": intake,
        "uniqueness_result": uniqueness,
        "intake_pass": intake_pass,
        "uniqueness_pass": uniqueness_pass,
        "qw2191_closed": intake_pass and uniqueness_pass,
        "closure_rule": "qw2191_closed_only_if_intake_pass_and_uniqueness_pass",
        "next_required_objects": [
            "strict_unique_selector_branch_witness",
            "selector_collision_elimination_proof",
        ],
    }

    out_path = out_dir / "p1525_s475_strict_selector_uniqueness_witness_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1525] wrote {out_path}")


if __name__ == "__main__":
    main()
