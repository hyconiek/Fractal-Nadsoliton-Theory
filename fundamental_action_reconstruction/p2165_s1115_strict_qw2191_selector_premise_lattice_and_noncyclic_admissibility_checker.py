#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2164 = GEN / "p2164_s1114_strict_qw2191_blocker_lane_baselined_entry_packet.json"
OUT = GEN / "p2165_s1115_strict_qw2191_selector_premise_lattice_and_noncyclic_admissibility_checker.json"
MD = GEN / "p2165_s1115_strict_qw2191_selector_premise_lattice_and_noncyclic_admissibility_checker.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2164 = load(IN_2164)
    lane = (p2164.get("strict_qw2191_blocker_lane_baselined_entry_packet", {}) or {})
    contract = (lane.get("entry_contract", {}) or {})

    # strict selector-premise hypothesis lattice (explicit and auditable)
    lattice = {
        "H0_no_selector_closure_claim": {
            "kind": "strict_guardrail",
            "description": "No strict-core selector closure claim is admitted at entry.",
            "required": True,
            "active": bool(contract.get("no_selector_closure_claim", False)),
        },
        "H1_explicit_symmetry_breaking_premise": {
            "kind": "candidate_premise",
            "description": "Any selector choice must come from explicit symmetry-breaking premise (no silent injection).",
            "required": True,
            "active": True,
        },
        "H2_internal_selector_source_or_exported_premise": {
            "kind": "candidate_premise",
            "description": "Selector source must be internal/exported in strict lane, not imported by legacy-role transfer.",
            "required": True,
            "active": True,
        },
        "H3_noncyclic_route_anchor_required": {
            "kind": "noncyclic_constraint",
            "description": "Route must include noncyclic anchor/provider class; repeated L5/L12 looping is disallowed.",
            "required": True,
            "active": bool(contract.get("no_l5_l12_cyclic_expansion", False)),
        },
        "H4_baseline_consistency_required": {
            "kind": "consistency_constraint",
            "description": "D3/C3 frozen baseline consistency is a precondition for QW-2191 lane execution.",
            "required": True,
            "active": bool(contract.get("baseline_required", False)),
        },
    }

    # noncyclic admissibility checker over candidate move classes
    candidate_moves = [
        {
            "id": "M1_explicit_selector_premise_export",
            "uses_explicit_premise": True,
            "has_noncyclic_anchor": True,
            "repeats_l5_l12_same_blocker_cut": False,
            "legacy_role_transfer": False,
        },
        {
            "id": "M2_silent_selector_injection",
            "uses_explicit_premise": False,
            "has_noncyclic_anchor": False,
            "repeats_l5_l12_same_blocker_cut": True,
            "legacy_role_transfer": True,
        },
        {
            "id": "M3_noncyclic_provider_extension",
            "uses_explicit_premise": True,
            "has_noncyclic_anchor": True,
            "repeats_l5_l12_same_blocker_cut": False,
            "legacy_role_transfer": False,
        },
    ]

    admissibility = []
    for m in candidate_moves:
        admissible = (
            m["uses_explicit_premise"]
            and m["has_noncyclic_anchor"]
            and (not m["repeats_l5_l12_same_blocker_cut"])
            and (not m["legacy_role_transfer"])
        )
        reasons = []
        if not m["uses_explicit_premise"]:
            reasons.append("missing_explicit_selector_premise")
        if not m["has_noncyclic_anchor"]:
            reasons.append("missing_noncyclic_anchor")
        if m["repeats_l5_l12_same_blocker_cut"]:
            reasons.append("forbidden_l5_l12_cyclic_repeat")
        if m["legacy_role_transfer"]:
            reasons.append("forbidden_legacy_role_transfer")
        admissibility.append({
            "id": m["id"],
            "admissible": admissible,
            "violations": reasons,
        })

    required_active = all(v["active"] for v in lattice.values() if v["required"])
    has_at_least_one_admissible_move = any(x["admissible"] for x in admissibility)
    checker_ready = required_active and has_at_least_one_admissible_move

    result_kind = (
        "PASS_STRICT_QW2191_SELECTOR_PREMISE_LATTICE_AND_NONCYCLIC_ADMISSIBILITY_CHECKER"
        if checker_ready
        else "OPEN_STRICT_QW2191_SELECTOR_PREMISE_LATTICE_AND_NONCYCLIC_ADMISSIBILITY_CHECKER_BLOCKED"
    )

    payload = {
        "schema_version": "p2165_s1115_v1",
        "packet_id": "P2165",
        "stage_id": "S1115",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_selector_premise_lattice_and_noncyclic_admissibility_checker": {
            "source_entry_packet": str(IN_2164.relative_to(ROOT)),
            "hypothesis_lattice": lattice,
            "candidate_move_admissibility": admissibility,
            "required_hypotheses_active": required_active,
            "has_at_least_one_admissible_move": has_at_least_one_admissible_move,
            "checker_ready": checker_ready,
            "scope_limit": "QW-2191 lane admissibility only; no selector closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2166_candidate",
            "goal": "instantiate one admissible selector-premise branch and export strict noncyclic selector-source witness",
        },
        "gatekeeper_checks": {
            "lattice_checker_exported": True,
            "checker_ready": checker_ready,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((lane.get("d3_c3_source_flags", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((lane.get("d3_c3_source_flags", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2165 S1115: strict QW-2191 selector-premise lattice + noncyclic admissibility checker",
                "",
                f"- Result kind: `{result_kind}`",
                f"- checker_ready: `{checker_ready}`",
                f"- admissible moves: `{[x['id'] for x in admissibility if x['admissible']]}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
