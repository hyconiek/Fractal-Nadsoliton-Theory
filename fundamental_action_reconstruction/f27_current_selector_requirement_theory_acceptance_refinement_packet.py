#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "f27_current_selector_requirement_theory_acceptance_refinement_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f27_current_selector_requirement_theory_acceptance_refinement_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def contains_any(text: str, needles: list[str]) -> bool:
    lowered = text.lower()
    return any(needle in lowered for needle in needles)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p108 = load_json(
        "fundamental_action_reconstruction/generated/p108_current_selector_symmetry_breaking_requirement_probe_summary.json"
    )
    n118 = load_json(
        "fundamental_action_reconstruction/generated/n118_current_selector_or_symmetry_breaking_requirement_theorem_for_qw2191_summary.json"
    )
    s2 = load_json(
        "fundamental_action_reconstruction/generated/s2_current_far_strategic_priority_reorientation_packet_summary.json"
    )

    authoritative_paths = [
        REPO / "README.md",
        ROOT / "README.md",
        ROOT / "N118_CURRENT_SELECTOR_OR_SYMMETRY_BREAKING_REQUIREMENT_THEOREM_FOR_QW2191.md",
        ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
        REPO
        / "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2192_MODE_INDEX_SELECTION_AXIOM_GATE.md",
        REPO
        / "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2193_SELECTION_AXIOM_FAMILY_ROBUSTNESS_GATE.md",
    ]
    authoritative_text = "\n".join(
        path.read_text(encoding="utf-8") for path in authoritative_paths
    )

    acceptance_needles = [
        "selector requirement is accepted as part of the theory",
        "symmetry-breaking requirement is accepted as part of the theory",
        "selector or symmetry-breaking requirement is now part of the theory",
        "theory-level acceptance of selector or symmetry-breaking requirement",
        "selection axiom is now accepted into the theory",
    ]
    deferral_needles = [
        "selector requirement decision is deferred",
        "symmetry-breaking requirement decision is deferred",
        "theory-level decision remains deferred",
        "selector requirement remains only an open option",
        "project declines theory-level acceptance of selector requirement",
    ]

    theory_level_acceptance_verdict_present = contains_any(
        authoritative_text, acceptance_needles
    )
    theory_level_deferral_verdict_present = contains_any(
        authoritative_text, deferral_needles
    )

    selector_requirement_supported = (
        p108["status"]
        == "CURRENT_REPO_SUPPORTS_THE_SELECTOR_OR_SYMMETRY_BREAKING_REQUIREMENT_CONCLUSION_FOR_THE_QW2191_UNIQUENESS_FRONTIER_AFTER_P108"
        and n118["theorem_result"][
            "selector_or_symmetry_breaking_requirement_supported_on_current_repo_state"
        ]
        is True
    )
    strategic_priority_formalized = (
        "explicit_selector_requirement_or_symmetry_breaking_after_QW2191"
        in s2["priority_order"]
    )

    checks_spec = [
        {
            "id": "selector_requirement_supported",
            "actual": selector_requirement_supported,
            "expected": True,
            "meaning": "P108/N118 already support the selector/symmetry-breaking requirement conclusion on the current repo state",
        },
        {
            "id": "strategic_priority_formalized",
            "actual": strategic_priority_formalized,
            "expected": True,
            "meaning": "S2 already elevates the selector/symmetry-breaking requirement question to top-level project priority",
        },
        {
            "id": "theory_level_acceptance_verdict_present",
            "actual": theory_level_acceptance_verdict_present,
            "expected": False,
            "meaning": "the current authoritative source set exports no explicit theory-level acceptance verdict",
        },
        {
            "id": "theory_level_deferral_verdict_present",
            "actual": theory_level_deferral_verdict_present,
            "expected": False,
            "meaning": "the current authoritative source set exports no explicit theory-level deferral verdict either",
        },
    ]

    checks: list[dict[str, Any]] = []
    for item in checks_spec:
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": item["actual"] == item["expected"],
                "meaning": item["meaning"],
            }
        )

    artifact = {
        "stage": "F27",
        "lane": "current_selector_requirement_theory_acceptance_refinement_qw2191_followup_only",
        "goal": "refine_the_remaining_theory_level_selector_requirement_gap_into_acceptance_vs_deferral_branches",
        "status": "F27_EXECUTED_CURRENT_SELECTOR_REQUIREMENT_THEORY_ACCEPTANCE_REFINEMENT_PACKET_NO_FALSE_PASS",
        "reason": "P108/N118 already support the selector requirement conclusion and S2 already elevates it to top-level priority, but the current authoritative source set exports neither an explicit theory-level acceptance verdict nor an explicit theory-level deferral verdict; therefore the narrowest honest refinement is to split the remaining gap into acceptance and deferral branches",
        "support_state": {
            "selector_requirement_supported": selector_requirement_supported,
            "strategic_priority_formalized": strategic_priority_formalized,
            "theory_level_acceptance_verdict_present": theory_level_acceptance_verdict_present,
            "theory_level_deferral_verdict_present": theory_level_deferral_verdict_present,
        },
        "remaining_missing_objects": [
            "explicit_theory_level_acceptance_verdict_adopting_the_selector_or_symmetry_breaking_requirement_into_axiom_augmented_theory_scope_if_no_internal_source_is_derived",
            "explicit_theory_level_deferral_verdict_keeping_the_selector_or_symmetry_breaking_requirement_as_an_active_boundary_while_internal_source_search_continues",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F27",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "support_state": artifact["support_state"],
        "remaining_missing_objects": artifact["remaining_missing_objects"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()
