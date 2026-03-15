#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

AS_OF = "2026-03-15"

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
F453_ASSIGNMENT = GENERATED / "mode_index_assignment_canonical_local_diagonal_strict_derived_v1.json"
F454_ASSIGNMENT = GENERATED / "mode_index_assignment_shannon_element_order_reference_strict_core_v1.json"


def load_if_exists(path: Path) -> dict | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def has_object(d: dict | None, obj: str) -> bool:
    return bool(d is not None and d.get("object") == obj)


diag_assignment = load_if_exists(F453_ASSIGNMENT)
shannon_assignment = load_if_exists(F454_ASSIGNMENT)

found_diagonal_axis_only = has_object(diag_assignment, "ModeIndexAssignment_canonical_local_diagonal_strict_derived_v1")
found_shannon_axis_only = has_object(shannon_assignment, "ModeIndexAssignment_shannon_element_order_reference_strict_core_v1")

strict_internal_selector_derivations_found = int(found_diagonal_axis_only) + int(found_shannon_axis_only)

status = "B2_EXECUTED_NO_STRICT_INTERNAL_ORIENTATION_DATUM_FOUND_AXIOM_FREE_UNIQUENESS_REMAINS_OPEN_NO_FALSE_PASS"
reduced_blocker = "no derived internal orientation datum exists in the current strict core to discharge QW-2191"
findings_strict_orientation_datum_status = "not_found_in_strict_core"
if strict_internal_selector_derivations_found > 0:
    status = "B2_UPDATED_STRICT_INTERNAL_ORIENTATION_DATUM_FOUND_ON_DIAGONAL_AND_SHANNON_LANES_AXIS_ONLY_RESIDUAL_Z2_REMAINS_NO_FALSE_PASS"
    findings_strict_orientation_datum_status = "found_axis_only_residual_z2"
    reduced_blocker = (
        "axis-only internal orientation datum exists on diagonal/local and/or Shannon ord-reference lanes (O(2)->Z2 cut), "
        "but residual Z2 sign remains and kernel-alone/global QW-2191 discharge remains open"
    )

summary = {
    "program_id": "fundamental_action_reconstruction",
    "as_of": AS_OF,
    "status": status,
    "source_policy": {
        "strict_admissible": [
            "QW-2190",
            "QW-2191",
            "QW-2192",
            "QW-2193",
            "A1",
            "A5",
            "A6",
            "A10",
        ],
        "ontology_context_only": [
            "TOE_FINAL_DOCUMENTATION.tex",
            "README.md",
        ],
        "heuristic_only": [
            "QW-1622",
            "QW-1210",
            "QW-1891",
        ],
    },
    "anti_overclaim": {
        "internal_orientation_datum_derived_claim": strict_internal_selector_derivations_found > 0,
        "axiom_free_uniqueness_closed_claim": False,
        "single_nadsoliton_ontology_discharge_claim": False,
        "fr_route_promoted_to_strict_core_claim": False,
    },
    "b2": {
        "audit_question": "Does the current strict core already contain an internal orientation datum that can replace the external selector?",
        "findings": [
            {
                "object": "single nadsoliton ontology",
                "status": "constructive_guidance_only",
                "source": "A1",
            },
            {
                "object": "strict internal orientation datum",
                "status": findings_strict_orientation_datum_status,
                "sources": [
                    "N487/F453/N492 (canonical local-diagonal lane)",
                    "N480/N488/N496/F454 (Shannon ord-reference lane)",
                ],
                "artifacts_present": {
                    "ModeIndexAssignment_canonical_local_diagonal_strict_derived_v1": found_diagonal_axis_only,
                    "ModeIndexAssignment_shannon_element_order_reference_strict_core_v1": found_shannon_axis_only,
                },
            },
            {
                "object": "kernel invariant selecting one O(2) point",
                "status": "not_found_in_strict_core",
                "source": "QW-2191",
            },
            {
                "object": "explicit selector axiom",
                "status": "control_route_only",
                "source": "QW-2192",
            },
            {
                "object": "robust selector family",
                "status": "control_family_only",
                "source": "QW-2193",
            },
            {
                "object": "FR/topological phase route",
                "status": "heuristically_plausible_unresolved",
                "sources": ["QW-1622", "QW-1210"],
            },
            {
                "object": "weak nadsoliton derivational constraints",
                "status": "insufficient",
                "source": "QW-1891",
            },
        ],
        "strict_internal_selector_derivations_found": strict_internal_selector_derivations_found,
        "heuristic_candidate_routes_found": 2,
        "reduced_blocker": reduced_blocker,
        "verdict": (
            "axis-only internal orientation datum found (lane-scoped); global kernel-alone uniqueness and residual Z2 sign remain open"
            if strict_internal_selector_derivations_found > 0
            else "blocker tightened further by eliminating hidden-source ambiguity, not discharged"
        ),
        "next_step": "B3",
    },
}

out = ROOT / "generated" / "b2_internal_orientation_datum_source_audit_summary.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
print(f"wrote {out}")
