#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_N674 = GENERATED / "n674_current_strict_t172_projective_closure_discharge_theorem_summary.json"
IN_N675 = GENERATED / "n675_current_strict_global_directed_selector_closure_obstruction_boundary_theorem_summary.json"
IN_N678 = GENERATED / "n678_current_strict_t172_directed_closure_discharge_theorem_summary.json"
IN_N676 = GENERATED / "n676_current_first_admissible_s_sel_int_source_object_discharge_theorem_summary.json"
IN_N680 = GENERATED / "n680_current_strict_t173_projective_strict_core_selector_closure_discharge_theorem_summary.json"

OUT = (
    GENERATED
    / "n679_current_strict_t172_strict_core_selector_closure_frontier_boundary_theorem_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    prereq = (IN_N674, IN_N675, IN_N678, IN_N676, IN_N680)
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        summary = {
            "step": "N679",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES_FOR_T172_FRONTIER_BOUNDARY_THEOREM",
            "scope": "current_strict_t172_remaining_frontier_boundary_only",
            "missing": missing,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    n674 = load_json(IN_N674)
    n675 = load_json(IN_N675)
    n678 = load_json(IN_N678)
    n676 = load_json(IN_N676)
    n680 = load_json(IN_N680)

    t172_projective_discharged = bool(
        n674.get("theorem_result", {}).get("discharged") is True
    )
    directed_raw_obstructed_without_sign_lift = bool(
        n675.get("theorem_result", {}).get(
            "directed_global_closure_obstructed_without_explicit_sign_lift"
        )
        is True
    )
    t172_directed_discharged = bool(
        n678.get("theorem_result", {}).get("discharged") is True
    )
    admissible_source_object_exported = bool(
        n676.get("theorem_result", {}).get(
            "admissible_S_sel_int_source_object_in_F34_sense"
        )
        is True
    )

    strict_core_selector_closure_projective_discharged = bool(
        n680.get("theorem_result", {}).get("discharged") is True
        and bool(n680.get("theorem_result", {}).get("strict_core_selector_closure")) is True
        and str(n680.get("theorem_result", {}).get("strict_core_selector_closure_scope") or "")
        == "projective_ray_state"
    )
    kernel_alone_qw2191_discharge_still_false = (
        bool(n675.get("theorem_result", {}).get("QW2191_kernel_alone_discharge"))
        is False
        and bool(n676.get("theorem_result", {}).get("QW2191_kernel_alone_discharge"))
        is False
    )

    checks_spec = [
        {
            "id": "t172_projective_closure_discharged",
            "actual": t172_projective_discharged,
            "expected": True,
            "meaning": "Projective (ray-level) global selector closure observable is discharged on C_v1 (N674).",
        },
        {
            "id": "directed_raw_output_obstruction_boundary_positive",
            "actual": directed_raw_obstructed_without_sign_lift,
            "expected": True,
            "meaning": "Raw directed closure outputs still require an explicit sign-lift/section choice (N675 boundary).",
        },
        {
            "id": "t172_directed_closure_discharged_premise_based",
            "actual": t172_directed_discharged,
            "expected": True,
            "meaning": "Directed (vector-level) global selector closure object is discharged in explicit premise-based scope (N678).",
        },
        {
            "id": "admissible_s_sel_int_source_object_exported_in_f34_sense",
            "actual": admissible_source_object_exported,
            "expected": True,
            "meaning": "An admissible strict-core source object for S_sel_int is exported in the sense of the F34 source-object contract (N676).",
        },
        {
            "id": "strict_core_selector_closure_projective_discharged",
            "actual": strict_core_selector_closure_projective_discharged,
            "expected": True,
            "meaning": "Projective strict-core selector closure is now discharged post-T172 (N680).",
        },
        {
            "id": "kernel_alone_global_qw2191_discharge_still_unclaimed",
            "actual": kernel_alone_qw2191_discharge_still_false,
            "expected": True,
            "meaning": "None of the above exports implies a kernel-alone/global discharge of QW-2191.",
        },
    ]

    checks: list[dict[str, Any]] = []
    mismatches: list[str] = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
                "meaning": item["meaning"],
            }
        )
        if not ok:
            mismatches.append(item["id"])

    discharged = len(mismatches) == 0
    status = (
        "N679_DISCHARGED_CURRENT_STRICT_T172_STRICT_CORE_SELECTOR_CLOSURE_FRONTIER_BOUNDARY_THEOREM_NO_FALSE_PASS"
        if discharged
        else "N679_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_T172_FRONTIER_STATE"
    )

    summary = {
        "step": "N679",
        "status": status,
        "scope": "current_strict_t172_remaining_frontier_boundary_only",
        "checks": checks,
        "blocking_mismatches": mismatches,
        "theorem_result": {
            "discharged": discharged,
            "t172_projective_closure_discharged": t172_projective_discharged,
            "t172_directed_closure_discharged_premise_based": t172_directed_discharged,
            "directed_raw_outputs_obstructed_without_explicit_sign_lift": directed_raw_obstructed_without_sign_lift,
            "admissible_s_sel_int_source_object_exported_in_f34_sense": admissible_source_object_exported,
            "strict_core_selector_closure": bool(strict_core_selector_closure_projective_discharged),
            "strict_core_selector_closure_scope": "projective_ray_state"
            if strict_core_selector_closure_projective_discharged
            else None,
            "QW2191_kernel_alone_discharge": False,
            "ToE_closure": False,
            "remaining_open_branch": "kernel_alone_global_QW2191_discharge_and_any_directed_physical_orientation_datum",
        },
        "hard_limits": [
            "no_directed_sign_sensitive_physical_orientation_claim",
            "no_global_kernel_alone_QW2191_discharge",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
