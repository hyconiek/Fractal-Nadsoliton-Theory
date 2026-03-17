#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_N680 = GENERATED / "n680_current_strict_t173_projective_strict_core_selector_closure_discharge_theorem_summary.json"
IN_N675 = GENERATED / "n675_current_strict_global_directed_selector_closure_obstruction_boundary_theorem_summary.json"
IN_P681 = GENERATED / "p681_current_strict_t173_reflection_breaking_sign_lift_from_s_sel_int_audit_probe_summary.json"

OUT = (
    GENERATED
    / "n681_current_strict_t173_directed_output_sign_lift_obstruction_boundary_theorem_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    prereq = (IN_N680, IN_N675, IN_P681)
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        summary = {
            "step": "N681",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES_FOR_T173_DIRECTED_OUTPUT_SIGN_LIFT_OBSTRUCTION_BOUNDARY",
            "scope": "current_strict_t173_directed_output_sign_lift_boundary_only",
            "missing": missing,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    n680 = load_json(IN_N680)
    n675 = load_json(IN_N675)
    p681 = load_json(IN_P681)

    projective_strict_core_closure_discharged = bool(
        (n680.get("theorem_result") or {}).get("discharged") is True
        and bool((n680.get("theorem_result") or {}).get("strict_core_selector_closure")) is True
        and str((n680.get("theorem_result") or {}).get("strict_core_selector_closure_scope") or "")
        == "projective_ray_state"
    )

    raw_directed_output_obstructed_without_sign_lift = bool(
        (n675.get("theorem_result") or {}).get("directed_global_closure_obstructed_without_explicit_sign_lift")
        is True
    )

    deterministic_sign_lift_candidate_available_on_all_charts = bool(
        (p681.get("all_pairs_supported") is True)
    )
    deterministic_sign_lift_candidate_output_sign_consistent = bool(
        (p681.get("directed_output_sign_consistent_if_all_pairs_supported") is True)
    )

    checks_spec = [
        {
            "id": "projective_strict_core_selector_closure_discharged",
            "actual": projective_strict_core_closure_discharged,
            "expected": True,
            "meaning": "Projective strict-core selector closure is discharged (N680).",
        },
        {
            "id": "raw_directed_output_obstruction_boundary_positive",
            "actual": raw_directed_output_obstructed_without_sign_lift,
            "expected": True,
            "meaning": "Raw directed closure outputs remain obstructed without explicit sign lift (N675 boundary).",
        },
        {
            "id": "deterministic_sign_lift_candidate_exists_from_exported_weights",
            "actual": deterministic_sign_lift_candidate_available_on_all_charts,
            "expected": True,
            "meaning": "A deterministic per-chart sign-lift candidate can be defined from exported S_sel_int seed payload weights (P681).",
        },
        {
            "id": "deterministic_sign_lift_candidate_fails_to_make_directed_output_chart_independent",
            "actual": deterministic_sign_lift_candidate_output_sign_consistent,
            "expected": False,
            "meaning": "Even under that deterministic sign-lift candidate, the directed output sign is not globally chart-independent (P681).",
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
        "N681_DISCHARGED_CURRENT_STRICT_T173_DIRECTED_OUTPUT_SIGN_LIFT_OBSTRUCTION_BOUNDARY_THEOREM_NO_FALSE_PASS"
        if discharged
        else "N681_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_T173_DIRECTED_OUTPUT_SIGN_LIFT_BOUNDARY_STATE"
    )

    summary = {
        "step": "N681",
        "status": status,
        "scope": "current_strict_t173_directed_output_sign_lift_boundary_only",
        "checks": checks,
        "blocking_mismatches": mismatches,
        "theorem_result": {
            "discharged": discharged,
            "projective_strict_core_selector_closure": True,
            "directed_output_sign_lift_determined_in_strict_core": False,
            "directed_sign_sensitive_physical_orientation_in_strict_core": False,
            "QW2191_kernel_alone_discharge": False,
            "ToE_closure": False,
            "evidence": {
                "N680": str(IN_N680.relative_to(REPO)),
                "N675": str(IN_N675.relative_to(REPO)),
                "P681": str(IN_P681.relative_to(REPO)),
            },
            "remaining_open_branch": "kernel_alone_global_QW2191_discharge_and_any_directed_physical_orientation_datum",
        },
        "hard_limits": [
            "no_kernel_alone_global_QW2191_discharge",
            "no_directed_sign_sensitive_physical_orientation_claim",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

