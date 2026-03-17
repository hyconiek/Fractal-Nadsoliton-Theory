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
IN_N676 = GENERATED / "n676_current_first_admissible_s_sel_int_source_object_discharge_theorem_summary.json"
IN_N546 = (
    GENERATED
    / "n546_current_exported_s_sel_int_strict_core_source_object_admissible_orientation_export_theorem_summary.json"
)
IN_P680 = (
    GENERATED
    / "p680_current_strict_t173_projective_strict_core_selector_closure_candidate_audit_probe_summary.json"
)

OUT = (
    GENERATED
    / "n680_current_strict_t173_projective_strict_core_selector_closure_discharge_theorem_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    prereq = (IN_N674, IN_N675, IN_N676, IN_N546, IN_P680)
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        summary = {
            "step": "N680",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES_FOR_T173_PROJECTIVE_STRICT_CORE_SELECTOR_CLOSURE_DISCHARGE_THEOREM",
            "scope": "current_strict_t173_projective_strict_core_selector_closure_only",
            "missing": missing,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    n674 = load_json(IN_N674)
    n675 = load_json(IN_N675)
    n676 = load_json(IN_N676)
    n546 = load_json(IN_N546)
    p680 = load_json(IN_P680)

    t172_projective_discharged = bool(
        (n674.get("theorem_result", {}) or {}).get("discharged") is True
    )
    admissible_s_sel_int_source_object = bool(
        (n676.get("theorem_result", {}) or {}).get("admissible_S_sel_int_source_object_in_F34_sense")
        is True
    )
    admissible_e_orient = bool(
        (n546.get("theorem_result", {}) or {}).get("admissible_E_orient") is True
    )
    directed_raw_outputs_obstructed_without_sign_lift = bool(
        (n675.get("theorem_result", {}) or {}).get(
            "directed_global_closure_obstructed_without_explicit_sign_lift"
        )
        is True
    )
    closure_candidate_supported = bool(p680.get("closure_candidate_supported") is True) and str(
        p680.get("status") or ""
    ).startswith("PASS_")

    checks_spec = [
        {
            "id": "t172_projective_closure_object_discharged",
            "actual": t172_projective_discharged,
            "expected": True,
            "meaning": "The global projective selector closure object on C_v1 is discharged (N674 packages F672/N672/N673).",
        },
        {
            "id": "admissible_s_sel_int_source_object_in_f34_sense",
            "actual": admissible_s_sel_int_source_object,
            "expected": True,
            "meaning": "An admissible strict-core source object for the S_sel_int step exists in the narrow F34 sense (N676).",
        },
        {
            "id": "admissible_e_orient_export_present",
            "actual": admissible_e_orient,
            "expected": True,
            "meaning": "An admissible strict-core orientation export exists from the S_sel_int source object (N546).",
        },
        {
            "id": "directed_raw_outputs_still_obstructed_without_sign_lift",
            "actual": directed_raw_outputs_obstructed_without_sign_lift,
            "expected": True,
            "meaning": "Directed raw outputs remain obstructed without explicit sign lift (N675 boundary), so no directed closure is promoted here.",
        },
        {
            "id": "p680_projective_closure_candidate_audit_pass",
            "actual": closure_candidate_supported,
            "expected": True,
            "meaning": "Probe-level audit supports a projective strict-core selector closure discharge attempt (P680).",
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
        "N680_DISCHARGED_CURRENT_STRICT_T173_PROJECTIVE_STRICT_CORE_SELECTOR_CLOSURE_DISCHARGE_THEOREM_NO_FALSE_PASS"
        if discharged
        else "N680_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_T173_PROJECTIVE_STRICT_CORE_CLOSURE_STATE"
    )

    summary = {
        "step": "N680",
        "status": status,
        "scope": "current_strict_t173_projective_strict_core_selector_closure_only",
        "checks": checks,
        "blocking_mismatches": mismatches,
        "theorem_result": {
            "discharged": discharged,
            "strict_core_selector_closure": bool(discharged),
            "strict_core_selector_closure_scope": "projective_ray_state",
            "QW2191_kernel_alone_discharge": False,
            "directed_sign_sensitive_physical_orientation": False,
            "ToE_closure": False,
            "evidence": {
                "N674": str(IN_N674.relative_to(REPO)),
                "N676": str(IN_N676.relative_to(REPO)),
                "N546": str(IN_N546.relative_to(REPO)),
                "N675": str(IN_N675.relative_to(REPO)),
                "P680": str(IN_P680.relative_to(REPO)),
            },
            "remaining_open_branch": "kernel_alone_global_QW2191_discharge_and_any_directed_physical_orientation_datum",
        },
        "hard_limits": [
            "no_directed_sign_sensitive_physical_orientation_claim",
            "no_kernel_alone_global_QW2191_discharge",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
