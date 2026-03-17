#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P692 = (
    GENERATED
    / "p692_current_strict_t173_projective_vs_directed_closure_output_compatibility_audit_probe_summary.json"
)
IN_F692 = GENERATED / "f692_current_strict_t175_sign_fixed_directed_closure_export_packet_summary.json"

OUT = (
    GENERATED
    / "n692_current_strict_t173_projective_directed_closure_output_ray_invariance_theorem_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    missing: list[str] = []
    if not IN_P692.exists():
        missing.append(str(IN_P692.relative_to(REPO)))
    if not IN_F692.exists():
        missing.append(str(IN_F692.relative_to(REPO)))

    if missing:
        summary = {
            "step": "N692",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES_FOR_PROJECTIVE_DIRECTED_CLOSURE_OUTPUT_RAY_INVARIANCE",
            "scope": "current_strict_t173_projective_output_ray_invariance_under_directed_lifts",
            "missing": missing,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    p692 = load_json(IN_P692)
    f692 = load_json(IN_F692)

    p692_pass = str(p692.get("status") or "").startswith("PASS_")
    f692_pass = str(f692.get("status") or "").startswith("PASS_")
    all_projectors_match = bool(p692.get("all_projectors_match_projective") is True)
    pairwise_consistent = bool(p692.get("pairwise_projector_consistent") is True)

    checks_spec = [
        {
            "id": "f692_export_passed",
            "actual": f692_pass,
            "expected": True,
            "meaning": "A sign-fixed directed closure object exists on C_v1 in strict_convention scope (F692).",
        },
        {
            "id": "p692_audit_passed",
            "actual": p692_pass,
            "expected": True,
            "meaning": "Projective vs directed closure output compatibility audit passed (P692).",
        },
        {
            "id": "all_directed_projectors_match_projective_output_projector",
            "actual": all_projectors_match,
            "expected": True,
            "meaning": "Each directed closure output vector induces the same rank‑1 output projector as the projective closure output projector (P692).",
        },
        {
            "id": "directed_projectors_pairwise_consistent",
            "actual": pairwise_consistent,
            "expected": True,
            "meaning": "Directed closure output projectors are pairwise consistent across exported directed scopes (P692).",
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
        "N692_DERIVABLE_CURRENT_STRICT_T173_PROJECTIVE_DIRECTED_CLOSURE_OUTPUT_RAY_INVARIANCE_THEOREM_NO_FALSE_PASS"
        if discharged
        else "N692_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_CLOSURE_OUTPUT_COMPATIBILITY"
    )

    summary = {
        "step": "N692",
        "status": status,
        "scope": "current_strict_t173_projective_output_ray_invariance_under_directed_lifts",
        "checks": checks,
        "blocking_mismatches": mismatches,
        "theorem_result": {
            "discharged": discharged,
            "projective_output_ray_invariant_across_exported_directed_closures": discharged,
            "directed_sign_sensitive_physical_orientation_in_strict_core": False,
            "QW2191_kernel_alone_discharge": False,
            "ToE_closure": False,
            "evidence": {
                "F692": str(IN_F692.relative_to(REPO)),
                "P692": str(IN_P692.relative_to(REPO)),
            },
            "note": (
                "This is a scope-limited invariance statement: directed closure outputs (in exported convention/premise scopes) "
                "collapse to the same projective output ray on Q_out. It does not upgrade any sign choice into strict physical orientation "
                "and does not claim kernel-alone/global QW-2191 discharge."
            ),
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

