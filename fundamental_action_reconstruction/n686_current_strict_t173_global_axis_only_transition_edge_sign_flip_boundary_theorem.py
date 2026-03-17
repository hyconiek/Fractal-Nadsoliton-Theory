#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P686 = (
    GENERATED
    / "p686_current_strict_t173_w_break_rooted_directed_state_full_transition_edge_compatibility_audit_probe_summary.json"
)

OUT = (
    GENERATED
    / "n686_current_strict_t173_global_axis_only_transition_edge_sign_flip_boundary_theorem_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    missing = [str(IN_P686.relative_to(REPO))] if not IN_P686.exists() else []
    if missing:
        summary = {
            "step": "N686",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES_FOR_T173_GLOBAL_EDGE_SIGN_FLIP_BOUNDARY",
            "scope": "current_strict_t173_global_axis_only_transition_edge_sign_boundary_only",
            "missing": missing,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    p686 = load_json(IN_P686)

    p686_pass = str(p686.get("status") or "").startswith("PASS_")
    all_edges_up_to_sign = bool(p686.get("all_edges_compatible_up_to_sign") is True)
    sign_flip_count = int(p686.get("sign_flip_count") or 0)
    sign_flips_present = sign_flip_count > 0

    checks_spec = [
        {
            "id": "p686_probe_passed",
            "actual": p686_pass,
            "expected": True,
            "meaning": "The full-edge compatibility probe ran and passed at least up to sign (P686).",
        },
        {
            "id": "all_edges_compatible_up_to_sign",
            "actual": all_edges_up_to_sign,
            "expected": True,
            "meaning": "Every exported overlap edge transports the directed state line correctly: O_ij u_i ≈ ± u_j (P686).",
        },
        {
            "id": "sign_flips_present_on_some_edges",
            "actual": sign_flips_present,
            "expected": True,
            "meaning": "Some overlap edges force a sign flip under axis-only (α mod π) transport representatives (P686).",
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
        "N686_DERIVABLE_CURRENT_STRICT_T173_GLOBAL_AXIS_ONLY_TRANSITION_EDGE_SIGN_FLIP_BOUNDARY_THEOREM_NO_FALSE_PASS"
        if discharged
        else "N686_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_T173_GLOBAL_EDGE_SIGN_FLIP_BOUNDARY_STATE"
    )

    summary = {
        "step": "N686",
        "status": status,
        "scope": "current_strict_t173_global_axis_only_transition_edge_sign_boundary_only",
        "checks": checks,
        "blocking_mismatches": mismatches,
        "theorem_result": {
            "discharged": discharged,
            "global_axis_only_transition_edge_sign_flips_present": sign_flips_present,
            "sign_flip_count": sign_flip_count,
            "directed_sign_sensitive_physical_orientation_in_strict_core": False,
            "QW2191_kernel_alone_discharge": False,
            "ToE_closure": False,
            "evidence": {
                "P686": str(IN_P686.relative_to(REPO)),
            },
            "note": "This boundary is about edgewise sign coherence under axis-only (α mod π) global transition representatives; it does not forbid convention-scoped oriented lifts.",
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

