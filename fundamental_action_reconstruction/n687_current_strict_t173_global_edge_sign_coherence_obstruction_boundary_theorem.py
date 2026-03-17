#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P687 = (
    GENERATED
    / "p687_current_strict_t173_global_edge_sign_coherence_solvability_audit_probe_summary.json"
)

OUT = (
    GENERATED
    / "n687_current_strict_t173_global_edge_sign_coherence_obstruction_boundary_theorem_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    missing = [str(IN_P687.relative_to(REPO))] if not IN_P687.exists() else []
    if missing:
        summary = {
            "step": "N687",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES_FOR_T173_GLOBAL_EDGE_SIGN_COHERENCE_OBSTRUCTION_BOUNDARY",
            "scope": "current_strict_t173_global_edge_sign_coherence_boundary_only",
            "missing": missing,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    p687 = load_json(IN_P687)

    p687_pass = str(p687.get("status") or "").startswith("PASS_")
    solvable = bool(p687.get("sign_system_solvable") is True)
    triangle_witness_present = bool(p687.get("triangle_witness_present") is True)

    checks_spec = [
        {
            "id": "p687_probe_passed",
            "actual": p687_pass,
            "expected": True,
            "meaning": "P687 computed a global edge sign-coherence solvability verdict.",
        },
        {
            "id": "sign_system_not_solvable",
            "actual": solvable,
            "expected": False,
            "meaning": "No per-chart Z2 sign relift solves the full-edge sign coherence system under fixed exported axis-only transition representatives (P687).",
        },
        {
            "id": "triangle_witness_present",
            "actual": triangle_witness_present,
            "expected": True,
            "meaning": "P687 provides an explicit negative 3-cycle witness (triangle product = -1).",
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
        "N687_DERIVABLE_CURRENT_STRICT_T173_GLOBAL_EDGE_SIGN_COHERENCE_OBSTRUCTION_BOUNDARY_THEOREM_NO_FALSE_PASS"
        if discharged
        else "N687_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_T173_GLOBAL_EDGE_SIGN_COHERENCE_BOUNDARY_STATE"
    )

    summary = {
        "step": "N687",
        "status": status,
        "scope": "current_strict_t173_global_edge_sign_coherence_boundary_only",
        "checks": checks,
        "blocking_mismatches": mismatches,
        "theorem_result": {
            "discharged": discharged,
            "global_edge_sign_coherence_solvable_by_chart_sign_relift": False,
            "directed_sign_sensitive_physical_orientation_in_strict_core": False,
            "QW2191_kernel_alone_discharge": False,
            "ToE_closure": False,
            "evidence": {
                "P687": str(IN_P687.relative_to(REPO)),
            },
            "note": "This boundary is about the nonexistence of a global per-chart sign relift solving the full-edge sign coherence constraints under fixed axis-only (α mod π) transition representatives.",
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

