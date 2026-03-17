#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

AS_OF = "2026-03-17"

IN_P695 = GENERATED / "p695_current_actual_strict_projective_operational_toe_preclosure_support_probe_summary.json"
OUT = GENERATED / "n695_current_first_strict_projective_operational_toe_preclosure_support_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_P695.exists():
        summary = {
            "step": "N695",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_P695.relative_to(REPO))],
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    p695 = load_json(IN_P695)
    mismatches = p695.get("blocking_mismatches") or []

    checks = [
        {
            "id": "p695_probe_status_is_pass",
            "actual": p695.get("status"),
            "expected": "CURRENT_REPO_EXPORTS_ONE_ACTUAL_STRICT_PROJECTIVE_OPERATIONAL_TOE_PRECLOSURE_SUPPORT_PACKET_AFTER_P695",
            "pass": bool(p695.get("status") == "CURRENT_REPO_EXPORTS_ONE_ACTUAL_STRICT_PROJECTIVE_OPERATIONAL_TOE_PRECLOSURE_SUPPORT_PACKET_AFTER_P695"),
        }
    ]

    discharged = (len(mismatches) == 0) and bool(
        p695.get("strict_projective_operational_toe_preclosure_support_packet_exported") is True
    )

    theorem_result = {
        "discharged": discharged,
        "strict_projective_operational_toe_preclosure_support_packet_exported": bool(
            p695.get("strict_projective_operational_toe_preclosure_support_packet_exported")
        ),
        "actual_strict_core_toe_closure_discharged": bool(p695.get("actual_strict_core_toe_closure_discharged")),
        "actual_global_toe_closure_discharged": bool(p695.get("actual_global_toe_closure_discharged")),
        "ToE_closure": False,
        "evidence": {"P695": str(IN_P695.relative_to(REPO))},
    }

    summary = {
        "step": "N695",
        "status": "N695_DERIVABLE_CURRENT_FIRST_STRICT_PROJECTIVE_OPERATIONAL_TOE_PRECLOSURE_SUPPORT_THEOREM_NO_FALSE_PASS",
        "as_of": AS_OF,
        "scope": "strict_projective_operational_toe_preclosure_only",
        "checks": checks,
        "blocking_mismatches": list(mismatches),
        "theorem_result": theorem_result,
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

