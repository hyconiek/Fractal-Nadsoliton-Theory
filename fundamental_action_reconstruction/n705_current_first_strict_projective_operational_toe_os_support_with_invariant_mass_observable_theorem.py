#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

AS_OF = "2026-03-17"

IN_P705 = GENERATED / "p705_current_actual_strict_projective_operational_toe_os_support_with_invariant_mass_observable_probe_summary.json"
OUT = (
    GENERATED
    / "n705_current_first_strict_projective_operational_toe_os_support_with_invariant_mass_observable_theorem_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_P705.exists():
        summary = {
            "step": "N705",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_P705.relative_to(REPO))],
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    p705 = load_json(IN_P705)

    ok = (
        p705.get("no_false_pass") is True
        and p705.get("status")
        == "CURRENT_REPO_EXPORTS_ONE_ACTUAL_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_V3_PACKET_AFTER_P705"
        and p705.get("strict_projective_operational_toe_os_support_packet_v3_exported") is True
        and p705.get("actual_strict_core_toe_closure_discharged") is False
        and p705.get("actual_global_toe_closure_discharged") is False
        and p705.get("actual_emergent_observer_closure") is False
    )

    status = (
        "N705_DERIVABLE_CURRENT_FIRST_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_V3_WITH_INVARIANT_MASS_OBSERVABLE_THEOREM_NO_FALSE_PASS"
    )
    if not ok:
        status = "N705_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_OPERATIONAL_TOE_OS_SUPPORT_V3_STATE"

    theorem_result = {
        "discharged": bool(ok),
        "strict_projective_operational_toe_os_support_packet_v3_exported": bool(
            p705.get("strict_projective_operational_toe_os_support_packet_v3_exported")
        ),
        "kernel_alone_qw2191_discharge": False,
        "ToE_closure": False,
        "actual_emergent_observer_closure": False,
        "evidence": {"P705": str(IN_P705.relative_to(REPO))},
    }

    summary = {
        "step": "N705",
        "status": status,
        "as_of": AS_OF,
        "scope": "strict_projective_operational_toe_os_support_v3_only",
        "theorem_result": theorem_result,
        "hard_limits": [
            "no_kernel_alone_global_QW2191_discharge",
            "no_directed_sign_sensitive_physical_orientation_claim",
            "no_ToE_closure",
            "no_actual_emergent_observer_closure",
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

