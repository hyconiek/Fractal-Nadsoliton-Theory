#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

AS_OF = "2026-03-17"

IN_P698 = GENERATED / "p698_current_actual_strict_projective_operational_toe_os_support_probe_summary.json"
OUT = GENERATED / "n698_current_first_strict_projective_operational_toe_os_support_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_P698.exists():
        summary = {
            "step": "N698",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_P698.relative_to(REPO))],
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    p698 = load_json(IN_P698)

    ok = (
        p698.get("no_false_pass") is True
        and p698.get("status")
        == "CURRENT_REPO_EXPORTS_ONE_ACTUAL_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_PACKET_AFTER_P698"
        and p698.get("strict_projective_operational_toe_os_support_packet_exported") is True
        and p698.get("actual_strict_core_toe_closure_discharged") is False
        and p698.get("actual_global_toe_closure_discharged") is False
    )

    status = "N698_DERIVABLE_CURRENT_FIRST_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_THEOREM_NO_FALSE_PASS"
    if not ok:
        status = "N698_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_OPERATIONAL_TOE_OS_SUPPORT_STATE"

    theorem_result = {
        "discharged": bool(ok),
        "strict_projective_operational_toe_os_support_packet_exported": bool(
            p698.get("strict_projective_operational_toe_os_support_packet_exported")
        ),
        "kernel_alone_qw2191_discharge": False,
        "ToE_closure": False,
        "actual_emergent_observer_closure": False,
        "evidence": {
            "P698": str(IN_P698.relative_to(REPO)),
        },
    }

    summary = {
        "step": "N698",
        "status": status,
        "as_of": AS_OF,
        "scope": "strict_projective_operational_toe_os_support_only",
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

