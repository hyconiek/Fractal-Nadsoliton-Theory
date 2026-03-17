#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

AS_OF = "2026-03-17"

IN_N700 = GENERATED / "n700_current_first_strict_projective_operational_toe_os_support_with_emergent_observer_chain_theorem_summary.json"
OUT = GENERATED / "n701_current_strict_projective_operational_toe_closure_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_N700.exists():
        summary = {
            "step": "N701",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_N700.relative_to(REPO))],
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    n700 = load_json(IN_N700)
    tr = (n700.get("theorem_result") or {}) if isinstance(n700, dict) else {}

    ok = (
        n700.get("no_false_pass") is True
        and n700.get("step") == "N700"
        and n700.get("status")
        == "N700_DERIVABLE_CURRENT_FIRST_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_V2_THEOREM_NO_FALSE_PASS"
        and tr.get("discharged") is True
        and tr.get("strict_projective_operational_toe_os_support_packet_v2_exported") is True
        and tr.get("kernel_alone_qw2191_discharge") is False
        and tr.get("ToE_closure") is False
        and tr.get("actual_emergent_observer_closure") is False
    )

    status = "N701_DERIVABLE_CURRENT_STRICT_PROJECTIVE_OPERATIONAL_TOE_CLOSURE_THEOREM_NO_FALSE_PASS"
    if not ok:
        status = "N701_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_OPERATIONAL_TOE_CLOSURE_STATE"

    theorem_result = {
        "discharged": bool(ok),
        "operational_toe_closure_strict_projective_os": bool(ok),
        "strict_projective_operational_toe_os_support_packet_v2_exported": bool(
            tr.get("strict_projective_operational_toe_os_support_packet_v2_exported")
        ),
        "kernel_alone_qw2191_discharge": False,
        "ToE_closure": False,
        "standard_model_host_matching": False,
        "actual_emergent_observer_closure": False,
        "evidence": {"N700": str(IN_N700.relative_to(REPO))},
    }

    summary = {
        "step": "N701",
        "status": status,
        "as_of": AS_OF,
        "scope": "strict_projective_operational_toe_closure_only",
        "theorem_result": theorem_result,
        "hard_limits": [
            "no_kernel_alone_global_QW2191_discharge",
            "no_directed_sign_sensitive_physical_orientation_claim",
            "no_standard_model_host_matching_claim",
            "no_ToE_closure",
            "no_actual_emergent_observer_closure",
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

