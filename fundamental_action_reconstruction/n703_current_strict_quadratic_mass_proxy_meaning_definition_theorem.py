#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

AS_OF = "2026-03-17"

IN_F459 = GENERATED / "psi_hessian_diagonal_local_strict_derived_value_instantiated_v1_summary.json"
IN_P694 = GENERATED / "p694_current_strict_physical_computability_mass_spectrum_proxy_from_projective_selector_closure_probe_summary.json"
IN_P696 = (
    GENERATED
    / "p696_current_strict_physical_computability_selector_aligned_channel_spectrum_proxy_from_projective_selector_closure_probe_summary.json"
)

OUT = GENERATED / "n703_current_strict_quadratic_mass_proxy_meaning_definition_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F459, IN_P694, IN_P696]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        summary = {
            "step": "N703",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    f459 = load_json(IN_F459)
    p694 = load_json(IN_P694)
    p696 = load_json(IN_P696)

    ok = (
        f459.get("no_false_pass") is True
        and f459.get("stage") == "F459"
        and isinstance(f459.get("outputs"), dict)
        and p694.get("no_false_pass") is True
        and p694.get("stage") == "P694"
        and str(p694.get("status", "")).startswith("PASS_")
        and p696.get("no_false_pass") is True
        and p696.get("stage") == "P696"
        and str(p696.get("status", "")).startswith("PASS_")
    )

    status = "N703_DERIVABLE_CURRENT_STRICT_QUADRATIC_MASS_PROXY_MEANING_DEFINITION_THEOREM_NO_FALSE_PASS"
    if not ok:
        status = "N703_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_MASS_PROXY_MEANING_STATE"

    theorem_result = {
        "discharged": bool(ok),
        "quadratic_mass_proxy_meaning_defined": bool(ok),
        "physical_unit_identification": False,
        "standard_model_host_matching": False,
        "kernel_alone_qw2191_discharge": False,
        "directed_sign_sensitive_physical_orientation_claim": False,
        "evidence": {
            "F459": str(IN_F459.relative_to(REPO)),
            "P694": str(IN_P694.relative_to(REPO)),
            "P696": str(IN_P696.relative_to(REPO)),
        },
    }

    summary = {
        "step": "N703",
        "status": status,
        "as_of": AS_OF,
        "scope": "strict_quadratic_proxy_meaning_only",
        "theorem_result": theorem_result,
        "hard_limits": [
            "no_physical_unit_map",
            "no_standard_model_host_matching_claim",
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

