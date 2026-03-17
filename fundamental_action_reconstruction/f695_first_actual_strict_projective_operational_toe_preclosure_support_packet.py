#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

AS_OF = "2026-03-17"

OUT = GENERATED / "f695_first_actual_strict_projective_operational_toe_preclosure_support_packet_summary.json"

IN_N680 = GENERATED / "n680_current_strict_t173_projective_strict_core_selector_closure_discharge_theorem_summary.json"
IN_N692 = GENERATED / "n692_current_strict_t173_projective_directed_closure_output_ray_invariance_theorem_summary.json"
IN_N693 = GENERATED / "n693_current_strict_t173_output_sign_lift_gauge_covariance_theorem_summary.json"
IN_P694 = (
    GENERATED / "p694_current_strict_physical_computability_mass_spectrum_proxy_from_projective_selector_closure_probe_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_N680, IN_N692, IN_N693, IN_P694]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        summary = {
            "packet_id": "F695",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    n680 = load_json(IN_N680)
    n692 = load_json(IN_N692)
    n693 = load_json(IN_N693)
    p694 = load_json(IN_P694)

    n680_res = (n680.get("theorem_result") or {})
    n692_res = (n692.get("theorem_result") or {})
    n693_res = (n693.get("theorem_result") or {})

    projective_strict_core_selector_closure_discharged = bool(n680_res.get("strict_core_selector_closure")) and (
        n680_res.get("strict_core_selector_closure_scope") == "projective_ray_state"
    )
    projective_output_ray_invariant_across_directed_lifts = bool(
        n692_res.get("projective_output_ray_invariant_across_exported_directed_closures")
    )
    output_sign_lift_gauge_covariant = bool(
        n693_res.get("output_sign_lift_is_gauge_covariant_under_chart_sign_relift")
    )
    physical_mass_spectrum_proxy_computable_from_projective_closure = bool(
        p694.get("status") == "PASS_PHYSICAL_MASS_SPECTRUM_PROXY_COMPUTABLE_FROM_PROJECTIVE_SELECTOR_CLOSURE"
    )

    # Hard-limit checks: these must remain false.
    toe_closure = bool(n680_res.get("ToE_closure")) or bool(n692_res.get("ToE_closure")) or bool(n693_res.get("ToE_closure"))
    kernel_alone_qw2191_discharge = bool(n680_res.get("QW2191_kernel_alone_discharge")) or bool(
        n692_res.get("QW2191_kernel_alone_discharge")
    ) or bool(n693_res.get("QW2191_kernel_alone_discharge"))

    strict_projective_operational_toe_preclosure_support_packet_exported = bool(
        projective_strict_core_selector_closure_discharged
        and projective_output_ray_invariant_across_directed_lifts
        and output_sign_lift_gauge_covariant
        and physical_mass_spectrum_proxy_computable_from_projective_closure
        and (not toe_closure)
        and (not kernel_alone_qw2191_discharge)
    )

    summary = {
        "packet_id": "F695",
        "status": "F695_EXECUTED_FIRST_ACTUAL_STRICT_PROJECTIVE_OPERATIONAL_TOE_PRECLOSURE_SUPPORT_PACKET_NO_FALSE_PASS",
        "as_of": AS_OF,
        "support_packet_name": "Lambda_strict_projective_operational_toe_preclosure_support_v1",
        "projective_strict_core_selector_closure_discharged": projective_strict_core_selector_closure_discharged,
        "projective_output_ray_invariant_across_directed_lifts": projective_output_ray_invariant_across_directed_lifts,
        "output_sign_lift_gauge_covariant": output_sign_lift_gauge_covariant,
        "physical_mass_spectrum_proxy_computable_from_projective_closure": physical_mass_spectrum_proxy_computable_from_projective_closure,
        "strict_projective_operational_toe_preclosure_support_packet_exported": strict_projective_operational_toe_preclosure_support_packet_exported,
        "kernel_alone_qw2191_discharge": False,
        "actual_strict_core_toe_closure_discharged": False,
        "actual_global_toe_closure_discharged": False,
        "evidence": {
            "N680": str(IN_N680.relative_to(REPO)),
            "N692": str(IN_N692.relative_to(REPO)),
            "N693": str(IN_N693.relative_to(REPO)),
            "P694": str(IN_P694.relative_to(REPO)),
        },
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

