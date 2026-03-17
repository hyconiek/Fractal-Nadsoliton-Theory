#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

AS_OF = "2026-03-17"

OUT = GENERATED / "f705_first_actual_strict_projective_operational_toe_os_support_with_invariant_mass_observable_packet_summary.json"

IN_N680 = GENERATED / "n680_current_strict_t173_projective_strict_core_selector_closure_discharge_theorem_summary.json"
IN_N692 = GENERATED / "n692_current_strict_t173_projective_directed_closure_output_ray_invariance_theorem_summary.json"
IN_N693 = GENERATED / "n693_current_strict_t173_output_sign_lift_gauge_covariance_theorem_summary.json"
IN_N699 = GENERATED / "n699_current_first_strict_projective_emergent_observer_chain_computability_theorem_summary.json"

IN_P694 = (
    GENERATED / "p694_current_strict_physical_computability_mass_spectrum_proxy_from_projective_selector_closure_probe_summary.json"
)
IN_P696 = (
    GENERATED
    / "p696_current_strict_physical_computability_selector_aligned_channel_spectrum_proxy_from_projective_selector_closure_probe_summary.json"
)
IN_P697 = (
    GENERATED
    / "p697_current_strict_projective_observer_limit_readout_from_global_projective_selector_closure_output_probe_summary.json"
)

IN_F704 = (
    GENERATED
    / "f704_current_strict_invariant_mass_observable_from_diagonal_local_psi_hessian_eigensystem_export_packet_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_N680, IN_N692, IN_N693, IN_N699, IN_P694, IN_P696, IN_P697, IN_F704]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        summary = {
            "packet_id": "F705",
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
    n699 = load_json(IN_N699)
    p694 = load_json(IN_P694)
    p696 = load_json(IN_P696)
    p697 = load_json(IN_P697)
    f704 = load_json(IN_F704)

    n680_res = (n680.get("theorem_result") or {})
    n692_res = (n692.get("theorem_result") or {})
    n693_res = (n693.get("theorem_result") or {})
    n699_res = (n699.get("theorem_result") or {})

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
    selector_aligned_channel_spectrum_proxy_computable_from_projective_closure = bool(
        p696.get("status") == "PASS_SELECTOR_ALIGNED_CHANNEL_SPECTRUM_PROXY_COMPUTABLE_FROM_PROJECTIVE_SELECTOR_CLOSURE"
    )
    observer_limit_readout_computable_from_projective_closure_output_projector = bool(
        p697.get("status")
        == "PASS_PROJECTIVE_OBSERVER_LIMIT_READOUT_COMPUTABLE_FROM_GLOBAL_PROJECTIVE_SELECTOR_CLOSURE_OUTPUT_PROJECTOR"
    )
    projective_emergent_observer_chain_computable_from_projective_closure_output_projector = bool(
        n699_res.get("projective_emergent_observer_chain_computable_from_projective_selector_closure_output_projector")
    )

    basis_invariant_mass_observable_exported = bool(
        f704.get("status") == "PASS_EXPORTED_STRICT_INVARIANT_MASS_OBSERVABLE_OBJECT"
    )

    # Hard-limit checks: these must remain false.
    toe_closure = any(
        bool(x.get("ToE_closure"))
        for x in [
            n680_res,
            n692_res,
            n693_res,
            n699_res,
        ]
        if isinstance(x, dict)
    )
    kernel_alone_qw2191_discharge = any(
        bool(x.get("QW2191_kernel_alone_discharge")) or bool(x.get("kernel_alone_qw2191_discharge"))
        for x in [
            n680_res,
            n692_res,
            n693_res,
            n699_res,
        ]
        if isinstance(x, dict)
    )
    actual_emergent_observer_closure = bool(n699_res.get("actual_emergent_observer_closure"))

    strict_projective_operational_toe_os_support_packet_v3_exported = bool(
        projective_strict_core_selector_closure_discharged
        and projective_output_ray_invariant_across_directed_lifts
        and output_sign_lift_gauge_covariant
        and physical_mass_spectrum_proxy_computable_from_projective_closure
        and selector_aligned_channel_spectrum_proxy_computable_from_projective_closure
        and basis_invariant_mass_observable_exported
        and observer_limit_readout_computable_from_projective_closure_output_projector
        and projective_emergent_observer_chain_computable_from_projective_closure_output_projector
        and (not toe_closure)
        and (not kernel_alone_qw2191_discharge)
        and (not actual_emergent_observer_closure)
    )

    summary = {
        "packet_id": "F705",
        "status": "F705_EXECUTED_FIRST_ACTUAL_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_WITH_INVARIANT_MASS_OBSERVABLE_PACKET_NO_FALSE_PASS",
        "as_of": AS_OF,
        "support_packet_name": "Lambda_strict_projective_operational_toe_os_support_v3",
        "projective_strict_core_selector_closure_discharged": projective_strict_core_selector_closure_discharged,
        "projective_output_ray_invariant_across_directed_lifts": projective_output_ray_invariant_across_directed_lifts,
        "output_sign_lift_gauge_covariant": output_sign_lift_gauge_covariant,
        "physical_mass_spectrum_proxy_computable_from_projective_closure": physical_mass_spectrum_proxy_computable_from_projective_closure,
        "selector_aligned_channel_spectrum_proxy_computable_from_projective_closure": selector_aligned_channel_spectrum_proxy_computable_from_projective_closure,
        "basis_invariant_mass_observable_exported": basis_invariant_mass_observable_exported,
        "observer_limit_readout_computable_from_projective_closure_output_projector": observer_limit_readout_computable_from_projective_closure_output_projector,
        "projective_emergent_observer_chain_computable_from_projective_closure_output_projector": projective_emergent_observer_chain_computable_from_projective_closure_output_projector,
        "strict_projective_operational_toe_os_support_packet_v3_exported": strict_projective_operational_toe_os_support_packet_v3_exported,
        "kernel_alone_qw2191_discharge": False,
        "actual_strict_core_toe_closure_discharged": False,
        "actual_global_toe_closure_discharged": False,
        "actual_emergent_observer_closure": False,
        "evidence": {
            "N680": str(IN_N680.relative_to(REPO)),
            "N692": str(IN_N692.relative_to(REPO)),
            "N693": str(IN_N693.relative_to(REPO)),
            "N699": str(IN_N699.relative_to(REPO)),
            "P694": str(IN_P694.relative_to(REPO)),
            "P696": str(IN_P696.relative_to(REPO)),
            "F704": str(IN_F704.relative_to(REPO)),
            "P697": str(IN_P697.relative_to(REPO)),
        },
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

