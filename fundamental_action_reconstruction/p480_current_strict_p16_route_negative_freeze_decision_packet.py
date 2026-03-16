#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

P475_SUMMARY = (
    GENERATED
    / "p475_current_strict_projective_only_continuation_decision_packet_summary.json"
)
P16_SUMMARY = (
    GENERATED
    / "p16_existing_kernel_feedback_legacy_chart_reduced_operator_export_probe_summary.json"
)

N520_SUMMARY = (
    GENERATED
    / "n520_current_first_strict_r18_declared_pair1_residual_zero_equations_value_instance_obstruction_theorem_summary.json"
)
N521_SUMMARY = (
    GENERATED
    / "n521_current_first_strict_t169_rordpow_sign_scan_r18_pair1_zero_equations_obstruction_theorem_summary.json"
)
N522_SUMMARY = (
    GENERATED
    / "n522_current_first_strict_reference_magnitude_family_sign_scan_r18_pair1_zero_equations_obstruction_theorem_summary.json"
)

DIRECT_FORMAL_FRONTIER_CANDIDATES: list[tuple[str, Path]] = [
    (
        "P629",
        GENERATED
        / "p629_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_t169_constrained_lift_g4_g6_family_shift_defect_zero_witness_packet_summary.json",
    ),
    (
        "P628",
        GENERATED
        / "p628_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi2_source_psi5_target_psi8_source_psi11_target_coherence_instances_summary.json",
    ),
    (
        "P627",
        GENERATED
        / "p627_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi2_psi5_and_psi8_psi11_slotwise_role_split_and_defect_polynomial_packets_summary.json",
    ),
    (
        "P626",
        GENERATED
        / "p626_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi7_source_eom_coherence_instance_summary.json",
    ),
    (
        "P625",
        GENERATED
        / "p625_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi10_target_coherence_instance_summary.json",
    ),
    (
        "P624",
        GENERATED
        / "p624_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi10_target_role_split_and_defect_polynomial_packets_summary.json",
    ),
    (
        "P623",
        GENERATED
        / "p623_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi2_psi5_and_psi8_psi11_common_plus3_assignment_slot_split_packets_summary.json",
    ),
    (
        "P622",
        GENERATED
        / "p622_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi2_psi5_and_psi8_psi11_role_matching_packets_summary.json",
    ),
    (
        "P61",
        GENERATED
        / "p61_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi7_source_eom_coefficient_defect_polynomial_packet_summary.json",
    ),
    (
        "P46",
        GENERATED
        / "p46_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_target_action_coefficient_defect_polynomial_packet_summary.json",
    ),
]

OUT = GENERATED / "p480_current_strict_p16_route_negative_freeze_decision_packet.json"
OUT_SUMMARY = (
    GENERATED / "p480_current_strict_p16_route_negative_freeze_decision_packet_summary.json"
)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _p(path: Path) -> str:
    try:
        return str(path.relative_to(REPO))
    except ValueError:
        return str(path)


def _pick_direct_formal_frontier() -> dict[str, Any]:
    for label, path in DIRECT_FORMAL_FRONTIER_CANDIDATES:
        if path.exists():
            try:
                obj = _read_json(path)
                missing = obj.get("remaining_missing_upstream_objects") or obj.get(
                    "remaining_missing_objects"
                )
                return {
                    "target": label,
                    "summary_path": _p(path),
                    "status": obj.get("status"),
                    "remaining_missing": missing,
                }
            except Exception:
                return {"target": label, "summary_path": _p(path), "parse_error": True}
    return {
        "target": "P46",
        "summary_path": None,
        "note": "No direct-formal summary file found; falling back to the root direct-formal frontier label.",
    }


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = {"P475_summary": P475_SUMMARY, "P16_summary": P16_SUMMARY}
    missing_required = [k for k, p in required.items() if not p.is_file()]
    if missing_required:
        payload = {
            "stage": "P480",
            "date": datetime.now(timezone.utc).date().isoformat(),
            "goal": "declare_professorial_negative_freeze_for_P16_route__proceed_on_direct_formal_route_under_projective_only_continuation",
            "status": "FAIL_MISSING_REQUIRED_INPUTS",
            "missing_required_inputs": missing_required,
            "required_paths": {k: _p(v) for k, v in required.items()},
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(
            json.dumps(
                {
                    "stage": payload["stage"],
                    "status": payload["status"],
                    "decision": None,
                    "recommended_next_strict_target": "P16",
                    "no_false_pass": True,
                },
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_SUMMARY)
        return

    p475 = _read_json(P475_SUMMARY)
    projective_only = (p475.get("selected_continuation") == "projective_only") or (
        str(p475.get("decision") or "") == "PROJECTIVE_ONLY_CONTINUATION_SELECTED"
    )
    if not projective_only:
        payload = {
            "stage": "P480",
            "date": datetime.now(timezone.utc).date().isoformat(),
            "goal": "declare_professorial_negative_freeze_for_P16_route__proceed_on_direct_formal_route_under_projective_only_continuation",
            "status": "FAIL_PROJECTIVE_ONLY_CONTINUATION_NOT_SELECTED",
            "p475_summary": _p(P475_SUMMARY),
            "p475_selected_continuation": p475.get("selected_continuation"),
            "decision": None,
            "no_false_pass": True,
            "hard_limits": [
                "Decision applies only under projective-only continuation; does not claim directed sign lift.",
            ],
        }
        OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(
            json.dumps(
                {
                    "stage": payload["stage"],
                    "status": payload["status"],
                    "decision": None,
                    "recommended_next_strict_target": "H37",
                    "no_false_pass": True,
                },
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_SUMMARY)
        return

    p16 = _read_json(P16_SUMMARY)
    direct_formal = _pick_direct_formal_frontier()

    evidence: dict[str, Any] = {"P16_summary": _p(P16_SUMMARY)}
    for name, path in (("N520", N520_SUMMARY), ("N521", N521_SUMMARY), ("N522", N522_SUMMARY)):
        if path.exists():
            try:
                obj = _read_json(path)
                evidence[name] = {"summary_path": _p(path), "status": obj.get("status")}
            except Exception:
                evidence[name] = {"summary_path": _p(path), "parse_error": True}
        else:
            evidence[name] = {"summary_path": _p(path), "missing": True}

    status = "PASS_DECISION_DECLARED_P16_ROUTE_NEGATIVE_FREEZE_SELECTED"
    decision = "P16_ROUTE_NEGATIVE_FREEZE_SELECTED"

    payload = {
        "stage": "P480",
        "date": datetime.now(timezone.utc).date().isoformat(),
        "goal": "declare_professorial_negative_freeze_for_P16_route__proceed_on_direct_formal_route_under_projective_only_continuation",
        "status": status,
        "decision": decision,
        "decision_basis": {
            "continuation_premise": {"projective_only_selected": True, "p475_summary": _p(P475_SUMMARY)},
            "p16_frontier_status": {
                "p16_status": p16.get("status"),
                "remaining_missing_upstream_objects": p16.get("remaining_missing_upstream_objects"),
                "required_next_step": p16.get("required_next_step"),
            },
            "evidence_chain": evidence,
            "strict_boundary_note": (
                "Freeze P16 as negative on current strict core: do not promote any coefficient-filled legacy chart-reduced operator export on pair1 "
                "while the residual zero/cancellation witness and QW-2191 canonicalization remain missing; proceed on strict-only closure lanes."
            ),
        },
        "continuation": {
            "selected": "direct_formal_c1s1_family_route",
            "meaning": "Proceed on the kernel-split-robust canonical-ontology-supported direct formal c1s1 family route (F3 priority), compatible with projective-only selector state semantics.",
        },
        "recommended_next_strict_target": direct_formal,
        "hard_limits": [
            "no theorem-level pass",
            "no full-closure pass",
            "does not discharge P16",
            "does not claim any strict zero/cancellation witness for R16–R18",
            "does not claim global QW-2191 discharge",
            "does not claim ToE closure",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": payload["stage"],
        "status": payload["status"],
        "decision": payload["decision"],
        "recommended_next_strict_target": direct_formal.get("target"),
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

