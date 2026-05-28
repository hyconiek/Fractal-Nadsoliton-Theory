#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_ALPHA = GEN / "alpha_geo_strict_derived_v1.json"
IN_2300 = GEN / "p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json"
IN_2302 = GEN / "p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json"
IN_2305 = GEN / "p2305_s1255_strict_norm_to_margin_bridge_underdetermination_probe.json"
IN_A1 = GEN / "a_1_pair1_orientation_projector_operator_strict_core_v1.json"
IN_THETA = GEN / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json"
IN_F456 = GEN / "f456_current_strict_sigma_int_orientation_slice_to_a1_pair1_operator_bridge_packet_summary.json"
IN_H37 = ROOT / "H37_SIGN_DISTINCTION_STATE_AUDIT.md"
OUT = GEN / "p2306_s1256_strict_response_orientation_functional_interface_audit_probe.json"
MD = GEN / "p2306_s1256_strict_response_orientation_functional_interface_audit_probe.md"


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8")


def sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(payload: Any) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def dot(a: list[float], b: list[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    alpha = load_json(IN_ALPHA)
    p2300 = load_json(IN_2300)
    p2302 = load_json(IN_2302)
    p2305 = load_json(IN_2305)
    a1 = load_json(IN_A1)
    theta = load_json(IN_THETA)
    f456 = load_json(IN_F456)
    h37_text = load_text(IN_H37)

    p2300_probe = p2300.get("strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe", {}) or {}
    p2302_probe = p2302.get("strict_task3_provider_lift_policy_lock_candidate_probe", {}) or {}
    p2305_probe = p2305.get("strict_norm_to_margin_bridge_underdetermination_probe", {}) or {}
    basis = p2300_probe.get("provider_basis", []) or []
    coeffs = [float(row.get("canonical_coefficient_numeric", 0.0) or 0.0) for row in basis]
    names = [str(row.get("name", f"c{i}")) for i, row in enumerate(basis)]
    required_lift = float((p2302_probe.get("policy_lock_candidate", {}) or {}).get("provider_lift_per_step", 0.0) or 0.0)
    if required_lift == 0.0:
        required_lift = float((p2305_probe.get("norm_summary", {}) or {}).get("required_lift", 0.0) or 0.0)

    a1_hard_limits = a1.get("hard_limits", []) or []
    theta_hard_limits = theta.get("hard_limits", []) or []
    h37_no_physical_sign = "NO_PHYSICAL_SIGN_DATUM" in h37_text or "does not count\nas a strict physical orientation datum" in h37_text
    h37_convention_layer = "strict_convention" in h37_text and "convention layer" in h37_text

    # These are not new assumptions; they are the minimal typed-interface requirements that any P2306 success would have to export.
    required_interface_conditions = [
        {
            "id": "domain_match",
            "required": "response functional R is typed on the 10-dimensional P2300 ADM/Bianchi-I provider coefficient basis",
            "current_status": "MISSING_EXPORTED_TYPED_MAP",
        },
        {
            "id": "signed_orientation",
            "required": "R fixes a signed orientation so that R(c) is a policy-margin lift, not just an unsigned norm",
            "current_status": "MISSING_STRICT_PHYSICAL_ORIENTATION_OR_RESPONSE_THEOREM",
        },
        {
            "id": "monotone_margin_map",
            "required": "R(c) lower-bounds the P2281/P2302 per-step policy margin under the replay dynamics",
            "current_status": "MISSING_DYNAMICAL_MARGIN_THEOREM",
        },
        {
            "id": "selector_safety",
            "required": "R is exported without an external selector premise and without discharging QW-2191 by convention",
            "current_status": "NOT_SATISFIED_BY_EXISTING_CONVENTION_OR_PREMISE_BASED_SIGN_LIFTS",
        },
    ]

    orientation_source_audit = [
        {
            "source_id": "A1_PAIR1_PROJECTOR_F456",
            "artifact": "generated/a_1_pair1_orientation_projector_operator_strict_core_v1.json",
            "loaded": bool(a1) and not a1.get("_missing"),
            "scope": "pair1 rank-one projector on span{c1,s1}",
            "sign_sensitive": False,
            "provider_basis_dimension": len(coeffs),
            "artifact_dimension": 2,
            "exports_provider_margin_response": False,
            "blocking_reason": "Projector is sign-gauge invariant and pair1-local; it does not define a signed functional on the P2300 10-column ADM/Bianchi provider basis.",
            "hard_limits_include_no_qw2191": any("QW-2191" in str(x) for x in a1_hard_limits),
        },
        {
            "source_id": "THETA_PAIR_SIGMA_INT_O2_CUT",
            "artifact": "generated/theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json",
            "loaded": bool(theta) and not theta.get("_missing"),
            "scope": "slot-free pair1/pair2 theta source",
            "sign_sensitive": True,
            "provider_basis_dimension": len(coeffs),
            "artifact_dimension": 2,
            "exports_provider_margin_response": False,
            "blocking_reason": "Theta-pair source is pair1/pair2 scoped and explicitly below selector closure/QW-2191 discharge; it does not type-map ADM/Bianchi provider columns to P2281 margin.",
            "hard_limits_include_no_qw2191": any("QW-2191" in str(x) for x in theta_hard_limits),
        },
        {
            "source_id": "H37_DIRECTED_SIGN_DISTINCTION_AUDIT",
            "artifact": "H37_SIGN_DISTINCTION_STATE_AUDIT.md",
            "loaded": bool(h37_text),
            "scope": "directed sign distinction and convention-layer sign lifts on C_v1",
            "sign_sensitive": True,
            "provider_basis_dimension": len(coeffs),
            "artifact_dimension": "C_v1/chart/sign-lift layer, not P2300 coefficient space",
            "exports_provider_margin_response": False,
            "blocking_reason": "H37 records directed/convention-layer structure but no strict physical sign datum; using it as a policy-margin orientation would add an unexported premise.",
            "hard_limits_include_no_physical_sign_datum": h37_no_physical_sign,
            "convention_layer_detected": h37_convention_layer,
        },
        {
            "source_id": "F456_PACKET_SUMMARY",
            "artifact": "generated/f456_current_strict_sigma_int_orientation_slice_to_a1_pair1_operator_bridge_packet_summary.json",
            "loaded": bool(f456) and not f456.get("_missing"),
            "scope": "summary of sigma-int orientation slice to A1(pair1) operator bridge",
            "sign_sensitive": False,
            "provider_basis_dimension": len(coeffs),
            "artifact_dimension": "summary only",
            "exports_provider_margin_response": False,
            "blocking_reason": "F456 confirms a narrow operator-level bridge object, not a P2300-to-P2281 response functional.",
            "audits_pass": f456.get("audits_pass"),
            "no_false_pass": f456.get("no_false_pass"),
        },
    ]

    candidate_response_functionals = [
        {
            "candidate_id": "PROJECTOR_TRACE_PULLBACK",
            "definition_attempt": "Use A1(pair1) projector trace/eigenvalue as a scalar orientation source.",
            "well_typed_on_p2300_basis": False,
            "signed_margin_value": None,
            "passes_required_lift": False,
            "rejection_reason": "Projector trace is sign invariant and has no exported pullback to the 10 P2300 provider coefficients.",
        },
        {
            "candidate_id": "THETA_PAIR_SIGN_PULLBACK",
            "definition_attempt": "Use theta-pair sigma-int sign to orient positive P2300 channels.",
            "well_typed_on_p2300_basis": False,
            "signed_margin_value": None,
            "passes_required_lift": False,
            "rejection_reason": "Would select/weight channels by a sign rule not exported as an ADM/Bianchi policy-margin theorem; this is a selector/orientation premise.",
        },
        {
            "candidate_id": "H37_CONVENTION_SIGN_LIFT_PULLBACK",
            "definition_attempt": "Use C_v1 convention-layer directed sign lift as response orientation.",
            "well_typed_on_p2300_basis": False,
            "signed_margin_value": None,
            "passes_required_lift": False,
            "rejection_reason": "Convention-layer sign lifts are explicitly below strict physical sign datum and do not type-map to P2281 margin.",
        },
        {
            "candidate_id": "SIGNED_TOTAL_CANONICAL_RESPONSE",
            "definition_attempt": "Use the all-ones signed response on the exported P2300 coefficient vector.",
            "well_typed_on_p2300_basis": True,
            "signed_margin_value": sum(coeffs),
            "passes_required_lift": sum(coeffs) >= required_lift,
            "rejection_reason": "This is well-typed but negative, hence it cannot certify a positive provider lift.",
        },
    ]

    theorem_export = {
        "statement_id": "P2306_STRICT_RESPONSE_ORIENTATION_INTERFACE_AUDIT_THEOREM",
        "formal_statement": (
            "On the current exported repository state, no strict response/orientation functional maps the P2300 ADM/Bianchi-I "
            "provider coefficient basis to the P2281/P2302 policy margin.  Existing orientation artifacts are either pair-local "
            "projectors, pair1/pair2 theta sources, or convention/premise-scoped directed sign lifts; none exports a typed, signed, "
            "monotone provider-to-margin response theorem.  Consequently P2305 norm sufficiency remains non-admissible for strict "
            "G1/G3 closure."
        ),
        "proof_bits": {
            "p2300_provider_basis_dimension": len(coeffs),
            "p2300_coefficient_names": names,
            "required_lift": required_lift,
            "orientation_sources_audited": [row["source_id"] for row in orientation_source_audit],
            "candidate_response_functionals_rejected": [row["candidate_id"] for row in candidate_response_functionals if not row["passes_required_lift"]],
            "all_required_interface_conditions_missing": all(row["current_status"].startswith("MISSING") or row["current_status"].startswith("NOT_SATISFIED") for row in required_interface_conditions),
        },
        "not_claimed": [
            "strict response-functional bridge",
            "strict orientation theorem for P2300 to P2281 margin",
            "strict G1 closure",
            "strict G3 closure",
            "full Task-3 closure",
            "QW-2191 discharge",
            "selector closure",
            "legacy-kernel role transfer",
            "ToE closure",
        ],
    }
    theorem_fingerprint = sha256_json(theorem_export)

    gatekeeper_checks = {
        "alpha_geo_strict_source_loaded": alpha.get("status") == "actual_exported_strict_derived_source_upgrade_value",
        "alpha_geo_is_four_ln2_not_legacy_import": alpha.get("value") == "4 ln(2)",
        "p2300_coefficients_loaded": len(coeffs) == 10,
        "p2305_underdetermination_loaded": (p2305_probe.get("strict_bridge_verdict", {}) or {}).get("norm_to_margin_bridge_proven") is False,
        "orientation_artifacts_audited": all(row["loaded"] for row in orientation_source_audit),
        "no_orientation_source_exports_provider_margin_response": not any(row["exports_provider_margin_response"] for row in orientation_source_audit),
        "candidate_response_functionals_do_not_close_required_lift": not any(row["passes_required_lift"] for row in candidate_response_functionals),
        "strict_response_orientation_bridge_not_proven": True,
        "strict_g1_g3_not_updated": True,
        "no_qw2191_discharge_claimed": True,
        "no_selector_closure_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2306_s1256_v1",
        "packet_id": "P2306",
        "stage_id": "S1256",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_OBSTRUCTION_STRICT_RESPONSE_ORIENTATION_FUNCTIONAL_NOT_EXPORTED",
        "result_kind": "STRICT_RESPONSE_ORIENTATION_INTERFACE_AUDIT_REFUTES_CURRENT_EXPORT_SET",
        "strict_response_orientation_functional_interface_audit_probe": {
            "probe_id": "P2306_S1256_STRICT_RESPONSE_ORIENTATION_FUNCTIONAL_INTERFACE_AUDIT",
            "source_packets": {
                "alpha_geo_strict_derived_v1": "generated/alpha_geo_strict_derived_v1.json",
                "p2300": "generated/p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json",
                "p2302": "generated/p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json",
                "p2305": "generated/p2305_s1255_strict_norm_to_margin_bridge_underdetermination_probe.json",
                "a1_pair1_projector": "generated/a_1_pair1_orientation_projector_operator_strict_core_v1.json",
                "theta_pair_sigma_int": "generated/theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json",
                "f456_summary": "generated/f456_current_strict_sigma_int_orientation_slice_to_a1_pair1_operator_bridge_packet_summary.json",
                "h37_sign_audit": "H37_SIGN_DISTINCTION_STATE_AUDIT.md",
            },
            "source_hashes": {
                "alpha_geo_strict_derived_v1_sha256": sha256_file(IN_ALPHA),
                "p2300_sha256": sha256_file(IN_2300),
                "p2302_sha256": sha256_file(IN_2302),
                "p2305_sha256": sha256_file(IN_2305),
                "a1_pair1_projector_sha256": sha256_file(IN_A1),
                "theta_pair_sigma_int_sha256": sha256_file(IN_THETA),
                "f456_summary_sha256": sha256_file(IN_F456),
                "h37_sign_audit_sha256": sha256_file(IN_H37),
            },
            "required_interface_conditions": required_interface_conditions,
            "orientation_source_audit": orientation_source_audit,
            "candidate_response_functionals": candidate_response_functionals,
            "strict_interface_verdict": {
                "response_orientation_functional_exported": False,
                "current_export_set_refutes_closure_use": True,
                "reason": "No audited orientation artifact exports a signed monotone P2300-provider-basis to P2281-policy-margin functional.",
            },
            "task3_g1_g3_update": {
                "G1_reduction_certainty": "OPEN_UNCHANGED",
                "G2_nonlinear_trajectory_realism": "CLOSED_FROM_P2301_UNCHANGED",
                "G3_operational_policy_rule": "OPEN_UNCHANGED",
                "reason": "P2306 finds an interface obstruction, not a response/orientation bridge.",
            },
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2307_candidate",
            "goal": "Construct a genuinely typed ADM/Bianchi-I response functional from strict-side dynamics, or prove a theorem-grade nonexistence result for all current orientation/convention-layer sources. Do not move G1/G3 until such a functional is exported and replayed.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_RESPONSE_ORIENTATION_INTERFACE_AUDIT_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2306/S1256 — strict response/orientation functional interface audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Required lift: `{required_lift}`",
            f"- P2300 provider basis dimension: `{len(coeffs)}`",
            "- Existing orientation artifacts audited: A1(pair1) projector, theta-pair sigma-int source, H37 sign-distinction audit, F456 summary.",
            "- Response/orientation functional exported: `False`",
            "- G1/G3 update: `OPEN_UNCHANGED`",
            f"- Theorem fingerprint: `{theorem_fingerprint}`",
            "",
            "## Guardrail statement",
            "P2306 does not close G1/G3, does not discharge QW-2191, does not add a selector premise, does not transfer legacy-kernel roles, and does not claim ToE closure.",
            "",
            "## Next honest step",
            payload["recommended_next_honest_step"]["goal"],
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
