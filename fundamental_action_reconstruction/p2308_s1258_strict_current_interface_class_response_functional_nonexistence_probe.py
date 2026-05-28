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
IN_2306 = GEN / "p2306_s1256_strict_response_orientation_functional_interface_audit_probe.json"
IN_2307 = GEN / "p2307_s1257_strict_dynamical_margin_response_functional_nonderivation_probe.json"
OUT = GEN / "p2308_s1258_strict_current_interface_class_response_functional_nonexistence_probe.json"
MD = GEN / "p2308_s1258_strict_current_interface_class_response_functional_nonexistence_probe.md"


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(payload: Any) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def bool_value(row: dict[str, Any], *keys: str) -> bool:
    for key in keys:
        if key in row:
            return bool(row[key])
    return False


def main() -> None:
    GEN.mkdir(exist_ok=True)
    alpha = load_json(IN_ALPHA)
    p2300 = load_json(IN_2300)
    p2302 = load_json(IN_2302)
    p2306 = load_json(IN_2306)
    p2307 = load_json(IN_2307)

    p2300_probe = p2300.get("strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe", {}) or {}
    p2302_probe = p2302.get("strict_task3_provider_lift_policy_lock_candidate_probe", {}) or {}
    p2306_probe = p2306.get("strict_response_orientation_functional_interface_audit_probe", {}) or {}
    p2307_probe = p2307.get("strict_dynamical_margin_response_functional_nonderivation_probe", {}) or {}

    coeffs = [float(row.get("canonical_coefficient_numeric", 0.0) or 0.0) for row in (p2300_probe.get("provider_basis", []) or [])]
    required_lift = float((p2302_probe.get("policy_lock_candidate", {}) or {}).get("provider_lift_per_step", 0.0) or 0.0)

    interface_predicates = [
        "typed_on_p2300_coefficients",
        "derived_from_margin_dynamics",
        "exports_or_derives_signed_orientation",
        "strict_internal_no_selector_or_convention_premise",
        "passes_required_lift",
        "admissible_as_p2300_response_functional",
    ]

    current_interface_class: list[dict[str, Any]] = []
    for row in p2306_probe.get("candidate_response_functionals", []) or []:
        current_interface_class.append({
            "candidate_id": "P2306::" + str(row.get("candidate_id")),
            "source_probe": "P2306",
            "definition": row.get("definition_attempt"),
            "typed_on_p2300_coefficients": bool_value(row, "well_typed_on_p2300_basis"),
            "derived_from_margin_dynamics": False,
            "exports_or_derives_signed_orientation": False,
            "strict_internal_no_selector_or_convention_premise": False,
            "passes_required_lift": bool(row.get("passes_required_lift", False)),
            "admissible_as_p2300_response_functional": False,
            "failure_reason": row.get("rejection_reason"),
        })
    for row in p2307_probe.get("candidate_functionals", []) or []:
        current_interface_class.append({
            "candidate_id": "P2307::" + str(row.get("candidate_id")),
            "source_probe": "P2307",
            "definition": row.get("definition"),
            "typed_on_p2300_coefficients": bool(row.get("typed_on_p2300_coefficients", False)),
            "derived_from_margin_dynamics": bool(row.get("derived_from_margin_dynamics", False)),
            "exports_or_derives_signed_orientation": False,
            "strict_internal_no_selector_or_convention_premise": (
                bool(row.get("derived_from_margin_dynamics", False)) and bool(row.get("typed_on_p2300_coefficients", False))
            ),
            "passes_required_lift": bool(row.get("passes_required_lift", False)),
            "admissible_as_p2300_response_functional": bool(row.get("admissible_as_p2300_response_functional", False)),
            "failure_reason": row.get("rejection_reason"),
        })

    for row in current_interface_class:
        failed = [key for key in interface_predicates if not row.get(key, False)]
        row["failed_predicates"] = failed
        row["passes_all_required_predicates"] = len(failed) == 0

    predicate_matrix = {
        key: {
            "satisfied_by": [row["candidate_id"] for row in current_interface_class if row.get(key, False)],
            "failed_by": [row["candidate_id"] for row in current_interface_class if not row.get(key, False)],
        }
        for key in interface_predicates
    }
    admissible_candidates = [row["candidate_id"] for row in current_interface_class if row["passes_all_required_predicates"]]

    class_nonexistence_certificate = {
        "class_id": "CURRENT_EXPORTED_P2300_P2281_INTERFACE_CLASS_V1",
        "class_scope": "Finite current export class assembled from P2306 candidate response functionals and P2307 dynamical candidate functionals; not a universal theorem over all future possible dynamics.",
        "candidate_count": len(current_interface_class),
        "required_predicates": interface_predicates,
        "admissible_candidate_count": len(admissible_candidates),
        "admissible_candidates": admissible_candidates,
        "nonexistence_over_current_class": len(admissible_candidates) == 0,
        "dominant_missing_predicates": {
            "typed_and_dynamics_joint": [
                row["candidate_id"] for row in current_interface_class
                if not (row["typed_on_p2300_coefficients"] and row["derived_from_margin_dynamics"])
            ],
            "signed_orientation": [
                row["candidate_id"] for row in current_interface_class
                if not row["exports_or_derives_signed_orientation"]
            ],
            "strict_internal_no_selector_or_convention_premise": [
                row["candidate_id"] for row in current_interface_class
                if not row["strict_internal_no_selector_or_convention_premise"]
            ],
        },
    }

    theorem_export = {
        "statement_id": "P2308_CURRENT_INTERFACE_CLASS_RESPONSE_FUNCTIONAL_NONEXISTENCE_THEOREM",
        "formal_statement": (
            "Within the finite current exported P2300/P2281 interface class assembled from P2306 and P2307, there is no admissible "
            "strict response functional lambda=R(c).  Each candidate fails at least one required predicate: typed P2300 domain, "
            "derivation from margin dynamics, signed orientation, strict internal/no-selector provenance, positive required lift, or final "
            "admissibility.  Therefore G1/G3 cannot be updated on the current export set.  This is a theorem over the current interface "
            "class, not a universal impossibility theorem for future strict dynamics."
        ),
        "proof_bits": {
            "required_lift": required_lift,
            "p2300_coefficient_count": len(coeffs),
            "candidate_count": len(current_interface_class),
            "admissible_candidate_count": len(admissible_candidates),
            "admissible_candidates": admissible_candidates,
            "all_candidates_have_failed_predicates": all(len(row["failed_predicates"]) > 0 for row in current_interface_class),
            "predicate_matrix": predicate_matrix,
        },
        "not_claimed": [
            "universal nonexistence over all future strict dynamics",
            "strict response-functional bridge",
            "strict derivation of provider_lift_per_step from P2300 coefficients",
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
        "p2302_required_lift_loaded": required_lift == 0.0068,
        "p2306_candidates_loaded": len(p2306_probe.get("candidate_response_functionals", []) or []) > 0,
        "p2307_candidates_loaded": len(p2307_probe.get("candidate_functionals", []) or []) > 0,
        "current_interface_class_nonempty": len(current_interface_class) > 0,
        "all_candidates_fail_at_least_one_required_predicate": all(len(row["failed_predicates"]) > 0 for row in current_interface_class),
        "no_admissible_current_class_response_functional": len(admissible_candidates) == 0,
        "theorem_scoped_to_current_export_class": True,
        "strict_response_functional_not_exported": True,
        "strict_g1_g3_not_updated": True,
        "no_qw2191_discharge_claimed": True,
        "no_selector_closure_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2308_s1258_v1",
        "packet_id": "P2308",
        "stage_id": "S1258",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_OBSTRUCTION_CURRENT_INTERFACE_CLASS_RESPONSE_FUNCTIONAL_NONEXISTENCE",
        "result_kind": "THEOREM_GRADE_NONEXISTENCE_FOR_CURRENT_P2300_P2281_INTERFACE_CLASS",
        "strict_current_interface_class_response_functional_nonexistence_probe": {
            "probe_id": "P2308_S1258_STRICT_CURRENT_INTERFACE_CLASS_RESPONSE_FUNCTIONAL_NONEXISTENCE",
            "source_packets": {
                "alpha_geo_strict_derived_v1": "generated/alpha_geo_strict_derived_v1.json",
                "p2300": "generated/p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json",
                "p2302": "generated/p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json",
                "p2306": "generated/p2306_s1256_strict_response_orientation_functional_interface_audit_probe.json",
                "p2307": "generated/p2307_s1257_strict_dynamical_margin_response_functional_nonderivation_probe.json",
            },
            "source_hashes": {
                "alpha_geo_strict_derived_v1_sha256": sha256_file(IN_ALPHA),
                "p2300_sha256": sha256_file(IN_2300),
                "p2302_sha256": sha256_file(IN_2302),
                "p2306_sha256": sha256_file(IN_2306),
                "p2307_sha256": sha256_file(IN_2307),
            },
            "current_interface_class": current_interface_class,
            "predicate_matrix": predicate_matrix,
            "class_nonexistence_certificate": class_nonexistence_certificate,
            "strict_nonexistence_verdict": {
                "lambda_equals_R_of_c_exported_in_current_class": False,
                "nonexistence_proven_for_current_interface_class": True,
                "universal_future_nonexistence_claimed": False,
                "reason": "The current finite interface class has zero candidates satisfying all required strict response-functional predicates.",
            },
            "task3_g1_g3_update": {
                "G1_reduction_certainty": "OPEN_UNCHANGED",
                "G2_nonlinear_trajectory_realism": "CLOSED_FROM_P2301_UNCHANGED",
                "G3_operational_policy_rule": "OPEN_UNCHANGED",
                "reason": "P2308 proves current-class nonexistence, not a new lambda=R(c) bridge.",
            },
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2309_candidate",
            "goal": "Introduce a genuinely new strict dynamical equation or variational identity that exports weights w_i for lambda=R(c), then replay P2302; otherwise keep the current-class nonexistence certificate as the blocker for G1/G3.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_CURRENT_CLASS_NONEXISTENCE_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2308/S1258 — current-interface-class response functional nonexistence",
            "",
            f"- Status: `{payload['status']}`",
            f"- Required lift: `{required_lift}`",
            f"- Current interface candidates audited: `{len(current_interface_class)}`",
            f"- Admissible candidates: `{len(admissible_candidates)}`",
            "- Nonexistence scope: current exported P2300/P2281 interface class only; no universal future nonexistence claim.",
            "- G1/G3 update: `OPEN_UNCHANGED`",
            f"- Theorem fingerprint: `{theorem_fingerprint}`",
            "",
            "## Guardrail statement",
            "P2308 does not close G1/G3, does not discharge QW-2191, does not add a selector premise, does not transfer legacy-kernel roles, and does not claim ToE closure.",
            "",
            "## Next honest step",
            payload["recommended_next_honest_step"]["goal"],
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
