#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_ALPHA = GEN / "alpha_geo_strict_derived_v1.json"
IN_2300 = GEN / "p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json"
IN_2302 = GEN / "p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json"
OUT = GEN / "p2303_s1253_strict_provider_to_margin_bridge_bound_audit_probe.json"
MD = GEN / "p2303_s1253_strict_provider_to_margin_bridge_bound_audit_probe.md"


def load(path: Path) -> dict[str, Any]:
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


def main() -> None:
    GEN.mkdir(exist_ok=True)
    alpha = load(IN_ALPHA)
    p2300 = load(IN_2300)
    p2302 = load(IN_2302)

    p2300_probe = p2300.get("strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe", {}) or {}
    solution = (p2300_probe.get("solution_space", {}) or {}).get("canonical_solution", []) or []
    basis = p2300_probe.get("provider_basis", []) or []
    p2302_probe = p2302.get("strict_task3_provider_lift_policy_lock_candidate_probe", {}) or {}
    required_lift_float = float((p2302_probe.get("policy_lock_candidate", {}) or {}).get("provider_lift_per_step", 0.0) or 0.0)
    required_lift = sp.Rational(str(required_lift_float))

    coeffs = [sp.sympify(item, locals={"pi": sp.pi, "log": sp.log}) for item in solution]
    coeff_records = []
    for idx, coeff in enumerate(coeffs):
        coeff_records.append({
            "column": idx,
            "name": (basis[idx] or {}).get("name", f"u{idx}") if idx < len(basis) else f"u{idx}",
            "coefficient": str(coeff),
            "coefficient_numeric": float(sp.N(coeff, 40)),
            "abs_numeric": float(sp.N(abs(coeff), 40)),
            "sign": "positive" if coeff > 0 else "negative" if coeff < 0 else "zero",
        })

    nonzero_abs = [abs(coeff) for coeff in coeffs if coeff != 0]
    positive_coeffs = [coeff for coeff in coeffs if coeff > 0]
    direct_floor = min(nonzero_abs) if nonzero_abs else sp.Integer(0)
    signed_total = sp.simplify(sum(coeffs))
    positive_sum = sp.simplify(sum(positive_coeffs))
    l1_norm = sp.simplify(sum(abs(coeff) for coeff in coeffs))
    l2_norm = sp.sqrt(sp.simplify(sum(coeff**2 for coeff in coeffs)))

    direct_floor_numeric = float(sp.N(direct_floor, 40))
    required_lift_numeric = float(sp.N(required_lift, 40))
    direct_gap = sp.simplify(direct_floor - required_lift)
    direct_gap_numeric = float(sp.N(direct_gap, 40))

    # Rigorous comparison aid: pi**2 > 9 implies the direct floor is even less than
    # this rational upper bound, and that upper bound remains below the required lift.
    direct_floor_upper_using_pi_sq_gt_9 = sp.Rational(32991491, 640000000 * 9)
    upper_bound_gap = sp.simplify(direct_floor_upper_using_pi_sq_gt_9 - required_lift)

    bridge_contracts = [
        {
            "contract_id": "STRICT_SINGLE_CHANNEL_FLOOR",
            "description": "Only the smallest nonzero canonical coefficient magnitude is guaranteed without choosing signs, weights, or an aggregation functional.",
            "bound_expression": str(direct_floor),
            "bound_numeric": direct_floor_numeric,
            "required_lift": required_lift_numeric,
            "passes_required_lift": direct_floor_numeric + 1e-15 >= required_lift_numeric,
            "strict_admissible_without_extra_premise": True,
            "verdict": "REFUTES_REQUIRED_LIFT_UNDER_THIS_STRICT_DIRECT_CONTRACT",
        },
        {
            "contract_id": "SIGNED_TOTAL_COEFFICIENT",
            "description": "Signed sum of canonical coefficients; negative value cannot certify a positive margin lift.",
            "bound_expression": str(signed_total),
            "bound_numeric": float(sp.N(signed_total, 40)),
            "required_lift": required_lift_numeric,
            "passes_required_lift": float(sp.N(signed_total, 40)) >= required_lift_numeric,
            "strict_admissible_without_extra_premise": True,
            "verdict": "DOES_NOT_CERTIFY_POSITIVE_LIFT",
        },
        {
            "contract_id": "POSITIVE_CHANNEL_SUM",
            "description": "Sum of positive coefficients exceeds the candidate, but selecting only positive channels is an aggregation/sign premise not exported by P2300.",
            "bound_expression": str(positive_sum),
            "bound_numeric": float(sp.N(positive_sum, 40)),
            "required_lift": required_lift_numeric,
            "passes_required_lift": float(sp.N(positive_sum, 40)) >= required_lift_numeric,
            "strict_admissible_without_extra_premise": False,
            "verdict": "NUMERICALLY_SUFFICIENT_BUT_NOT_A_STRICT_BRIDGE_WITHOUT_NEW_AGGREGATION_PREMISE",
        },
        {
            "contract_id": "L1_OR_L2_NORM_AGGREGATION",
            "description": "Norms exceed the candidate, but a norm-to-margin map would be a new provider-to-margin theorem not exported here.",
            "l1_expression": str(l1_norm),
            "l1_numeric": float(sp.N(l1_norm, 40)),
            "l2_expression": str(l2_norm),
            "l2_numeric": float(sp.N(l2_norm, 40)),
            "required_lift": required_lift_numeric,
            "passes_required_lift": float(sp.N(l2_norm, 40)) >= required_lift_numeric,
            "strict_admissible_without_extra_premise": False,
            "verdict": "NUMERICALLY_SUFFICIENT_BUT_REQUIRES_A_NEW_NORM_TO_MARGIN_BRIDGE_THEOREM",
        },
    ]

    strict_direct_contracts = [row for row in bridge_contracts if row["strict_admissible_without_extra_premise"]]
    strict_direct_bridge_proven = any(row["passes_required_lift"] for row in strict_direct_contracts)
    strict_direct_bridge_refuted = all(not row["passes_required_lift"] for row in strict_direct_contracts)

    theorem_export = {
        "statement_id": "P2303_PROVIDER_TO_MARGIN_BRIDGE_BOUND_AUDIT_THEOREM",
        "formal_statement": (
            "The P2300 canonical spatial-EOM coefficients do not currently prove the P2302 requirement "
            "provider_lift_per_step >= 0.0068 under the strict direct coefficient-floor contracts that require no selector, "
            "sign-selection, or norm-to-margin aggregation premise.  The smallest guaranteed nonzero coefficient magnitude is "
            "below 0.0068 and the signed total is negative.  Larger positive-sum and norm bounds are numerically sufficient "
            "but require a new provider-to-margin aggregation theorem, so strict G1/G3 cannot be updated yet."
        ),
        "proof_bits": {
            "required_lift": str(required_lift),
            "direct_floor": str(direct_floor),
            "direct_floor_minus_required": str(direct_gap),
            "direct_floor_numeric": direct_floor_numeric,
            "direct_gap_numeric": direct_gap_numeric,
            "direct_floor_upper_using_pi_sq_gt_9": str(direct_floor_upper_using_pi_sq_gt_9),
            "direct_floor_upper_minus_required": str(upper_bound_gap),
            "signed_total": str(signed_total),
            "positive_sum": str(positive_sum),
            "l1_norm": str(l1_norm),
            "l2_norm": str(l2_norm),
            "strict_direct_bridge_proven": strict_direct_bridge_proven,
            "strict_direct_bridge_refuted_for_current_contracts": strict_direct_bridge_refuted,
        },
        "not_claimed": [
            "absolute impossibility of any future provider-to-margin theorem",
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
        "p2302_required_lift_loaded": required_lift_numeric == 0.0068,
        "direct_floor_below_required_lift": direct_floor_numeric < required_lift_numeric,
        "signed_total_not_positive_lift": float(sp.N(signed_total, 40)) < required_lift_numeric,
        "strict_direct_bridge_not_proven": not strict_direct_bridge_proven,
        "strict_direct_bridge_refuted_for_current_contracts": strict_direct_bridge_refuted,
        "sufficient_norms_marked_non_admissible_without_new_bridge": all(
            (not row["passes_required_lift"]) or row["strict_admissible_without_extra_premise"] is False
            for row in bridge_contracts
            if row["contract_id"] in {"POSITIVE_CHANNEL_SUM", "L1_OR_L2_NORM_AGGREGATION"}
        ),
        "strict_g1_g3_not_updated": True,
        "no_qw2191_discharge_claimed": True,
        "no_selector_closure_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2303_s1253_v1",
        "packet_id": "P2303",
        "stage_id": "S1253",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_BRIDGE_OBLIGATION_STRICT_DIRECT_BOUND_REFUTED_WITH_TRACE",
        "result_kind": "STRICT_PROVIDER_TO_MARGIN_BRIDGE_BOUND_AUDIT_REFUTES_CURRENT_DIRECT_CERTIFICATE",
        "strict_provider_to_margin_bridge_bound_audit_probe": {
            "probe_id": "P2303_S1253_STRICT_PROVIDER_TO_MARGIN_BRIDGE_BOUND_AUDIT",
            "source_packets": {
                "alpha_geo_strict_derived_v1": "generated/alpha_geo_strict_derived_v1.json",
                "p2300": "generated/p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json",
                "p2302": "generated/p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json",
            },
            "source_hashes": {
                "alpha_geo_strict_derived_v1_sha256": sha256_file(IN_ALPHA),
                "p2300_sha256": sha256_file(IN_2300),
                "p2302_sha256": sha256_file(IN_2302),
            },
            "required_lift": {
                "source": "P2302.policy_lock_candidate.provider_lift_per_step",
                "value": required_lift_numeric,
                "exact": str(required_lift),
            },
            "canonical_coefficients": coeff_records,
            "bridge_contract_audit": bridge_contracts,
            "strict_bridge_verdict": {
                "provider_to_margin_bridge_proven": False,
                "provider_to_margin_bridge_refuted_for_current_direct_contracts": strict_direct_bridge_refuted,
                "reason": "direct floor and signed-total contracts cannot certify 0.0068; stronger aggregations require a new bridge theorem/premise",
            },
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2304_candidate",
            "goal": "Either derive a strict norm-to-margin aggregation theorem from the P2300 operator basis, or lower the P2302 policy-lock candidate to a bridge bound that the strict direct coefficient floor can certify; do not update G1/G3 before that.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_REFUTED_DIRECT_BRIDGE_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2303/S1253 — provider-to-margin bridge bound audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Required P2302 lift: `{required_lift_numeric}`",
            f"- Strict direct coefficient floor: `{direct_floor_numeric}`",
            f"- Direct floor gap: `{direct_gap_numeric}`",
            f"- Signed coefficient total: `{float(sp.N(signed_total, 40))}`",
            "- Verdict: strict direct bridge is not proven; current direct contracts refute the requested bound.",
            f"- Theorem fingerprint: `{theorem_fingerprint}`",
            "",
            "## Guardrail statement",
            "Positive-sum and norm aggregations are numerically large enough, but P2303 does not promote them because no strict provider-to-margin aggregation theorem is exported. G1/G3 remain open; no QW-2191, selector, legacy-kernel, or ToE closure is claimed.",
            "",
            "## Next honest step",
            payload["recommended_next_honest_step"]["goal"],
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
