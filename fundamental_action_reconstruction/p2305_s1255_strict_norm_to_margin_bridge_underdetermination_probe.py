#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_ALPHA = GEN / "alpha_geo_strict_derived_v1.json"
IN_2300 = GEN / "p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json"
IN_2302 = GEN / "p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json"
IN_2303 = GEN / "p2303_s1253_strict_provider_to_margin_bridge_bound_audit_probe.json"
IN_2304 = GEN / "p2304_s1254_strict_direct_floor_policy_lock_replay_probe.json"
OUT = GEN / "p2305_s1255_strict_norm_to_margin_bridge_underdetermination_probe.json"
MD = GEN / "p2305_s1255_strict_norm_to_margin_bridge_underdetermination_probe.md"


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


def dot(a: list[float], b: list[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def norm2(v: list[float]) -> float:
    return math.sqrt(dot(v, v))


def scale(v: list[float], s: float) -> list[float]:
    return [s * x for x in v]


def unit_orthogonal(v: list[float]) -> list[float]:
    n2 = dot(v, v)
    if n2 <= 0.0:
        return [1.0] + [0.0] * (len(v) - 1)
    # Start from a coordinate vector least aligned with v, then subtract the projection.
    idx = min(range(len(v)), key=lambda i: abs(v[i]))
    e = [0.0] * len(v)
    e[idx] = 1.0
    proj = dot(e, v) / n2
    w = [ei - proj * vi for ei, vi in zip(e, v)]
    nw = norm2(w)
    if nw <= 1e-15:
        # Degenerate fallback: use another coordinate; P2300 has dimension 10, so this path should not be used.
        j = (idx + 1) % len(v)
        e = [0.0] * len(v)
        e[j] = 1.0
        proj = dot(e, v) / n2
        w = [ei - proj * vi for ei, vi in zip(e, v)]
        nw = norm2(w)
    return scale(w, 1.0 / nw)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    alpha = load(IN_ALPHA)
    p2300 = load(IN_2300)
    p2302 = load(IN_2302)
    p2303 = load(IN_2303)
    p2304 = load(IN_2304)

    p2300_probe = p2300.get("strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe", {}) or {}
    p2302_probe = p2302.get("strict_task3_provider_lift_policy_lock_candidate_probe", {}) or {}
    p2303_probe = p2303.get("strict_provider_to_margin_bridge_bound_audit_probe", {}) or {}
    p2304_probe = p2304.get("strict_direct_floor_policy_lock_replay_probe", {}) or {}

    basis = p2300_probe.get("provider_basis", []) or []
    coeffs = [float(row.get("canonical_coefficient_numeric", 0.0) or 0.0) for row in basis]
    names = [str(row.get("name", f"c{i}")) for i, row in enumerate(basis)]
    required_lift = float((p2302_probe.get("policy_lock_candidate", {}) or {}).get("provider_lift_per_step", 0.0) or 0.0)
    if required_lift == 0.0:
        required_lift = float((p2303_probe.get("required_lift", {}) or {}).get("value", 0.0) or 0.0)

    l1 = sum(abs(c) for c in coeffs)
    l2 = norm2(coeffs)
    linf = max((abs(c) for c in coeffs), default=0.0)
    positive_sum = sum(c for c in coeffs if c > 0)
    signed_sum = sum(coeffs)

    if l2 > 0.0:
        aligned_response = scale(coeffs, 1.0 / l2)
        anti_aligned_response = scale(coeffs, -1.0 / l2)
    else:
        aligned_response = [0.0 for _ in coeffs]
        anti_aligned_response = [0.0 for _ in coeffs]
    orthogonal_response = unit_orthogonal(coeffs)

    response_witnesses = [
        {
            "witness_id": "UNIT_ALIGNED_POSITIVE_RESPONSE",
            "response_vector": aligned_response,
            "response_norm_l2": norm2(aligned_response),
            "lift_value": dot(aligned_response, coeffs),
            "passes_required_lift": dot(aligned_response, coeffs) >= required_lift,
            "strict_admissible_without_exported_response_functional": False,
            "meaning": "Numerically sufficient only if an aligned provider-to-margin response functional is added.",
        },
        {
            "witness_id": "UNIT_ANTI_ALIGNED_NEGATIVE_RESPONSE",
            "response_vector": anti_aligned_response,
            "response_norm_l2": norm2(anti_aligned_response),
            "lift_value": dot(anti_aligned_response, coeffs),
            "passes_required_lift": dot(anti_aligned_response, coeffs) >= required_lift,
            "strict_admissible_without_exported_response_functional": True,
            "meaning": "Same coefficient vector and norm can produce a negative margin response under a unit functional.",
        },
        {
            "witness_id": "UNIT_ORTHOGONAL_ZERO_RESPONSE",
            "response_vector": orthogonal_response,
            "response_norm_l2": norm2(orthogonal_response),
            "lift_value": dot(orthogonal_response, coeffs),
            "passes_required_lift": dot(orthogonal_response, coeffs) >= required_lift,
            "strict_admissible_without_exported_response_functional": True,
            "meaning": "Same coefficient vector can be invisible to a unit response functional if orientation is not fixed.",
        },
    ]

    aggregation_contracts = [
        {
            "contract_id": "L2_NORM_AS_MARGIN",
            "bound_numeric": l2,
            "required_lift": required_lift,
            "passes_required_lift": l2 >= required_lift,
            "strict_admissible_without_orientation_theorem": False,
            "obstruction": "A norm is unsigned; P2300/P2303/P2304 do not export a theorem mapping unsigned operator energy to signed policy margin.",
        },
        {
            "contract_id": "L1_NORM_AS_MARGIN",
            "bound_numeric": l1,
            "required_lift": required_lift,
            "passes_required_lift": l1 >= required_lift,
            "strict_admissible_without_orientation_theorem": False,
            "obstruction": "Absolute aggregation erases cancellations and requires a sign/monotonicity premise not exported by current strict artifacts.",
        },
        {
            "contract_id": "POSITIVE_CHANNEL_SUM_AS_MARGIN",
            "bound_numeric": positive_sum,
            "required_lift": required_lift,
            "passes_required_lift": positive_sum >= required_lift,
            "strict_admissible_without_orientation_theorem": False,
            "obstruction": "Selecting only positive channels is a selector/orientation premise unless independently derived.",
        },
        {
            "contract_id": "SIGNED_TOTAL_AS_MARGIN",
            "bound_numeric": signed_sum,
            "required_lift": required_lift,
            "passes_required_lift": signed_sum >= required_lift,
            "strict_admissible_without_orientation_theorem": True,
            "obstruction": "The exported signed total is negative and therefore cannot certify a positive policy lift.",
        },
    ]

    theorem_export = {
        "statement_id": "P2305_STRICT_NORM_TO_MARGIN_UNDERDETERMINATION_THEOREM",
        "formal_statement": (
            "For the P2300 canonical coefficient vector c, the numerical L1/L2/positive-channel aggregations exceed the P2302 "
            "candidate lift, but this does not prove a strict provider-to-margin bridge.  In the absence of an exported response "
            "functional or monotone orientation theorem, unit response witnesses can map the same c to +||c||_2, 0, or -||c||_2.  "
            "Therefore norm aggregation remains underdetermined and cannot move strict G1/G3 after the P2304 direct-floor refutation."
        ),
        "proof_bits": {
            "coefficient_names": names,
            "coefficient_l2_norm": l2,
            "coefficient_l1_norm": l1,
            "coefficient_linf_norm": linf,
            "positive_channel_sum": positive_sum,
            "signed_total": signed_sum,
            "required_lift": required_lift,
            "aligned_lift": response_witnesses[0]["lift_value"],
            "orthogonal_lift_abs": abs(response_witnesses[2]["lift_value"]),
            "anti_aligned_lift": response_witnesses[1]["lift_value"],
            "underdetermination_witness_count": len(response_witnesses),
        },
        "not_claimed": [
            "strict norm-to-margin bridge",
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
        "p2300_coefficients_loaded": len(coeffs) == 10 and l2 > 0.0,
        "p2302_required_lift_loaded": required_lift == 0.0068,
        "p2304_direct_floor_refutation_loaded": p2304_probe.get("strict_closure_status") == "HELD_OPEN_DIRECT_FLOOR_ROUTE_REFUTED",
        "norms_numerically_exceed_required_lift": l1 >= required_lift and l2 >= required_lift,
        "aligned_witness_can_pass_only_with_extra_functional": response_witnesses[0]["passes_required_lift"] and not response_witnesses[0]["strict_admissible_without_exported_response_functional"],
        "anti_aligned_witness_refutes_norm_only_implication": response_witnesses[1]["lift_value"] < 0.0,
        "orthogonal_witness_refutes_norm_only_implication": abs(response_witnesses[2]["lift_value"]) < 1e-12,
        "strict_norm_to_margin_bridge_not_proven": True,
        "strict_g1_g3_not_updated": True,
        "no_qw2191_discharge_claimed": True,
        "no_selector_closure_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2305_s1255_v1",
        "packet_id": "P2305",
        "stage_id": "S1255",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_OBSTRUCTION_STRICT_NORM_TO_MARGIN_BRIDGE_UNDERDETERMINED",
        "result_kind": "STRICT_NORM_AGGREGATION_NUMERICALLY_SUFFICIENT_BUT_THEOREM_UNDERDETERMINED",
        "strict_norm_to_margin_bridge_underdetermination_probe": {
            "probe_id": "P2305_S1255_STRICT_NORM_TO_MARGIN_BRIDGE_UNDERDETERMINATION",
            "source_packets": {
                "alpha_geo_strict_derived_v1": "generated/alpha_geo_strict_derived_v1.json",
                "p2300": "generated/p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json",
                "p2302": "generated/p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json",
                "p2303": "generated/p2303_s1253_strict_provider_to_margin_bridge_bound_audit_probe.json",
                "p2304": "generated/p2304_s1254_strict_direct_floor_policy_lock_replay_probe.json",
            },
            "source_hashes": {
                "alpha_geo_strict_derived_v1_sha256": sha256_file(IN_ALPHA),
                "p2300_sha256": sha256_file(IN_2300),
                "p2302_sha256": sha256_file(IN_2302),
                "p2303_sha256": sha256_file(IN_2303),
                "p2304_sha256": sha256_file(IN_2304),
            },
            "coefficient_vector": [
                {"name": name, "coefficient_numeric": coeff} for name, coeff in zip(names, coeffs)
            ],
            "norm_summary": {
                "l1_norm": l1,
                "l2_norm": l2,
                "linf_norm": linf,
                "positive_channel_sum": positive_sum,
                "signed_total": signed_sum,
                "required_lift": required_lift,
            },
            "aggregation_contracts": aggregation_contracts,
            "response_functional_underdetermination_witnesses": response_witnesses,
            "strict_bridge_verdict": {
                "norm_to_margin_bridge_proven": False,
                "norm_to_margin_bridge_refuted_as_norm_only_theorem": True,
                "reason": "Norm magnitude alone does not determine signed policy-margin lift; current strict artifacts lack the required response/orientation theorem.",
            },
            "task3_g1_g3_update": {
                "G1_reduction_certainty": "OPEN_UNCHANGED",
                "G2_nonlinear_trajectory_realism": "CLOSED_FROM_P2301_UNCHANGED",
                "G3_operational_policy_rule": "OPEN_UNCHANGED",
                "reason": "P2305 supplies an underdetermination theorem, not a strict policy-lock bridge.",
            },
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2306_candidate",
            "goal": "Export or refute a strict response-functional/orientation theorem mapping the P2300 operator basis into the P2281 policy margin. Without that theorem, keep P2302/P2303 norm and positive-channel bounds non-strict and do not move G1/G3.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_NORM_TO_MARGIN_UNDERDETERMINATION_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2305/S1255 — strict norm-to-margin bridge underdetermination",
            "",
            f"- Status: `{payload['status']}`",
            f"- Required lift: `{required_lift}`",
            f"- L2 norm: `{l2}`",
            f"- L1 norm: `{l1}`",
            f"- Positive-channel sum: `{positive_sum}`",
            f"- Signed total: `{signed_sum}`",
            "- Norms numerically exceed the required lift, but no strict response-functional/orientation theorem is exported.",
            f"- Strict bridge proven: `{payload['strict_norm_to_margin_bridge_underdetermination_probe']['strict_bridge_verdict']['norm_to_margin_bridge_proven']}`",
            f"- Theorem fingerprint: `{theorem_fingerprint}`",
            "",
            "## Guardrail statement",
            "P2305 does not close G1/G3, does not discharge QW-2191, does not add a selector premise, does not transfer legacy-kernel roles, and does not claim ToE closure.",
            "",
            "## Next honest step",
            payload["recommended_next_honest_step"]["goal"],
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
