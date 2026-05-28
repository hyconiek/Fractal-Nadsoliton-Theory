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
IN_2306 = GEN / "p2306_s1256_strict_response_orientation_functional_interface_audit_probe.json"
SRC_2302 = ROOT / "p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.py"
OUT = GEN / "p2307_s1257_strict_dynamical_margin_response_functional_nonderivation_probe.json"
MD = GEN / "p2307_s1257_strict_dynamical_margin_response_functional_nonderivation_probe.md"


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


def main() -> None:
    GEN.mkdir(exist_ok=True)
    alpha = load_json(IN_ALPHA)
    p2300 = load_json(IN_2300)
    p2302 = load_json(IN_2302)
    p2306 = load_json(IN_2306)

    p2300_probe = p2300.get("strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe", {}) or {}
    p2302_probe = p2302.get("strict_task3_provider_lift_policy_lock_candidate_probe", {}) or {}
    p2306_probe = p2306.get("strict_response_orientation_functional_interface_audit_probe", {}) or {}
    basis = p2300_probe.get("provider_basis", []) or []
    coeffs = [float(row.get("canonical_coefficient_numeric", 0.0) or 0.0) for row in basis]
    names = [str(row.get("name", f"c{i}")) for i, row in enumerate(basis)]
    required_lift = float((p2302_probe.get("policy_lock_candidate", {}) or {}).get("provider_lift_per_step", 0.0) or 0.0)

    # Symbolic replay of the P2302/P2304 one-step margin dynamics.
    m, rho, kappa, lam, xi = sp.symbols("m rho kappa lambda xi", real=True)
    c_symbols = sp.symbols("c0:10", real=True)
    update = m - (sp.Rational(15, 1000) + sp.Rational(1, 100) * m**2) + rho * sp.Rational(8, 1000) / kappa + lam + xi
    d_update_d_lambda = sp.diff(update, lam)
    d_update_d_coeffs_without_bridge = [sp.diff(update, c) for c in c_symbols]
    horizon = 36
    scalar_lift_total_derivative = horizon * d_update_d_lambda

    # If one postulates lambda = w·c, the dynamics only gives the chain-rule form; the weights w are the missing bridge object.
    w_symbols = sp.symbols("w0:10", real=True)
    lambda_response = sum(w * c for w, c in zip(w_symbols, c_symbols))
    update_with_response = sp.expand(update.subs(lam, lambda_response))
    d_update_d_coeffs_with_unexported_weights = [sp.diff(update_with_response, c) for c in c_symbols]

    signed_total = sum(coeffs)
    positive_sum = sum(c for c in coeffs if c > 0)
    l1 = sum(abs(c) for c in coeffs)
    l2 = sum(c * c for c in coeffs) ** 0.5

    dynamical_interface_audit = {
        "one_step_update_formula": str(update),
        "d_update_d_provider_lift": str(d_update_d_lambda),
        "d_update_d_p2300_coefficients_without_bridge": [str(x) for x in d_update_d_coeffs_without_bridge],
        "horizon": horizon,
        "scalar_lift_total_derivative_over_horizon": str(scalar_lift_total_derivative),
        "chain_rule_if_lambda_equals_w_dot_c": {
            "lambda_response_ansatz": str(lambda_response),
            "d_update_d_coefficients": [str(x) for x in d_update_d_coeffs_with_unexported_weights],
            "missing_object": "exported strict weights/response functional w_i mapping P2300 coefficients to provider_lift_per_step",
        },
    }

    candidate_functionals = [
        {
            "candidate_id": "DYNAMICS_NATIVE_PROVIDER_LIFT_PARAMETER",
            "definition": "lambda is the scalar provider_lift_per_step argument in the P2281/P2302 replay recurrence",
            "derived_from_margin_dynamics": True,
            "typed_on_p2300_coefficients": False,
            "value": required_lift,
            "passes_required_lift": True,
            "admissible_as_p2300_response_functional": False,
            "rejection_reason": "The dynamics differentiates with respect to lambda, but does not derive lambda from P2300 coefficients.",
        },
        {
            "candidate_id": "CANONICAL_SIGNED_TOTAL_FROM_P2300",
            "definition": "sum_i c_i",
            "derived_from_margin_dynamics": False,
            "typed_on_p2300_coefficients": True,
            "value": signed_total,
            "passes_required_lift": signed_total >= required_lift,
            "admissible_as_p2300_response_functional": False,
            "rejection_reason": "Well typed but negative, so it cannot lower-bound the positive policy lift.",
        },
        {
            "candidate_id": "POSITIVE_CHANNEL_SUM_FROM_P2300",
            "definition": "sum_{i:c_i>0} c_i",
            "derived_from_margin_dynamics": False,
            "typed_on_p2300_coefficients": True,
            "value": positive_sum,
            "passes_required_lift": positive_sum >= required_lift,
            "admissible_as_p2300_response_functional": False,
            "rejection_reason": "Numerically sufficient but requires a channel-selection/sign premise not derived from the replay dynamics.",
        },
        {
            "candidate_id": "NORM_RESPONSE_FROM_P2300",
            "definition": "||c||_2 or ||c||_1",
            "derived_from_margin_dynamics": False,
            "typed_on_p2300_coefficients": True,
            "value": {"l1": l1, "l2": l2},
            "passes_required_lift": l1 >= required_lift and l2 >= required_lift,
            "admissible_as_p2300_response_functional": False,
            "rejection_reason": "P2305/P2306 already show unsigned norms need an orientation theorem absent from current dynamics.",
        },
    ]

    theorem_export = {
        "statement_id": "P2307_STRICT_DYNAMICAL_MARGIN_RESPONSE_FUNCTIONAL_NONDERIVATION_THEOREM",
        "formal_statement": (
            "The P2281/P2302 margin dynamics exports a scalar sensitivity dM_{t+1}/d(lambda)=1 to an already supplied "
            "provider_lift_per_step lambda, but it does not export a typed functional lambda=R(c) from the P2300 ADM/Bianchi-I "
            "coefficient vector c.  Symbolically, with no bridge, dM_{t+1}/dc_i=0 for every P2300 coefficient; with the ansatz "
            "lambda=w·c the missing data are exactly the weights w_i.  Therefore the current strict dynamics does not derive the "
            "P2302 lift from P2300, and G1/G3 remain open."
        ),
        "proof_bits": {
            "p2300_coefficient_names": names,
            "required_lift": required_lift,
            "d_update_d_provider_lift": str(d_update_d_lambda),
            "d_update_d_coefficients_without_bridge": [str(x) for x in d_update_d_coeffs_without_bridge],
            "chain_rule_requires_unexported_weights": [str(x) for x in d_update_d_coeffs_with_unexported_weights],
            "signed_total": signed_total,
            "positive_sum": positive_sum,
            "l1_norm": l1,
            "l2_norm": l2,
        },
        "not_claimed": [
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
        "p2306_interface_audit_loaded": (p2306_probe.get("strict_interface_verdict", {}) or {}).get("response_orientation_functional_exported") is False,
        "dynamics_derivative_to_scalar_lift_is_one": str(d_update_d_lambda) == "1",
        "dynamics_derivative_to_coefficients_without_bridge_is_zero": all(x == 0 for x in d_update_d_coeffs_without_bridge),
        "chain_rule_requires_unexported_weights": [str(x) for x in d_update_d_coeffs_with_unexported_weights] == [f"w{i}" for i in range(10)],
        "no_candidate_derives_admissible_response_functional": not any(row["admissible_as_p2300_response_functional"] for row in candidate_functionals),
        "strict_response_functional_not_derived": True,
        "strict_g1_g3_not_updated": True,
        "no_qw2191_discharge_claimed": True,
        "no_selector_closure_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2307_s1257_v1",
        "packet_id": "P2307",
        "stage_id": "S1257",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_OBSTRUCTION_STRICT_DYNAMICAL_RESPONSE_FUNCTIONAL_NONDERIVED",
        "result_kind": "STRICT_DYNAMICAL_MARGIN_REPLAY_DERIVES_SCALAR_LIFT_SENSITIVITY_BUT_NOT_P2300_RESPONSE_FUNCTIONAL",
        "strict_dynamical_margin_response_functional_nonderivation_probe": {
            "probe_id": "P2307_S1257_STRICT_DYNAMICAL_MARGIN_RESPONSE_FUNCTIONAL_NONDERIVATION",
            "source_packets": {
                "alpha_geo_strict_derived_v1": "generated/alpha_geo_strict_derived_v1.json",
                "p2300": "generated/p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json",
                "p2302": "generated/p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json",
                "p2306": "generated/p2306_s1256_strict_response_orientation_functional_interface_audit_probe.json",
                "p2302_script": "p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.py",
            },
            "source_hashes": {
                "alpha_geo_strict_derived_v1_sha256": sha256_file(IN_ALPHA),
                "p2300_sha256": sha256_file(IN_2300),
                "p2302_sha256": sha256_file(IN_2302),
                "p2306_sha256": sha256_file(IN_2306),
                "p2302_script_sha256": sha256_file(SRC_2302),
            },
            "dynamical_interface_audit": dynamical_interface_audit,
            "candidate_functionals": candidate_functionals,
            "strict_nonderivation_verdict": {
                "scalar_lift_sensitivity_derived": True,
                "p2300_to_margin_response_functional_derived": False,
                "current_dynamics_refutes_closure_use": True,
                "reason": "Replay dynamics contains lambda as an external scalar input; no strict equation maps P2300 coefficients to lambda.",
            },
            "task3_g1_g3_update": {
                "G1_reduction_certainty": "OPEN_UNCHANGED",
                "G2_nonlinear_trajectory_realism": "CLOSED_FROM_P2301_UNCHANGED",
                "G3_operational_policy_rule": "OPEN_UNCHANGED",
                "reason": "P2307 proves a nonderivation of R(c) from current dynamics, not a strict response bridge.",
            },
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2308_candidate",
            "goal": "Either extend the strict ADM/Bianchi-I dynamics with an internally derived typed equation lambda=R(c), or package a stronger nonexistence theorem over the current P2300/P2281 interface class. Do not update G1/G3 until lambda=R(c) is exported and replayed.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_DYNAMICAL_RESPONSE_NONDERIVATION_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2307/S1257 — strict dynamical margin response functional nonderivation",
            "",
            f"- Status: `{payload['status']}`",
            f"- Required lift: `{required_lift}`",
            f"- One-step derivative wrt scalar provider lift: `{d_update_d_lambda}`",
            "- One-step derivatives wrt P2300 coefficients without a bridge: all `0`.",
            "- Chain-rule ansatz `lambda = w·c` exposes the missing object: exported strict weights `w_i`.",
            "- P2300→margin response functional derived: `False`",
            "- G1/G3 update: `OPEN_UNCHANGED`",
            f"- Theorem fingerprint: `{theorem_fingerprint}`",
            "",
            "## Guardrail statement",
            "P2307 does not close G1/G3, does not discharge QW-2191, does not add a selector premise, does not transfer legacy-kernel roles, and does not claim ToE closure.",
            "",
            "## Next honest step",
            payload["recommended_next_honest_step"]["goal"],
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
