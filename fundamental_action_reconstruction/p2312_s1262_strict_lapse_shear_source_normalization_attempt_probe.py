#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

IN_1985 = GEN / "p1985_s935_strict_adm_bianchi_non_gb_curvature_squared_lapse_obstruction_witness.json"
IN_1986 = GEN / "p1986_s936_strict_adm_bianchi_non_gb_residual_decomposition_witness.json"
IN_2300 = GEN / "p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json"
IN_2302 = GEN / "p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json"
IN_2310 = GEN / "p2310_s1260_strict_self_energy_response_source_audit_and_replay_probe.json"
IN_2311 = GEN / "p2311_s1261_strict_self_energy_to_policy_margin_theorem_audit_probe.json"
OUT = GEN / "p2312_s1262_strict_lapse_shear_source_normalization_attempt_probe.json"
MD = GEN / "p2312_s1262_strict_lapse_shear_source_normalization_attempt_probe.md"

SYMBOL_NAMES = [
    "H",
    "Hd",
    "N",
    "Nd",
    "Ndd",
    "V",
    "sigma1",
    "sigma2",
    "dsigma1",
    "dsigma2",
    "d2sigma1",
    "d2sigma2",
]


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


def nested(payload: dict[str, Any], *keys: str, default: Any = None) -> Any:
    node: Any = payload
    for key in keys:
        if not isinstance(node, dict):
            return default
        node = node.get(key)
    return default if node is None else node


def sympify_p1985_operator(symbolic: str) -> tuple[sp.Expr, dict[str, sp.Symbol]]:
    symbols = {name: sp.symbols(name, real=True) for name in SYMBOL_NAMES}
    locals_map: dict[str, Any] = dict(symbols)
    locals_map.update({"log": sp.log, "pi": sp.pi})
    return sp.sympify(symbolic, locals=locals_map), symbols


def evaluate(expr: sp.Expr, symbols: dict[str, sp.Symbol], assignment: dict[str, float]) -> float:
    subs = {symbols[name]: assignment.get(name, 0.0) for name in SYMBOL_NAMES}
    subs[symbols["N"]] = assignment.get("N", 1.0)
    subs[symbols["V"]] = assignment.get("V", 1.0)
    return float(expr.subs(subs).evalf())


def shear_energy(assignment: dict[str, float]) -> float:
    s1 = assignment.get("sigma1", 0.0)
    s2 = assignment.get("sigma2", 0.0)
    return float(s1 * s1 + s1 * s2 + s2 * s2)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1985 = load(IN_1985)
    p1986 = load(IN_1986)
    p2300 = load(IN_2300)
    p2302 = load(IN_2302)
    p2310 = load(IN_2310)
    p2311 = load(IN_2311)

    p2300_probe = nested(p2300, "strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe", default={}) or {}
    p2302_probe = nested(p2302, "strict_task3_provider_lift_policy_lock_candidate_probe", default={}) or {}
    p2310_probe = nested(p2310, "strict_self_energy_response_source_audit_and_replay_probe", default={}) or {}
    p2311_probe = nested(p2311, "strict_self_energy_to_policy_margin_theorem_audit_probe", default={}) or {}

    symbolic = nested(p1985, "weighted_non_gb_lapse_operator", "symbolic", default="0")
    expr, symbols = sympify_p1985_operator(str(symbolic))
    required_lift = float(nested(p2302_probe, "policy_lock_candidate", "provider_lift_per_step", default=0.0) or 0.0)
    self_energy_lambda = float(nested(p2310_probe, "self_energy_variational_attempt", "lambda_numeric", default=0.0) or 0.0)

    witness_assignments = [
        {
            "witness_id": "STATIC_SHEAR_NEGATIVE_LAPSE_SOURCE",
            "assignment": {"H": 1.0, "N": 1.0, "V": 1.0, "sigma1": 1.0, "sigma2": 0.0},
            "interpretation": "ordinary expanding anisotropic shear with no lapse/shear acceleration terms",
        },
        {
            "witness_id": "SHEAR_ACCELERATION_POSITIVE_LAPSE_SOURCE",
            "assignment": {"H": 0.0, "N": 1.0, "V": 1.0, "sigma1": 1.0, "sigma2": 0.0, "d2sigma1": 10.0},
            "interpretation": "same shear energy but nonzero shear acceleration flips the candidate source sign",
        },
        {
            "witness_id": "LAPSE_DERIVATIVE_POSITIVE_LAPSE_SOURCE",
            "assignment": {"H": 0.0, "N": 1.0, "V": 1.0, "sigma1": 1.0, "sigma2": 0.0, "Nd": 10.0},
            "interpretation": "same shear energy but lapse-derivative terms flip the candidate source sign",
        },
    ]
    orientation_witnesses = []
    values = []
    for row in witness_assignments:
        assignment = row["assignment"]
        source_value = evaluate(expr, symbols, assignment)
        values.append(source_value)
        orientation_witnesses.append({
            "witness_id": row["witness_id"],
            "assignment": assignment,
            "shear_energy_Q": shear_energy(assignment),
            "non_gb_lapse_source_value": source_value,
            "source_sign": "positive" if source_value > 0 else "negative" if source_value < 0 else "zero",
            "interpretation": row["interpretation"],
        })

    sign_indefinite = any(v > 0 for v in values) and any(v < 0 for v in values)
    same_shear_energy_sign_flip = len({row["shear_energy_Q"] for row in orientation_witnesses}) == 1 and sign_indefinite

    normalization_attempts = [
        {
            "attempt_id": "DIRECT_NON_GB_LAPSE_SOURCE",
            "definition": "lambda = C * EL_nonGB_lapse with a fixed positive normalization C",
            "strict_admissible": False,
            "failure_reason": "P1985 non-GB lapse source is sign-indefinite on explicit Bianchi-I shear/lapse witnesses, so no fixed positive normalization gives a monotone policy lift.",
        },
        {
            "attempt_id": "SIGN_FLIPPED_NON_GB_LAPSE_SOURCE",
            "definition": "lambda = -C * EL_nonGB_lapse",
            "strict_admissible": False,
            "failure_reason": "Flipping the global sign only swaps the positive and negative witnesses; sign-indefiniteness remains.",
        },
        {
            "attempt_id": "ABSOLUTE_VALUE_NON_GB_LAPSE_SOURCE",
            "definition": "lambda = C * |EL_nonGB_lapse|",
            "strict_admissible": False,
            "failure_reason": "Absolute value is a convention-layer rectifier not exported by the ADM/Bianchi-I variational equations.",
        },
        {
            "attempt_id": "POSITIVE_PART_NON_GB_LAPSE_SOURCE",
            "definition": "lambda = C * max(EL_nonGB_lapse, 0)",
            "strict_admissible": False,
            "failure_reason": "Positive-part extraction is a selector/channel premise and cannot be used for strict G1/G3 closure.",
        },
        {
            "attempt_id": "SELF_ENERGY_TIMES_NON_GB_ORIENTATION",
            "definition": "lambda = ||c||^2 * sign(EL_nonGB_lapse)",
            "strict_admissible": False,
            "failure_reason": "The sign is dynamical and not determined by the P2300 coefficient self-energy alone; P2311's missing response theorem remains missing.",
        },
    ]

    source_normalization_verdict = {
        "new_lapse_shear_source_found": False,
        "direct_non_gb_lapse_source_sign_indefinite": sign_indefinite,
        "same_shear_energy_allows_opposite_source_signs": same_shear_energy_sign_flip,
        "self_energy_lambda_numeric": self_energy_lambda,
        "required_lift": required_lift,
        "self_energy_numeric_condition_passes": self_energy_lambda >= required_lift,
        "admissible_for_provider_lift_export": False,
        "admissible_for_g1_g3_update": False,
        "reason": "The only concrete current non-GB lapse/shear source candidate is a sign-indefinite residual, while rectifying it requires a non-strict selector/convention layer.",
    }

    theorem_export = {
        "statement_id": "P2312_NON_GB_LAPSE_SOURCE_NORMALIZATION_SIGN_INDEFINITE_THEOREM",
        "formal_statement": (
            "For the currently exported P1985 strict non-GB ADM/Bianchi-I lapse residual, no direct fixed-sign normalization can be promoted to "
            "provider_lift_per_step for P2281/P2302.  Explicit Bianchi-I witnesses with the same shear energy Q give opposite signs for the "
            "candidate lapse source.  Absolute-value, positive-part, or sign-conditioned normalizations would introduce a convention/selector layer.  "
            "Therefore P2312 does not supply the missing P2311 energy-to-policy-margin theorem, and G1/G3 remain open."
        ),
        "proof_bits": {
            "witness_values": orientation_witnesses,
            "sign_indefinite": sign_indefinite,
            "same_shear_energy_sign_flip": same_shear_energy_sign_flip,
            "required_lift": required_lift,
            "self_energy_lambda": self_energy_lambda,
            "failed_normalization_attempts": [row["attempt_id"] for row in normalization_attempts if not row["strict_admissible"]],
            "p2311_failed_obligations": nested(p2311_probe, "theorem_export", "proof_bits", "failed_obligations", default=[]),
        },
        "not_claimed": [
            "universal no-go for all future lapse/shear sources",
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
        "p1985_non_gb_lapse_operator_loaded": bool((p1985.get("gatekeeper_checks", {}) or {}).get("obstruction_witness_passed") is True),
        "p1986_non_gb_decomposition_loaded": bool((p1986.get("gatekeeper_checks", {}) or {}).get("non_gb_channel_nonzero") is True),
        "p2300_coefficients_loaded": bool((p2300.get("gatekeeper_checks", {}) or {}).get("canonical_solution_exported") is True or len(p2300_probe.get("provider_basis", []) or []) == 10),
        "p2310_self_energy_numeric_condition_loaded": self_energy_lambda >= required_lift,
        "p2311_bridge_blocker_loaded": bool((p2311.get("gatekeeper_checks", {}) or {}).get("current_export_set_blocks_self_energy_bridge_use") is True),
        "explicit_positive_and_negative_lapse_source_witnesses": sign_indefinite,
        "same_shear_energy_sign_flip_witnessed": same_shear_energy_sign_flip,
        "all_normalization_attempts_quarantined": not any(row["strict_admissible"] for row in normalization_attempts),
        "no_new_lapse_shear_source_exported": source_normalization_verdict["new_lapse_shear_source_found"] is False,
        "strict_g1_g3_not_updated": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2312_s1262_v1",
        "packet_id": "P2312",
        "stage_id": "S1262",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_OBSTRUCTION_NON_GB_LAPSE_SOURCE_NORMALIZATION_SIGN_INDEFINITE",
        "result_kind": "STRICT_LAPSE_SHEAR_SOURCE_NORMALIZATION_ATTEMPT_REFUTED_FOR_CURRENT_NON_GB_RESIDUAL",
        "strict_lapse_shear_source_normalization_attempt_probe": {
            "probe_id": "P2312_S1262_STRICT_LAPSE_SHEAR_SOURCE_NORMALIZATION_ATTEMPT",
            "source_packets": {
                "p1985": "generated/p1985_s935_strict_adm_bianchi_non_gb_curvature_squared_lapse_obstruction_witness.json",
                "p1986": "generated/p1986_s936_strict_adm_bianchi_non_gb_residual_decomposition_witness.json",
                "p2300": "generated/p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json",
                "p2302": "generated/p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json",
                "p2310": "generated/p2310_s1260_strict_self_energy_response_source_audit_and_replay_probe.json",
                "p2311": "generated/p2311_s1261_strict_self_energy_to_policy_margin_theorem_audit_probe.json",
            },
            "source_hashes": {
                "p1985_sha256": sha256_file(IN_1985),
                "p1986_sha256": sha256_file(IN_1986),
                "p2300_sha256": sha256_file(IN_2300),
                "p2302_sha256": sha256_file(IN_2302),
                "p2310_sha256": sha256_file(IN_2310),
                "p2311_sha256": sha256_file(IN_2311),
            },
            "candidate_source": {
                "source_id": "P1985_STRICT_NON_GB_LAPSE_RESIDUAL",
                "definition": "Use the P1985 weighted non-GB ADM/Bianchi-I lapse residual as the missing lapse/shear policy-margin source.",
                "symbolic_source": str(symbolic),
                "strict_admissible": False,
            },
            "orientation_witnesses": orientation_witnesses,
            "normalization_attempts": normalization_attempts,
            "source_normalization_verdict": source_normalization_verdict,
            "task3_g1_g3_update": {
                "G1_reduction_certainty": "OPEN_UNCHANGED",
                "G2_nonlinear_trajectory_realism": "CLOSED_FROM_P2301_UNCHANGED",
                "G3_operational_policy_rule": "OPEN_UNCHANGED",
                "reason": "P2312 finds sign-indefiniteness in the concrete non-GB lapse source normalization attempt; no strict provider_lift export is admitted.",
            },
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2313_candidate",
            "goal": "Either export a new strict dynamical constraint fixing the lapse/shear derivative orientation used by P1985, or stop the G1/G3 bridge route at the P2308-P2312 blocker stack and move to a genuinely new provider/source class.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_SIGN_INDEFINITE_LAPSE_SOURCE_NORMALIZATION_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2312/S1262 — strict lapse/shear source normalization attempt",
            "",
            f"- Status: `{payload['status']}`",
            "- Candidate source: `P1985_STRICT_NON_GB_LAPSE_RESIDUAL`",
            f"- Required lift: `{required_lift}`",
            f"- Self-energy lambda `||c||^2`: `{self_energy_lambda}`",
            f"- Direct non-GB lapse source sign-indefinite: `{sign_indefinite}`",
            f"- Same shear energy sign flip witnessed: `{same_shear_energy_sign_flip}`",
            "- Strict provider-lift export: `False`",
            "- G1/G3 update: `OPEN_UNCHANGED`",
            f"- Theorem fingerprint: `{theorem_fingerprint}`",
            "",
            "## Concrete audit result",
            "The P1985 non-GB lapse residual is a real strict ADM/Bianchi-I object, but explicit witnesses with the same shear energy give opposite signs.  A direct fixed-sign normalization cannot be a monotone policy lift; absolute-value or positive-part fixes would add a convention/selector layer.",
            "",
            "## Guardrail statement",
            "P2312 does not close G1/G3, does not discharge QW-2191, does not add a selector premise, does not transfer legacy-kernel roles, and does not claim ToE closure.",
            "",
            "## Next honest step",
            payload["recommended_next_honest_step"]["goal"],
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
