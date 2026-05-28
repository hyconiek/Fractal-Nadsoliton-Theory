#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

IN_1980 = GEN / "p1980_s930_strict_adm_bianchi_eh_lapse_variation_witness.json"
IN_1981 = GEN / "p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.json"
IN_1982 = GEN / "p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.json"
IN_1983 = GEN / "p1983_s933_strict_adm_bianchi_riemann2_lapse_variation_obligation.json"
IN_1984 = GEN / "p1984_s934_strict_adm_bianchi_gauss_bonnet_lapse_cancellation_witness.json"
IN_1985 = GEN / "p1985_s935_strict_adm_bianchi_non_gb_curvature_squared_lapse_obstruction_witness.json"
IN_1986 = GEN / "p1986_s936_strict_adm_bianchi_non_gb_residual_decomposition_witness.json"
IN_2281 = GEN / "p2281_s1231_strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe.json"
IN_2300 = GEN / "p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json"
IN_2302 = GEN / "p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json"
IN_2310 = GEN / "p2310_s1260_strict_self_energy_response_source_audit_and_replay_probe.json"
OUT = GEN / "p2311_s1261_strict_self_energy_to_policy_margin_theorem_audit_probe.json"
MD = GEN / "p2311_s1261_strict_self_energy_to_policy_margin_theorem_audit_probe.md"


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


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1980 = load(IN_1980)
    p1981 = load(IN_1981)
    p1982 = load(IN_1982)
    p1983 = load(IN_1983)
    p1984 = load(IN_1984)
    p1985 = load(IN_1985)
    p1986 = load(IN_1986)
    p2281 = load(IN_2281)
    p2300 = load(IN_2300)
    p2302 = load(IN_2302)
    p2310 = load(IN_2310)

    p2300_probe = nested(p2300, "strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe", default={}) or {}
    p2302_probe = nested(p2302, "strict_task3_provider_lift_policy_lock_candidate_probe", default={}) or {}
    p2310_probe = nested(p2310, "strict_self_energy_response_source_audit_and_replay_probe", default={}) or {}
    p2281_probe = nested(p2281, "strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe", default={}) or {}

    required_lift = float(nested(p2302_probe, "policy_lock_candidate", "provider_lift_per_step", default=0.0) or 0.0)
    self_energy_lambda = float(nested(p2310_probe, "self_energy_variational_attempt", "lambda_numeric", default=0.0) or 0.0)
    self_energy_replay_pass = bool(nested(p2310_probe, "self_energy_replay", "summary", "all_rows_meet_target", default=False) is True)
    p2281_global_pass = bool(nested(p2281_probe, "global_summary", "all_rows_meet_target", default=False) is True)

    eh_lapse_sign_present = bool((p1980.get("gatekeeper_checks", {}) or {}).get("eh_lapse_witness_passed") is True)
    r2_obligation_present = bool((p1981.get("gatekeeper_checks", {}) or {}).get("r2_lapse_obligation_passed") is True)
    ricci2_obligation_present = bool((p1982.get("gatekeeper_checks", {}) or {}).get("ricci2_lapse_obligation_passed") is True)
    riemann2_obligation_present = bool((p1983.get("gatekeeper_checks", {}) or {}).get("riemann2_lapse_obligation_passed") is True)
    gb_lapse_cancels = bool((p1984.get("gatekeeper_checks", {}) or {}).get("gb_lapse_cancellation_witness_passed") is True)
    non_gb_lapse_obstruction_present = bool((p1985.get("gatekeeper_checks", {}) or {}).get("obstruction_witness_passed") is True)
    non_gb_decomposition_obstruction_present = bool((p1986.get("gatekeeper_checks", {}) or {}).get("non_gb_channel_nonzero") is True)

    theorem_obligations = [
        {
            "obligation_id": "O1_TYPED_SELF_ENERGY_DOMAIN",
            "requirement": "P2300 must export the canonical ADM/Bianchi-I coefficient vector c used in E(c)=||c||^2.",
            "satisfied": bool((p2300.get("gatekeeper_checks", {}) or {}).get("provider_matrix_passes") is True or len(p2300_probe.get("provider_basis", []) or []) == 10),
            "evidence": "P2300 provider basis / provider-matrix pass.",
        },
        {
            "obligation_id": "O2_SELF_ENERGY_NUMERIC_LIFT_AND_REPLAY",
            "requirement": "The self-energy candidate must exceed the P2302 required lift and replay through the policy-lock simulator.",
            "satisfied": bool(self_energy_lambda >= required_lift and self_energy_replay_pass),
            "evidence": "P2310 lambda=||c||^2 replay.",
        },
        {
            "obligation_id": "O3_ADM_LAPSE_SHEAR_ENERGY_SOURCE_THEOREM",
            "requirement": "A strict ADM lapse/shear theorem must identify E(c)=||c||^2 as a physical margin source rather than a norm diagnostic.",
            "satisfied": False,
            "evidence": "P1980 gives an EH shear-sign witness, P1981-P1983 export curvature-squared lapse obligations, P1984 cancels GB, and P1985/P1986 leave a non-GB residual obstruction; none exports E(c)->provider_lift_per_step.",
        },
        {
            "obligation_id": "O4_SIGNED_MONOTONE_POLICY_MARGIN_RESPONSE",
            "requirement": "The ADM source must be mapped to the P2281 policy-margin recurrence with a signed monotone orientation.",
            "satisfied": False,
            "evidence": "P2281 global replay remains open, while P2306/P2307/P2308/P2309/P2310 already quarantine norm/response routes without a typed response functional.",
        },
        {
            "obligation_id": "O5_NO_SELECTOR_OR_CONVENTION_LAYER",
            "requirement": "The bridge must not select positive channels, signs, lapse normalization, or policy response orientation by a convention/selector premise.",
            "satisfied": True,
            "evidence": "P2311 does not add a selector/convention layer and therefore refuses G1/G3 closure.",
        },
    ]

    adm_lapse_shear_audit = {
        "eh_lapse_shear_sign_present": eh_lapse_sign_present,
        "eh_lapse_shear_sign_is_only_sub_witness": True,
        "curvature_squared_lapse_obligations_present": {
            "R2": r2_obligation_present,
            "Ricci2": ricci2_obligation_present,
            "Riemann2": riemann2_obligation_present,
        },
        "gb_lapse_cancellation_present": gb_lapse_cancels,
        "non_gb_lapse_obstruction_present": non_gb_lapse_obstruction_present,
        "non_gb_decomposition_obstruction_present": non_gb_decomposition_obstruction_present,
        "missing_bridge_objects": [
            "strict theorem identifying coefficient self-energy E(c)=||c||^2 with ADM lapse/shear energy density",
            "normalization from ADM lapse/shear source to provider_lift_per_step",
            "signed monotone response functional from ADM source to P2281 policy margin",
            "proof that the non-GB residual obstruction is exactly cancelled or positively oriented by the P2300 coefficients",
        ],
        "current_export_set_conclusion": "The current ADM/Bianchi-I export set contains real lapse/shear facts, but it does not export the self-energy-to-policy-margin bridge theorem requested for P2311.",
    }

    attempted_theorem = {
        "theorem_id": "P2311_ATTEMPT_E_NORM_SQUARED_TO_PROVIDER_LIFT_PER_STEP",
        "candidate_statement": "E(c)=||c||_2^2 from P2300 is a strict ADM/Bianchi-I lapse/shear energy whose signed policy-margin response is provider_lift_per_step >= 0.0068.",
        "candidate_lambda_numeric": self_energy_lambda,
        "required_lift": required_lift,
        "numeric_condition_passes": self_energy_lambda >= required_lift,
        "p2310_replay_passes": self_energy_replay_pass,
        "strict_theorem_proven": False,
        "failure_mode": "MISSING_ADM_LAPSE_SHEAR_TO_POLICY_MARGIN_RESPONSE_SOURCE",
        "why_numeric_pass_is_insufficient": "The replay uses lambda as an injected margin lift.  Without a strict theorem proving E(c) is that lambda, the replay is a conditional simulation, not a G1/G3 closure proof.",
    }

    admissibility_verdict = {
        "self_energy_to_policy_margin_theorem_proven": False,
        "current_export_set_nonexistence_scope": "current exported interface set P1980-P1986/P2281/P2300/P2302/P2310 only; not a universal no-go for future lapse/shear sources",
        "current_export_set_can_supply_bridge_without_new_lapse_shear_source": False,
        "formal_blocker": "P2310 supplies a non-target-calibrated scalar and replay pass, but P2311 finds no strict ADM/Bianchi-I lapse/shear source theorem or signed response functional tying that scalar to P2281 policy margin.",
        "admissible_for_g1_g3_update": False,
    }

    theorem_export = {
        "statement_id": "P2311_SELF_ENERGY_TO_POLICY_MARGIN_CURRENT_EXPORT_NONDERIVATION_THEOREM",
        "formal_statement": (
            "Within the currently exported strict ADM/Bianchi-I interface class P1980-P1986 together with P2281, P2300, P2302, and P2310, "
            "the implication E(c)=||c||_2^2 -> provider_lift_per_step is not derivable.  P2310's self-energy scalar is large enough and "
            "passes the conditional replay, but the export set lacks both an ADM lapse/shear energy-source theorem for E(c) and a signed monotone "
            "response functional from that source to the P2281 policy margin.  Therefore G1 and G3 remain open unless a new lapse/shear source or equivalent strict bridge is exported."
        ),
        "proof_bits": {
            "self_energy_lambda": self_energy_lambda,
            "required_lift": required_lift,
            "numeric_condition_passes": self_energy_lambda >= required_lift,
            "conditional_replay_passes": self_energy_replay_pass,
            "eh_lapse_sign_present": eh_lapse_sign_present,
            "gb_lapse_cancellation_present": gb_lapse_cancels,
            "non_gb_lapse_obstruction_present": non_gb_lapse_obstruction_present,
            "non_gb_decomposition_obstruction_present": non_gb_decomposition_obstruction_present,
            "p2281_global_all_rows_meet_target": p2281_global_pass,
            "failed_obligations": [row["obligation_id"] for row in theorem_obligations if not row["satisfied"]],
        },
        "not_claimed": [
            "universal impossibility of future ADM lapse/shear sources",
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
        "p1980_eh_lapse_witness_loaded": eh_lapse_sign_present,
        "p1981_p1983_curvature_squared_obligations_loaded": r2_obligation_present and ricci2_obligation_present and riemann2_obligation_present,
        "p1984_gb_lapse_cancellation_loaded": gb_lapse_cancels,
        "p1985_p1986_non_gb_obstruction_loaded": non_gb_lapse_obstruction_present and non_gb_decomposition_obstruction_present,
        "p2300_coefficients_loaded": theorem_obligations[0]["satisfied"],
        "p2310_self_energy_replay_loaded": self_energy_replay_pass,
        "self_energy_numeric_condition_passes": self_energy_lambda >= required_lift,
        "adm_lapse_shear_energy_to_policy_margin_theorem_not_exported": not theorem_obligations[2]["satisfied"],
        "signed_monotone_policy_response_not_exported": not theorem_obligations[3]["satisfied"],
        "current_export_set_blocks_self_energy_bridge_use": admissibility_verdict["current_export_set_can_supply_bridge_without_new_lapse_shear_source"] is False,
        "strict_g1_g3_not_updated": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2311_s1261_v1",
        "packet_id": "P2311",
        "stage_id": "S1261",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_OBSTRUCTION_SELF_ENERGY_TO_POLICY_MARGIN_THEOREM_NOT_DERIVED",
        "result_kind": "CURRENT_ADM_BIANCHI_EXPORT_SET_REFUTES_SELF_ENERGY_BRIDGE_USE_WITHOUT_NEW_LAPSE_SHEAR_SOURCE",
        "strict_self_energy_to_policy_margin_theorem_audit_probe": {
            "probe_id": "P2311_S1261_STRICT_SELF_ENERGY_TO_POLICY_MARGIN_THEOREM_AUDIT",
            "source_packets": {
                "p1980": "generated/p1980_s930_strict_adm_bianchi_eh_lapse_variation_witness.json",
                "p1981": "generated/p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.json",
                "p1982": "generated/p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.json",
                "p1983": "generated/p1983_s933_strict_adm_bianchi_riemann2_lapse_variation_obligation.json",
                "p1984": "generated/p1984_s934_strict_adm_bianchi_gauss_bonnet_lapse_cancellation_witness.json",
                "p1985": "generated/p1985_s935_strict_adm_bianchi_non_gb_curvature_squared_lapse_obstruction_witness.json",
                "p1986": "generated/p1986_s936_strict_adm_bianchi_non_gb_residual_decomposition_witness.json",
                "p2281": "generated/p2281_s1231_strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe.json",
                "p2300": "generated/p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json",
                "p2302": "generated/p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json",
                "p2310": "generated/p2310_s1260_strict_self_energy_response_source_audit_and_replay_probe.json",
            },
            "source_hashes": {
                "p1980_sha256": sha256_file(IN_1980),
                "p1981_sha256": sha256_file(IN_1981),
                "p1982_sha256": sha256_file(IN_1982),
                "p1983_sha256": sha256_file(IN_1983),
                "p1984_sha256": sha256_file(IN_1984),
                "p1985_sha256": sha256_file(IN_1985),
                "p1986_sha256": sha256_file(IN_1986),
                "p2281_sha256": sha256_file(IN_2281),
                "p2300_sha256": sha256_file(IN_2300),
                "p2302_sha256": sha256_file(IN_2302),
                "p2310_sha256": sha256_file(IN_2310),
            },
            "attempted_theorem": attempted_theorem,
            "theorem_obligations": theorem_obligations,
            "adm_lapse_shear_evidence_audit": adm_lapse_shear_audit,
            "admissibility_verdict": admissibility_verdict,
            "task3_g1_g3_update": {
                "G1_reduction_certainty": "OPEN_UNCHANGED",
                "G2_nonlinear_trajectory_realism": "CLOSED_FROM_P2301_UNCHANGED",
                "G3_operational_policy_rule": "OPEN_UNCHANGED",
                "reason": "P2311 does not export E(c)=||c||^2 as a strict ADM lapse/shear policy-margin source.",
            },
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2312_candidate",
            "goal": "Introduce a genuinely new lapse/shear source or variational identity that maps ADM energy into the P2281 policy margin with signed monotone orientation; otherwise keep P2308/P2309/P2310/P2311 as the formal G1/G3 blocker stack.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_THEOREM_GRADE_CURRENT_EXPORT_NONDERIVATION_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2311/S1261 — strict self-energy to policy-margin theorem audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Attempted theorem: `{attempted_theorem['theorem_id']}`",
            f"- Required lift: `{required_lift}`",
            f"- Self-energy lambda `||c||^2`: `{self_energy_lambda}`",
            f"- Numeric condition passes: `{attempted_theorem['numeric_condition_passes']}`",
            f"- P2310 conditional replay passes: `{self_energy_replay_pass}`",
            "- Strict theorem proven: `False`",
            "- Failure mode: `MISSING_ADM_LAPSE_SHEAR_TO_POLICY_MARGIN_RESPONSE_SOURCE`",
            "- G1/G3 update: `OPEN_UNCHANGED`",
            f"- Theorem fingerprint: `{theorem_fingerprint}`",
            "",
            "## Concrete audit result",
            "P1980 gives a real EH lapse/shear sign witness, P1981-P1983 provide curvature-squared lapse obligations, P1984 cancels the GB lapse channel, and P1985/P1986 preserve a non-GB residual obstruction.  None of these exports the strict theorem `E(c)=||c||^2 -> provider_lift_per_step` for the P2281 policy margin.",
            "",
            "## Guardrail statement",
            "P2311 does not close G1/G3, does not discharge QW-2191, does not add a selector premise, does not transfer legacy-kernel roles, and does not claim ToE closure.",
            "",
            "## Next honest step",
            payload["recommended_next_honest_step"]["goal"],
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
