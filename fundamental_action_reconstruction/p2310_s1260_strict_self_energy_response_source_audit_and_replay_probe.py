#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import random
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_ALPHA = GEN / "alpha_geo_strict_derived_v1.json"
IN_1979 = GEN / "p1979_s929_strict_lapse_shear_provider_export_audit.json"
IN_1988 = GEN / "p1988_s938_strict_non_gb_to_spatial_eom_family_projection_witness.json"
IN_2274 = GEN / "p2274_s1224_strict_nu_branch_group_policy_robustness_region_certificate_probe.json"
IN_2300 = GEN / "p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json"
IN_2302 = GEN / "p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json"
IN_2309 = GEN / "p2309_s1259_strict_min_norm_response_weights_replay_quarantine_probe.json"
OUT = GEN / "p2310_s1260_strict_self_energy_response_source_audit_and_replay_probe.json"
MD = GEN / "p2310_s1260_strict_self_energy_response_source_audit_and_replay_probe.md"


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


def simulate_pass_rate(rho: float, kappa_scale: float, seed: int, *, trials: int, horizon: int, noise_amp: float, provider_lift: float) -> float:
    rng = random.Random(seed)
    passed = 0
    for _ in range(trials):
        margin = 0.20
        for _ in range(horizon):
            noise = (rng.random() - 0.5) * noise_amp
            drift = 0.015 + 0.01 * (margin ** 2)
            control_relief = rho * 0.008 / max(kappa_scale, 1e-12)
            margin = margin - drift + control_relief + provider_lift + noise
            if margin < -1e-12:
                break
        if margin >= -1e-12:
            passed += 1
    return passed / trials


def replay_rows_for_lift(cert_rows: list[dict[str, Any]], provider_lift: float, *, delta: float, trial_multiplier: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    seeds = [12311, 12312, 12313, 12314, 12315, 12316]
    base_trials = 700
    trials = base_trials * trial_multiplier
    horizon = 36
    noise_amp = 0.02
    n_eff = max(1, len(seeds) * trials * 4)
    adaptive_margin = math.sqrt(math.log(2.0 / delta) / (2.0 * n_eff))
    rows: list[dict[str, Any]] = []
    for cert in cert_rows:
        risk = float(cert.get("risk_tolerance", 0.05) or 0.05)
        target = max(0.0, 1.0 - risk)
        box = cert.get("certified_box", {}) or {}
        corners = [
            (float(box.get("rho_min", 0.8) or 0.8), float(box.get("kappa_scale_min", 1.0) or 1.0)),
            (float(box.get("rho_min", 0.8) or 0.8), float(box.get("kappa_scale_max", 1.0) or 1.0)),
            (float(box.get("rho_max", 0.8) or 0.8), float(box.get("kappa_scale_min", 1.0) or 1.0)),
            (float(box.get("rho_max", 0.8) or 0.8), float(box.get("kappa_scale_max", 1.0) or 1.0)),
        ]
        corner_rows = []
        corner_mins = []
        for rho, kappa in corners:
            samples = [
                simulate_pass_rate(rho, kappa, seed, trials=trials, horizon=horizon, noise_amp=noise_amp, provider_lift=provider_lift)
                for seed in seeds
            ]
            sample_min = min(samples)
            corner_mins.append(sample_min)
            corner_rows.append({
                "rho": rho,
                "kappa_scale": kappa,
                "samples": samples,
                "sample_min": sample_min,
                "sample_mean": sum(samples) / len(samples),
            })
        empirical_floor = min(corner_mins) if corner_mins else 0.0
        certified_floor = max(0.0, empirical_floor - adaptive_margin)
        rows.append({
            "risk_tolerance": risk,
            "target_floor": target,
            "provider_lift_per_step": provider_lift,
            "empirical_corner_floor": empirical_floor,
            "adaptive_margin": adaptive_margin,
            "predicted_certified_floor": certified_floor,
            "meets_target": certified_floor + 1e-12 >= target,
            "margin_to_target": certified_floor - target,
            "corner_rows": corner_rows,
        })
    summary = {
        "all_rows_meet_target": all(row["meets_target"] for row in rows),
        "worst_margin_to_target": min((row["margin_to_target"] for row in rows), default=0.0),
        "delta": delta,
        "trial_multiplier": trial_multiplier,
        "trials": trials,
        "n_eff": n_eff,
        "adaptive_margin": adaptive_margin,
    }
    return rows, summary


def main() -> None:
    GEN.mkdir(exist_ok=True)
    alpha = load(IN_ALPHA)
    p1979 = load(IN_1979)
    p1988 = load(IN_1988)
    p2274 = load(IN_2274)
    p2300 = load(IN_2300)
    p2302 = load(IN_2302)
    p2309 = load(IN_2309)

    p2300_probe = p2300.get("strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe", {}) or {}
    p2302_probe = p2302.get("strict_task3_provider_lift_policy_lock_candidate_probe", {}) or {}
    p2309_probe = p2309.get("strict_min_norm_response_weights_replay_quarantine_probe", {}) or {}
    basis = p2300_probe.get("provider_basis", []) or []
    coeffs = [float(row.get("canonical_coefficient_numeric", 0.0) or 0.0) for row in basis]
    names = [str(row.get("name", f"c{i}")) for i, row in enumerate(basis)]
    required_lift = float((p2302_probe.get("policy_lock_candidate", {}) or {}).get("provider_lift_per_step", 0.0) or 0.0)
    cert_rows = (p2274.get("strict_nu_branch_group_policy_robustness_region_certificate_probe", {}) or {}).get("certified_rows", []) or []

    coeff_norm_sq = dot(coeffs, coeffs)
    signed_total = sum(coeffs)
    positive_sum = sum(c for c in coeffs if c > 0)
    p1988_severity = float((p1988.get("numeric_severity", {}) or {}).get("p1987_probe_l2_norm", 0.0) or 0.0)

    c_symbols = sp.symbols("c0:10", real=True)
    self_energy_expr = sp.Rational(1, 2) * sum(c * c for c in c_symbols)
    self_energy_gradient = [sp.diff(self_energy_expr, c) for c in c_symbols]

    delta = float((p2302_probe.get("policy_lock_candidate", {}) or {}).get("delta", 0.1) or 0.1)
    trial_multiplier = int((p2302_probe.get("policy_lock_candidate", {}) or {}).get("trial_multiplier", 1) or 1)
    replay_rows, replay_summary = replay_rows_for_lift(cert_rows, coeff_norm_sq, delta=delta, trial_multiplier=trial_multiplier)

    source_candidates = [
        {
            "candidate_id": "ADM_COEFFICIENT_SELF_ENERGY_NORM_SQUARED",
            "definition": "lambda_candidate = ||c||_2^2 using the P2300 coefficient vector",
            "lambda_numeric": coeff_norm_sq,
            "passes_required_lift": coeff_norm_sq >= required_lift,
            "target_calibrated": False,
            "weights_or_gradient": [float(c) for c in coeffs],
            "replay_all_rows_meet_target": replay_summary["all_rows_meet_target"],
            "strict_admissible": False,
            "failure_reason": "The self-energy scalar is generated from P2300 coefficients, but no ADM/Bianchi-I lapse/shear or policy-margin theorem exports ||c||^2 as provider_lift_per_step.",
        },
        {
            "candidate_id": "P1988_PROJECTION_SEVERITY_NORM",
            "definition": "lambda_candidate = P1988 numeric non-GB projection severity",
            "lambda_numeric": p1988_severity,
            "passes_required_lift": p1988_severity >= required_lift,
            "target_calibrated": False,
            "weights_or_gradient": None,
            "replay_all_rows_meet_target": None,
            "strict_admissible": False,
            "failure_reason": "Projection severity is a diagnostic residual norm, not a signed monotone provider lift.",
        },
        {
            "candidate_id": "SIGNED_TOTAL_COEFFICIENT_SOURCE",
            "definition": "lambda_candidate = sum_i c_i",
            "lambda_numeric": signed_total,
            "passes_required_lift": signed_total >= required_lift,
            "target_calibrated": False,
            "weights_or_gradient": [1.0 for _ in coeffs],
            "replay_all_rows_meet_target": None,
            "strict_admissible": False,
            "failure_reason": "Signed total is negative and cannot certify a positive policy lift.",
        },
        {
            "candidate_id": "POSITIVE_CHANNEL_SOURCE",
            "definition": "lambda_candidate = sum_{i:c_i>0} c_i",
            "lambda_numeric": positive_sum,
            "passes_required_lift": positive_sum >= required_lift,
            "target_calibrated": False,
            "weights_or_gradient": [1.0 if c > 0 else 0.0 for c in coeffs],
            "replay_all_rows_meet_target": None,
            "strict_admissible": False,
            "failure_reason": "Positive-channel selection is a sign/channel-selection premise, not derived from strict ADM/Bianchi-I equations.",
        },
    ]

    adm_source_audit = {
        "p1979_status": p1979.get("status"),
        "p1979_no_exported_lapse_shear_provider_certificate": (p1979.get("gatekeeper_checks", {}) or {}).get("no_exported_lapse_shear_provider_certificate_in_current_registries"),
        "p1979_false_pass_guard": p1979.get("false_pass_guard"),
        "p1988_status": p1988.get("status"),
        "p1988_obstruction_tags": p1988.get("obstruction_tags", []),
        "p2300_scope_reason": (p2300_probe.get("task3_g1_g2_g3_impact", {}) or {}).get("reason"),
        "conclusion": "Current strict ADM/Bianchi-I exports a spatial-EOM provider coefficient vector, but not a lapse/shear or policy-margin identity mapping coefficient self-energy to provider_lift_per_step.",
    }

    theorem_export = {
        "statement_id": "P2310_SELF_ENERGY_RESPONSE_SOURCE_AUDIT_THEOREM",
        "formal_statement": (
            "A non-target-calibrated self-energy candidate lambda=||c||_2^2 can be formed from the P2300 coefficients and numerically passes "
            "the P2302 replay.  However the current ADM/Bianchi-I export set does not prove that this self-energy is a signed monotone "
            "policy-margin source: P1979 still reports no exported lapse/shear provider certificate, P1988 is a residual projection witness, "
            "and P2300 is scoped to a spatial-EOM provider matrix pass.  Therefore the candidate is quarantined and G1/G3 remain open."
        ),
        "proof_bits": {
            "required_lift": required_lift,
            "coefficient_names": names,
            "self_energy_lambda": coeff_norm_sq,
            "self_energy_replay_all_rows_meet_target": replay_summary["all_rows_meet_target"],
            "self_energy_replay_worst_margin_to_target": replay_summary["worst_margin_to_target"],
            "self_energy_gradient_symbolic": [str(x) for x in self_energy_gradient],
            "p1979_no_exported_lapse_shear_provider_certificate": adm_source_audit["p1979_no_exported_lapse_shear_provider_certificate"],
            "strict_admissible_candidates": [row["candidate_id"] for row in source_candidates if row["strict_admissible"]],
        },
        "not_claimed": [
            "strict ADM/Bianchi-I derivation of lambda=R(c)",
            "strict admissible self-energy-to-margin theorem",
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
        "p2300_coefficients_loaded": len(coeffs) == 10 and coeff_norm_sq > 0.0,
        "p2302_required_lift_loaded": required_lift == 0.0068,
        "p2309_quarantine_loaded": (p2309_probe.get("strict_admissibility_verdict", {}) or {}).get("admissible_for_g1_g3_update") is False,
        "self_energy_candidate_not_target_calibrated": source_candidates[0]["target_calibrated"] is False,
        "self_energy_candidate_passes_replay": replay_summary["all_rows_meet_target"] is True,
        "p1979_blocks_adm_lapse_shear_source_export": adm_source_audit["p1979_no_exported_lapse_shear_provider_certificate"] is True,
        "no_source_candidate_strict_admissible": not any(row["strict_admissible"] for row in source_candidates),
        "strict_self_energy_bridge_not_claimed": True,
        "strict_g1_g3_not_updated": True,
        "no_qw2191_discharge_claimed": True,
        "no_selector_closure_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2310_s1260_v1",
        "packet_id": "P2310",
        "stage_id": "S1260",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_OBSTRUCTION_SELF_ENERGY_SOURCE_CANDIDATE_REPLAY_PASS_QUARANTINED",
        "result_kind": "STRICT_SELF_ENERGY_LAMBDA_CANDIDATE_NUMERICALLY_PASSES_BUT_ADM_MARGIN_SOURCE_NOT_EXPORTED",
        "strict_self_energy_response_source_audit_and_replay_probe": {
            "probe_id": "P2310_S1260_STRICT_SELF_ENERGY_RESPONSE_SOURCE_AUDIT_AND_REPLAY",
            "source_packets": {
                "alpha_geo_strict_derived_v1": "generated/alpha_geo_strict_derived_v1.json",
                "p1979": "generated/p1979_s929_strict_lapse_shear_provider_export_audit.json",
                "p1988": "generated/p1988_s938_strict_non_gb_to_spatial_eom_family_projection_witness.json",
                "p2274": "generated/p2274_s1224_strict_nu_branch_group_policy_robustness_region_certificate_probe.json",
                "p2300": "generated/p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json",
                "p2302": "generated/p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json",
                "p2309": "generated/p2309_s1259_strict_min_norm_response_weights_replay_quarantine_probe.json",
            },
            "source_hashes": {
                "alpha_geo_strict_derived_v1_sha256": sha256_file(IN_ALPHA),
                "p1979_sha256": sha256_file(IN_1979),
                "p1988_sha256": sha256_file(IN_1988),
                "p2274_sha256": sha256_file(IN_2274),
                "p2300_sha256": sha256_file(IN_2300),
                "p2302_sha256": sha256_file(IN_2302),
                "p2309_sha256": sha256_file(IN_2309),
            },
            "self_energy_variational_attempt": {
                "functional": "E(c)=1/2 ||c||_2^2",
                "gradient_symbolic": [str(x) for x in self_energy_gradient],
                "lambda_candidate": "||c||_2^2",
                "lambda_numeric": coeff_norm_sq,
                "target_calibrated": False,
                "strict_admissible": False,
            },
            "source_candidates": source_candidates,
            "self_energy_replay": {
                "summary": replay_summary,
                "rows": replay_rows,
            },
            "adm_source_audit": adm_source_audit,
            "strict_admissibility_verdict": {
                "non_target_lambda_candidate_found": True,
                "non_target_lambda_candidate_replay_passes": replay_summary["all_rows_meet_target"],
                "adm_bianchi_margin_source_exported": False,
                "admissible_for_g1_g3_update": False,
                "reason": "Self-energy is a useful candidate scalar, but the current strict ADM/Bianchi-I artifacts do not export it as a signed monotone policy-margin source.",
            },
            "task3_g1_g3_update": {
                "G1_reduction_certainty": "OPEN_UNCHANGED",
                "G2_nonlinear_trajectory_realism": "CLOSED_FROM_P2301_UNCHANGED",
                "G3_operational_policy_rule": "OPEN_UNCHANGED",
                "reason": "P2310 finds a replay-passing self-energy candidate, but no strict ADM margin-source theorem.",
            },
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2311_candidate",
            "goal": "Derive or refute an ADM lapse/shear energy-to-policy-margin theorem for E(c)=||c||^2. If this theorem is proven, replay P2302 and only then update G1/G3; otherwise keep P2310 quarantined.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_SELF_ENERGY_REPLAY_PASS_QUARANTINE_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2310/S1260 — strict self-energy response source audit and replay",
            "",
            f"- Status: `{payload['status']}`",
            f"- Required lift: `{required_lift}`",
            f"- Self-energy lambda candidate `||c||^2`: `{coeff_norm_sq}`",
            f"- Self-energy replay all rows meet target: `{replay_summary['all_rows_meet_target']}`",
            f"- Self-energy replay worst margin: `{replay_summary['worst_margin_to_target']}`",
            "- Strict admissible: `False` because no ADM/Bianchi-I energy-to-policy-margin theorem is exported.",
            "- G1/G3 update: `OPEN_UNCHANGED`",
            f"- Theorem fingerprint: `{theorem_fingerprint}`",
            "",
            "## Guardrail statement",
            "P2310 does not close G1/G3, does not discharge QW-2191, does not add a selector premise, does not transfer legacy-kernel roles, and does not claim ToE closure.",
            "",
            "## Next honest step",
            payload["recommended_next_honest_step"]["goal"],
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
