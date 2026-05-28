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
IN_2274 = GEN / "p2274_s1224_strict_nu_branch_group_policy_robustness_region_certificate_probe.json"
IN_2300 = GEN / "p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json"
IN_2302 = GEN / "p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json"
IN_2308 = GEN / "p2308_s1258_strict_current_interface_class_response_functional_nonexistence_probe.json"
OUT = GEN / "p2309_s1259_strict_min_norm_response_weights_replay_quarantine_probe.json"
MD = GEN / "p2309_s1259_strict_min_norm_response_weights_replay_quarantine_probe.md"


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
    p2274 = load(IN_2274)
    p2300 = load(IN_2300)
    p2302 = load(IN_2302)
    p2308 = load(IN_2308)

    p2300_probe = p2300.get("strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe", {}) or {}
    p2302_probe = p2302.get("strict_task3_provider_lift_policy_lock_candidate_probe", {}) or {}
    p2308_probe = p2308.get("strict_current_interface_class_response_functional_nonexistence_probe", {}) or {}
    basis = p2300_probe.get("provider_basis", []) or []
    coeffs = [float(row.get("canonical_coefficient_numeric", 0.0) or 0.0) for row in basis]
    names = [str(row.get("name", f"c{i}")) for i, row in enumerate(basis)]
    required_lift = float((p2302_probe.get("policy_lock_candidate", {}) or {}).get("provider_lift_per_step", 0.0) or 0.0)
    cert_rows = (p2274.get("strict_nu_branch_group_policy_robustness_region_certificate_probe", {}) or {}).get("certified_rows", []) or []

    c_norm_sq = dot(coeffs, coeffs)
    weights = [required_lift * c / c_norm_sq for c in coeffs] if c_norm_sq > 0 else [0.0 for _ in coeffs]
    induced_lift = dot(weights, coeffs)
    weight_l2 = math.sqrt(dot(weights, weights))

    # Lagrange multiplier derivation for min ||w||^2 subject to w·c = lambda_star.
    mu = sp.symbols("mu", real=True)
    lambda_star = sp.symbols("lambda_star", real=True)
    c_symbols = sp.symbols("c0:10", real=True)
    norm_sq_symbolic = sum(ci**2 for ci in c_symbols)
    lagrange_solution = [sp.simplify(lambda_star * ci / norm_sq_symbolic) for ci in c_symbols]

    delta = float((p2302_probe.get("policy_lock_candidate", {}) or {}).get("delta", 0.1) or 0.1)
    trial_multiplier = int((p2302_probe.get("policy_lock_candidate", {}) or {}).get("trial_multiplier", 1) or 1)
    replay_rows, replay_summary = replay_rows_for_lift(cert_rows, induced_lift, delta=delta, trial_multiplier=trial_multiplier)

    variational_identity_attempt = {
        "identity_id": "TARGET_CALIBRATED_MIN_NORM_WEIGHTS",
        "optimization_problem": "minimize ||w||_2^2 subject to w·c = lambda_star",
        "symbolic_solution": [str(x) for x in lagrange_solution],
        "lambda_star_source": "P2302 required policy-lock lift, not independently derived from strict ADM/Bianchi-I dynamics",
        "coefficient_norm_sq": c_norm_sq,
        "required_lift": required_lift,
        "weights": [
            {"name": name, "coefficient_numeric": coeff, "weight_numeric": weight}
            for name, coeff, weight in zip(names, coeffs, weights)
        ],
        "induced_lift": induced_lift,
        "weight_l2_norm": weight_l2,
        "weights_exported": True,
        "strict_admissible": False,
        "nonadmissibility_reason": "Weights are target-calibrated using the already-required lift; this is a tautological projection, not an internally derived variational identity from the strict dynamics.",
    }

    theorem_export = {
        "statement_id": "P2309_TARGET_CALIBRATED_RESPONSE_WEIGHTS_QUARANTINE_THEOREM",
        "formal_statement": (
            "The minimum-norm constrained identity min ||w||_2^2 subject to w·c=lambda_star exports concrete weights for the P2300 "
            "coefficient vector and, when lambda_star is set to the P2302 required lift, the replay passes.  However lambda_star is "
            "imported from the target policy-lock requirement rather than derived from strict ADM/Bianchi-I dynamics.  Therefore the "
            "weights are quarantined as target-calibrated and cannot close G1/G3."
        ),
        "proof_bits": {
            "required_lift": required_lift,
            "induced_lift": induced_lift,
            "lift_error_abs": abs(induced_lift - required_lift),
            "replay_all_rows_meet_target": replay_summary["all_rows_meet_target"],
            "replay_worst_margin_to_target": replay_summary["worst_margin_to_target"],
            "strict_admissible": variational_identity_attempt["strict_admissible"],
            "target_calibrated": True,
        },
        "not_claimed": [
            "strict derivation of lambda=R(c)",
            "strict admissible response-functional bridge",
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
        "p2300_coefficients_loaded": len(coeffs) == 10 and c_norm_sq > 0.0,
        "p2302_required_lift_loaded": required_lift == 0.0068,
        "p2308_current_class_nonexistence_loaded": (p2308_probe.get("strict_nonexistence_verdict", {}) or {}).get("nonexistence_proven_for_current_interface_class") is True,
        "min_norm_weights_exported": len(weights) == len(coeffs) and len(weights) > 0,
        "induced_lift_matches_required_lift": abs(induced_lift - required_lift) < 1e-12,
        "replay_passes_if_target_calibrated_weights_are_admitted": replay_summary["all_rows_meet_target"] is True,
        "weights_are_target_calibrated_not_strictly_derived": True,
        "strict_admissible_response_bridge_not_claimed": variational_identity_attempt["strict_admissible"] is False,
        "strict_g1_g3_not_updated": True,
        "no_qw2191_discharge_claimed": True,
        "no_selector_closure_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2309_s1259_v1",
        "packet_id": "P2309",
        "stage_id": "S1259",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_OBSTRUCTION_TARGET_CALIBRATED_WEIGHTS_REPLAY_QUARANTINED",
        "result_kind": "MIN_NORM_RESPONSE_WEIGHTS_EXPORT_REPLAY_PASS_BUT_STRICT_ADMISSIBILITY_FAILS",
        "strict_min_norm_response_weights_replay_quarantine_probe": {
            "probe_id": "P2309_S1259_STRICT_MIN_NORM_RESPONSE_WEIGHTS_REPLAY_QUARANTINE",
            "source_packets": {
                "alpha_geo_strict_derived_v1": "generated/alpha_geo_strict_derived_v1.json",
                "p2274": "generated/p2274_s1224_strict_nu_branch_group_policy_robustness_region_certificate_probe.json",
                "p2300": "generated/p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json",
                "p2302": "generated/p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json",
                "p2308": "generated/p2308_s1258_strict_current_interface_class_response_functional_nonexistence_probe.json",
            },
            "source_hashes": {
                "alpha_geo_strict_derived_v1_sha256": sha256_file(IN_ALPHA),
                "p2274_sha256": sha256_file(IN_2274),
                "p2300_sha256": sha256_file(IN_2300),
                "p2302_sha256": sha256_file(IN_2302),
                "p2308_sha256": sha256_file(IN_2308),
            },
            "variational_identity_attempt": variational_identity_attempt,
            "target_calibrated_replay": {
                "summary": replay_summary,
                "rows": replay_rows,
            },
            "strict_admissibility_verdict": {
                "lambda_equals_R_of_c_exported": True,
                "lambda_equals_R_of_c_strictly_derived": False,
                "replay_passes": replay_summary["all_rows_meet_target"],
                "admissible_for_g1_g3_update": False,
                "reason": "The only exported weights are calibrated to the P2302 target lift, so the replay pass is circular evidence rather than a strict derivation.",
            },
            "task3_g1_g3_update": {
                "G1_reduction_certainty": "OPEN_UNCHANGED",
                "G2_nonlinear_trajectory_realism": "CLOSED_FROM_P2301_UNCHANGED",
                "G3_operational_policy_rule": "OPEN_UNCHANGED",
                "reason": "P2309 replay passes only under quarantined target-calibrated weights.",
            },
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2310_candidate",
            "goal": "Replace the target-calibrated constraint w·c=lambda_star with a strict variational source for lambda_star or w_i from ADM/Bianchi-I equations; otherwise keep P2308/P2309 as blockers for G1/G3.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_TARGET_CALIBRATED_REPLAY_QUARANTINE_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2309/S1259 — min-norm response weights replay quarantine",
            "",
            f"- Status: `{payload['status']}`",
            f"- Required lift: `{required_lift}`",
            f"- Induced lift from exported weights: `{induced_lift}`",
            f"- Replay all rows meet target: `{replay_summary['all_rows_meet_target']}`",
            f"- Replay worst margin: `{replay_summary['worst_margin_to_target']}`",
            "- Strict admissible: `False` because weights are calibrated to the target lift rather than derived from strict dynamics.",
            "- G1/G3 update: `OPEN_UNCHANGED`",
            f"- Theorem fingerprint: `{theorem_fingerprint}`",
            "",
            "## Guardrail statement",
            "P2309 does not close G1/G3, does not discharge QW-2191, does not add a selector premise, does not transfer legacy-kernel roles, and does not claim ToE closure.",
            "",
            "## Next honest step",
            payload["recommended_next_honest_step"]["goal"],
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
