#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import random
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_ALPHA = GEN / "alpha_geo_strict_derived_v1.json"
IN_2274 = GEN / "p2274_s1224_strict_nu_branch_group_policy_robustness_region_certificate_probe.json"
IN_2280 = GEN / "p2280_s1230_strict_nu_branch_group_policy_minimal_lock_criterion_grid_probe.json"
IN_2281 = GEN / "p2281_s1231_strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe.json"
IN_2300 = GEN / "p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json"
IN_2301 = GEN / "p2301_s1251_strict_task3_provider_corrected_g1_g2_g3_replay_probe.json"
OUT = GEN / "p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json"
MD = GEN / "p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.md"


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
    p2280 = load(IN_2280)
    p2281 = load(IN_2281)
    p2300 = load(IN_2300)
    p2301 = load(IN_2301)

    cert_rows = (p2274.get("strict_nu_branch_group_policy_robustness_region_certificate_probe", {}) or {}).get("certified_rows", []) or []
    p2300_probe = p2300.get("strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe", {}) or {}
    p2300_solution = p2300_probe.get("solution_space", {}) or {}
    p2301_probe = p2301.get("strict_task3_provider_corrected_g1_g2_g3_replay_probe", {}) or {}
    p2301_statuses = [row.get("status") for row in p2301_probe.get("recomputed_gap_rows", [])]

    delta = 0.1
    trial_multiplier = 1
    # Coarse-to-refined deterministic search around the first viable lift window.
    # Earlier P2281/P2301 rows show no lift and the 0.004/0.006 sentinels fail;
    # the 0.0002 refinement then resolves the first candidate at 0.0068.
    lift_grid = [0.0, 0.004, 0.006, 0.0062, 0.0064, 0.0066, 0.0068, 0.0070]
    scan_rows = []
    feasible_replay: dict[str, Any] | None = None
    feasible_rows: list[dict[str, Any]] = []
    for provider_lift in lift_grid:
        rows, summary = replay_rows_for_lift(cert_rows, provider_lift, delta=delta, trial_multiplier=trial_multiplier)
        scan_row = {
            "provider_lift_per_step": provider_lift,
            "all_rows_meet_target": summary["all_rows_meet_target"],
            "worst_margin_to_target": summary["worst_margin_to_target"],
            "min_empirical_corner_floor": min((row["empirical_corner_floor"] for row in rows), default=0.0),
        }
        scan_rows.append(scan_row)
        if summary["all_rows_meet_target"] and feasible_replay is None:
            feasible_replay = {"rows": rows, "summary": summary}
            feasible_rows = rows
            break

    minimal_lift = feasible_replay["summary"]["delta"] and feasible_rows[0]["provider_lift_per_step"] if feasible_replay else None
    p2280_probe = p2280.get("strict_nu_branch_group_policy_minimal_lock_criterion_grid_probe", {}) or {}
    p2281_probe = p2281.get("strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe", {}) or {}
    baseline_worst_margin = (p2281_probe.get("global_summary", {}) or {}).get("worst_margin_to_target")
    baseline_feasible_count = int(p2280_probe.get("feasible_count", 0) or 0)

    conditional_gap_rows = [
        {
            "id": "G1_reduction_certainty",
            "status": "CONDITIONALLY_CLOSED_UNDER_PROVIDER_MARGIN_LIFT" if feasible_replay else "OPEN",
            "metric": feasible_replay["summary"]["worst_margin_to_target"] if feasible_replay else baseline_worst_margin,
            "criterion": "provider-lift replay rows meet all P2281 target floors",
        },
        {
            "id": "G2_nonlinear_trajectory_realism",
            "status": "CLOSED_FROM_P2301",
            "metric": (p2301_probe.get("provider_corrected_transport_replay", {}) or {}).get("provider_corrected_max_transport_residual_l1"),
            "criterion": "P2301 provider-corrected transport residual below threshold",
        },
        {
            "id": "G3_operational_policy_rule",
            "status": "CONDITIONALLY_CLOSED_UNDER_PROVIDER_MARGIN_LIFT" if feasible_replay else "OPEN",
            "metric": minimal_lift if feasible_replay else -1.0,
            "criterion": "minimal provider_lift_per_step exists on the P2280/P2281 replay contract",
        },
    ]

    policy_lock_candidate = {
        "candidate_id": "P2302_PROVIDER_MARGIN_LIFT_POLICY_LOCK_CANDIDATE_V1",
        "status": "CONDITIONAL_CANDIDATE_REQUIRES_PROVIDER_TO_MARGIN_BRIDGE" if feasible_replay else "NO_FEASIBLE_CANDIDATE_FOUND",
        "delta": delta,
        "trial_multiplier": trial_multiplier,
        "provider_lift_per_step": minimal_lift,
        "cost_proxy": (trial_multiplier * math.log(2.0 / delta) + float(minimal_lift or 0.0) * 1000.0) if feasible_replay else None,
        "bridge_obligation": "prove that P2300 canonical spatial-EOM provider coefficients induce at least this provider_lift_per_step in the P2281 margin dynamics without adding a selector premise",
    }

    theorem_export = {
        "statement_id": "P2302_PROVIDER_MARGIN_LIFT_POLICY_LOCK_CANDIDATE_THEOREM",
        "formal_statement": (
            "A deterministic search over a provider-lift augmentation of the P2281 margin replay finds a minimal candidate "
            "provider_lift_per_step that would make all P2281 target floors pass and would supply a conditional policy-lock "
            "candidate for G3.  This is not yet strict Task-3 closure: the required provider-to-margin bridge from the P2300 "
            "spatial-EOM coefficients to the replay lift is still an explicit proof obligation."
        ),
        "proof_bits": {
            "baseline_p2280_feasible_count": baseline_feasible_count,
            "baseline_p2281_worst_margin": baseline_worst_margin,
            "p2301_gap_pattern": p2301_statuses,
            "conditional_candidate_found": feasible_replay is not None,
            "minimal_provider_lift_per_step": minimal_lift,
            "conditional_worst_margin": feasible_replay["summary"]["worst_margin_to_target"] if feasible_replay else None,
            "provider_to_margin_bridge_proven": False,
        },
        "not_claimed": [
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
        "p2300_canonical_coefficients_loaded": len(p2300_solution.get("canonical_solution", []) or []) == 10,
        "p2301_only_g2_closed_input": p2301_statuses == ["OPEN", "CLOSED", "OPEN"],
        "provider_lift_candidate_found": feasible_replay is not None,
        "conditional_replay_all_rows_meet_target": bool(feasible_replay and feasible_replay["summary"]["all_rows_meet_target"]),
        "baseline_p2280_had_no_feasible_lock": baseline_feasible_count == 0,
        "provider_to_margin_bridge_still_open": True,
        "strict_task3_closure_not_claimed": True,
        "no_qw2191_discharge_claimed": True,
        "no_selector_closure_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2302_s1252_v1",
        "packet_id": "P2302",
        "stage_id": "S1252",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_BRIDGE_OBLIGATION_WITH_CONDITIONAL_G1_G3_POLICY_LOCK_CANDIDATE_TRACE",
        "result_kind": "STRICT_TASK3_PROVIDER_MARGIN_LIFT_POLICY_LOCK_CANDIDATE_WITH_OPEN_BRIDGE_OBLIGATION",
        "strict_task3_provider_lift_policy_lock_candidate_probe": {
            "probe_id": "P2302_S1252_STRICT_TASK3_PROVIDER_LIFT_POLICY_LOCK_CANDIDATE",
            "source_packets": {
                "alpha_geo_strict_derived_v1": "generated/alpha_geo_strict_derived_v1.json",
                "p2274": "generated/p2274_s1224_strict_nu_branch_group_policy_robustness_region_certificate_probe.json",
                "p2280": "generated/p2280_s1230_strict_nu_branch_group_policy_minimal_lock_criterion_grid_probe.json",
                "p2281": "generated/p2281_s1231_strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe.json",
                "p2300": "generated/p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json",
                "p2301": "generated/p2301_s1251_strict_task3_provider_corrected_g1_g2_g3_replay_probe.json",
            },
            "source_hashes": {
                "alpha_geo_strict_derived_v1_sha256": sha256_file(IN_ALPHA),
                "p2274_sha256": sha256_file(IN_2274),
                "p2280_sha256": sha256_file(IN_2280),
                "p2281_sha256": sha256_file(IN_2281),
                "p2300_sha256": sha256_file(IN_2300),
                "p2301_sha256": sha256_file(IN_2301),
            },
            "policy_lock_candidate": policy_lock_candidate,
            "provider_lift_scan_rows": scan_rows,
            "conditional_fresh_replay": feasible_replay,
            "conditional_gap_rows": conditional_gap_rows,
            "strict_closure_status": "HELD_OPEN_UNTIL_PROVIDER_TO_MARGIN_BRIDGE_PROVEN",
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2303_candidate",
            "goal": "Prove or refute the provider-to-margin bridge: derive the P2302 provider_lift_per_step bound directly from the P2300 canonical spatial-EOM coefficients without adding a selector premise; only then update strict G1/G3.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_CONDITIONAL_POLICY_LOCK_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2302/S1252 — provider-lift policy-lock candidate",
            "",
            f"- Status: `{payload['status']}`",
            f"- Conditional candidate found: `{feasible_replay is not None}`",
            f"- Minimal provider_lift_per_step: `{minimal_lift}`",
            f"- Conditional gap rows: `{[(row['id'], row['status']) for row in conditional_gap_rows]}`",
            "- Strict closure status: `HELD_OPEN_UNTIL_PROVIDER_TO_MARGIN_BRIDGE_PROVEN`",
            f"- Theorem fingerprint: `{theorem_fingerprint}`",
            "",
            "## Guardrail statement",
            "P2302 finds a numerical conditional G1/G3 policy-lock candidate, but it does not prove the provider-to-margin bridge from P2300 coefficients. Therefore strict Task-3 closure remains open; no QW-2191, selector, legacy-kernel, or ToE closure is claimed.",
            "",
            "## Next honest step",
            payload["recommended_next_honest_step"]["goal"],
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
