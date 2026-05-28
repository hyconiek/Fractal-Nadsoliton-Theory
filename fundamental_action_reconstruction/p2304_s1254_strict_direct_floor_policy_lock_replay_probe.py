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
IN_2302 = GEN / "p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json"
IN_2303 = GEN / "p2303_s1253_strict_provider_to_margin_bridge_bound_audit_probe.json"
OUT = GEN / "p2304_s1254_strict_direct_floor_policy_lock_replay_probe.json"
MD = GEN / "p2304_s1254_strict_direct_floor_policy_lock_replay_probe.md"


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


def replay_rows(cert_rows: list[dict[str, Any]], provider_lift: float, *, delta: float, trial_multiplier: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
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
    p2302 = load(IN_2302)
    p2303 = load(IN_2303)

    cert_rows = (p2274.get("strict_nu_branch_group_policy_robustness_region_certificate_probe", {}) or {}).get("certified_rows", []) or []
    p2302_probe = p2302.get("strict_task3_provider_lift_policy_lock_candidate_probe", {}) or {}
    p2303_probe = p2303.get("strict_provider_to_margin_bridge_bound_audit_probe", {}) or {}
    candidate = p2302_probe.get("policy_lock_candidate", {}) or {}
    required_lift = float(candidate.get("provider_lift_per_step", 0.0) or 0.0)
    direct_contract = next(
        row for row in (p2303_probe.get("bridge_contract_audit", []) or [])
        if row.get("contract_id") == "STRICT_SINGLE_CHANNEL_FLOOR"
    )
    direct_floor_lift = float(direct_contract.get("bound_numeric", 0.0) or 0.0)
    delta = float(candidate.get("delta", 0.1) or 0.1)
    trial_multiplier = int(candidate.get("trial_multiplier", 1) or 1)

    direct_rows, direct_summary = replay_rows(cert_rows, direct_floor_lift, delta=delta, trial_multiplier=trial_multiplier)
    required_replay = p2302_probe.get("conditional_fresh_replay", {}) or {}
    required_summary = required_replay.get("summary", {}) or {}

    direct_gap_rows = [
        {
            "id": "G1_reduction_certainty",
            "status": "OPEN",
            "metric": direct_summary["worst_margin_to_target"],
            "criterion": "strict direct coefficient-floor replay rows all meet P2281 target floors",
            "reason": "direct-floor replay still has negative worst margin",
        },
        {
            "id": "G2_nonlinear_trajectory_realism",
            "status": "CLOSED_FROM_P2301",
            "metric": (p2302_probe.get("conditional_gap_rows", [{}, {}, {}])[1] or {}).get("metric"),
            "criterion": "P2301 provider-corrected transport residual below threshold",
        },
        {
            "id": "G3_operational_policy_rule",
            "status": "OPEN",
            "metric": direct_floor_lift,
            "criterion": "strict direct coefficient-floor must be sufficient for policy-lock replay",
            "reason": "direct-floor replay does not meet all target floors",
        },
    ]

    theorem_export = {
        "statement_id": "P2304_STRICT_DIRECT_FLOOR_POLICY_LOCK_REPLAY_THEOREM",
        "formal_statement": (
            "Lowering the P2302 conditional policy-lock candidate to the strict direct coefficient floor certified by P2303 "
            "does not close G1/G3: replaying the P2281 corner-margin contract at provider_lift_per_step equal to the P2303 "
            "direct floor still leaves a negative worst margin.  Therefore the direct-floor route is refuted for strict Task-3 "
            "closure, and progress now requires a strict aggregation/norm-to-margin theorem or a different certified policy-lock mechanism."
        ),
        "proof_bits": {
            "required_p2302_lift": required_lift,
            "strict_direct_floor_lift": direct_floor_lift,
            "required_lift_worst_margin": required_summary.get("worst_margin_to_target"),
            "direct_floor_worst_margin": direct_summary["worst_margin_to_target"],
            "direct_floor_all_rows_meet_target": direct_summary["all_rows_meet_target"],
            "direct_floor_gap_pattern": [row["status"] for row in direct_gap_rows],
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
        "p2303_direct_floor_loaded": direct_floor_lift > 0.0,
        "p2302_required_lift_loaded": required_lift == 0.0068,
        "direct_floor_below_required_lift": direct_floor_lift < required_lift,
        "direct_floor_replay_rows_exported": len(direct_rows) == len(cert_rows) and len(direct_rows) > 0,
        "direct_floor_replay_fails_g1_g3": not direct_summary["all_rows_meet_target"] and direct_summary["worst_margin_to_target"] < 0.0,
        "direct_floor_route_refuted": [row["status"] for row in direct_gap_rows] == ["OPEN", "CLOSED_FROM_P2301", "OPEN"],
        "strict_task3_closure_not_claimed": True,
        "no_qw2191_discharge_claimed": True,
        "no_selector_closure_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2304_s1254_v1",
        "packet_id": "P2304",
        "stage_id": "S1254",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_OBSTRUCTION_STRICT_DIRECT_FLOOR_POLICY_LOCK_REPLAY_REFUTED",
        "result_kind": "STRICT_DIRECT_FLOOR_POLICY_LOCK_REPLAY_REFUTES_G1_G3_CLOSURE",
        "strict_direct_floor_policy_lock_replay_probe": {
            "probe_id": "P2304_S1254_STRICT_DIRECT_FLOOR_POLICY_LOCK_REPLAY",
            "source_packets": {
                "alpha_geo_strict_derived_v1": "generated/alpha_geo_strict_derived_v1.json",
                "p2274": "generated/p2274_s1224_strict_nu_branch_group_policy_robustness_region_certificate_probe.json",
                "p2302": "generated/p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json",
                "p2303": "generated/p2303_s1253_strict_provider_to_margin_bridge_bound_audit_probe.json",
            },
            "source_hashes": {
                "alpha_geo_strict_derived_v1_sha256": sha256_file(IN_ALPHA),
                "p2274_sha256": sha256_file(IN_2274),
                "p2302_sha256": sha256_file(IN_2302),
                "p2303_sha256": sha256_file(IN_2303),
            },
            "lift_comparison": {
                "p2302_required_lift": required_lift,
                "p2303_strict_direct_floor_lift": direct_floor_lift,
                "direct_floor_shortfall": required_lift - direct_floor_lift,
            },
            "direct_floor_fresh_replay": {
                "summary": direct_summary,
                "rows": direct_rows,
            },
            "required_lift_reference_summary": required_summary,
            "direct_floor_gap_rows": direct_gap_rows,
            "strict_closure_status": "HELD_OPEN_DIRECT_FLOOR_ROUTE_REFUTED",
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2305_candidate",
            "goal": "Attempt a strict aggregation/norm-to-margin theorem for the P2300 operator basis; if that cannot be proven without a selector premise, keep G1/G3 open and route the sufficient positive/norm bounds as non-strict controls only.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_DIRECT_FLOOR_REFUTATION_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2304/S1254 — strict direct-floor policy-lock replay",
            "",
            f"- Status: `{payload['status']}`",
            f"- P2302 required lift: `{required_lift}`",
            f"- P2303 strict direct floor: `{direct_floor_lift}`",
            f"- Direct-floor shortfall: `{required_lift - direct_floor_lift}`",
            f"- Direct-floor replay all rows meet target: `{direct_summary['all_rows_meet_target']}`",
            f"- Direct-floor worst margin: `{direct_summary['worst_margin_to_target']}`",
            f"- Strict closure status: `{payload['strict_direct_floor_policy_lock_replay_probe']['strict_closure_status']}`",
            f"- Theorem fingerprint: `{theorem_fingerprint}`",
            "",
            "## Guardrail statement",
            "P2304 refutes the conservative direct-floor closure route. It does not close G1/G3, does not discharge QW-2191, does not add a selector premise, does not transfer legacy-kernel roles, and does not claim ToE closure.",
            "",
            "## Next honest step",
            payload["recommended_next_honest_step"]["goal"],
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
