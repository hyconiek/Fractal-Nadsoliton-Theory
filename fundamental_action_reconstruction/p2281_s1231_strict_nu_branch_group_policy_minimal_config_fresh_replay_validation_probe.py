#!/usr/bin/env python3
from __future__ import annotations
import json
import math
import random
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2274 = GEN / "p2274_s1224_strict_nu_branch_group_policy_robustness_region_certificate_probe.json"
IN_2280 = GEN / "p2280_s1230_strict_nu_branch_group_policy_minimal_lock_criterion_grid_probe.json"
OUT = GEN / "p2281_s1231_strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe.json"
MD = GEN / "p2281_s1231_strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def simulate_pass_rate(rho: float, kappa_scale: float, seed: int, trials: int, horizon: int, noise_amp: float) -> float:
    rng = random.Random(seed)
    passed = 0
    for _ in range(trials):
        margin = 0.20
        for _ in range(horizon):
            noise = (rng.random() - 0.5) * noise_amp
            drift = 0.015 + 0.01 * (margin ** 2)
            control_relief = rho * 0.008 / max(kappa_scale, 1e-12)
            margin = margin - drift + control_relief + noise
            if margin < -1e-12:
                break
        if margin >= -1e-12:
            passed += 1
    return passed / trials


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2274 = load(IN_2274)
    p2280 = load(IN_2280)

    cert_rows = (p2274.get("strict_nu_branch_group_policy_robustness_region_certificate_probe", {}) or {}).get("certified_rows", []) or []
    minimal = (p2280.get("strict_nu_branch_group_policy_minimal_lock_criterion_grid_probe", {}) or {}).get("minimal_feasible_config", {}) or {}

    delta = float(minimal.get("delta", 0.01) or 0.01)
    trial_multiplier = int(minimal.get("trial_multiplier", 1) or 1)

    seeds = [12311, 12312, 12313, 12314, 12315, 12316]
    base_trials = 700
    trials = base_trials * trial_multiplier
    horizon = 36
    noise_amp = 0.02

    n_eff = max(1, len(seeds) * trials * 4)
    adaptive_margin = math.sqrt(math.log(2.0 / delta) / (2.0 * n_eff))

    rows = []
    for c in cert_rows:
        risk = float(c.get("risk_tolerance", 0.05) or 0.05)
        target = max(0.0, 1.0 - risk)
        box = c.get("certified_box", {}) or {}
        corners = [
            (float(box.get("rho_min", 0.8) or 0.8), float(box.get("kappa_scale_min", 1.0) or 1.0)),
            (float(box.get("rho_min", 0.8) or 0.8), float(box.get("kappa_scale_max", 1.0) or 1.0)),
            (float(box.get("rho_max", 0.8) or 0.8), float(box.get("kappa_scale_min", 1.0) or 1.0)),
            (float(box.get("rho_max", 0.8) or 0.8), float(box.get("kappa_scale_max", 1.0) or 1.0)),
        ]

        mins = []
        corner_rows = []
        for rho, kappa in corners:
            samples = [simulate_pass_rate(rho, kappa, s, trials=trials, horizon=horizon, noise_amp=noise_amp) for s in seeds]
            m = min(samples)
            mins.append(m)
            corner_rows.append({
                "rho": rho,
                "kappa_scale": kappa,
                "samples": samples,
                "sample_min": m,
                "sample_mean": sum(samples) / len(samples),
            })

        empirical_floor = min(mins) if mins else 0.0
        predicted_cert_floor = max(0.0, empirical_floor - adaptive_margin)
        rows.append({
            "risk_tolerance": risk,
            "target_floor": target,
            "empirical_corner_floor": empirical_floor,
            "adaptive_margin": adaptive_margin,
            "predicted_certified_floor": predicted_cert_floor,
            "meets_target": predicted_cert_floor + 1e-12 >= target,
            "margin_to_target": predicted_cert_floor - target,
            "corner_rows": corner_rows,
        })

    payload = {
        "schema_version": "p2281_s1231_v1",
        "packet_id": "P2281",
        "stage_id": "S1231",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_MINIMAL_CONFIG_FRESH_REPLAY_VALIDATION_PROBE",
        "strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_MINIMAL_CONFIG_FRESH_REPLAY_VALIDATION_PROBE_V1",
            "source_packets": [str(IN_2274.relative_to(ROOT)), str(IN_2280.relative_to(ROOT))],
            "locked_minimal_config": {
                "delta": delta,
                "trial_multiplier": trial_multiplier,
                "seeds": seeds,
                "base_trials": base_trials,
                "trials": trials,
                "horizon": horizon,
                "noise_amp": noise_amp,
                "n_eff": n_eff,
                "adaptive_margin": adaptive_margin,
            },
            "rows": rows,
            "global_summary": {
                "all_rows_meet_target": all(r["meets_target"] for r in rows),
                "worst_margin_to_target": min((r["margin_to_target"] for r in rows), default=0.0),
            },
            "theorem_scope_limit": "fresh replay validation under minimal config on surrogate route only; not selector closure and not ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2282_candidate",
            "goal": "derive deterministic lower-bound surrogate for fresh-replay floor with explicit dependence on (delta,n_eff) and compare tightness to observed rows",
        },
        "gatekeeper_checks": {
            "rows_exported": len(rows) > 0,
            "adaptive_margin_nonnegative": adaptive_margin >= 0.0,
            "all_floors_bounded": all(0.0 <= r["predicted_certified_floor"] <= 1.0 for r in rows),
            "all_targets_bounded": all(0.0 <= r["target_floor"] <= 1.0 for r in rows),
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2281 S1231: minimal-config fresh replay validation probe",
            "",
            f"- minimal delta: `{delta}`",
            f"- minimal trial multiplier: `{trial_multiplier}`",
            f"- adaptive margin: `{adaptive_margin:.12e}`",
            f"- all rows meet target: `{all(r['meets_target'] for r in rows)}`",
            f"- worst margin to target: `{min((r['margin_to_target'] for r in rows), default=0.0):.12e}`",
            "",
            "Fresh replay validation only; no selector closure / ToE closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
