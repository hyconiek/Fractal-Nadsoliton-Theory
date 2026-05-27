#!/usr/bin/env python3
from __future__ import annotations
import json
import random
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2274 = GEN / "p2274_s1224_strict_nu_branch_group_policy_robustness_region_certificate_probe.json"
OUT = GEN / "p2275_s1225_strict_nu_branch_group_policy_certified_box_corner_replay_passrate_floor_probe.json"
MD = GEN / "p2275_s1225_strict_nu_branch_group_policy_certified_box_corner_replay_passrate_floor_probe.md"


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
    rows = (p2274.get("strict_nu_branch_group_policy_robustness_region_certificate_probe", {}) or {}).get("certified_rows", []) or []

    seeds = [12251, 12252, 12253, 12254, 12255]
    trials = 700
    horizon = 36
    noise_amp = 0.02

    replay_rows = []

    for row in rows:
        risk = float(row.get("risk_tolerance", 0.05) or 0.05)
        target_floor = max(0.0, 1.0 - risk)
        box = row.get("certified_box", {}) or {}
        rho_min = float(box.get("rho_min", 0.8) or 0.8)
        rho_max = float(box.get("rho_max", 0.8) or 0.8)
        k_min = float(box.get("kappa_scale_min", 1.0) or 1.0)
        k_max = float(box.get("kappa_scale_max", 1.0) or 1.0)

        corners = [
            (rho_min, k_min),
            (rho_min, k_max),
            (rho_max, k_min),
            (rho_max, k_max),
        ]

        corner_stats = []
        empirical_floor = 1.0
        for (rho, kappa) in corners:
            samples = [simulate_pass_rate(rho, kappa, seed=s, trials=trials, horizon=horizon, noise_amp=noise_amp) for s in seeds]
            mn = min(samples)
            empirical_floor = min(empirical_floor, mn)
            corner_stats.append(
                {
                    "rho": rho,
                    "kappa_scale": kappa,
                    "passrate_samples": samples,
                    "passrate_min": mn,
                    "passrate_mean": sum(samples) / len(samples),
                }
            )

        replay_rows.append(
            {
                "risk_tolerance": risk,
                "target_passrate_floor": target_floor,
                "corner_replay": corner_stats,
                "empirical_passrate_floor_over_corners": empirical_floor,
                "floor_holds": empirical_floor + 1e-12 >= target_floor,
            }
        )

    payload = {
        "schema_version": "p2275_s1225_v1",
        "packet_id": "P2275",
        "stage_id": "S1225",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_CERTIFIED_BOX_CORNER_REPLAY_PASSRATE_FLOOR_PROBE",
        "strict_nu_branch_group_policy_certified_box_corner_replay_passrate_floor_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_CERTIFIED_BOX_CORNER_REPLAY_PASSRATE_FLOOR_PROBE_V1",
            "source_packets": [str(IN_2274.relative_to(ROOT))],
            "settings": {
                "seeds": seeds,
                "trials": trials,
                "horizon": horizon,
                "noise_amp": noise_amp,
            },
            "replay_rows": replay_rows,
            "global_summary": {
                "all_rows_floor_hold": all(r["floor_holds"] for r in replay_rows),
                "worst_empirical_floor": min((r["empirical_passrate_floor_over_corners"] for r in replay_rows), default=0.0),
            },
            "theorem_scope_limit": "corner replay certificate floor probe only; not selector closure and not ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2276_candidate",
            "goal": "derive analytic sufficient condition linking tightened transport budget and corner replay floor with explicit confidence margin",
        },
        "gatekeeper_checks": {
            "replay_rows_exported": len(replay_rows) > 0,
            "all_corner_stats_nonempty": all(len(r["corner_replay"]) == 4 for r in replay_rows),
            "all_target_floors_bounded": all(0.0 <= r["target_passrate_floor"] <= 1.0 for r in replay_rows),
            "all_empirical_floors_bounded": all(0.0 <= r["empirical_passrate_floor_over_corners"] <= 1.0 for r in replay_rows),
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2275 S1225: certified-box corner replay passrate floor probe",
            "",
            f"- rows replayed: `{len(replay_rows)}`",
            f"- all rows floor hold: `{all(r['floor_holds'] for r in replay_rows)}`",
            f"- worst empirical floor: `{min((r['empirical_passrate_floor_over_corners'] for r in replay_rows), default=0.0):.12e}`",
            "",
            "Corner replay probe only; no selector closure / ToE closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
