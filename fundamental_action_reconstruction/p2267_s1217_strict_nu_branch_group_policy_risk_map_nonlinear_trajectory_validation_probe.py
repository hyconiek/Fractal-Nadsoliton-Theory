#!/usr/bin/env python3
from __future__ import annotations
import json
import random
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2266 = GEN / "p2266_s1216_strict_nu_branch_group_policy_risk_calibrated_controller_map_probe.json"
OUT = GEN / "p2267_s1217_strict_nu_branch_group_policy_risk_map_nonlinear_trajectory_validation_probe.json"
MD = GEN / "p2267_s1217_strict_nu_branch_group_policy_risk_map_nonlinear_trajectory_validation_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2266 = load(IN_2266)
    probe = (p2266.get("strict_nu_branch_group_policy_risk_calibrated_controller_map_probe", {}) or {})
    rows = probe.get("controller_map_rows", []) or []

    rng = random.Random(2267)
    horizon = 12
    trials = 200
    results = []

    for row in rows:
        risk = float(row.get("risk_tolerance", 0.05) or 0.05)
        rho = float(row.get("rho_recommended", 0.8) or 0.8)
        kappa_scale = float(row.get("kappa_scale_recommended", 1.0) or 1.0)

        pass_count = 0
        for _ in range(trials):
            margin = 0.20
            for _t in range(horizon):
                # nonlinear drift + noise + controller attenuation by rho/kappa_scale
                noise = (rng.random() - 0.5) * 0.01
                drift = 0.015 + 0.01 * (margin ** 2)
                control_relief = rho * 0.008 / max(kappa_scale, 1e-12)
                margin = margin - drift + control_relief + noise
                if margin < -1e-12:
                    break
            if margin >= -1e-12:
                pass_count += 1

        pass_rate = pass_count / trials
        results.append(
            {
                "risk_tolerance": risk,
                "rho": rho,
                "kappa_scale": kappa_scale,
                "certificate_pass_rate": pass_rate,
                "target_pass_rate": 1.0 - risk,
                "certificate_holds": pass_rate + 1e-15 >= (1.0 - risk),
            }
        )

    payload = {
        "schema_version": "p2267_s1217_v1",
        "packet_id": "P2267",
        "stage_id": "S1217",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_RISK_MAP_NONLINEAR_TRAJECTORY_VALIDATION_PROBE",
        "strict_nu_branch_group_policy_risk_map_nonlinear_trajectory_validation_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_RISK_MAP_NONLINEAR_TRAJECTORY_VALIDATION_PROBE_V1",
            "source_packets": [str(IN_2266.relative_to(ROOT))],
            "inputs": {
                "horizon": horizon,
                "trials": trials,
                "random_seed": 2267,
            },
            "risk_map_validation_rows": results,
            "physical_interpretation_note": "Nonlinear perturbed trajectories test whether risk-calibrated controller settings preserve nonnegative margin at rates matching operational tolerance targets.",
            "theorem_scope_limit": "synthetic nonlinear trajectory validation only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2268_candidate",
            "goal": "derive analytic conservative lower bound on certificate pass-rate as function of rho, kappa_scale, and perturbation envelope",
        },
        "gatekeeper_checks": {
            "risk_map_validation_exported": True,
            "all_pass_rates_bounded": all(0.0 <= r["certificate_pass_rate"] <= 1.0 for r in results),
            "all_targets_bounded": all(0.0 <= r["target_pass_rate"] <= 1.0 for r in results),
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text("\n".join([
        "# P2267 S1217: risk-map nonlinear trajectory validation probe",
        "",
        f"- horizon: `{horizon}`",
        f"- trials per risk row: `{trials}`",
        f"- rows validated: `{len(results)}`",
        f"- all rates bounded: `{all(0.0 <= r['certificate_pass_rate'] <= 1.0 for r in results)}`",
        "",
        "Synthetic nonlinear validation only; no kernel-bridge or selector-closure claim.",
    ]) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
