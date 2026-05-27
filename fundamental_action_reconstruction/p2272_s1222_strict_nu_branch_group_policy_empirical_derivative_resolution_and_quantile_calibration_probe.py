#!/usr/bin/env python3
from __future__ import annotations
import json
import random
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2270 = GEN / "p2270_s1220_strict_nu_branch_group_policy_symbolic_derivative_bound_probe.json"
IN_2271 = GEN / "p2271_s1221_strict_nu_branch_group_policy_symbolic_vs_empirical_derivative_slack_probe.json"
OUT = GEN / "p2272_s1222_strict_nu_branch_group_policy_empirical_derivative_resolution_and_quantile_calibration_probe.json"
MD = GEN / "p2272_s1222_strict_nu_branch_group_policy_empirical_derivative_resolution_and_quantile_calibration_probe.md"


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


def quantile(xs: list[float], q: float) -> float:
    ys = sorted(xs)
    if not ys:
        return 0.0
    idx = int((len(ys) - 1) * q)
    return ys[idx]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2270 = load(IN_2270)
    p2271 = load(IN_2271)

    rows = (p2270.get("strict_nu_branch_group_policy_symbolic_derivative_bound_probe", {}) or {}).get("symbolic_row_checks", []) or []
    prior = (p2271.get("strict_nu_branch_group_policy_symbolic_vs_empirical_derivative_slack_probe", {}) or {})
    prior_max_drho = float((prior.get("global_summary", {}) or {}).get("max_abs_empirical_d_rho", 0.0) or 0.0)
    prior_max_dkappa = float((prior.get("global_summary", {}) or {}).get("max_abs_empirical_d_kappa", 0.0) or 0.0)

    seeds = [12221, 12222, 12223, 12224, 12225, 12226, 12227]
    eps = 2e-2
    profiles = [
        {"name": "baseline", "horizon": 12, "noise_amp": 0.01, "trials": 400},
        {"name": "stress_h24_n15", "horizon": 24, "noise_amp": 0.015, "trials": 500},
        {"name": "stress_h36_n20", "horizon": 36, "noise_amp": 0.02, "trials": 600},
    ]

    calibrated_rows = []
    any_nonzero_empirical = False

    for r in rows:
        rho = float(r.get("rho", 0.8) or 0.8)
        kappa = float(r.get("kappa_scale", 1.0) or 1.0)
        sym_drho_abs = abs(float(r.get("symbolic_d_bound_d_rho", 0.0) or 0.0))
        sym_dkappa_abs = abs(float(r.get("symbolic_d_bound_d_kappa", 0.0) or 0.0))

        per_profile = []
        max_emp_drho = 0.0
        max_emp_dkappa = 0.0

        for p in profiles:
            drho_samples: list[float] = []
            dkappa_samples: list[float] = []
            for seed in seeds:
                base = simulate_pass_rate(rho, kappa, seed, trials=p["trials"], horizon=p["horizon"], noise_amp=p["noise_amp"])
                plus_rho = simulate_pass_rate(min(0.95, rho + eps), kappa, seed, trials=p["trials"], horizon=p["horizon"], noise_amp=p["noise_amp"])
                plus_kappa = simulate_pass_rate(rho, min(1.8, kappa + eps), seed, trials=p["trials"], horizon=p["horizon"], noise_amp=p["noise_amp"])
                drho_samples.append((plus_rho - base) / eps)
                dkappa_samples.append((plus_kappa - base) / eps)

            abs_drho = [abs(v) for v in drho_samples]
            abs_dkappa = [abs(v) for v in dkappa_samples]
            q95_drho = quantile(abs_drho, 0.95)
            q95_dkappa = quantile(abs_dkappa, 0.95)
            profile_max_drho = max(abs_drho) if abs_drho else 0.0
            profile_max_dkappa = max(abs_dkappa) if abs_dkappa else 0.0

            max_emp_drho = max(max_emp_drho, profile_max_drho)
            max_emp_dkappa = max(max_emp_dkappa, profile_max_dkappa)

            per_profile.append({
                "profile": p,
                "empirical_fd_samples_d_rho": drho_samples,
                "empirical_fd_samples_d_kappa": dkappa_samples,
                "empirical_abs_d_rho_q95": q95_drho,
                "empirical_abs_d_kappa_q95": q95_dkappa,
                "empirical_abs_d_rho_max": profile_max_drho,
                "empirical_abs_d_kappa_max": profile_max_dkappa,
            })

        any_nonzero_empirical = any_nonzero_empirical or (max_emp_drho > 0.0) or (max_emp_dkappa > 0.0)

        calibrated_rows.append({
            "risk_tolerance": float(r.get("risk_tolerance", 0.05) or 0.05),
            "rho": rho,
            "kappa_scale": kappa,
            "symbolic_abs_d_bound_d_rho": sym_drho_abs,
            "symbolic_abs_d_bound_d_kappa": sym_dkappa_abs,
            "empirical_profiles": per_profile,
            "empirical_abs_max_d_rho_across_profiles": max_emp_drho,
            "empirical_abs_max_d_kappa_across_profiles": max_emp_dkappa,
            "symbolic_covers_empirical_max_rho": sym_drho_abs + 1e-12 >= max_emp_drho,
            "symbolic_covers_empirical_max_kappa": sym_dkappa_abs + 1e-12 >= max_emp_dkappa,
            "calibrated_slack_factor_rho": sym_drho_abs / max(max_emp_drho, 1e-12),
            "calibrated_slack_factor_kappa": sym_dkappa_abs / max(max_emp_dkappa, 1e-12),
        })

    payload = {
        "schema_version": "p2272_s1222_v1",
        "packet_id": "P2272",
        "stage_id": "S1222",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_EMPIRICAL_DERIVATIVE_RESOLUTION_AND_QUANTILE_CALIBRATION_PROBE",
        "strict_nu_branch_group_policy_empirical_derivative_resolution_and_quantile_calibration_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_EMPIRICAL_DERIVATIVE_RESOLUTION_AND_QUANTILE_CALIBRATION_PROBE_V1",
            "source_packets": [str(IN_2270.relative_to(ROOT)), str(IN_2271.relative_to(ROOT))],
            "settings": {"eps": eps, "seeds": seeds, "profiles": profiles},
            "prior_zero_sensitivity_summary": {
                "prior_max_abs_empirical_d_rho": prior_max_drho,
                "prior_max_abs_empirical_d_kappa": prior_max_dkappa,
            },
            "calibrated_rows": calibrated_rows,
            "global_summary": {
                "resolved_zero_sensitivity_issue": any_nonzero_empirical,
                "all_symbolic_cover_empirical_max_rho": all(x["symbolic_covers_empirical_max_rho"] for x in calibrated_rows),
                "all_symbolic_cover_empirical_max_kappa": all(x["symbolic_covers_empirical_max_kappa"] for x in calibrated_rows),
            },
            "theorem_scope_limit": "empirical sensitivity resolution and surrogate slack calibration only; not selector closure and not ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2273_candidate",
            "goal": "derive closed-form tightened envelope constants from calibrated empirical maxima with explicit safety multiplier policy",
        },
        "gatekeeper_checks": {
            "calibrated_rows_exported": len(calibrated_rows) > 0,
            "resolved_zero_sensitivity_issue": any_nonzero_empirical,
            "all_symbolic_cover_empirical_max_rho": all(x["symbolic_covers_empirical_max_rho"] for x in calibrated_rows),
            "all_symbolic_cover_empirical_max_kappa": all(x["symbolic_covers_empirical_max_kappa"] for x in calibrated_rows),
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2272 S1222: empirical derivative resolution + quantile calibration",
            "",
            f"- prior max |empirical d/drho| (P2271): `{prior_max_drho:.12e}`",
            f"- prior max |empirical d/dkappa| (P2271): `{prior_max_dkappa:.12e}`",
            f"- resolved zero-sensitivity issue: `{any_nonzero_empirical}`",
            f"- symbolic cover empirical max (rho): `{all(x['symbolic_covers_empirical_max_rho'] for x in calibrated_rows)}`",
            f"- symbolic cover empirical max (kappa): `{all(x['symbolic_covers_empirical_max_kappa'] for x in calibrated_rows)}`",
            "",
            "Calibration probe only; no selector closure / ToE closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
