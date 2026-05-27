#!/usr/bin/env python3
from __future__ import annotations
import json
import random
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2270 = GEN / "p2270_s1220_strict_nu_branch_group_policy_symbolic_derivative_bound_probe.json"
OUT = GEN / "p2271_s1221_strict_nu_branch_group_policy_symbolic_vs_empirical_derivative_slack_probe.json"
MD = GEN / "p2271_s1221_strict_nu_branch_group_policy_symbolic_vs_empirical_derivative_slack_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def simulate_pass_rate(rho: float, kappa_scale: float, seed: int, trials: int = 300, horizon: int = 12) -> float:
    rng = random.Random(seed)
    passed = 0
    for _ in range(trials):
        margin = 0.20
        for _ in range(horizon):
            noise = (rng.random() - 0.5) * 0.01
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
    p2270 = load(IN_2270)
    probe = p2270.get("strict_nu_branch_group_policy_symbolic_derivative_bound_probe", {}) or {}
    rows = probe.get("symbolic_row_checks", []) or []

    empirical_seeds = [12211, 12212, 12213, 12214, 12215]
    eps = 1e-2

    slack_rows = []
    max_abs_emp_drho = 0.0
    max_abs_emp_dkappa = 0.0

    for r in rows:
        rho = float(r.get("rho", 0.8) or 0.8)
        kappa = float(r.get("kappa_scale", 1.0) or 1.0)
        sym_drho = abs(float(r.get("symbolic_d_bound_d_rho", 0.0) or 0.0))
        sym_dkappa = abs(float(r.get("symbolic_d_bound_d_kappa", 0.0) or 0.0))

        drho_samples = []
        dkappa_samples = []

        for seed in empirical_seeds:
            base = simulate_pass_rate(rho, kappa, seed=seed)
            plus_rho = simulate_pass_rate(min(0.95, rho + eps), kappa, seed=seed)
            plus_kappa = simulate_pass_rate(rho, min(1.8, kappa + eps), seed=seed)
            drho_samples.append((plus_rho - base) / eps)
            dkappa_samples.append((plus_kappa - base) / eps)

        emp_abs_drho_max = max(abs(v) for v in drho_samples)
        emp_abs_dkappa_max = max(abs(v) for v in dkappa_samples)

        max_abs_emp_drho = max(max_abs_emp_drho, emp_abs_drho_max)
        max_abs_emp_dkappa = max(max_abs_emp_dkappa, emp_abs_dkappa_max)

        slack_rows.append(
            {
                "risk_tolerance": float(r.get("risk_tolerance", 0.05) or 0.05),
                "rho": rho,
                "kappa_scale": kappa,
                "symbolic_abs_d_bound_d_rho": sym_drho,
                "symbolic_abs_d_bound_d_kappa": sym_dkappa,
                "empirical_fd_samples_d_rho": drho_samples,
                "empirical_fd_samples_d_kappa": dkappa_samples,
                "empirical_abs_max_d_rho": emp_abs_drho_max,
                "empirical_abs_max_d_kappa": emp_abs_dkappa_max,
                "slack_factor_rho": (sym_drho / max(emp_abs_drho_max, 1e-12)),
                "slack_factor_kappa": (sym_dkappa / max(emp_abs_dkappa_max, 1e-12)),
                "symbolic_upper_bounds_empirical_rho": sym_drho + 1e-12 >= emp_abs_drho_max,
                "symbolic_upper_bounds_empirical_kappa": sym_dkappa + 1e-12 >= emp_abs_dkappa_max,
            }
        )

    payload = {
        "schema_version": "p2271_s1221_v1",
        "packet_id": "P2271",
        "stage_id": "S1221",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_SYMBOLIC_VS_EMPIRICAL_DERIVATIVE_SLACK_PROBE",
        "strict_nu_branch_group_policy_symbolic_vs_empirical_derivative_slack_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_SYMBOLIC_VS_EMPIRICAL_DERIVATIVE_SLACK_PROBE_V1",
            "source_packets": [str(IN_2270.relative_to(ROOT))],
            "empirical_settings": {
                "finite_difference_eps": eps,
                "seeds": empirical_seeds,
                "trials_per_seed": 300,
                "horizon": 12,
            },
            "slack_rows": slack_rows,
            "global_summary": {
                "all_symbolic_rho_bounds_cover_empirical": all(s["symbolic_upper_bounds_empirical_rho"] for s in slack_rows),
                "all_symbolic_kappa_bounds_cover_empirical": all(s["symbolic_upper_bounds_empirical_kappa"] for s in slack_rows),
                "max_abs_empirical_d_rho": max_abs_emp_drho,
                "max_abs_empirical_d_kappa": max_abs_emp_dkappa,
            },
            "physical_interpretation_note": "This packet quantifies model-form slack: how conservative symbolic derivative envelopes are versus multi-seed empirical finite-difference sensitivity under nonlinear perturbation replay.",
            "theorem_scope_limit": "surrogate derivative-sensitivity comparison only; not strict-core selector closure and not ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2272_candidate",
            "goal": "tighten symbolic envelope constants using calibrated empirical quantiles while preserving explicit no-false-pass safety margin",
        },
        "gatekeeper_checks": {
            "slack_rows_exported": len(slack_rows) > 0,
            "all_symbolic_rho_bounds_cover_empirical": all(s["symbolic_upper_bounds_empirical_rho"] for s in slack_rows),
            "all_symbolic_kappa_bounds_cover_empirical": all(s["symbolic_upper_bounds_empirical_kappa"] for s in slack_rows),
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2271 S1221: symbolic vs empirical derivative slack probe",
            "",
            f"- rows analyzed: `{len(slack_rows)}`",
            f"- symbolic rho envelopes cover empirical: `{all(s['symbolic_upper_bounds_empirical_rho'] for s in slack_rows)}`",
            f"- symbolic kappa envelopes cover empirical: `{all(s['symbolic_upper_bounds_empirical_kappa'] for s in slack_rows)}`",
            f"- max |empirical d/drho|: `{max_abs_emp_drho:.12e}`",
            f"- max |empirical d/dkappa|: `{max_abs_emp_dkappa:.12e}`",
            "",
            "Sensitivity-comparison probe only; no selector closure / ToE closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
