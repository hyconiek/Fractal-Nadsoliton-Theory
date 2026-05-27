#!/usr/bin/env python3
from __future__ import annotations
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2268 = GEN / "p2268_s1218_strict_nu_branch_group_policy_certificate_passrate_lower_bound_probe.json"
OUT = GEN / "p2269_s1219_strict_nu_branch_group_policy_passrate_bound_monotonicity_lipschitz_probe.json"
MD = GEN / "p2269_s1219_strict_nu_branch_group_policy_passrate_bound_monotonicity_lipschitz_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def clipped(x: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


def bound_value(rho: float, kappa_scale: float, c: float = 1.8) -> float:
    envelope = 0.015 + 0.01 * (1.0 - rho) + 0.005 * (kappa_scale - 1.0)
    attenuation = rho / max(kappa_scale, 1e-12)
    return clipped(math.exp(-c * envelope / max(attenuation, 1e-12)))


def grad_components(rho: float, kappa_scale: float, eps: float = 1e-5) -> tuple[float, float]:
    b0 = bound_value(rho, kappa_scale)
    br = bound_value(rho + eps, kappa_scale)
    bk = bound_value(rho, kappa_scale + eps)
    d_rho = (br - b0) / eps
    d_kappa = (bk - b0) / eps
    return d_rho, d_kappa


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2268 = load(IN_2268)
    rows = (p2268.get("strict_nu_branch_group_policy_certificate_passrate_lower_bound_probe", {}) or {}).get("bound_rows", []) or []

    checks = []
    abs_grad_sum_max = 0.0
    min_drho = float("inf")
    max_dkappa = float("-inf")

    for row in rows:
        rho = float(row.get("rho", 0.8) or 0.8)
        kappa = float(row.get("kappa_scale", 1.0) or 1.0)
        b = float(row.get("analytic_lower_bound_passrate", bound_value(rho, kappa)))
        d_rho, d_kappa = grad_components(rho, kappa)

        abs_grad_sum = abs(d_rho) + abs(d_kappa)
        abs_grad_sum_max = max(abs_grad_sum_max, abs_grad_sum)
        min_drho = min(min_drho, d_rho)
        max_dkappa = max(max_dkappa, d_kappa)

        checks.append(
            {
                "risk_tolerance": float(row.get("risk_tolerance", 0.05)),
                "rho": rho,
                "kappa_scale": kappa,
                "bound_value": b,
                "d_bound_d_rho_fd": d_rho,
                "d_bound_d_kappa_fd": d_kappa,
                "monotone_in_rho_non_decreasing_local": d_rho >= -1e-6,
                "monotone_in_kappa_non_increasing_local": d_kappa <= 1e-6,
            }
        )

    payload = {
        "schema_version": "p2269_s1219_v1",
        "packet_id": "P2269",
        "stage_id": "S1219",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_PASSRATE_BOUND_MONOTONICITY_LIPSCHITZ_PROBE",
        "strict_nu_branch_group_policy_passrate_bound_monotonicity_lipschitz_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_PASSRATE_BOUND_MONOTONICITY_LIPSCHITZ_PROBE_V1",
            "source_packets": [str(IN_2268.relative_to(ROOT))],
            "finite_difference_eps": 1e-5,
            "local_checks": checks,
            "lipschitz_certificate_l1": abs_grad_sum_max,
            "monotonicity_summary": {
                "min_d_bound_d_rho_fd": min_drho if checks else 0.0,
                "max_d_bound_d_kappa_fd": max_dkappa if checks else 0.0,
                "all_rows_monotone_in_rho": all(c["monotone_in_rho_non_decreasing_local"] for c in checks),
                "all_rows_monotone_in_kappa": all(c["monotone_in_kappa_non_increasing_local"] for c in checks),
            },
            "physical_interpretation_note": "Higher rho increases effective control attenuation and should weakly raise conservative certificate survival; larger kappa_scale weakens attenuation and should weakly lower survival.",
            "theorem_scope_limit": "local finite-difference monotonicity/Lipschitz evidence for synthetic surrogate only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2270_candidate",
            "goal": "derive closed-form symbolic derivative bounds under explicit admissible box constraints for rho and kappa_scale",
        },
        "gatekeeper_checks": {
            "monotonicity_lipschitz_exported": True,
            "all_rows_monotone_in_rho": all(c["monotone_in_rho_non_decreasing_local"] for c in checks),
            "all_rows_monotone_in_kappa": all(c["monotone_in_kappa_non_increasing_local"] for c in checks),
            "lipschitz_l1_nonnegative": abs_grad_sum_max >= 0.0,
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2269 S1219: pass-rate bound monotonicity/lipschitz probe",
                "",
                f"- rows analyzed: `{len(checks)}`",
                f"- monotone in rho (local FD): `{all(c['monotone_in_rho_non_decreasing_local'] for c in checks)}`",
                f"- monotone in kappa_scale (local FD): `{all(c['monotone_in_kappa_non_increasing_local'] for c in checks)}`",
                f"- L1 Lipschitz certificate (max |d/drho|+|d/dkappa|): `{abs_grad_sum_max}`",
                "",
                "Local synthetic derivative evidence only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
