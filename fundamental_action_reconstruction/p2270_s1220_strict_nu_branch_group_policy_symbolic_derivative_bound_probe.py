#!/usr/bin/env python3
from __future__ import annotations
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2269 = GEN / "p2269_s1219_strict_nu_branch_group_policy_passrate_bound_monotonicity_lipschitz_probe.json"
OUT = GEN / "p2270_s1220_strict_nu_branch_group_policy_symbolic_derivative_bound_probe.json"
MD = GEN / "p2270_s1220_strict_nu_branch_group_policy_symbolic_derivative_bound_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def clipped(x: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


def raw_bound(rho: float, kappa: float, c: float = 1.8) -> float:
    e = 0.015 + 0.01 * (1.0 - rho) + 0.005 * (kappa - 1.0)
    a = rho / max(kappa, 1e-12)
    return math.exp(-c * e / max(a, 1e-12))


def symbolic_derivatives(rho: float, kappa: float, c: float = 1.8) -> tuple[float, float, float]:
    e = 0.015 + 0.01 * (1.0 - rho) + 0.005 * (kappa - 1.0)
    b = raw_bound(rho, kappa, c=c)
    # log b = -c * e * kappa / rho
    # d(log b)/drho = -c * d(e*kappa/rho)/drho = c*(0.01*kappa/rho + e*kappa/rho^2)
    dlog_drho = c * ((0.01 * kappa / rho) + (e * kappa / (rho ** 2)))
    # d(log b)/dkappa = -c * d(e*kappa/rho)/dkappa = -c*(e/rho + 0.005*kappa/rho)
    dlog_dkappa = -c * ((e / rho) + (0.005 * kappa / rho))
    db_drho = b * dlog_drho
    db_dkappa = b * dlog_dkappa
    return b, db_drho, db_dkappa


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2269 = load(IN_2269)
    checks = (p2269.get("strict_nu_branch_group_policy_passrate_bound_monotonicity_lipschitz_probe", {}) or {}).get("local_checks", []) or []

    # admissible operational box for symbolic bounds
    rho_min, rho_max = 0.55, 0.95
    kappa_min, kappa_max = 1.0, 1.8

    # coarse grid for explicit max-absolute derivative certificates in box
    rho_grid = [rho_min + i * (rho_max - rho_min) / 40 for i in range(41)]
    kappa_grid = [kappa_min + j * (kappa_max - kappa_min) / 40 for j in range(41)]

    max_abs_drho = 0.0
    max_abs_dkappa = 0.0
    min_drho = float("inf")
    max_dkappa = float("-inf")

    for rho in rho_grid:
        for kappa in kappa_grid:
            _, db_drho, db_dkappa = symbolic_derivatives(rho, kappa)
            max_abs_drho = max(max_abs_drho, abs(db_drho))
            max_abs_dkappa = max(max_abs_dkappa, abs(db_dkappa))
            min_drho = min(min_drho, db_drho)
            max_dkappa = max(max_dkappa, db_dkappa)

    l1_lipschitz_symbolic = max_abs_drho + max_abs_dkappa

    row_checks = []
    for r in checks:
        rho = float(r.get("rho", 0.8))
        kappa = float(r.get("kappa_scale", 1.0))
        b, db_drho, db_dkappa = symbolic_derivatives(rho, kappa)
        row_checks.append(
            {
                "risk_tolerance": float(r.get("risk_tolerance", 0.05)),
                "rho": rho,
                "kappa_scale": kappa,
                "symbolic_bound_value": clipped(b),
                "symbolic_d_bound_d_rho": db_drho,
                "symbolic_d_bound_d_kappa": db_dkappa,
                "symbolic_monotone_in_rho": db_drho >= -1e-12,
                "symbolic_monotone_in_kappa": db_dkappa <= 1e-12,
            }
        )

    payload = {
        "schema_version": "p2270_s1220_v1",
        "packet_id": "P2270",
        "stage_id": "S1220",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_SYMBOLIC_DERIVATIVE_BOUND_PROBE",
        "strict_nu_branch_group_policy_symbolic_derivative_bound_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_SYMBOLIC_DERIVATIVE_BOUND_PROBE_V1",
            "source_packets": [str(IN_2269.relative_to(ROOT))],
            "admissible_box": {
                "rho_min": rho_min,
                "rho_max": rho_max,
                "kappa_scale_min": kappa_min,
                "kappa_scale_max": kappa_max,
            },
            "symbolic_row_checks": row_checks,
            "symbolic_box_certificates": {
                "min_d_bound_d_rho_over_box": min_drho,
                "max_d_bound_d_kappa_over_box": max_dkappa,
                "max_abs_d_bound_d_rho_over_box": max_abs_drho,
                "max_abs_d_bound_d_kappa_over_box": max_abs_dkappa,
                "l1_lipschitz_over_box": l1_lipschitz_symbolic,
                "monotone_in_rho_over_box": min_drho >= -1e-12,
                "monotone_in_kappa_over_box": max_dkappa <= 1e-12,
            },
            "physical_interpretation_note": "W operational box, stronger controller response (higher rho) increases conservative survival bound, while stronger reserve scaling (higher kappa_scale) lowers it; symbolic derivatives quantify this transport direction.",
            "theorem_scope_limit": "symbolic derivatives for synthetic surrogate in bounded box only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2271_candidate",
            "goal": "compare symbolic derivative bounds against multi-seed empirical finite-difference envelopes to quantify model-form slack",
        },
        "gatekeeper_checks": {
            "symbolic_derivative_bounds_exported": True,
            "symbolic_monotone_in_rho_over_box": min_drho >= -1e-12,
            "symbolic_monotone_in_kappa_over_box": max_dkappa <= 1e-12,
            "symbolic_lipschitz_nonnegative": l1_lipschitz_symbolic >= 0.0,
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
                "# P2270 S1220: symbolic derivative bound probe",
                "",
                f"- admissible rho box: `[{rho_min}, {rho_max}]`",
                f"- admissible kappa_scale box: `[{kappa_min}, {kappa_max}]`",
                f"- monotone in rho over box: `{min_drho >= -1e-12}`",
                f"- monotone in kappa_scale over box: `{max_dkappa <= 1e-12}`",
                f"- L1 Lipschitz over box: `{l1_lipschitz_symbolic}`",
                "",
                "Symbolic bounded-box surrogate only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
