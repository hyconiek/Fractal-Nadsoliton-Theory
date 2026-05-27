#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2252 = GEN / "p2252_s1202_strict_nu_branch_group_policy_adaptive_rho_curvature_proxy_probe.json"
OUT = GEN / "p2253_s1203_strict_nu_branch_group_policy_adaptive_vs_fixed_rho_trajectory_probe.json"
MD = GEN / "p2253_s1203_strict_nu_branch_group_policy_adaptive_vs_fixed_rho_trajectory_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2252 = load(IN_2252)
    probe = (p2252.get("strict_nu_branch_group_policy_adaptive_rho_curvature_proxy_probe", {}) or {})
    inp = (probe.get("inputs", {}) or {})
    cmodel = (probe.get("curvature_proxy_model", {}) or {})

    n0 = float(inp.get("initial_group_count", 1.0) or 1.0)
    l0 = float(inp.get("initial_load_ratio", 1.0) or 1.0)
    total = float(inp.get("total_coverage_mass", 1.0) or 1.0)
    kappa = float(inp.get("kappa_reference", 1.0) or 1.0)

    rho_fixed = float(cmodel.get("rho_baseline", 0.8) or 0.8)
    rho_min = float(cmodel.get("rho_min", 0.55) or 0.55)
    rho_max = float(cmodel.get("rho_max", 0.9) or 0.9)

    dir_dn = 0.10
    dir_dl = 0.05 * max(l0, 1e-12)
    horizon = 10

    def margin(nv: float, lv: float) -> float:
        return total - ((nv - 1.0) + kappa * lv)

    def rho_adaptive(m: float) -> float:
        norm = max(abs(total), 1.0)
        proxy = max(0.0, min(1.0, 1.0 - m / norm))
        return rho_max - (rho_max - rho_min) * proxy

    def run(mode: str) -> dict[str, Any]:
        n, l = n0, l0
        rows = []
        stop_reason = "max_steps_reached"
        for t in range(horizon):
            m = margin(n, l)
            rho = rho_fixed if mode == "fixed" else rho_adaptive(m)
            cap = rho * max(m, 0.0)
            raw = dir_dn + kappa * dir_dl
            scale = 0.0 if raw <= 1e-15 else min(1.0, cap / raw)
            dn = scale * dir_dn
            dl = scale * dir_dl
            use = dn + kappa * dl
            n += dn
            l += dl
            m2 = margin(n, l)
            rows.append({"step": t, "rho": rho, "scale": scale, "budget_use": use, "margin_after": m2})
            if scale <= 1e-12:
                stop_reason = "controller_stopped_near_boundary"
                break
        min_margin = min((r["margin_after"] for r in rows), default=margin(n, l))
        cum_use = sum(r["budget_use"] for r in rows)
        return {
            "mode": mode,
            "rows": rows,
            "executed_steps": len(rows),
            "stop_reason": stop_reason,
            "minimum_margin": min_margin,
            "cumulative_budget_use": cum_use,
            "nonnegative_margin": min_margin >= -1e-15,
        }

    fixed = run("fixed")
    adaptive = run("adaptive")

    payload = {
        "schema_version": "p2253_s1203_v1",
        "packet_id": "P2253",
        "stage_id": "S1203",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_ADAPTIVE_VS_FIXED_RHO_TRAJECTORY_PROBE",
        "strict_nu_branch_group_policy_adaptive_vs_fixed_rho_trajectory_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_ADAPTIVE_VS_FIXED_RHO_TRAJECTORY_PROBE_V1",
            "source_packets": [str(IN_2252.relative_to(ROOT))],
            "inputs": {
                "initial_group_count": n0,
                "initial_load_ratio": l0,
                "total_coverage_mass": total,
                "kappa_reference": kappa,
                "horizon": horizon,
                "direction": {"delta_group_count": dir_dn, "delta_load_ratio": dir_dl},
            },
            "fixed_rho_result": fixed,
            "adaptive_rho_result": adaptive,
            "comparison": {
                "adaptive_minus_fixed_cumulative_progress": adaptive["cumulative_budget_use"] - fixed["cumulative_budget_use"],
                "adaptive_minus_fixed_min_margin": adaptive["minimum_margin"] - fixed["minimum_margin"],
            },
            "physical_interpretation_note": "Direct trajectory comparison quantifies trade-off between stability reserve and policy progress: adaptive rho modulates local aggressiveness near boundary while preserving nonnegative margin invariants.",
            "theorem_scope_limit": "finite-horizon comparative controller diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2254_candidate",
            "goal": "derive sufficient dominance condition when adaptive-rho outperforms fixed-rho in cumulative progress under minimum-margin floor constraint",
        },
        "gatekeeper_checks": {
            "adaptive_vs_fixed_exported": True,
            "fixed_nonnegative_margin": fixed["nonnegative_margin"],
            "adaptive_nonnegative_margin": adaptive["nonnegative_margin"],
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
                "# P2253 S1203: adaptive-vs-fixed rho trajectory probe",
                "",
                f"- fixed cumulative progress: `{fixed['cumulative_budget_use']:.12e}`",
                f"- adaptive cumulative progress: `{adaptive['cumulative_budget_use']:.12e}`",
                f"- delta progress (adaptive-fixed): `{(adaptive['cumulative_budget_use']-fixed['cumulative_budget_use']):.12e}`",
                f"- fixed minimum margin: `{fixed['minimum_margin']:.12e}`",
                f"- adaptive minimum margin: `{adaptive['minimum_margin']:.12e}`",
                "",
                "Finite-horizon comparative diagnostic only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
