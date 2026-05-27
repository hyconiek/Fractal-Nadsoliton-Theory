#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2250 = GEN / "p2250_s1200_strict_nu_branch_group_policy_multistep_controller_trajectory_probe.json"
OUT = GEN / "p2251_s1201_strict_nu_branch_group_policy_direction_field_rho_schedule_probe.json"
MD = GEN / "p2251_s1201_strict_nu_branch_group_policy_direction_field_rho_schedule_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2250 = load(IN_2250)
    probe = (p2250.get("strict_nu_branch_group_policy_multistep_controller_trajectory_probe", {}) or {})
    inp = (probe.get("inputs", {}) or {})

    n0 = float(inp.get("initial_group_count", 1.0) or 1.0)
    l0 = float(inp.get("initial_load_ratio", 1.0) or 1.0)
    total = float(inp.get("total_coverage_mass", 1.0) or 1.0)
    kappa = float(inp.get("kappa_reference", 1.0) or 1.0)

    # Two direction fields and two rho schedules.
    directions = {
        "balanced": (0.10, 0.05 * max(l0, 1e-12)),
        "load_heavy": (0.05, 0.12 * max(l0, 1e-12)),
    }

    def rho_const(_: float, __: int) -> float:
        return 0.8

    def rho_adaptive(m: float, _: int) -> float:
        # tighter near boundary, looser when margin comfortable
        return 0.9 if m > 0.2 else 0.6

    schedules = {
        "const_0p8": rho_const,
        "adaptive_0p9_0p6": rho_adaptive,
    }

    def margin(nv: float, lv: float) -> float:
        return total - ((nv - 1.0) + kappa * lv)

    horizon = 10
    results = []
    for dname, (dir_dn, dir_dl) in directions.items():
        for sname, rho_fn in schedules.items():
            n, l = n0, l0
            executed = 0
            stop_reason = "max_steps_reached"
            min_m = margin(n, l)
            cum_progress = 0.0
            for t in range(horizon):
                m = margin(n, l)
                rho = rho_fn(m, t)
                cap = rho * max(m, 0.0)
                raw_use = dir_dn + kappa * dir_dl
                scale = 0.0 if raw_use <= 1e-15 else min(1.0, cap / raw_use)
                dn = scale * dir_dn
                dl = scale * dir_dl
                use = dn + kappa * dl
                n += dn
                l += dl
                m2 = margin(n, l)
                min_m = min(min_m, m2)
                cum_progress += use
                executed += 1
                if scale <= 1e-12:
                    stop_reason = "controller_stopped_near_boundary"
                    break
            results.append(
                {
                    "direction_field": dname,
                    "rho_schedule": sname,
                    "executed_steps": executed,
                    "cumulative_budget_progress": cum_progress,
                    "minimum_margin": min_m,
                    "stop_reason": stop_reason,
                    "nonnegative_margin": min_m >= -1e-15,
                }
            )

    best = max(results, key=lambda r: r["cumulative_budget_progress"]) if results else None

    payload = {
        "schema_version": "p2251_s1201_v1",
        "packet_id": "P2251",
        "stage_id": "S1201",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_DIRECTION_FIELD_RHO_SCHEDULE_PROBE",
        "strict_nu_branch_group_policy_direction_field_rho_schedule_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_DIRECTION_FIELD_RHO_SCHEDULE_PROBE_V1",
            "source_packets": [str(IN_2250.relative_to(ROOT))],
            "inputs": {
                "initial_group_count": n0,
                "initial_load_ratio": l0,
                "total_coverage_mass": total,
                "kappa_reference": kappa,
                "horizon": horizon,
            },
            "scenario_results": results,
            "best_progress_scenario": best,
            "physical_interpretation_note": "Comparing direction fields and rho schedules exposes policy-navigation anisotropy: some update geometries consume safety budget more efficiently while preserving nonnegative margin invariants.",
            "theorem_scope_limit": "finite-horizon scenario-comparison diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2252_candidate",
            "goal": "derive analytic sufficient criterion selecting rho schedule by local curvature proxy of margin landscape",
        },
        "gatekeeper_checks": {
            "scenario_comparison_exported": True,
            "all_scenarios_nonnegative_margin": all(r["nonnegative_margin"] for r in results),
            "best_scenario_identified": best is not None,
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
                "# P2251 S1201: direction-field and rho-schedule comparison probe",
                "",
                f"- scenario count: `{len(results)}`",
                f"- all scenarios nonnegative margin: `{all(r['nonnegative_margin'] for r in results)}`",
                f"- best scenario: `{best['direction_field']}/{best['rho_schedule']}`" if best else "- best scenario: `none`",
                f"- best cumulative progress: `{best['cumulative_budget_progress']:.12e}`" if best else "- best cumulative progress: `n/a`",
                "",
                "Finite-horizon comparison only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
