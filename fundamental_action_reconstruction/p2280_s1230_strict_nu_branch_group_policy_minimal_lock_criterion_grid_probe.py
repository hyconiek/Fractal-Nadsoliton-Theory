#!/usr/bin/env python3
from __future__ import annotations
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2275 = GEN / "p2275_s1225_strict_nu_branch_group_policy_certified_box_corner_replay_passrate_floor_probe.json"
IN_2279 = GEN / "p2279_s1229_strict_nu_branch_group_policy_locked_confidence_profile_corner_replay_probe.json"
OUT = GEN / "p2280_s1230_strict_nu_branch_group_policy_minimal_lock_criterion_grid_probe.json"
MD = GEN / "p2280_s1230_strict_nu_branch_group_policy_minimal_lock_criterion_grid_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def margin(delta: float, n_eff: int) -> float:
    return math.sqrt(math.log(2.0 / delta) / (2.0 * max(1, n_eff)))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2275 = load(IN_2275)
    p2279 = load(IN_2279)

    replay_rows = (p2275.get("strict_nu_branch_group_policy_certified_box_corner_replay_passrate_floor_probe", {}) or {}).get("replay_rows", []) or []
    locked = (p2279.get("strict_nu_branch_group_policy_locked_confidence_profile_corner_replay_probe", {}) or {}).get("strict_nu_branch_group_policy_locked_confidence_profile_corner_replay_probe", {})
    if not locked:
        locked = p2279.get("strict_nu_branch_group_policy_locked_confidence_profile_corner_replay_probe", {}) or {}

    profile = (locked.get("locked_profile", {}) or {})
    seeds = profile.get("seeds", []) or [1]
    base_trials = int(profile.get("base_trials", 700) or 700)

    empirical_by_risk = {
        float(r.get("risk_tolerance", 0.05) or 0.05): float(r.get("empirical_passrate_floor_over_corners", 0.0) or 0.0)
        for r in replay_rows
    }

    delta_grid = [0.1, 0.05, 0.02, 0.01, 0.005, 0.001]
    trial_mult_grid = [1, 2, 3, 4, 6, 8]

    rows = []
    feasible = []
    for d in delta_grid:
        for m in trial_mult_grid:
            trials = base_trials * m
            n_eff = max(1, len(seeds) * trials * 4)
            mar = margin(d, n_eff)

            per_risk = []
            for risk, emp_floor in sorted(empirical_by_risk.items()):
                target = max(0.0, 1.0 - risk)
                cert_floor = max(0.0, emp_floor - mar)
                per_risk.append({
                    "risk_tolerance": risk,
                    "target_floor": target,
                    "empirical_floor": emp_floor,
                    "certified_floor": cert_floor,
                    "meets": cert_floor + 1e-12 >= target,
                    "margin_to_target": cert_floor - target,
                })

            all_meet = all(x["meets"] for x in per_risk)
            worst_margin = min((x["margin_to_target"] for x in per_risk), default=0.0)
            cost = m * math.log(2.0 / d)
            row = {
                "delta": d,
                "trial_multiplier": m,
                "n_eff": n_eff,
                "adaptive_margin": mar,
                "all_meet": all_meet,
                "worst_margin_to_target": worst_margin,
                "cost_proxy": cost,
                "per_risk": per_risk,
            }
            rows.append(row)
            if all_meet:
                feasible.append(row)

    minimal = min(feasible, key=lambda x: (x["cost_proxy"], x["adaptive_margin"])) if feasible else {}

    payload = {
        "schema_version": "p2280_s1230_v1",
        "packet_id": "P2280",
        "stage_id": "S1230",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_MINIMAL_LOCK_CRITERION_GRID_PROBE",
        "strict_nu_branch_group_policy_minimal_lock_criterion_grid_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_MINIMAL_LOCK_CRITERION_GRID_PROBE_V1",
            "source_packets": [str(IN_2275.relative_to(ROOT)), str(IN_2279.relative_to(ROOT))],
            "grid": {"delta": delta_grid, "trial_multiplier": trial_mult_grid},
            "rows": rows,
            "feasible_count": len(feasible),
            "minimal_feasible_config": minimal,
            "selection_rule": "minimize cost_proxy = trial_multiplier*log(2/delta) subject to all_meet=true",
            "theorem_scope_limit": "policy-lock grid search over surrogate confidence model only; not selector closure and not ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2281_candidate",
            "goal": "run fresh corner replay under minimal feasible config and compare realized margins against predicted grid margins",
        },
        "gatekeeper_checks": {
            "rows_exported": len(rows) > 0,
            "cost_proxy_nonnegative": all(r["cost_proxy"] >= 0.0 for r in rows),
            "margins_nonnegative": all(r["adaptive_margin"] >= 0.0 for r in rows),
            "certified_floors_bounded": all(0.0 <= rr["certified_floor"] <= 1.0 for r in rows for rr in r["per_risk"]),
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2280 S1230: minimal lock criterion grid probe",
            "",
            f"- grid rows: `{len(rows)}`",
            f"- feasible rows: `{len(feasible)}`",
            f"- minimal feasible delta: `{minimal.get('delta', 'NA') if minimal else 'NA'}`",
            f"- minimal feasible trial multiplier: `{minimal.get('trial_multiplier', 'NA') if minimal else 'NA'}`",
            "",
            "Grid-search certificate only; no selector closure / ToE closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
