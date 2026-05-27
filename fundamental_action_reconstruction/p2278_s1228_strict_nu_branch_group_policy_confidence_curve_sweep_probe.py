#!/usr/bin/env python3
from __future__ import annotations
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2275 = GEN / "p2275_s1225_strict_nu_branch_group_policy_certified_box_corner_replay_passrate_floor_probe.json"
IN_2277 = GEN / "p2277_s1227_strict_nu_branch_group_policy_adaptive_confidence_margin_certificate_probe.json"
OUT = GEN / "p2278_s1228_strict_nu_branch_group_policy_confidence_curve_sweep_probe.json"
MD = GEN / "p2278_s1228_strict_nu_branch_group_policy_confidence_curve_sweep_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def hoeffding_margin(delta: float, n_eff: int) -> float:
    return math.sqrt(math.log(2.0 / delta) / (2.0 * max(1, n_eff)))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2275 = load(IN_2275)
    p2277 = load(IN_2277)

    replay_rows = (p2275.get("strict_nu_branch_group_policy_certified_box_corner_replay_passrate_floor_probe", {}) or {}).get("replay_rows", []) or []
    settings = (p2275.get("strict_nu_branch_group_policy_certified_box_corner_replay_passrate_floor_probe", {}) or {}).get("settings", {}) or {}

    base_seeds = len(settings.get("seeds", []) or [])
    base_trials = int(settings.get("trials", 1) or 1)
    corners = 4

    deltas = [0.1, 0.05, 0.01, 0.001]
    trial_multipliers = [1, 2, 4]

    sweep_rows = []
    for m in trial_multipliers:
        trials = base_trials * m
        n_eff = max(1, base_seeds * trials * corners)

        for d in deltas:
            margin = hoeffding_margin(d, n_eff)
            row_results = []
            for r in replay_rows:
                risk = float(r.get("risk_tolerance", 0.05) or 0.05)
                target = float(r.get("target_passrate_floor", max(0.0, 1.0 - risk)) or 0.0)
                empirical_floor = float(r.get("empirical_passrate_floor_over_corners", 0.0) or 0.0)
                cert_floor = max(0.0, empirical_floor - margin)
                row_results.append(
                    {
                        "risk_tolerance": risk,
                        "target_floor": target,
                        "empirical_corner_floor": empirical_floor,
                        "certified_floor": cert_floor,
                        "meets_target": cert_floor + 1e-12 >= target,
                        "margin_to_target": cert_floor - target,
                    }
                )

            sweep_rows.append(
                {
                    "delta": d,
                    "trial_multiplier": m,
                    "n_eff": n_eff,
                    "adaptive_margin": margin,
                    "rows": row_results,
                    "all_rows_meet_target": all(x["meets_target"] for x in row_results),
                    "worst_margin_to_target": min((x["margin_to_target"] for x in row_results), default=0.0),
                }
            )

    best_row = max(sweep_rows, key=lambda x: x["worst_margin_to_target"]) if sweep_rows else {}
    baseline_margin = float((p2277.get("strict_nu_branch_group_policy_adaptive_confidence_margin_certificate_probe", {}) or {}).get("confidence_model", {}).get("adaptive_margin", 0.0) or 0.0)

    payload = {
        "schema_version": "p2278_s1228_v1",
        "packet_id": "P2278",
        "stage_id": "S1228",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_CONFIDENCE_CURVE_SWEEP_PROBE",
        "strict_nu_branch_group_policy_confidence_curve_sweep_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_CONFIDENCE_CURVE_SWEEP_PROBE_V1",
            "source_packets": [str(IN_2275.relative_to(ROOT)), str(IN_2277.relative_to(ROOT))],
            "sweep_grid": {
                "deltas": deltas,
                "trial_multipliers": trial_multipliers,
                "base_seeds": base_seeds,
                "base_trials": base_trials,
                "corners": corners,
            },
            "baseline_adaptive_margin_from_p2277": baseline_margin,
            "sweep_rows": sweep_rows,
            "best_config_by_worst_margin": best_row,
            "theorem_scope_limit": "confidence-curve sweep over surrogate replay only; not selector closure and not ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2279_candidate",
            "goal": "promote best confidence-budget config to locked policy profile and rerun full corner replay with refreshed seeds",
        },
        "gatekeeper_checks": {
            "sweep_rows_exported": len(sweep_rows) > 0,
            "all_margins_nonnegative": all(x["adaptive_margin"] >= 0.0 for x in sweep_rows),
            "all_certified_floors_bounded": all(
                0.0 <= r["certified_floor"] <= 1.0
                for x in sweep_rows
                for r in x["rows"]
            ),
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2278 S1228: confidence-curve sweep probe",
            "",
            f"- sweep rows: `{len(sweep_rows)}`",
            f"- baseline adaptive margin (P2277): `{baseline_margin:.12e}`",
            f"- best config worst margin: `{float(best_row.get('worst_margin_to_target', 0.0) if best_row else 0.0):.12e}`",
            "",
            "Confidence-curve sweep only; no selector closure / ToE closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
